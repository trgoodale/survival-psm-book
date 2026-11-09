// Copyright (c) 2026 Tom Goodale. SPDX-License-Identifier: MIT
// extras.js — simulation and estimation helpers for the apps in Chapters 4,
// 9 and 10. Like survtools.js, these are teaching implementations, written
// for clarity. Use R for real analyses.

import { makeRandom, kmFit, coxFit, logisticFit, nnMatch, smd, pooledSD, mean, expit, invert } from "./survtools.js";

// ---------------------------------------------------------------------------
// Cox model with start-stop (counting-process) data and one covariate.
// Risk set at an event time t: rows with start < t <= stop.
// ---------------------------------------------------------------------------
export function coxTV(start, stop, event, x, maxIter = 60) {
  const n = stop.length;
  const evTimes = [...new Set(stop.filter((_, i) => event[i]))].sort((a, b) => a - b);
  let beta = 0, info = 0;
  for (let it = 0; it < maxIter; it++) {
    let U = 0, I = 0;
    for (const t of evTimes) {
      let S0 = 0, S1 = 0, S2 = 0, d = 0, dx = 0;
      for (let i = 0; i < n; i++) {
        if (start[i] < t && stop[i] >= t) {
          const w = Math.exp(beta * x[i]);
          S0 += w; S1 += w * x[i]; S2 += w * x[i] * x[i];
        }
        if (stop[i] === t && event[i]) { d++; dx += x[i]; }
      }
      if (S0 === 0) continue;
      U += dx - d * S1 / S0;
      I += d * (S2 / S0 - (S1 / S0) ** 2);
    }
    if (!(I > 0)) break;
    // cap the Newton step, so near-separated data can't make it overshoot
    const step = Math.max(-1, Math.min(1, U / I));
    beta += step; info = I;
    if (Math.abs(step) < 1e-9) break;
  }
  const se = 1 / Math.sqrt(info);
  return { hr: Math.exp(beta), lower: Math.exp(beta - 1.96 * se), upper: Math.exp(beta + 1.96 * se), coef: beta, se };
}

// ---------------------------------------------------------------------------
// Immortal time: a proportion pTreat of patients would start treatment at a
// random time in [0, maxWait], if still alive; the rest never would. Before treatment the hazard is `rate`; after
// starting, it is rate * hr. Follow-up ends at `fu` years.
// ---------------------------------------------------------------------------
export function simImmortal({ n, rate, hr, maxWait, fu, seed, pTreat = 1 }) {
  const r = makeRandom(seed);
  return Array.from({ length: n }, (_, i) => {
    const willTreat = r.unif() < pTreat;    // would this patient ever be offered treatment?
    const s = willTreat ? r.unif() * maxWait : Infinity;   // when treatment would start
    const E = -Math.log(r.unif());          // cumulative hazard at death
    let death;
    if (E < rate * s) death = E / rate;     // dies before treatment is due
    else death = s + (E - rate * s) / (rate * hr);
    const time = Math.min(death, fu);
    const status = death <= fu ? 1 : 0;
    const treated = s < time ? 1 : 0;       // started treatment during follow-up
    return { id: i + 1, start: s, time, status, treated };
  });
}

export function immortalEstimates(d, landmark) {
  const zeros = d.map(() => 0);
  const naive = coxTV(zeros, d.map(p => p.time), d.map(p => p.status), d.map(p => p.treated));
  // landmark: alive and followed at the landmark; group = started by then
  const L = d.filter(p => p.time > landmark);
  let lm = { hr: NaN, lower: NaN, upper: NaN }, nL = L.length, nLt = 0;
  if (L.length > 10) {
    const g = L.map(p => (p.start < landmark ? 1 : 0));
    nLt = g.reduce((a, b) => a + b, 0);
    if (nLt > 0 && nLt < L.length)
      lm = coxTV(L.map(() => 0), L.map(p => p.time - landmark), L.map(p => p.status), g);
  }
  // time-varying: split treated patients at their start time
  const st = [], sp = [], ev = [], x = [];
  for (const p of d) {
    if (p.treated) {
      if (p.start > 0) { st.push(0); sp.push(p.start); ev.push(0); x.push(0); }
      st.push(p.start); sp.push(p.time); ev.push(p.status); x.push(1);
    } else { st.push(0); sp.push(p.time); ev.push(p.status); x.push(0); }
  }
  const tv = coxTV(st, sp, ev, x);
  return { naive, lm, tv, nL, nLt };
}

// ---------------------------------------------------------------------------
// Competing risks: event of interest at rate1, competing event at rate2,
// censoring uniform on [fu/3, fu].
// ---------------------------------------------------------------------------
export function simCompeting({ n, rate1, rate2, fu, seed }) {
  const r = makeRandom(seed);
  return Array.from({ length: n }, () => {
    const t1 = r.exp(rate1);
    const t2 = rate2 > 0 ? r.exp(rate2) : Infinity;
    const c = fu / 3 + r.unif() * (2 * fu / 3);
    const time = Math.min(t1, t2, c);
    const cause = time === t1 ? 1 : time === t2 ? 2 : 0;
    return { time, cause };
  });
}

// 1 - KM for cause 1 (cause 2 treated as censored), and the Aalen–Johansen
// cumulative incidence of cause 1. Returns step-function rows.
export function competingCurves(d) {
  const ord = [...d].sort((a, b) => a.time - b.time);
  let atRisk = ord.length, S_km = 1, S_all = 1, cif = 0, i = 0;
  const rows = [{ time: 0, oneMinusKM: 0, cif: 0 }];
  while (i < ord.length) {
    const t = ord[i].time;
    let d1 = 0, d2 = 0, c = 0;
    while (i < ord.length && ord[i].time === t) {
      if (ord[i].cause === 1) d1++; else if (ord[i].cause === 2) d2++; else c++;
      i++;
    }
    if (d1 + d2 > 0) {
      cif += S_all * d1 / atRisk;              // overall survival just before t
      S_km *= 1 - d1 / atRisk;
      S_all *= 1 - (d1 + d2) / atRisk;
      rows.push({ time: t, oneMinusKM: 1 - S_km, cif });
    }
    atRisk -= d1 + d2 + c;
  }
  return rows;
}

export function stepAt(rows, key, t) {
  let v = rows[0][key];
  for (const r of rows) { if (r.time <= t) v = r[key]; else break; }
  return v;
}

// ---------------------------------------------------------------------------
// Weighting app: two confounders (age, severity), a continuous outcome, as
// in the Chapter 5 matching app.
// ---------------------------------------------------------------------------
export function simConfounded({ n, gamma, beta, tau, seed }) {
  const r = makeRandom(seed);
  return Array.from({ length: n }, (_, i) => {
    const age = r.norm(), sev = r.norm();
    const t = r.bern(expit(-0.4 + gamma * (age + sev)));
    const y = 50 + tau * t - 4 * beta * (age + sev) + r.norm(0, 5);
    return { id: i, age, sev, treat: t, y };
  });
}

const wmean = (x, w) => x.reduce((s, v, i) => s + v * w[i], 0) / w.reduce((a, b) => a + b, 0);
export const ess = (w) => (w.reduce((a, b) => a + b, 0) ** 2) / w.reduce((a, b) => a + b * b, 0);

export function weightingResults(d, { sevMeasured, cap, caliper }) {
  const ps = logisticFit(d.map(p => p.treat), d.map(p => sevMeasured ? [p.age, p.sev] : [p.age])).fitted;
  let w = d.map((p, i) => p.treat ? 1 / ps[i] : 1 / (1 - ps[i]));
  if (cap < Infinity) w = w.map(v => Math.min(v, cap));
  const T = d.map((p, i) => i).filter(i => d[i].treat === 1);
  const C = d.map((p, i) => i).filter(i => d[i].treat === 0);
  const naive = mean(T.map(i => d[i].y)) - mean(C.map(i => d[i].y));
  const iptw = wmean(T.map(i => d[i].y), T.map(i => w[i])) - wmean(C.map(i => d[i].y), C.map(i => w[i]));
  const pairs = nnMatch(ps, d.map(p => p.treat), caliper);
  const matched = pairs.length ? mean(pairs.map(([t, c]) => d[t].y - d[c].y)) : NaN;
  const balance = [];
  for (const [v, label] of [["age", "Age"], ["sev", "Severity"]]) {
    const x1 = T.map(i => d[i][v]), x0 = C.map(i => d[i][v]);
    const sd = pooledSD(x1, x0);
    balance.push({ variable: label, sample: "Before", smd: smd(x1, x0, sd) });
    if (pairs.length > 1) balance.push({ variable: label, sample: "Matched", smd: smd(pairs.map(p => d[p[0]][v]), pairs.map(p => d[p[1]][v]), sd) });
    balance.push({ variable: label, sample: "Weighted", smd: (wmean(x1, T.map(i => w[i])) - wmean(x0, C.map(i => w[i]))) / sd });
  }
  return {
    ps, w, naive, iptw, matched, nPairs: pairs.length, nT: T.length, nC: C.length,
    maxW: Math.max(...w), essT: ess(T.map(i => w[i])), essC: ess(C.map(i => w[i])), balance,
  };
}

// ---------------------------------------------------------------------------
// RMST app: a two-arm trial whose effect switches at `cut` years.
// ---------------------------------------------------------------------------
export function simSwitch({ n, h0, hrEarly, hrLate, cut, seed }) {
  const r = makeRandom(seed);
  return Array.from({ length: n }, (_, i) => {
    const arm = i < n / 2 ? 0 : 1;
    const E = -Math.log(r.unif());
    let t;
    if (arm === 0) t = E / h0;
    else t = E < h0 * hrEarly * cut ? E / (h0 * hrEarly) : cut + (E - h0 * hrEarly * cut) / (h0 * hrLate);
    const c = 2 + 3 * r.unif();
    return { arm, time: Math.min(t, c), status: t <= c ? 1 : 0 };
  });
}

// Area under a KM step function from 0 to tau
export function rmst(fit, tau) {
  let area = 0, prevT = 0, prevS = 1;
  for (const r of fit) {
    if (r.time > tau) break;
    area += prevS * (r.time - prevT);
    prevT = r.time; prevS = r.surv;
  }
  return area + prevS * (tau - prevT);
}
