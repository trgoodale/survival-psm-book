// Copyright (c) 2026 Tom Goodale. SPDX-License-Identifier: MIT
// survtools.js — small, dependency-free helpers used by the interactive
// (Observable JS) apps in this book. Everything here runs in the reader's
// browser, so the rendered book needs no R server.
//
// These are teaching implementations: correct, but written for clarity
// rather than speed or edge-case robustness. Use R for real analyses.

// ---------------------------------------------------------------------------
// Random numbers (seeded, so a given seed always gives the same "study")
// ---------------------------------------------------------------------------
function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export function makeRandom(seed = 1) {
  const u = mulberry32(seed);
  const unif = () => { let x; do { x = u(); } while (x === 0); return x; };
  return {
    unif,
    norm: (m = 0, s = 1) =>
      m + s * Math.sqrt(-2 * Math.log(unif())) * Math.cos(2 * Math.PI * unif()),
    exp: (rate) => -Math.log(unif()) / rate,
    bern: (p) => (u() < p ? 1 : 0),
  };
}

// ---------------------------------------------------------------------------
// Distribution helpers
// ---------------------------------------------------------------------------
// Standard normal CDF (Abramowitz & Stegun 7.1.26, |error| < 1.5e-7)
export function normCdf(z) {
  const t = 1 / (1 + 0.3275911 * Math.abs(z) / Math.SQRT2);
  const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t
    - 0.284496736) * t + 0.254829592) * t * Math.exp(-(z * z) / 2);
  return z >= 0 ? 0.5 * (1 + y) : 0.5 * (1 - y);
}
// Upper tail p-value of a chi-square with 1 df
export function chisq1P(x) { return 2 * (1 - normCdf(Math.sqrt(Math.max(x, 0)))); }

export const expit = (x) => 1 / (1 + Math.exp(-x));
export const logit = (p) => Math.log(p / (1 - p));

export function mean(a) { return a.reduce((s, v) => s + v, 0) / a.length; }
export function variance(a) {
  const m = mean(a);
  return a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1);
}

// ---------------------------------------------------------------------------
// Kaplan–Meier estimator with Greenwood variance and log(-log) CI
// (the same CI type that R's survfit() uses by default since survival 3.x)
// ---------------------------------------------------------------------------
export function kmFit(time, event) {
  const n = time.length;
  const z = 1.959964; // 95% CI
  const order = [...time.keys()].sort((a, b) => time[a] - time[b]);
  const rows = [{ time: 0, nrisk: n, nevent: 0, ncens: 0, surv: 1, lower: 1, upper: 1 }];
  let S = 1, gw = 0, atRisk = n, i = 0;
  while (i < n) {
    const t = time[order[i]];
    let d = 0, c = 0;
    while (i < n && time[order[i]] === t) { event[order[i]] ? d++ : c++; i++; }
    if (d > 0) {
      S *= 1 - d / atRisk;
      if (atRisk - d > 0) gw += d / (atRisk * (atRisk - d));
    }
    let lower = S, upper = S;
    if (S > 0 && S < 1) {
      const se = Math.sqrt(gw) / Math.abs(Math.log(S));
      lower = Math.pow(S, Math.exp(z * se));
      upper = Math.pow(S, Math.exp(-z * se));
    }
    rows.push({ time: t, nrisk: atRisk, nevent: d, ncens: c, surv: S, lower, upper });
    atRisk -= d + c;
  }
  return rows;
}

// Median survival: first time at which S(t) <= 0.5 (NaN if never reached)
export function kmMedian(fit) {
  const r = fit.find((r) => r.surv <= 0.5);
  return r ? r.time : NaN;
}

// Number at risk at time t
export function nAtRisk(time, t) { return time.filter((x) => x >= t).length; }

// ---------------------------------------------------------------------------
// Log-rank test for two groups (group coded 0/1)
// ---------------------------------------------------------------------------
export function logrank(time, event, group) {
  const ev = [...new Set(time.filter((_, i) => event[i]))].sort((a, b) => a - b);
  let O1 = 0, E1 = 0, V = 0, O0 = 0, E0 = 0;
  for (const t of ev) {
    let n = 0, n1 = 0, d = 0, d1 = 0;
    for (let i = 0; i < time.length; i++) {
      if (time[i] >= t) { n++; if (group[i]) n1++; }
      if (time[i] === t && event[i]) { d++; if (group[i]) d1++; }
    }
    const e1 = (d * n1) / n;
    O1 += d1; E1 += e1; O0 += d - d1; E0 += d - e1;
    if (n > 1) V += (d * (n1 / n) * (1 - n1 / n) * (n - d)) / (n - 1);
  }
  const chisq = V > 0 ? ((O1 - E1) ** 2) / V : 0;
  return { O1, E1, O0, E0, chisq, p: chisq1P(chisq) };
}

// ---------------------------------------------------------------------------
// Small linear-algebra helpers (Gaussian elimination)
// ---------------------------------------------------------------------------
export function invert(A) {
  const n = A.length;
  const M = A.map((row, i) => [...row, ...Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))]);
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
    [M[c], M[p]] = [M[p], M[c]];
    const piv = M[c][c];
    if (Math.abs(piv) < 1e-12) return null;
    for (let j = 0; j < 2 * n; j++) M[c][j] /= piv;
    for (let r = 0; r < n; r++) {
      if (r === c) continue;
      const f = M[r][c];
      for (let j = 0; j < 2 * n; j++) M[r][j] -= f * M[c][j];
    }
  }
  return M.map((row) => row.slice(n));
}
const matVec = (A, v) => A.map((row) => row.reduce((s, a, j) => s + a * v[j], 0));

// ---------------------------------------------------------------------------
// Cox proportional hazards model (Breslow ties) by Newton–Raphson.
// X is an array of rows (n x p). Returns coefficients, SEs, HRs and 95% CIs.
// ---------------------------------------------------------------------------
export function coxFit(time, event, X, maxIter = 25) {
  const n = time.length, p = X[0].length;
  const order = [...time.keys()].sort((a, b) => time[b] - time[a]); // descending
  let beta = new Array(p).fill(0), info = null, loglik = 0;
  for (let iter = 0; iter < maxIter; iter++) {
    const U = new Array(p).fill(0);
    const I = Array.from({ length: p }, () => new Array(p).fill(0));
    let S0 = 0; const S1 = new Array(p).fill(0);
    const S2 = Array.from({ length: p }, () => new Array(p).fill(0));
    loglik = 0;
    let k = 0;
    while (k < n) {
      const t = time[order[k]];
      const tied = [];
      // add everyone with this time to the risk set first (risk set = time >= t)
      while (k < n && time[order[k]] === t) {
        const i = order[k], x = X[i];
        const w = Math.exp(x.reduce((s, xj, j) => s + xj * beta[j], 0));
        S0 += w;
        for (let a = 0; a < p; a++) {
          S1[a] += w * x[a];
          for (let b = 0; b < p; b++) S2[a][b] += w * x[a] * x[b];
        }
        if (event[i]) tied.push(i);
        k++;
      }
      for (const i of tied) {
        const x = X[i];
        loglik += x.reduce((s, xj, j) => s + xj * beta[j], 0) - Math.log(S0);
        for (let a = 0; a < p; a++) {
          U[a] += x[a] - S1[a] / S0;
          for (let b = 0; b < p; b++)
            I[a][b] += S2[a][b] / S0 - (S1[a] * S1[b]) / (S0 * S0);
        }
      }
    }
    info = invert(I);
    if (!info) break;
    const step = matVec(info, U);
    beta = beta.map((b, j) => b + step[j]);
    if (Math.max(...step.map(Math.abs)) < 1e-9) break;
  }
  const se = info ? info.map((row, j) => Math.sqrt(row[j])) : beta.map(() => NaN);
  return beta.map((b, j) => ({
    coef: b, se: se[j], hr: Math.exp(b),
    lower: Math.exp(b - 1.96 * se[j]), upper: Math.exp(b + 1.96 * se[j]),
    p: chisq1P((b / se[j]) ** 2),
  }));
}

// ---------------------------------------------------------------------------
// Logistic regression by IRLS (an intercept is added automatically)
// ---------------------------------------------------------------------------
export function logisticFit(y, X, maxIter = 30) {
  const Z = X.map((row) => [1, ...row]);
  const p = Z[0].length;
  let beta = new Array(p).fill(0);
  for (let iter = 0; iter < maxIter; iter++) {
    const H = Array.from({ length: p }, () => new Array(p).fill(0));
    const g = new Array(p).fill(0);
    for (let i = 0; i < Z.length; i++) {
      const mu = expit(Z[i].reduce((s, z, j) => s + z * beta[j], 0));
      const w = mu * (1 - mu);
      for (let a = 0; a < p; a++) {
        g[a] += (y[i] - mu) * Z[i][a];
        for (let b = 0; b < p; b++) H[a][b] += w * Z[i][a] * Z[i][b];
      }
    }
    const Hi = invert(H);
    if (!Hi) break;
    const step = matVec(Hi, g);
    beta = beta.map((b, j) => b + step[j]);
    if (Math.max(...step.map(Math.abs)) < 1e-9) break;
  }
  const fitted = Z.map((z) => expit(z.reduce((s, zj, j) => s + zj * beta[j], 0)));
  return { beta, fitted };
}

// ---------------------------------------------------------------------------
// Greedy 1:1 nearest-neighbour matching on the logit of the propensity score,
// without replacement, with an optional caliper (in SDs of the logit PS).
// Treated units are matched in descending order of PS, as MatchIt's default
// ("largest") does. Returns an array of [treatedIndex, controlIndex] pairs.
// ---------------------------------------------------------------------------
export function nnMatch(ps, treat, caliperSD = Infinity) {
  const lp = ps.map(logit);
  const sd = Math.sqrt(variance(lp));
  const cal = caliperSD * sd;
  const treated = [...ps.keys()].filter((i) => treat[i] === 1).sort((a, b) => ps[b] - ps[a]);
  const used = new Set();
  const controls = [...ps.keys()].filter((i) => treat[i] === 0);
  const pairs = [];
  for (const t of treated) {
    let best = -1, bestD = Infinity;
    for (const c of controls) {
      if (used.has(c)) continue;
      const d = Math.abs(lp[t] - lp[c]);
      if (d < bestD) { bestD = d; best = c; }
    }
    if (best >= 0 && bestD <= cal) { pairs.push([t, best]); used.add(best); }
  }
  return pairs;
}

// Standardised mean difference, using a fixed denominator (the pooled SD in
// the *unmatched* sample) so that before/after values are comparable.
export function smd(x1, x0, sdPooled) {
  return (mean(x1) - mean(x0)) / sdPooled;
}
export function pooledSD(x1, x0) { return Math.sqrt((variance(x1) + variance(x0)) / 2); }

// ---------------------------------------------------------------------------
// Self-check quiz. questions = [{q, options: [...], answer: index, why}]
// Returns a DOM element. Feedback appears as soon as an option is chosen,
// and readers can change their answer as often as they like.
// ---------------------------------------------------------------------------
let quizCounter = 0;
export function quiz(questions) {
  const id = `quiz${++quizCounter}`;
  const root = document.createElement("div");
  root.className = "quiz";
  const score = document.createElement("p");
  score.className = "quiz-score";
  const correct = new Array(questions.length).fill(null);
  const updateScore = () => {
    const done = correct.filter((c) => c !== null).length;
    const right = correct.filter((c) => c === true).length;
    score.textContent = done === 0
      ? `${questions.length} questions. Choose an answer to see feedback.`
      : `Score: ${right} of ${questions.length} correct` +
        (done < questions.length ? ` (${questions.length - done} still to answer)` : right === questions.length ? ". All correct!" : ".");
  };
  questions.forEach((raw, qi) => {
    // Shuffle the options in a fixed (seeded) order so the correct answer
    // isn't always in the same position, but the order is stable on reload.
    let h = 0;
    for (const ch of raw.q) h = (Math.imul(h, 31) + ch.charCodeAt(0)) | 0;
    const r = mulberry32(h);
    const order = raw.options.map((_, i) => i);
    for (let i = order.length - 1; i > 0; i--) {
      const j = Math.floor(r() * (i + 1));
      [order[i], order[j]] = [order[j], order[i]];
    }
    const item = {
      ...raw,
      options: order.map((i) => raw.options[i]),
      hints: raw.hints ? order.map((i) => raw.hints[i]) : undefined,
      answer: order.indexOf(raw.answer),
    };
    const box = document.createElement("fieldset");
    box.className = "quiz-q";
    const legend = document.createElement("legend");
    legend.innerHTML = `<span class="quiz-num">${qi + 1}</span> ${item.q}`;
    box.appendChild(legend);
    const feedback = document.createElement("div");
    feedback.className = "quiz-feedback";
    feedback.setAttribute("aria-live", "polite");
    const labels = item.options.map((opt, oi) => {
      const label = document.createElement("label");
      label.className = "quiz-opt";
      const input = document.createElement("input");
      input.type = "radio";
      input.name = `${id}-${qi}`;
      input.addEventListener("change", (e) => {
        e.stopPropagation(); // don't trigger an OJS re-evaluation
        const ok = oi === item.answer;
        correct[qi] = ok;
        labels.forEach((l, li) => {
          l.classList.remove("is-right", "is-wrong");
          if (li === item.answer && ok) l.classList.add("is-right");
          if (li === oi && !ok) l.classList.add("is-wrong");
        });
        feedback.className = "quiz-feedback " + (ok ? "ok" : "no");
        feedback.innerHTML = (ok ? "<b>✓ Correct.</b> " : "<b>✗ Not quite. Try again.</b> ") +
          (ok ? item.why : (item.hints && item.hints[oi]) || "");
        updateScore();
      });
      label.appendChild(input);
      label.insertAdjacentHTML("beforeend", ` <span>${opt}</span>`);
      box.appendChild(label);
      return label;
    });
    box.appendChild(feedback);
    root.appendChild(box);
  });
  updateScore();
  root.appendChild(score);
  return root;
}
