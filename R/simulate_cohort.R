# Copyright (c) 2026 Tom Goodale. SPDX-License-Identifier: MIT
# simulate_cohort.R ----------------------------------------------------------
#
# Simulates an *observational* cohort in which patients choose (or are given)
# a treatment non-randomly, and are then followed up until death or censoring.
#
# Because we simulate the data ourselves we KNOW the truth, so we can check
# whether each method recovers it. This is the best way to learn these tools.
#
#   * Older patients, diabetics, people with higher BMI and smokers are MORE
#     likely to be treated  -> confounding by indication.
#   * The same factors also raise the hazard of death.
#   * The true conditional hazard ratio for treatment is 0.70 (a 30% lower
#     hazard of death, holding the covariates fixed).
#
# Time is measured in years.

simulate_cohort <- function(n = 2000, hr_treat = 0.70, seed = 2026) {
  set.seed(seed)

  age      <- round(pmin(pmax(rnorm(n, 65, 10), 35), 95))
  female   <- rbinom(n, 1, 0.45)
  diabetes <- rbinom(n, 1, plogis(-1.4 + 0.04 * (age - 65)))
  bmi      <- round(rnorm(n, 28, 4.5), 1)
  smoker   <- rbinom(n, 1, 0.25)

  # --- Treatment assignment (the "propensity score" model) ------------------
  lp_treat <- -0.6 + 0.06 * (age - 65) + 1.0 * diabetes + 0.08 * (bmi - 28) -
    0.3 * female + 0.4 * smoker
  treated <- rbinom(n, 1, plogis(lp_treat))

  # --- Survival time: Weibull proportional hazards --------------------------
  lp_death <- log(hr_treat) * treated + 0.06 * (age - 65) + 0.6 * diabetes +
    0.03 * (bmi - 28) - 0.25 * female + 0.55 * smoker
  shape  <- 1.2
  lambda <- 0.05                       # baseline hazard scale (per year)
  t_event <- (-log(runif(n)) / (lambda * exp(lp_death)))^(1 / shape)

  # --- Censoring: staggered entry (administrative) + random loss to follow-up
  t_admin   <- runif(n, 3, 8)          # study ends 3-8 years after entry
  t_dropout <- rexp(n, rate = 0.04)
  t_cens    <- pmin(t_admin, t_dropout)

  time   <- pmin(t_event, t_cens)
  status <- as.integer(t_event <= t_cens)

  data.frame(
    id       = seq_len(n),
    age, female, diabetes, bmi, smoker,
    treated,
    time     = round(time, 3),
    status
  )
}
