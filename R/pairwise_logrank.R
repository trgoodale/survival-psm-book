# Copyright (c) 2026 Tom Goodale. SPDX-License-Identifier: MIT
# pairwise_logrank() ---------------------------------------------------------
#
# Follow-up comparisons after an omnibus log-rank test with 3+ groups.
# Runs a two-group log-rank test (survdiff) for every pair of groups, using
# only the patients in those two groups, and adjusts the p-values for
# multiple testing (Holm by default; see ?p.adjust for other methods).
#
#   formula : Surv(time, status) ~ group   (one grouping variable)
#   data    : a data frame
#
# Returns one row per pair: the two groups, n, chi-square, raw p, adjusted p.

pairwise_logrank <- function(formula, data, p.adjust.method = "holm") {
  grp <- all.vars(formula)[length(all.vars(formula))]   # the grouping variable
  levs <- levels(factor(data[[grp]]))
  pairs <- combn(levs, 2, simplify = FALSE)

  res <- do.call(rbind, lapply(pairs, function(pr) {
    sub <- data[!is.na(data[[grp]]) & data[[grp]] %in% pr, ]
    sub[[grp]] <- droplevels(factor(sub[[grp]]))
    lr <- survival::survdiff(formula, data = sub)
    data.frame(group1 = pr[1], group2 = pr[2], n = sum(lr$n),
               chisq  = lr$chisq,
               p      = pchisq(lr$chisq, df = 1, lower.tail = FALSE))
  }))

  res$p_adj <- p.adjust(res$p, method = p.adjust.method)
  res
}
