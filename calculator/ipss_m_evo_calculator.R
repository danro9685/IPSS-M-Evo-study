# Minimal IPSS-M-Evo calculator with missing-data handling
#
# Required columns:
#   ipss_m, age, asxl1_kras, srsf2_nras, nras_runx1, atrx, jak2
#
# Binary variables must be coded as:
#   1 = present
#   0 = absent
#   NA = missing
#
# Missing binary variables are handled by deterministic imputation:
#   average: missing = 0.5
#   best: missing = value that minimizes the score contribution
#   worst: missing = value that maximizes the score contribution
#
# The output includes:
#   ipss_m_evo_score: average-imputed score
#   ipss_m_evo_best: best-case score
#   ipss_m_evo_worst: worst-case score
#   ipss_m_evo_range: worst - best
#   ipss_m_evo_reliable: TRUE if range <= 1.0
#   ipss_m_evo_risk: risk category from the average-imputed score
calculate_ipss_m_evo <- function(x) {

  x <- as.data.frame(x)

  coef <- c(
    ipss_m = 0.62,
    age = 0.02,
    asxl1_kras = 0.55,
    srsf2_nras = 0.65,
    nras_runx1 = -0.805,
    atrx = 0.41,
    jak2 = 0.46
  )

  binary_vars <- c("asxl1_kras", "srsf2_nras", "nras_runx1", "atrx", "jak2")

  avg <- x
  best <- x
  worst <- x

  for (v in binary_vars) {
    avg[[v]][is.na(avg[[v]])] <- 0.5
    if (coef[v] >= 0) {
      best[[v]][is.na(best[[v]])] <- 0
      worst[[v]][is.na(worst[[v]])] <- 1
    } else {
      best[[v]][is.na(best[[v]])] <- 1
      worst[[v]][is.na(worst[[v]])] <- 0
    }
  }

  score <- function(d) {
    coef["ipss_m"] * d$ipss_m +
      coef["age"] * d$age +
      coef["asxl1_kras"] * d$asxl1_kras +
      coef["srsf2_nras"] * d$srsf2_nras +
      coef["nras_runx1"] * d$nras_runx1 +
      coef["atrx"] * d$atrx +
      coef["jak2"] * d$jak2
  }

  x$ipss_m_evo_score <- score(avg)
  x$ipss_m_evo_best <- score(best)
  x$ipss_m_evo_worst <- score(worst)
  x$ipss_m_evo_range <- x$ipss_m_evo_worst - x$ipss_m_evo_best
  x$ipss_m_evo_reliable <- x$ipss_m_evo_range <= 1.0

  x$ipss_m_evo_risk <- cut(
    x$ipss_m_evo_score,
    breaks = c(-Inf, 0.50, 1.00, 1.50, 2.00, 3.00, Inf),
    labels = c("Very Low", "Low", "Moderate-Low", "Moderate-High", "High", "Very High"),
    right = TRUE
  )

  x
}
