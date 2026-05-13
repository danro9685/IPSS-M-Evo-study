# IPSS-M-Evo calculator

This folder contains a minimal R implementation of the IPSS-M-Evo score calculator.

## Files

- `ipss_m_evo_calculator.R`: defines `calculate_ipss_m_evo()`, a single function that computes the IPSS-M-Evo score, risk category, and missing-data uncertainty.

## Input

The function expects a `data.frame` or matrix with one row per patient and the following columns:

| Column | Meaning |
|---|---|
| `ipss_m` | Continuous IPSS-M score |
| `age` | Age in years |
| `asxl1_kras` | ASXL1→KRAS evolutionary route |
| `srsf2_nras` | SRSF2→NRAS evolutionary route |
| `nras_runx1` | NRAS/RUNX1 co-occurrence |
| `atrx` | ATRX mutation |
| `jak2` | JAK2 mutation |

The five molecular/evolutionary variables must be coded as:

- `1`: present
- `0`: absent
- `NA`: missing

## Usage

```r
source("ipss_m_evo_calculator.R")

patients <- data.frame(
  ipss_m = c(1.2, 2.4),
  age = c(70, 63),
  asxl1_kras = c(1, 0),
  srsf2_nras = c(0, NA),
  nras_runx1 = c(0, 1),
  atrx = c(NA, 0),
  jak2 = c(0, 1)
)

result <- calculate_ipss_m_evo(patients)
print(result)
```

## Output

The input table is returned with these additional columns:

| Column | Meaning |
|---|---|
| `ipss_m_evo_score` | Score using average imputation for missing molecular variables |
| `ipss_m_evo_best` | Best-case score under deterministic imputation |
| `ipss_m_evo_worst` | Worst-case score under deterministic imputation |
| `ipss_m_evo_range` | Difference between worst and best scores |
| `ipss_m_evo_reliable` | `TRUE` if `ipss_m_evo_range <= 1.0` |
| `ipss_m_evo_risk` | Risk category from the average-imputed score |

## Notes

The script assumes that `ipss_m` and `age` are available. Missingness is handled only for the five molecular/evolutionary IPSS-M-Evo variables. The evolutionary-route variables (`asxl1_kras`, `srsf2_nras`) and the co-occurrence variable (`nras_runx1`) must already be encoded in the input table.
