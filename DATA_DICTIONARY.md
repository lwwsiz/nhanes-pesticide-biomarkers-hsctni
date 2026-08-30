# Data dictionary

The derived analytic dataset contains one row per participant and 50 columns. Original NHANES codebooks remain the authoritative source for laboratory units, category labels, and survey-design definitions.

## Identifiers and survey design

| Variable | Meaning |
|---|---|
| `SEQN` | Public NHANES participant identifier, stored as an integer. |
| `cycle` | NHANES survey cycle: `2001-2002` or `1999-2000`. |
| `WTSPP2YR` | Two-year urinary-pesticide subsample weight used in the primary models. |
| `WTSSCB2Y` | Two-year stored-specimen weight used in sensitivity analysis. |
| `WT4YR` | Combined four-year weight used for the cycle-interaction analysis. |
| `SDMVSTRA` | Masked variance pseudo-stratum. |
| `SDMVPSU` | Masked variance pseudo-primary sampling unit. |

## Demographic and clinical covariates

| Variable | Meaning or coding |
|---|---|
| `RIDAGEYR` | Age in years. |
| `RIAGENDR` | NHANES sex code (`1` male, `2` female). |
| `RIDRETH1` | NHANES race/ethnicity category code. |
| `INDFMPIR` | Family poverty-income ratio. |
| `BMXBMI` | Body mass index, kg/m2. |
| `MEAN_SBP` | Mean available systolic blood-pressure measurements, mm Hg. |
| `LOG_COT` | Natural logarithm of serum cotinine. |
| `LBDSCR` | Serum creatinine, mg/dL. |
| `EGFR` | 2021 CKD-EPI creatinine-based estimated glomerular filtration rate. |
| `LOG_UACR` | Natural logarithm of urinary albumin-to-creatinine ratio. |
| `PREVALENT_CVD` | Indicator for self-reported physician-diagnosed heart failure, coronary disease, angina, myocardial infarction, or stroke. |
| `URXUCR` | Urinary creatinine concentration, mg/dL. |
| `LOG_UCR` | Natural logarithm of urinary creatinine. |

## Urinary pesticide-related biomarkers

| Concentration | Detection flag | Interpretation |
|---|---|---|
| `URX24D` | `URD24DLC` | Urinary 2,4-dichlorophenoxyacetic acid (2,4-D) and its NHANES lower-limit flag. |
| `URXCPM` | `URDCPMLC` | Urinary 3,5,6-trichloro-2-pyridinol (TCPy) and its lower-limit flag. |
| `URXPAR` | `URDPARLC` | Urinary para-nitrophenol/4-nitrophenol and its lower-limit flag. |
| `URXOPM` | `URDOPMLC` | Urinary 3-phenoxybenzoic acid (3-PBA) and its lower-limit flag. |

The concentration units and flag definitions follow the corresponding NHANES laboratory documentation. Exact SAS transport zero sentinels were normalized to numeric zero in the released analytic dataset.

## Cardiac biomarkers

| Concentration | Lower-measurement flag | Interpretation |
|---|---|---|
| `SSTNIA` | `SSTNIALC` | Abbott high-sensitivity cardiac troponin I. |
| `SSTNIS` | `SSTNISLC` | Siemens high-sensitivity cardiac troponin I. |
| `SSTNIO` | `SSTNIOLC` | Ortho high-sensitivity cardiac troponin I. |
| `SSTNT` | `SSTNTLC` | Roche high-sensitivity cardiac troponin T comparator. |
| `SSBNP` | Not applicable | B-type natriuretic peptide comparator. |

`RZ_SSTNIA`, `RZ_SSTNIS`, `RZ_SSTNIO`, `RZ_SSTNT`, and `RZ_BNP` are within-cycle inverse-normal rank scores for the corresponding assays.

## Derived exposure and outcome scores

| Variable | Derivation |
|---|---|
| `PESTICIDE_MIX` | Equal-weight mean of the four within-cycle inverse-normal pesticide-biomarker rank scores. |
| `MIX_PRIMARY_Z` | Standardized primary mixture score. |
| `MIX_QUARTILE` | Cycle-specific survey-weighted quartile category (`1`-`4`). |
| `MIX_QUARTILE_Z` | Standardized quartile-coded mixture score. |
| `MIX_CREAT_RATIO_Z` | Standardized mixture score based on creatinine-ratio encodings. |
| `MIX_PCA1_Z` | Standardized first principal-component exposure score. |
| `CTNI_COMPOSITE` | Equal-weight mean of `RZ_SSTNIA`, `RZ_SSTNIS`, and `RZ_SSTNIO`. |
| `CTNI_COMPOSITE_Z` | Standardized three-assay hs-cTnI composite used as the primary outcome. |

## Missing values and zero handling

- Blank CSV fields represent missing values.
- The SAS/XPT numeric-zero sentinel `5.397605346934028e-79` was replaced by exact numeric zero for public release.
- The analysis code removes rows with missing variables required by a given model.
- No restricted-use geography or linked-mortality fields are included.
