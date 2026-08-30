# Data notice

The released analytic dataset was derived from deidentified NHANES public-use files. `SEQN` is the public NHANES participant identifier and is retained to support traceability to the source records. No restricted-use variables or linked-mortality records are included.

This repository does not redistribute original CDC transport files. Users should consult the NHANES documentation for survey design, laboratory methods, analytic weights, disclosure rules, and citation requirements.

The original SAS transport reader represents numeric zero as approximately `5.397605346934028e-79` in several fields. Exact occurrences of this transport sentinel were normalized to numeric zero in the released analytic dataset. This operation does not change the rank-based analysis. `SEQN` is stored as an integer.

The derived dataset is provided without warranty for reproducibility. It must not be used for participant re-identification or clinical decision-making.
