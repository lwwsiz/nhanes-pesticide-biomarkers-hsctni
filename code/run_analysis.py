from __future__ import annotations

import json
import platform
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from analysis_utils import BASE_COVARIATES, fixed_effect_summary, survey_linear_regression


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "analysis_dataset.csv.gz"
OUT = ROOT / "results" / "reproduced"
CYCLES = ["2001-2002", "1999-2000"]
ASSAYS = {
    "Abbott hs-cTnI": "RZ_SSTNIA",
    "Siemens hs-cTnI": "RZ_SSTNIS",
    "Ortho hs-cTnI": "RZ_SSTNIO",
}


def cycle_models(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    primary = []
    weight_checks = []
    progression = []
    stages = [
        ("Model 1", []),
        ("Model 2", ["INDFMPIR", "LOG_COT"]),
        ("Model 3", BASE_COVARIATES),
    ]
    for cycle in CYCLES:
        frame = data.loc[data["cycle"].eq(cycle)].copy()
        result = survey_linear_regression(frame, "CTNI_COMPOSITE_Z")
        primary.append({"cycle": cycle, **result})
        for weight, label in [
            ("WTSPP2YR", "pesticide subsample weight"),
            ("WTSSCB2Y", "stored-specimen weight"),
        ]:
            check = survey_linear_regression(frame, "CTNI_COMPOSITE_Z", weight=weight)
            weight_checks.append(
                {"cycle": cycle, "weight_definition": label, "weight_variable": weight, **check}
            )
        for label, covariates in stages:
            model = survey_linear_regression(
                frame, "CTNI_COMPOSITE_Z", covariates=covariates
            )
            progression.append({"cycle": cycle, "model": label, **model})
    pooled = fixed_effect_summary(primary)
    primary.append(
        {
            "cycle": "descriptive pooled estimate",
            "n": int(sum(row["n"] for row in primary)),
            "design_df": np.nan,
            **pooled,
        }
    )
    return pd.DataFrame(primary), pd.DataFrame(weight_checks), pd.DataFrame(progression)


def assay_models(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cycle in CYCLES:
        frame = data.loc[data["cycle"].eq(cycle)].copy()
        for assay, outcome in ASSAYS.items():
            result = survey_linear_regression(frame, outcome)
            rows.append({"cycle": cycle, "assay": assay, **result})
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(DATA)
    primary, weights, progression = cycle_models(data)
    assays = assay_models(data)

    primary.to_csv(OUT / "primary_cycle_results.csv", index=False)
    weights.to_csv(OUT / "official_weight_sensitivity.csv", index=False)
    progression.to_csv(OUT / "sequential_adjustment_models.csv", index=False)
    assays.to_csv(OUT / "assay_specific_results.csv", index=False)

    summary = {
        "dataset_rows": int(len(data)),
        "cycles": data["cycle"].value_counts().to_dict(),
        "primary_results": primary.to_dict(orient="records"),
        "notes": [
            "Cycle-specific models use Taylor-linearized PSU-level sandwich variance.",
            "The pooled row is a descriptive inverse-variance summary, not a meta-analysis.",
            "The two official weights are evaluated separately and are never multiplied.",
        ],
    }
    (OUT / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    environment = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }
    (OUT / "environment_versions.json").write_text(
        json.dumps(environment, indent=2), encoding="utf-8"
    )
    print(primary.to_string(index=False))


if __name__ == "__main__":
    main()
