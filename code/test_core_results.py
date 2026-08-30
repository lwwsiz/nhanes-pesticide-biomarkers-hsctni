from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "reproduced" / "primary_cycle_results.csv"

EXPECTED = {
    "2001-2002": {"n": 740, "beta": 0.138693, "p": 0.002429},
    "1999-2000": {"n": 556, "beta": 0.097072, "p": 0.036228},
    "descriptive pooled estimate": {"n": 1296, "beta": 0.119849, "p": 0.000021},
}


def main() -> None:
    frame = pd.read_csv(RESULTS).set_index("cycle")
    for label, expected in EXPECTED.items():
        observed = frame.loc[label]
        assert int(observed["n"]) == expected["n"]
        assert abs(float(observed["beta"]) - expected["beta"]) < 5e-7
        assert abs(float(observed["p"]) - expected["p"]) < 5e-7
    print("Core result checks passed.")


if __name__ == "__main__":
    main()
