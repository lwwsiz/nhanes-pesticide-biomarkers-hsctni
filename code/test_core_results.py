from __future__ import annotations

from pathlib import Path
import json
import math

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
    expected_all=json.loads((Path(__file__).parent/'expected_core_results.json').read_text())
    specifications={
        'primary_cycle_results.csv':lambda r:'descriptive pooled' if r['cycle']=='descriptive pooled estimate' else 'primary '+r['cycle'],
        'official_weight_sensitivity.csv':lambda r:'weight '+r['cycle']+' '+r['weight_variable'],
        'sequential_adjustment_models.csv':lambda r:'stage '+r['cycle']+' '+r['model'],
        'assay_specific_results.csv':lambda r:'assay '+r['cycle']+' '+r['assay'],
    }
    observed_keys=set()
    for filename,key_fn in specifications.items():
        table=pd.read_csv(RESULTS.parent/filename)
        for _,row in table.iterrows():
            key=key_fn(row)
            assert key not in observed_keys,'duplicate result '+key
            observed_keys.add(key)
            for field,value in expected_all[key].items():
                observed=float(row[field])
                assert math.isfinite(observed),(key,field)
                if field in ('n','design_df'):
                    assert observed==value,(key,field)
                else:
                    assert abs(observed-value)<1e-8,(key,field,observed,value)
    assert observed_keys==set(expected_all),'missing or extra core result rows'
    print('Extended checks passed: 19 rows, n/df/beta/SE/P/CI.')


if __name__ == "__main__":
    main()
