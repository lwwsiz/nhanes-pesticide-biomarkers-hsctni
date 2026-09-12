# Tested environment

The original release recorded Windows, Python 3.11.15, NumPy 2.4.6, and pandas 2.3.3. Its dependency pins remain in `requirements.txt`.

The current local candidate was rerun on Windows with Python 3.12.14, NumPy 2.3.5, and pandas 3.0.1. All four core result CSV files were unchanged. `results/reproduced/environment_versions.json` records this latest execution environment; it is regenerated on each run.

The latest check did not perform a fresh installation of the original pinned environment. Successful execution in the environment above should not be interpreted as testing every dependency combination. SciPy was used separately for independent numerical cross-checking, but is not required by the released analysis or tests.

Install compatible dependencies with:

```powershell
python -m pip install -r requirements.txt
```
