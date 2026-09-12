"""Input-failure checks; run after the core analysis."""
from pathlib import Path
import math
import numpy as np
import pandas as pd
from analysis_utils import survey_linear_regression, _invert_full_rank, _t_critical

root=Path(__file__).resolve().parents[1]
d=pd.read_csv(root/'data/analysis_dataset.csv.gz')
d=d.loc[d.cycle.eq('2001-2002')].copy()
def fails(label, call, expected):
    try: call()
    except expected: print('PASS',label)
    else: raise AssertionError(label+' did not reject invalid input')
for v in [0,-1,np.nan,np.inf]:
    bad=d.copy();bad.loc[bad.index[0],'WTSPP2YR']=v
    fails('invalid weight '+str(v),lambda:survey_linear_regression(bad,'CTNI_COMPOSITE_Z'),ValueError)
bad=d.copy();bad['SDMVPSU']=1
fails('singleton strata',lambda:survey_linear_regression(bad,'CTNI_COMPOSITE_Z'),ValueError)
fails('empty data',lambda:survey_linear_regression(d.iloc[:0],'CTNI_COMPOSITE_Z'),ValueError)
fails('rank deficiency',lambda:_invert_full_rank(np.ones((3,3))),np.linalg.LinAlgError)
fails('invalid df',lambda:_t_critical(0),ValueError)
fails('invalid alpha',lambda:_t_critical(15,0),ValueError)
fails('outcome equals exposure',lambda:survey_linear_regression(d,'MIX_PRIMARY_Z'),ValueError)
fails('design name collision',lambda:survey_linear_regression(d,'CTNI_COMPOSITE_Z',covariates=['INDFMPIR','INDFMPIR']),ValueError)
small=d.groupby(['SDMVSTRA','SDMVPSU']).head(1)
first_stratum=small.SDMVSTRA.iloc[0]
small=small[small.SDMVSTRA.eq(first_stratum)]
fails('insufficient observations',lambda:survey_linear_regression(small,'CTNI_COMPOSITE_Z'),ValueError)
bad=d.copy();bad['MIX_PRIMARY_Z']=1.0
fails('constant exposure',lambda:survey_linear_regression(bad,'CTNI_COMPOSITE_Z'),np.linalg.LinAlgError)
assert abs(_t_critical(1)-1/math.tan(math.pi*.025))<1e-9
base=survey_linear_regression(d,'CTNI_COMPOSITE_Z')
scaled=d.copy();scaled['WTSPP2YR']*=7
alt=survey_linear_regression(scaled,'CTNI_COMPOSITE_Z')
for key in ['beta','se','p','ci_low','ci_high']:assert abs(base[key]-alt[key])<1e-10,key
print('PASS df=1 critical value and constant-weight-scale invariance')
