from pathlib import Path
import json,numpy as np
from statistics_core import boot,metric,interval
R=Path(__file__).resolve().parents[1];checks=[]
a=np.arange(1,21,dtype=float)
for sign in [1,-1]:
 z=metric(a,sign*a,20);assert abs(z['rho']-sign)<1e-12;assert z['proportion']==(1 if sign==1 else 0);assert abs(z['rho_CI95'][0]-sign)<1e-12;checks.append({'case':'analytic aligned/anti-aligned '+str(sign),'status':'PASS'})
z=metric([1,0,np.nan,-1],[1,1,2,1],4);assert z['finite_N']==3 and z['signed_N']==2 and z['concordant_N']==1;checks.append({'case':'finite and nonzero denominators separate','status':'PASS'})
b=np.sin(a)*3+a/10;x=boot(a,b,B=1000,seed=20261011);y=boot(a,b,B=1000,seed=20261012);assert x['rho']!=y['rho'];assert x==boot(a,b,B=1000,seed=20261011);checks.append({'case':'independent deterministic gap streams','status':'PASS'})
x=boot(a,b,B=100,cluster=np.repeat(np.arange(10),2));assert len(x['rho'])==100 and np.isfinite(x['rho']).all();checks.append({'case':'cluster bootstrap valid variable-count resampling','status':'PASS'})
z=metric([np.nan],[np.nan],1);assert z['finite_N']==0 and z['rho'] is None and z['rho_CI95']==[None,None];checks.append({'case':'empty scores not represented as zero performance','status':'PASS'})
(R/'11_statistics/ANALYTIC_STATISTICAL_QA.json').write_text(json.dumps({'status':'PASS','synthetic_inputs':'analytic scorer/statistics QA only;never experimental observations or figure source','checks':checks},indent=2));print('STATISTICAL_QA_PASS',flush=True)
