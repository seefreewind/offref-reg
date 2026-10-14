from pathlib import Path
import sys,json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];sys.path[:0]=[str(R/'00_protocol'),str(R/'01_model_metadata')]
from audit_common import *
from borzoi_runtime import unsquash
assert_freezes();full=pd.read_csv(R/'02_track_harmonization/TRACK_METADATA_FULL.tsv',sep='\t').fillna('');targets=pd.read_csv(P/'16_phase1r_rescue/03_borzoi/targets_human.txt',sep='\t',index_col=0).fillna('');truth=pd.read_csv(P/'17_phase1r2/03_gnrs_truth/GNRS_CLASSB_PRIMARY33.tsv',sep='\t').set_index('OFFREF_ID');jobs=[j for j in readjson(OLD/'00_protocol/JOBS_FROZEN.json') if j['role']=='SENTINEL' and j['panel']=='BIOLOGICAL_TRUTH'];rows=[]
for j in jobs:
 for m in ['AlphaGenome','Borzoi']:
  old=readjson(OLD/'10_statistics'/m/(j['job_ID']+'.json'));W=1048576 if m=='AlphaGenome' else 524288;start=j['anchor0']-W//2;q=dict(j,W=W,start0=start,strict=False);fixed=domain(j['anchor0'],len(j['insertion']),j['anchor0']-262144,524288,163840,32,False);iv=regions(q,m,fixed)
  for method,key in [('original_strict_common','common_RNA'),('original_native','native_RNA')]:rows.append({'OFFREF_ID':j['job_ID'],'gene':j['gene'],'beta':float(truth.loc[j['job_ID'],'beta_presence']),'model':m,'method':method,'delta':old[key].get('delta'),'mapped_bp':old[key]['mapped_bp'],'status':old[key]['status']})
  row={'OFFREF_ID':j['job_ID'],'gene':j['gene'],'beta':float(truth.loc[j['job_ID'],'beta_presence']),'model':m,'method':'fractional_common_homology','delta':None,'mapped_bp':sum(b-a for a,b,o in iv),'status':'NO_COMMON_EXON'}
  if iv:
   results={'status':'PASS','strand':j['strand']};md=None
   for allele in ['REF','ALT']:
    cid=old['cache_'+allele];meta=readjson(OLD/'09_predictions'/m/(cid+'.json'))['metadata'];md=[]
    for tr in meta:
     idx=int(tr['track_index']) if m=='Borzoi' else int(full[(full.model==m)&(full.track_name==tr['name'])&(full.strand==tr['strand'])].iloc[0].track_ID);t=full[(full.model==m)&(full.track_ID.astype(int)==idx)].iloc[0].to_dict()
     if m=='Borzoi':t.update(targets.loc[idx].to_dict())
     md.append(t)
    with np.load(OLD/'09_predictions'/m/(cid+'.npz')) as z:v=z['rna'];res=int(z['resolution'])
    if m=='Borzoi':v=unsquash(v,md)
    results[allele]=summarize(v,iv,start,163840 if m=='Borzoi' else 0,res,allele);del v
   results['metadata']=md;row['delta']=delta_score(results);row['status']='PASS'
  rows.append(row)
pd.DataFrame(rows).to_csv(R/'08_truth_coverage/GNRS_SCORING_METHOD_SENSITIVITIES.tsv',sep='\t',index=False);print('CACHED_GNRS_SCORING_SENSITIVITIES_COMPLETE',flush=True)
