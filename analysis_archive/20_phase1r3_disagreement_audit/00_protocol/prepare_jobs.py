from pathlib import Path
import sys,json,hashlib
import pandas as pd
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'00_protocol'));from audit_common import *
assert readjson(R/'15_freeze/PREVIOUS_FREEZES_RECHECK.json')['overall']=='PASS'
oldjobs=[j for j in readjson(OLD/'00_protocol/JOBS_FROZEN.json') if j['role']=='SENTINEL'];core=[j for j in oldjobs if j['panel']=='UNBIASED_CLASSB_CORE'];truth=[j for j in oldjobs if j['panel']=='BIOLOGICAL_TRUTH'];subset=set(pd.read_csv(R/'04_context_robustness/SUBSET100_FROZEN.tsv',sep='\t').OFFREF_ID);ref=readjson(R/'07_reference_benchmark/REFERENCE_PANEL_FROZEN.json');jobs=[]
for model in ['AlphaGenome','Borzoi']:
 W=1048576 if model=='AlphaGenome' else 524288;crop=0 if model=='AlphaGenome' else 163840
 for j in core:
  q=dict(j,model=model,W=W,start0=j['anchor0']-W//2,member='ALL_FOLDS' if model=='AlphaGenome' else '0',stage='GLOBAL_CORE',track_scope='GLOBAL',strict=True)
  q['fixed_domains']=domain(j['anchor0'],len(j['insertion']),j['anchor0']-262144,524288,163840,32);jobs.append(q)
 for j in core:
  if j['job_ID'] not in subset:continue
  A=j['anchor0'];L=len(j['insertion']);native=domain(A,L,A-262144,524288,163840,32);shared=intersect_domains(native,[domain(A,L,A-262144+k,524288,163840,32) for k in [-32768,32768]])
  for k in [0,-32768,32768]:jobs.append(dict(j,model=model,W=W,start0=A-W//2+k,member='ALL_FOLDS' if model=='AlphaGenome' else '0',stage='PLACEMENT',track_scope='LCL',strict=True,fixed_domains=shared,placement_shift=k))
  # Fixed homologous target support across all valid short/native contexts.
  context=intersect_domains(native,[domain(A,L,A-8192,16384,0,1)])
  for length in [16384,131072,W]:
   jobs.append(dict(j,model=model,W=length if model=='AlphaGenome' else W,start0=A-(length if model=='AlphaGenome' else W)//2,mask_length=length if model=='Borzoi' and length!=W else None,context_length=length,member='ALL_FOLDS' if model=='AlphaGenome' else '0',stage='CONTEXT',track_scope='LCL',strict=True,fixed_domains=context))
  for member in (['FOLD_0','FOLD_1','FOLD_2','FOLD_3'] if model=='AlphaGenome' else ['0','1','2','3']):jobs.append(dict(j,model=model,W=W,start0=A-W//2,member=member,stage='ENSEMBLE',track_scope='LCL',strict=True,fixed_domains=native))
 for j in truth:
  if j['exons']:
   center=int((min(j['anchor0'],j['exons'][0][0])+max(j['anchor0']+len(j['insertion']),j['exons'][-1][1]))//2)
  else:center=j['anchor0']
  start=center-W//2;fixed=domain(j['anchor0'],len(j['insertion']),center-262144,524288,163840,32)
  jobs.append(dict(j,model=model,W=W,start0=start,member='ALL_FOLDS' if model=='AlphaGenome' else '0',stage='GNRS_GENE_CENTER',track_scope='LCL',strict=True,fixed_domains=fixed,gene_center0=center))
 for j in ref:
  jobs.append(dict(j,model=model,W=W,start0=j['anchor0']-W//2,member='ALL_FOLDS' if model=='AlphaGenome' else '0',stage='REFERENCE',track_scope='LCL',strict=True,fixed_domains=domain(j['anchor0'],0,j['anchor0']-262144,524288,163840,32)))
for j in jobs:
 j['audit_job_ID']=hashlib.sha256(json.dumps(j,sort_keys=True).encode()).hexdigest()
p=R/'00_protocol/PREDICTION_JOBS_FROZEN.json';p.write_text(json.dumps(jobs,indent=2))
files=['00_protocol/PHASE1R3_PROTOCOL.md','00_protocol/audit_common.py','00_protocol/prepare_jobs.py','00_protocol/PREDICTION_JOBS_FROZEN.json','02_track_harmonization/TRACK_METADATA_FULL.tsv','02_track_harmonization/TRACK_MAPPING_FROZEN.json','04_context_robustness/SUBSET100_FROZEN.tsv','07_reference_benchmark/REFERENCE_PANEL_FROZEN.json']
(R/'15_freeze/PREDICTION_INPUT_SHA256.json').write_text(json.dumps([{'path':p,'sha256':sha(R/p)} for p in files],indent=2));print('ALL_PREDICTION_JOBS_FROZEN',len(jobs),flush=True)
