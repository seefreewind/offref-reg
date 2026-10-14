from pathlib import Path
import sys,json,datetime
R=Path(__file__).resolve().parents[1]
def archive_for_one_retry(model,stages=None,shards=1,shard_index=0,members=None):
 jobs=[j for j in json.loads((R/'00_protocol/PREDICTION_JOBS_FROZEN.json').read_text()) if j['model']==model and (not stages or j['stage'] in stages)]
 order={'GLOBAL_CORE':0,'REFERENCE':1,'GNRS_GENE_CENTER':2,'PLACEMENT':3,'CONTEXT':4,'ENSEMBLE':5};jobs.sort(key=lambda j:(str(j['member']),order[j['stage']],j['audit_job_ID']));jobs=[j for i,j in enumerate(jobs) if i%shards==shard_index];folders={'GLOBAL_CORE':'02_track_harmonization','ENSEMBLE':'03_within_model_ensembles','PLACEMENT':'04_context_robustness','CONTEXT':'04_context_robustness','GNRS_GENE_CENTER':'08_truth_coverage','REFERENCE':'07_reference_benchmark'};records=[]
 for j in jobs:
  if members is not None and j['member'] not in members:continue
  p=R/folders[j['stage']]/model/(j['audit_job_ID']+'.json')
  if not p.exists():continue
  z=json.loads(p.read_text());error=z.get('error','');kind=z.get('error_type','')
  transient=(kind=='UnboundLocalError' and "local variable 'net'" in error and model=='Borzoi' and j['stage']=='ENSEMBLE' and (R/'01_model_metadata'/f"borzoi_official_rep{j['member']}_model0_best.h5").exists()) or (kind=='AssertionError' and error=='OFFICIAL_WEIGHT_NOT_READY' and model=='Borzoi' and j['stage']=='ENSEMBLE' and (R/'01_model_metadata'/f"borzoi_official_rep{j['member']}_model0_best.h5").exists()) or kind in ['_MultiThreadedRendezvous','_InactiveRpcError','TimeoutError','ConnectionError'] or (kind=='FileNotFoundError' and 'borzoi_official_rep' in error)
  if z.get('status')!='UNAVAILABLE' or not transient:continue
  # Missing official weights are not retried until the exact corresponding artifact is present.
  if kind=='FileNotFoundError' and not (R/'01_model_metadata'/f"borzoi_official_rep{j['member']}_model0_best.h5").exists():continue
  hist=p.parent/'bounded_retry_history';hist.mkdir(exist_ok=True);target=hist/(p.stem+'_attempt1.json')
  if target.exists():continue
  p.rename(target);records.append({'ID':j['audit_job_ID'],'model':model,'stage':j['stage'],'initial_error_type':kind,'action':'one additional controller resumption;original result retained'})
 if records:
  out=R/'00_protocol'/f'RESUMPTION_AUDIT_{model}_{shard_index}.json';out.write_text(json.dumps({'UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'records':records},indent=2))
 return records
