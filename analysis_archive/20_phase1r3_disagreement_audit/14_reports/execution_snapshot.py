from pathlib import Path
import json,collections,subprocess,datetime
R=Path(__file__).resolve().parents[1];jobs=json.loads((R/'00_protocol/PREDICTION_JOBS_FROZEN.json').read_text());folders={'GLOBAL_CORE':'02_track_harmonization','ENSEMBLE':'03_within_model_ensembles','PLACEMENT':'04_context_robustness','CONTEXT':'04_context_robustness','GNRS_GENE_CENTER':'08_truth_coverage','REFERENCE':'07_reference_benchmark'};counts=collections.defaultdict(collections.Counter)
for j in jobs:
 p=R/folders[j['stage']]/j['model']/(j['audit_job_ID']+'.json');counts[j['model']+' '+j['stage']][json.loads(p.read_text())['status'] if p.exists() else 'PENDING']+=1
records=[]
for p in [R/'00_protocol/CONTROLLER_PROCESS.json',R/'00_protocol/BORZOI_FORMAL_PROCESS.json',R/'02_track_harmonization/GLOBAL_CHECKPOINT_PILOT_PROCESS.json']:
 if p.exists():
  z=json.loads(p.read_text());records.append(z)
checkpoints=[]
for p in (R/'02_track_harmonization/AlphaGenome/ontology_checkpoints').rglob('[!._]*.json'):
 if p.name.endswith('.transport_failure.json'):continue
 z=json.loads(p.read_text());checkpoints.append({'RNA_tracks_N':len(z['metadata']),'seconds':z['elapsed_seconds']})
output={'status':'RUNNING_NOT_FINAL','UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'Phase2':'NOT_ENTERED','frozen_jobs_N':len(jobs),'job_status_counts':{k:dict(v) for k,v in counts.items()},'process_records':records,'global_pilot_completed_ontology_chunks_N':len(checkpoints),'global_pilot_chunks':checkpoints,'runtime_estimates':'Borzoi approximately5–8hours total;AlphaGenome matched-tissue approximately1–2days;global all667RNAtrack sensitivity may take weeks at observed transport rate','pending_preference':'Full400all-track sensitivity retained in plan until user explicitly defers it;main andpilot jobs continue','completion_guards':'No final report/scientific verdict while predictions missing;tenfigurefinalsizeQA andfinalhash required'};(R/'14_reports/RUN_SNAPSHOT.json').write_text(json.dumps(output,indent=2));print(json.dumps(output['job_status_counts']))
