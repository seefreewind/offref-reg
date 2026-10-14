from pathlib import Path
import json,subprocess,time,os,datetime,hashlib,sys,signal
R=Path(__file__).resolve().parents[1];AG='/Users/zy/.codex/tools/alphagenome-api-venv/bin/python';B='/Users/zy/.codex/tools/offref-reg-venv/bin/python'
def state(status,**more):(R/'14_reports/EXECUTION_STATUS.json').write_text(json.dumps({'status':status,'UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'Phase2':'NOT_ENTERED',**more},indent=2))
def alive(pid):
 try:os.kill(pid,0);z=subprocess.run(['ps','-p',str(pid),'-o','state='],capture_output=True,text=True);return z.returncode==0 and 'Z' not in z.stdout
 except ProcessLookupError:return False
def launch(label,args,py=AG):
 f=(R/'01_model_metadata'/f'{label}.log').open('a');p=subprocess.Popen([py,str(R/'00_protocol/run_predictions.py'),*args],stdout=f,stderr=subprocess.STDOUT,start_new_session=True);return p.pid

from retry_transient_failures import archive_for_one_retry
def retry(w):
 args=w['args'];model=args[args.index('--model')+1];stages=None
 if '--stages' in args:
  index=args.index('--stages')+1;stages=[]
  while index<len(args) and not args[index].startswith('--'):stages.append(args[index]);index+=1
 shards=int(args[args.index('--shards')+1]) if '--shards' in args else 1;shard_index=int(args[args.index('--shard-index')+1]) if '--shard-index' in args else 0
 archive_for_one_retry(model,['ENSEMBLE'] if model=='Borzoi' else stages,shards,shard_index,members=['2','3'] if model=='Borzoi' else None)
def main(existing_b):
 workers={'Borzoi':{'pid':existing_b,'args':['--model','Borzoi'],'python':B,'resumptions':0},'AG_LCL_0':{'args':['--model','AlphaGenome','--stages','REFERENCE','GNRS_GENE_CENTER','PLACEMENT','CONTEXT','ENSEMBLE','--shards','2','--shard-index','0'],'python':AG,'resumptions':0},'AG_LCL_1':{'args':['--model','AlphaGenome','--stages','REFERENCE','GNRS_GENE_CENTER','PLACEMENT','CONTEXT','ENSEMBLE','--shards','2','--shard-index','1'],'python':AG,'resumptions':0},'AG_GLOBAL':{'args':['--model','AlphaGenome','--stages','GLOBAL_CORE'],'python':AG,'resumptions':0}}
 for label in ['AG_LCL_0','AG_LCL_1']:
  workers[label]['pid']=int(sys.argv[2 if label=='AG_LCL_0' else 3]) if len(sys.argv)>3 else launch(label,workers[label]['args'],workers[label]['python'])
 # Await successful official global batching smoke, so a technical limitation cannot generate400 repeated failures.
 state('FORMAL_BORZOI_RUNNING_GLOBAL_AG_SMOKE_PENDING',workers=workers)
 while True:
  valid=[]
  for p in (R/'02_track_harmonization/AlphaGenome').glob('*.json'):
   if p.name.startswith('._'):continue
   z=json.loads(p.read_text())
   if z['status']=='PASS':valid.append(z)
  if valid:
   pilot=R/'02_track_harmonization/GLOBAL_CHECKPOINT_PILOT_PROCESS.json'
   if pilot.exists() and alive(json.loads(pilot.read_text())['PID']):os.kill(json.loads(pilot.read_text())['PID'],signal.SIGTERM)
   break
  if (R/'02_track_harmonization/GLOBAL_TECHNICAL_LIMITATION_FINAL.json').exists():break
  time.sleep(30)
 for label,w in workers.items():
  if label=='AG_GLOBAL':w['pid']=int(sys.argv[4]) if len(sys.argv)>4 else launch(label,w['args'],w['python'])
 state('FORMAL_PREDICTIONS_RUNNING',workers=workers)
 # Reinvoke resumable script once after any exiting worker; unchanged per-job JSON prevents repeated completed inference.
 while True:
  statuses=[]
  for label,w in workers.items():
   if alive(w['pid']):statuses.append(True)
   elif w['resumptions']<1:
    if label=='Borzoi' and not all((R/'01_model_metadata'/f'borzoi_official_rep{k}_model0_best.h5').exists() for k in [1,2,3]):
     download=json.loads((R/'01_model_metadata/DOWNLOAD_JOB.json').read_text())
     if alive(download['PID']):statuses.append(True);continue
     if download.get('controller_resume_N',0)<1:
      f=(R/'01_model_metadata/weight_download_final_resume.log').open('a');dp=subprocess.Popen([AG,str(R/'01_model_metadata/download_official_weights.py')],stdout=f,stderr=subprocess.STDOUT,start_new_session=True);download.update(PID=dp.pid,controller_resume_N=1,status='RUNNING_RESUMED');(R/'01_model_metadata/DOWNLOAD_JOB.json').write_text(json.dumps(download,indent=2));statuses.append(True);continue
    retry(w);w['resumptions']+=1;w['pid']=launch(label+'_resume',w['args'],w['python']);statuses.append(True)
   else:statuses.append(False)
  if not any(statuses):break
  time.sleep(45)
 state('PREDICTIONS_COMPLETE_ANALYSIS_PENDING',workers=workers)
 for file in ['11_statistics/analyze.py','09_disagreement_models/fit_models.py','11_statistics/finalize_endpoints.py']:
  if not (R/file).exists():state('PREDICTIONS_COMPLETE_ANALYSIS_CODE_REQUIRED',required=file);return
  with (R/'11_statistics/analysis.log').open('a') as f:subprocess.run([AG,str(R/file)],stdout=f,stderr=subprocess.STDOUT,check=True)
 state('STATISTICS_COMPLETE_FIGURES_PENDING')
 for file in ['12_figures/render_figures.py','14_reports/build_report.py']:
  if not (R/file).exists():state('STATISTICS_COMPLETE_DELIVERABLE_CODE_REQUIRED',required=file);return
  with (R/'14_reports/completion.log').open('a') as f:subprocess.run([AG,str(R/file)],stdout=f,stderr=subprocess.STDOUT,check=True)
 state('FIGURES_AND_DRAFT_REPORT_COMPLETE_MANUAL_QA_REQUIRED')
if __name__=='__main__':
 try:main(int(sys.argv[1]))
 except Exception as e:state('CONTROLLER_FAILED',error_type=type(e).__name__,error=str(e)[:500]);raise
