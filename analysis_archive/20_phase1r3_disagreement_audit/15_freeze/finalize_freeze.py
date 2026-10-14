from pathlib import Path
import json,hashlib,sys,subprocess,datetime
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'00_protocol'));from audit_common import assert_freezes
assert_freezes();status=json.loads((R/'11_statistics/ANALYSIS_STATUS.json').read_text());assert status['counts'].get('MISSING',0)==0
qa=json.loads((R/'12_figures/FINAL_SIZE_VISUAL_QA.json').read_text());assert qa['status']=='PASS' and sorted(qa['inspected_figures'])==list(range(1,11)),'All ten figures require final-size visual QA before freeze'
for figure in qa['figures']:
 for filename,expected in figure['files_sha256'].items():
  assert hashlib.sha256((R/'12_figures'/filename).read_bytes()).hexdigest()==expected,'Figure changed after visual QA: '+filename
assert (R/'14_reports/PHASE1R3_FINAL_REPORT.md').exists();route=json.loads((R/'13_manuscript_gate/FINAL_SCIENTIFIC_ROUTE.json').read_text());assert route['Phase2']=='NOT_ENTERED'
for i in range(1,11):
 for suffix in ['.pdf','.svg','.png','_source.tsv']:assert (R/'12_figures'/f'Figure_{i:02d}{suffix}').is_file()
# Hash prior freezes again at close; no old artifact is edited.
with (R/'15_freeze/final_prior_recheck.log').open('w') as f:subprocess.run([sys.executable,str(R/'00_protocol/verify_previous.py')],stdout=f,stderr=subprocess.STDOUT,check=True)
assert json.loads((R/'15_freeze/PREVIOUS_FREEZES_RECHECK.json').read_text())['overall']=='PASS'
# Rebuild the report with the current close-time integrity evidence.
subprocess.run([sys.executable,str(R/'14_reports/build_report.py')],check=True)
report=R/'14_reports/PHASE1R3_FINAL_REPORT.md';text=report.read_text().replace('Figures are rendered separately; final visual review and output hashing remain required until their audit records confirm completion.','All ten figures passed final-size visual review and automated font/source audits. Output hashes are recorded in the final manifest.');report.write_text(text)
records=[];excluded=[]
for p in sorted(R.rglob('*')):
 if not p.is_file():continue
 relative=str(p.relative_to(R));skip=(p.name.startswith('._') or '__pycache__' in p.parts or p.suffix=='.log' or p.name.endswith('.partial') or p.name.endswith('_PROCESS.json') or p.name in ['SHA256_MANIFEST.json','SHA256_MANIFEST.tsv','FINAL_FREEZE_STATUS.json','EXECUTION_STATUS.json'])
 if skip:excluded.append(relative);continue
 h=hashlib.sha256()
 with p.open('rb') as f:
  while chunk:=f.read(8*1024*1024):h.update(chunk)
 records.append({'path':relative,'bytes':p.stat().st_size,'sha256':h.hexdigest()})
(R/'15_freeze/SHA256_MANIFEST.json').write_text(json.dumps(records,indent=2));(R/'15_freeze/SHA256_MANIFEST.tsv').write_text('path\tbytes\tsha256\n'+''.join(f"{z['path']}\t{z['bytes']}\t{z['sha256']}\n" for z in records));(R/'15_freeze/FINAL_FREEZE_STATUS.json').write_text(json.dumps({'status':'PASS','UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scientific_route':route['FINAL_VERDICT'],'files_N':len(records),'excluded_mutable_operational_files':excluded,'Phase2':'NOT_ENTERED','stop':True},indent=2));print('FINAL_FREEZE_PASS_STOP_BEFORE_PHASE2',flush=True)
