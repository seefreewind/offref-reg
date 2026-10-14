from pathlib import Path
import os,csv,json,hashlib
R=Path(__file__).resolve().parents[1];prior=json.loads((R/'14_freeze/PRIOR_FREEZE_INTEGRITY.json').read_text());checks=[]
for key,n in [('historical',590),('Phase1R',221)]:checks.append([key,'PASS' if prior[key]['pass_N']==prior[key]['N']==n and prior[key]['status']=='PASS' else 'FAIL',f'Complete prior byte rehash {n}/{n}'])
for name in ['AMENDMENT_SHA256.tsv','IMPLEMENTATION_SHA256.tsv']:
 r=next(csv.DictReader((R/'14_freeze'/name).open(),delimiter='\t'));sha=hashlib.sha256((R/r['path']).read_bytes()).hexdigest();checks.append([r['path'],'PASS' if sha==r['sha256'] else 'FAIL','Unchanged before output freeze'])
for key,p in [('historical',R.parent/'15_freeze/SHA256_MANIFEST.tsv'),('Phase1R',R.parent/'16_phase1r_rescue/13_freeze/SHA256_MANIFEST.tsv')]:checks.append([key+'_manifest','PASS' if hashlib.sha256(p.read_bytes()).hexdigest()==prior[key]['manifest_sha256'] else 'FAIL','Prior manifest unchanged'])
s=json.loads((R/'07_sentinel/EXECUTION_STATUS.json').read_text());checks.append(['No_formal_sentinel_or_Phase2','PASS' if s['sentinel_frozen_N']==0 and s['formal_prediction_N']==0 and s['Phase2']=='NOT_ENTERED' and not (R/'07_sentinel/SENTINEL500_FROZEN.tsv').exists() else 'FAIL','S1 HOLD enforced'])
assert all(r[1]=='PASS' for r in checks),checks
with (R/'14_freeze/FREEZE_VALIDATION.tsv').open('w') as f:w=csv.writer(f,delimiter='\t');w.writerow(['check','status','evidence']);w.writerows(checks)
rows=[]
for root,dirs,files in os.walk(R):
 dirs[:]=[d for d in dirs if not d.startswith('._') and d!='__pycache__']
 for name in sorted(files):
  p=Path(root)/name
  if name.startswith('._') or (p.parent==R/'14_freeze' and name in ['SHA256_MANIFEST.tsv','FILE_INVENTORY.tsv']):continue
  with p.open('rb') as f:sha=hashlib.file_digest(f,'sha256').hexdigest()
  rows.append([str(p.relative_to(R)),p.stat().st_size,sha])
rows.sort()
for name in ['SHA256_MANIFEST.tsv','FILE_INVENTORY.tsv']:
 with (R/'14_freeze'/name).open('w') as f:w=csv.writer(f,delimiter='\t');w.writerow(['path','bytes','sha256']);w.writerows(rows)
print('Phase1R2 output freeze PASS:',len(rows),'files;all checks PASS;no Phase2',flush=True)
