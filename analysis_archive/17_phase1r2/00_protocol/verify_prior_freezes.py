from pathlib import Path
import csv,hashlib,json,time
R=Path(__file__).resolve().parents[1];P=R.parent;results={}
for label,base,manifest in [('historical',P,P/'15_freeze/SHA256_MANIFEST.tsv'),('Phase1R',P/'16_phase1r_rescue',P/'16_phase1r_rescue/13_freeze/SHA256_MANIFEST.tsv')]:
 rows=[]
 for r in csv.DictReader(manifest.open(),delimiter='\t'):
  p=base/r['path']
  try:
   with p.open('rb') as f:sha=hashlib.file_digest(f,'sha256').hexdigest()
   status='PASS' if sha==r['sha256'] and p.stat().st_size==int(r['bytes']) else 'FAIL'
  except Exception as e:sha=type(e).__name__;status='FAIL'
  rows.append([r['path'],r['sha256'],sha,status])
 with (R/'14_freeze'/f'{label}_RECHECK.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['path','expected_sha256','actual_sha256','status']);w.writerows(rows)
 results[label]={'N':len(rows),'pass_N':sum(r[3]=='PASS' for r in rows),'status':'PASS' if all(r[3]=='PASS' for r in rows) else 'FAIL','manifest_sha256':hashlib.sha256(manifest.read_bytes()).hexdigest()}
 print(label,results[label],flush=True)
(R/'14_freeze/PRIOR_FREEZE_INTEGRITY.json').write_text(json.dumps(results,indent=2))
assert results['historical']['N']==590 and all(v['status']=='PASS' for v in results.values()),results
