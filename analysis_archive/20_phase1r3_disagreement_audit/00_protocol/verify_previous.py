from pathlib import Path
import csv,json,hashlib,datetime
R=Path(__file__).resolve().parents[1];P=R.parent
spec=[('Phase0',P,P/'15_freeze/SHA256_MANIFEST.tsv'),('Phase1R',P/'16_phase1r_rescue',P/'16_phase1r_rescue/13_freeze/SHA256_MANIFEST.tsv'),('Phase1R2_preparation',P/'17_phase1r2',P/'17_phase1r2/14_freeze/SHA256_MANIFEST.tsv'),('Phase1R2_execution',P/'19_phase1r2_execution',P/'19_phase1r2_execution/14_freeze/SHA256_MANIFEST.tsv')]
out=[];cache={}
for name,base,manifest in spec:
 rows=list(csv.DictReader(manifest.open(),delimiter='\t'));results=[]
 for row in rows:
  p=base/row['path'];key=str(p.resolve());h=hashlib.sha256()
  if key not in cache:
   with p.open('rb') as f:
    for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
   cache[key]=h.hexdigest()
  results.append({'path':row['path'],'expected':row['sha256'],'observed':cache[key],'status':'PASS' if cache[key]==row['sha256'] else 'FAIL'})
 with (R/'15_freeze'/f'{name}_RECHECK.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=results[0].keys(),delimiter='\t');w.writeheader();w.writerows(results)
 out.append({'phase':name,'N':len(rows),'PASS':sum(x['status']=='PASS' for x in results)});print(name,out[-1],flush=True)
status={'overall':'PASS' if all(x['N']==x['PASS'] for x in out) else 'FAIL','checks':out,'time_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()};(R/'15_freeze/PREVIOUS_FREEZES_RECHECK.json').write_text(json.dumps(status,indent=2));assert status['overall']=='PASS';print('ALL_PREVIOUS_FREEZES_PASS',flush=True)
