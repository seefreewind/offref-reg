from pathlib import Path
import subprocess,hashlib,json
R=Path(__file__).resolve().parents[1];rows=[]
for i in [1,2,3]:
 url=f'https://storage.googleapis.com/seqnn-share/borzoi/f{i}/model0_best.h5';p=R/'01_model_metadata'/f'borzoi_official_rep{i}_model0_best.h5';temp=p.with_suffix('.partial')
 if not p.exists():
  subprocess.run(['curl','-L','--fail','--retry','4','--retry-all-errors','--connect-timeout','30','--max-time','1800','-C','-','-o',str(temp),url],check=True);temp.rename(p)
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 rows.append({'replicate':i,'URL':url,'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':h.hexdigest(),'source':'official Calico download_models.sh;same data splits, not independent folds'});print('OFFICIAL_REPLICATE_DOWNLOADED',i,flush=True)
(R/'01_model_metadata/BORZOI_OFFICIAL_ENSEMBLE_WEIGHTS.json').write_text(json.dumps(rows,indent=2))
