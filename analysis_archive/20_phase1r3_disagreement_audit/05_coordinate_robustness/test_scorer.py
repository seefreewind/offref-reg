from pathlib import Path
import sys,json,numpy as np
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'00_protocol'));from audit_common import *
checks=[]
for L in [0,13,32,96,3701]:
 anchor=1000000;W=524288;start=anchor-W//2;d=domain(anchor,L,start,W,163840,32);j={'anchor0':anchor,'insertion':'A'*L,'W':W,'start0':start,'exons':[[anchor-1000,anchor-200],[anchor+4000,anchor+6000]],'strict':True};ints=regions(j,'Borzoi',d);a=np.ones((6144,2));z=a*3
 ra=summarize(a,ints,start,163840,32,'REF');rz=summarize(z,ints,start,163840,32,'ALT');assert np.allclose(np.array(rz['track_means'])-ra['track_means'],2)
 assert sum(b-a for a,b,_ in ints)==2800;checks.append({'insertion_length':L,'expected_delta':2,'status':'PASS'})
assert merge([[1,5],[3,8],[8,9]])==[[1,9]]
# Independently rescore cached arrays for one locus per model; compare original native common score.
import pandas as pd
jobs=[j for j in readjson(OLD/'00_protocol/JOBS_FROZEN.json') if j['role']=='SENTINEL' and j['panel']=='UNBIASED_CLASSB_CORE'];j=next(j for j in jobs if readjson(OLD/'10_statistics/AlphaGenome'/f"{j['job_ID']}.json")['common_RNA']['status']=='PASS');md=pd.read_csv(R/'02_track_harmonization/TRACK_METADATA_FULL.tsv',sep='\t').fillna('')
for model in ['AlphaGenome','Borzoi']:
 W=1048576 if model=='AlphaGenome' else 524288;start=j['anchor0']-W//2;score=readjson(OLD/'10_statistics'/model/f"{j['job_ID']}.json");domainB=domain(j['anchor0'],len(j['insertion']),j['anchor0']-262144,524288,163840,32);local=dict(j,W=W,start0=start);ints=regions(local,model,domainB);means=[]
 for allele in ['REF','ALT']:
  c=score['cache_'+allele];meta=readjson(OLD/'09_predictions'/model/(c+'.json'));data=np.load(OLD/'09_predictions'/model/(c+'.npz'));v=data['rna'];res=int(data['resolution']);data.close()
  if model=='Borzoi':
   x=v.astype(float)
   for k,t in enumerate(meta['metadata']):
    q=x[:,k]/float(t['scale']);clip=float(t['clip_soft']);q=np.where(q>clip,(q-clip)**2+clip,q);x[:,k]=q**(1/.75)/32
   v=x
  out=summarize(v,ints,start,163840 if model=='Borzoi' else 0,res,allele);idx=[k for k,t in enumerate(meta['metadata']) if t['strand']==j['strand']];means.append(np.mean(np.array(out['track_means'])[idx]))
 delta=means[1]-means[0];expected=score['common_RNA']['delta'];assert np.isclose(delta,expected,rtol=1e-7,atol=1e-10),(model,delta,expected);checks.append({'model':model,'cached_locus':j['job_ID'],'independent_delta':delta,'frozen_delta':expected,'status':'PASS'})
(R/'05_coordinate_robustness/SCORER_QA.json').write_text(json.dumps({'status':'PASS','analytic_and_independent_cache_checks':checks,'overlapping_exon_union':'PASS','fractional_overlap_and_downstream_nonmultiple32':'PASS'},indent=2));print('COORDINATE_AND_ORIGINAL_SCORE_QA_PASS')
