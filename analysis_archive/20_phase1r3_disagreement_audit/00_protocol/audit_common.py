from pathlib import Path
import json,hashlib,math,csv
import numpy as np
R=Path(__file__).resolve().parents[1];P=R.parent;OLD=P/'19_phase1r2_execution'
def readjson(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def assert_freezes():
 assert readjson(R/'15_freeze/PREVIOUS_FREEZES_RECHECK.json')['overall']=='PASS'
 z=readjson(R/'15_freeze/PROTOCOL_SHA256.json');assert sha(R/z['path'])==z['sha256']
 for x in readjson(R/'15_freeze/PREDICTION_INPUT_SHA256.json'):assert sha(R/x['path'])==x['sha256']
def merge(xs):
 z=[]
 for a,b in sorted(xs):
  if b<=a:continue
  if z and a<=z[-1][1]:z[-1][1]=max(b,z[-1][1])
  else:z.append([int(a),int(b)])
 return z

def domain(anchor,L,start,W,crop,res,strict=True):
 idx=anchor-start;out=[]
 for lo,hi,off in [(crop,min(idx,W-crop),0),(max(crop,idx),W-crop-L,L)]:
  if strict:
   lo=max(lo,math.ceil(lo/res)*res,math.ceil((lo+off)/res)*res-off);hi=min(hi,math.floor(hi/res)*res,math.floor((hi+off)/res)*res-off)
  if hi>lo:out.append((start+lo,start+hi,off))
 return out

def intersect_domains(ds,others):
 for od in others:
  ds=[(max(a,c),min(b,d),o) for a,b,o in ds for c,d,_ in od if min(b,d)>max(a,c)]
 return ds

def regions(job,model='Borzoi',fixed=None):
 anchor=job['anchor0'];L=len(job.get('insertion',''));W=job['W'];start=job['start0'];crop=163840 if model=='Borzoi' else 0;res=32 if model=='Borzoi' else 1
 ds=domain(anchor,L,start,W,crop,res,job.get('strict',True))
 if fixed is not None:ds=intersect_domains(ds,[fixed])
 return [(max(a,x),min(b,y),o) for a,b,o in ds for x,y in merge(job['exons']) if min(b,y)>max(a,x)]

def integrate(values,lo,hi,res):
 # Piecewise-constant bin integration; only target slices, no full-array cumsum/copy.
 n=len(values);l=lo/res;h=hi/res;assert l>=-1e-8 and h<=n+1e-8;(l,h)=(max(0,l),min(n,h));i=int(math.floor(l));j=int(math.floor(h))
 if i==j:return values[i].astype(np.float64)*(h-l)*res if i<n else np.zeros(values.shape[1])
 v=values[i].astype(np.float64)*(1-(l-i))
 if j>i+1:v+=np.sum(values[i+1:j],axis=0,dtype=np.float64)
 if j<n and h>j:v+=values[j].astype(np.float64)*(h-j)
 return v*res

def summarize(values,intervals,start,crop,res,allele):
 bp=sum(b-a for a,b,_ in intervals)
 if not bp:return {'status':'NO_COMMON_EXON','mapped_bp':0,'track_means':None}
 mass=np.zeros(values.shape[1])
 for a,b,off in intervals:
  shift=off if allele=='ALT' else 0;mass+=integrate(values,a-start-crop+shift,b-start-crop+shift,res)
 return {'status':'PASS','mapped_bp':bp,'track_means':(mass/bp).tolist()}

def sequence(genome,job,model,allele):
 start=job['start0'];W=job['W'];chrom=job['chrom'];anchor=job['anchor0'];L=len(job.get('insertion',''));i=anchor-start
 assert 0<=start and start+W<=genome.chroms(chrom),'REFERENCE_BOUNDARY'
 assert 0<=i and i+L<W,'COMPLETE_INSERTION_NOT_CONTAINED'
 s=genome.sequence(chrom,start,start+W).upper();assert len(s)==W
 if allele=='ALT':
  if job.get('allele_type')=='SNV':
   q=job['pos1']-1-start;assert s[q]==job['REF'],'REFERENCE_ALLELE_MISMATCH';s=s[:q]+job['ALT']+s[q+1:]
  else:s=s[:i]+job['insertion']+s[i:W-L];assert s[i:i+L]==job['insertion']
 if job.get('mask_length'):
  radius=job['mask_length']//2;lo=max(0,i-radius);hi=min(W,i+radius);assert i+L<=hi;s='N'*lo+s[lo:hi]+'N'*(W-hi)
 assert len(s)==W;return s

def delta_score(result,families=None):
 if result.get('status')!='PASS':return None
 a=np.array(result['REF']['track_means']);b=np.array(result['ALT']['track_means']);meta=result['metadata'];strand=result['strand'];ix=[i for i,t in enumerate(meta) if t['strand']==strand]
 if families is not None:
  means=[]
  for fam,ids in families.items():
   select=[i for i in ix if int(meta[i]['track_ID']) in ids]
   if not select:return None
   means.append(float((b-a)[select].mean()))
  return float(np.mean(means)) if means else None
 return float((b-a)[ix].mean()) if ix else None
