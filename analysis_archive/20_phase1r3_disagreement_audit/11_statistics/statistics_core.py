import numpy as np
from scipy.stats import spearmanr,pearsonr,rankdata,binomtest
SEED=20261011

def interval(v):
 v=np.asarray(v);v=v[np.isfinite(v)];return np.quantile(v,[.025,.975]).tolist() if len(v) else [None,None]
def rowcorr(a,b):
 a=a-a.mean(axis=1,keepdims=True);b=b-b.mean(axis=1,keepdims=True);den=np.sqrt(np.sum(a*a,axis=1)*np.sum(b*b,axis=1));return np.divide(np.sum(a*b,axis=1),den,out=np.full(len(a),np.nan),where=den>0)
def boot(a,b,B=10000,cluster=None,seed=SEED):
 rng=np.random.default_rng(seed);a=np.array(a,float);b=np.array(b,float);n=len(a);r=[];p=[];s=[]
 if cluster is None:
  for _ in range(0,B,500):
   ix=rng.integers(0,n,(min(500,B-len(r)),n));x=a[ix];y=b[ix];r.extend(rowcorr(rankdata(x,axis=1),rankdata(y,axis=1)));p.extend(rowcorr(x,y));sg=(x!=0)&(y!=0);s.extend(np.divide(np.sum((np.sign(x)==np.sign(y))&sg,axis=1),sg.sum(axis=1),out=np.full(len(ix),np.nan),where=sg.sum(axis=1)>0))
 else:
  labels=np.array(cluster);unique=np.unique(labels);groups=[np.where(labels==z)[0] for z in unique]
  for _ in range(B):
   ix=np.concatenate([groups[i] for i in rng.integers(0,len(groups),len(groups))]);x=a[ix];y=b[ix];r.append(spearmanr(x,y).statistic if np.ptp(x) and np.ptp(y) else np.nan);p.append(pearsonr(x,y).statistic if np.ptp(x) and np.ptp(y) else np.nan);sg=(x!=0)&(y!=0);s.append(np.mean(np.sign(x[sg])==np.sign(y[sg])) if sg.any() else np.nan)
 return {'rho':r,'pearson':p,'sign':s}
def metric(a,b,frozen_N=None,cluster=None):
 a=np.array(a,float);b=np.array(b,float);good=np.isfinite(a)&np.isfinite(b);a=a[good];b=b[good];cluster=np.array(cluster)[good] if cluster is not None else None;n=len(a);signed=(a!=0)&(b!=0);N=int(signed.sum());K=int(np.sum(np.sign(a[signed])==np.sign(b[signed])));bt=binomtest(K,N,.5) if N else None;ci=bt.proportion_ci() if bt else None
 out={'frozen_N':frozen_N,'finite_N':n,'signed_N':N,'concordant_N':K,'proportion':K/N if N else None,'exact_CI95':[ci.low,ci.high] if ci else [None,None],'binomial_P':bt.pvalue if bt else None,'rho':float(spearmanr(a,b).statistic) if n>=3 and np.ptp(a) and np.ptp(b) else None,'pearson':float(pearsonr(a,b).statistic) if n>=3 and np.ptp(a) and np.ptp(b) else None,'median_abs_AG':float(np.median(abs(a))) if n else None,'median_abs_Borzoi':float(np.median(abs(b))) if n else None}
 if n>=3:
  v=boot(a,b,cluster=cluster);out.update(rho_CI95=interval(v['rho']),pearson_CI95=interval(v['pearson']),sign_bootstrap_CI95=interval(v['sign']),bootstrap='10000 gene/locus-cluster' if cluster is not None else '10000 paired locus')
 else:out.update(rho_CI95=[None,None],pearson_CI95=[None,None],sign_bootstrap_CI95=[None,None])
 return out
