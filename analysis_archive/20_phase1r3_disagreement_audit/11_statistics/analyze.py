from pathlib import Path
import json,sys,collections
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];sys.path[:0]=[str(R/'00_protocol'),str(R/'11_statistics')];from audit_common import *
from statistics_core import metric,boot,interval

def clean(x):
 if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
 if isinstance(x,list):return [clean(v) for v in x]
 if isinstance(x,(np.integer,)):return int(x)
 if isinstance(x,(float,np.floating)):return float(x) if np.isfinite(x) else None
 return x

def save(name,x):(R/'11_statistics'/name).write_text(json.dumps(clean(x),indent=2))
def values(z,tier='TIER2_NEAR'):
 if z.get('status')!='PASS':return None
 if tier=='TIER1_EXACT':return z.get('tier1_delta')
 if tier=='TIER3_GLOBAL':return z.get('delta')
 mask=[i for i,t in enumerate(z['metadata']) if t['strand']==z['strand'] and t['cell_tissue']=='GM12878'];d=np.array(z['ALT']['track_means'])-z['REF']['track_means'];return float(d[mask].mean()) if mask else None

def main():
 assert_freezes();labels=readjson(R/'15_freeze/CORRECTED_LABEL_INPUT_SHA256.json');assert sha(R/labels['path'])==labels['sha256'];s=pd.read_csv(R/'04_context_robustness/SUBSET100_FROZEN.tsv',sep='\t');core=pd.read_csv(R/'06_offref_core/CORE400_DERIVED_STRATA_CORRECTED.tsv',sep='\t');cached=pd.read_csv(R/'11_statistics/CACHED_ORIGINAL_AND_TIER1_SCORES.tsv',sep='\t');core=core.merge(cached.drop(columns=['panel','gene','strand','anchor0','length']),on='OFFREF_ID',validate='one_to_one');joblist=readjson(R/'00_protocol/PREDICTION_JOBS_FROZEN.json');results={};statuses=[]
 for j in joblist:
  folder={'GLOBAL_CORE':'02_track_harmonization','ENSEMBLE':'03_within_model_ensembles','PLACEMENT':'04_context_robustness','CONTEXT':'04_context_robustness','GNRS_GENE_CENTER':'08_truth_coverage','REFERENCE':'07_reference_benchmark'}[j['stage']];p=R/folder/j['model']/(j['audit_job_ID']+'.json');z=readjson(p) if p.exists() else {'status':'MISSING'};results[j['audit_job_ID']]=z;statuses.append({'ID':j['audit_job_ID'],'model':j['model'],'stage':j['stage'],'status':z['status'],'error':z.get('error')})
 assert not any(x['status']=='MISSING' for x in statuses),'Unfinished jobs: statistical finalization prohibited'
 pd.DataFrame(statuses).to_csv(R/'11_statistics/JOB_OUTCOME_AUDIT.tsv',sep='\t',index=False)
 for j in joblist:
  if j['stage']=='GLOBAL_CORE':
   z=results[j['audit_job_ID']];idx=core.index[core.OFFREF_ID==j['job_ID']]
   for tier in ['TIER1_EXACT','TIER2_NEAR','TIER3_GLOBAL']:core.loc[idx,j['model']+'_'+tier]=values(z,tier)
 # Cached exact-tier is valid for old strand-corrected16 Borzoi/4 AG tracks; all other GM tracks are unstranded. Report separately, never replace failed new requests.
 summary={}
 for tier in ['original','tier1_cached','log1p','TIER1_EXACT','TIER2_NEAR','TIER3_GLOBAL']:
  for m in ['AlphaGenome','Borzoi']:
   if m+'_'+tier not in core:core[m+'_'+tier]=np.nan
  summary[tier]=metric(core['AlphaGenome_'+tier],core['Borzoi_'+tier],400);print('AGREEMENT',tier,summary[tier]['finite_N'],summary[tier]['rho'],flush=True)
 core.to_csv(R/'11_statistics/CORE400_ALL_SCORING_DEFINITIONS.tsv',sep='\t',index=False);save('HARMONIZATION_AGREEMENT.json',summary)
 strat=[]
 for variable in ['repeat_stratum','length_bin']:
  cats=['NON_REPEAT','SINE/Alu','LINE','LTR','VNTR/STR','OTHER_REPEAT'] if variable=='repeat_stratum' else ['50–199','200–499','500–999','1–2kb','2–5kb','>5kb']
  for cat in cats:
   z=core[core[variable]==cat];a=z.AlphaGenome_original;b=z.Borzoi_original;r=metric(a,b,len(z));r.update(variable=variable,stratum=cat,coverage_failure_fraction=1-r['finite_N']/len(z) if len(z) else None);strat.append(r)
 for rich in [True,False]:
  z=core[(core.repeat_fraction>=.8)==rich];r=metric(z.AlphaGenome_original,z.Borzoi_original,len(z));r.update(variable='repeat_fraction',stratum='repeat-rich' if rich else 'low-repeat');strat.append(r)
 save('OOD_STRATA.json',strat)
 robustness=[];ensemble=[]
 for model in ['AlphaGenome','Borzoi']:
  for stage,dimension,baseline,alternatives in [('PLACEMENT','placement_shift',0,[-32768,32768]),('CONTEXT','context_length',1048576 if model=='AlphaGenome' else 524288,[16384,131072])]:
   table=[]
   for j in joblist:
    if j['stage']==stage and j['model']==model:
     z=results[j['audit_job_ID']];table.append({'OFFREF_ID':j['job_ID'],'setting':j[dimension],'delta':values(z),'status':z['status']})
   table=pd.DataFrame(table);table.to_csv(R/'11_statistics'/f'{model}_{stage}_ALL_IDENTITIES.tsv',sep='\t',index=False);wide=table.pivot(index='OFFREF_ID',columns='setting',values='delta')
   for alt in alternatives:
    a=wide[baseline] if baseline in wide else [np.nan]*100;b=wide[alt] if alt in wide else [np.nan]*100;r=metric(a,b,100);r.update(model=model,analysis=stage,baseline=baseline,alternative=alt,interpretation='same fixed homologous exon intersection;Borzoi short=mask ablation' if stage=='CONTEXT' else 'same fixed homologous exon intersection');robustness.append(r)
   stability=[]
   for ID,z in wide.iterrows():
    x=z.dropna().astype(float);nz=x[x!=0];stability.append({'OFFREF_ID':ID,'model':model,'analysis':stage,'valid_settings_N':len(x),'nonzero_settings_N':len(nz),'dominant_sign_fraction':max(int((nz>0).sum()),int((nz<0).sum()))/len(nz) if len(nz) else None,'any_sign_flip':bool((nz>0).any() and (nz<0).any()) if len(nz)>=2 else None})
   pd.DataFrame(stability).to_csv(R/'11_statistics'/f'{model}_{stage}_LOCUS_STABILITY.tsv',sep='\t',index=False)
  table=[]
  for j in joblist:
   if j['stage']=='ENSEMBLE' and j['model']==model:table.append({'OFFREF_ID':j['job_ID'],'member':j['member'],'delta':values(results[j['audit_job_ID']]),'status':results[j['audit_job_ID']]['status']})
  table=pd.DataFrame(table);table.to_csv(R/'11_statistics'/f'{model}_ENSEMBLE_ALL_IDENTITIES.tsv',sep='\t',index=False);wide=table.pivot(index='OFFREF_ID',columns='member',values='delta');loci=[]
  for ID,z in wide.iterrows():
   x=z.dropna().astype(float);nz=x[x!=0];mean=float(x.mean()) if len(x) else None;sd=float(x.std(ddof=1)) if len(x)>1 else None;loci.append({'OFFREF_ID':ID,'model':model,'valid_members_N':len(x),'nonzero_members_N':len(nz),'mean_delta':mean,'SD_delta':sd,'variance':sd**2 if sd is not None else None,'CV':sd/abs(mean) if sd is not None and mean else None,'dominant_sign_fraction':max(int((nz>0).sum()),int((nz<0).sum()))/len(nz) if len(nz) else None,'fraction_same_as_mean_sign':float((np.sign(nz)==np.sign(mean)).mean()) if len(nz) and mean else None})
  pd.DataFrame(loci).to_csv(R/'11_statistics'/f'{model}_ENSEMBLE_LOCUS_STABILITY.tsv',sep='\t',index=False)
  members=wide.columns.tolist()
  for i,a in enumerate(members):
   for b in members[i+1:]:r=metric(wide[a],wide[b],100);r.update(model=model,memberA=a,memberB=b);ensemble.append(r)
 save('PLACEMENT_AND_CONTEXT_AGREEMENT.json',robustness);save('ENSEMBLE_PAIRWISE_AGREEMENT.json',ensemble)
 ref=readjson(R/'07_reference_benchmark/REFERENCE_PANEL_FROZEN.json');refrows=[]
 for j in ref:
  row={'ID':j['job_ID'],'gene':j['gene'],'beta':j['beta'],'PIP':j['PIP']}
  for model in ['AlphaGenome','Borzoi']:
   q=next(x for x in joblist if x['stage']=='REFERENCE' and x['model']==model and x['job_ID']==j['job_ID']);row[model]=values(results[q['audit_job_ID']]);row[model+'_TIER1']=values(results[q['audit_job_ID']],'TIER1_EXACT')
  refrows.append(row)
 refdf=pd.DataFrame(refrows);refdf.to_csv(R/'11_statistics/REFERENCE_VARIANT_SCORES.tsv',sep='\t',index=False);refsummary={'cross_model':metric(refdf.AlphaGenome,refdf.Borzoi,len(ref),refdf.gene),'strict_cross_model':metric(refdf.AlphaGenome_TIER1,refdf.Borzoi_TIER1,len(ref),refdf.gene)}
 for m in ['AlphaGenome','Borzoi']:refsummary[m+'_measured_beta']=metric(refdf[m+'_TIER1'],refdf.beta,len(ref),refdf.gene)
 refsummary['genes_N']=refdf.gene.nunique();save('REFERENCE_PERFORMANCE.json',refsummary)
 save('ANALYSIS_STATUS.json',{'status':'COMPLETE_ALL_FROZEN_IDENTITIES_REPORTED','jobs_N':len(joblist),'counts':dict(collections.Counter(z['status'] for z in results.values())),'Phase2':'NOT_ENTERED'})
if __name__=='__main__':main()
