from pathlib import Path
import sys,json,math,collections
import numpy as np,pandas as pd
from scipy.stats import binomtest,rankdata
from statistics_core import rowcorr
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parents[1];sys.path[:0]=[str(R/'00_protocol'),str(R/'11_statistics')];from audit_common import *
from statistics_core import metric,boot,interval

def write(folder,name,data):(R/folder/name).write_text(json.dumps(data,indent=2))
def contrast(a,b,variable,label1,label2):
 rng=np.random.default_rng(20261011);a=a[['AlphaGenome_original','Borzoi_original']].dropna();b=b[['AlphaGenome_original','Borzoi_original']].dropna();av=a.to_numpy();bv=b.to_numpy();out={'variable':variable,'group1':label1,'group2':label2,'N1':len(a),'N2':len(b)}
 if len(a)<3 or len(b)<3:out['status']='TOO_SPARSE';return out
 A=boot(av[:,0],av[:,1]);B=boot(bv[:,0],bv[:,1],seed=20261012);rho=np.array(A['rho'])-B['rho'];sign=np.array(A['sign'])-B['sign'];out.update(status='DESCRIPTIVE_PAIRED_LOCUS_BOOTSTRAP',rho_difference_CI95=interval(rho),sign_difference_CI95=interval(sign),sign_difference=float(np.mean(np.sign(av[:,0])==np.sign(av[:,1]))-np.mean(np.sign(bv[:,0])==np.sign(bv[:,1]))),P_bootstrap_sign=min(1,2*min((np.sum(sign<=0)+1)/10001,(np.sum(sign>=0)+1)/10001)));return out

def main():
 assert_freezes();harm=readjson(R/'11_statistics/HARMONIZATION_AGREEMENT.json');refsum=readjson(R/'11_statistics/REFERENCE_PERFORMANCE.json');rob=readjson(R/'11_statistics/PLACEMENT_AND_CONTEXT_AGREEMENT.json');core=pd.read_csv(R/'11_statistics/CORE400_ALL_SCORING_DEFINITIONS.tsv',sep='\t');refs=pd.read_csv(R/'11_statistics/REFERENCE_VARIANT_SCORES.tsv',sep='\t');cross=harm['tier1_cached']
 # Strict exact-tier cached arrays contain every strand-compatible exact-family source track, avoiding full-global transport availability becoming a scoring exclusion.
 write('02_track_harmonization','STRICT_TIER1_PROVENANCE.json',{'primary_phase1r3_diagnostic':'cached exact-family scoring on all400; same frozen arrays and unchanged common regions','scored_N':cross['finite_N'],'new_full_global_requests':'additional model-wide sensitivity; trackwise equal-family recomputation cross-check wherever both obtained','original_Phase1R2_unchanged':True,'tier2_rule':'same original GM12878 tracks; not a new independent experiment'})
 gap={'status':'INSUFFICIENT_REFERENCE_OUTPUT','rho_difference':None,'sign_difference':None,'rho_difference_CI95':[None,None],'sign_difference_CI95':[None,None],'reference_genes_usable_N':0}
 good=refs[['AlphaGenome_TIER1','Borzoi_TIER1']].notna().all(axis=1);rg=refs.loc[good];cg=core[['AlphaGenome_tier1_cached','Borzoi_tier1_cached']].dropna();gap['reference_genes_usable_N']=rg.gene.nunique()
 if len(rg)>=3:
  rv=boot(rg.AlphaGenome_TIER1,rg.Borzoi_TIER1,cluster=rg.gene);cv=boot(cg.AlphaGenome_tier1_cached,cg.Borzoi_tier1_cached,seed=20261012);rd=np.array(rv['rho'])-cv['rho'];sd=np.array(rv['sign'])-cv['sign'];rm=refsum['strict_cross_model'];gap.update(status='GENE_CLUSTER_REFERENCE_AND_LOCUS_OOD_BOOTSTRAP',rho_difference=rm['rho']-cross['rho'] if rm['rho'] is not None else None,sign_difference=rm['proportion']-cross['proportion'] if rm['proportion'] is not None else None,rho_difference_CI95=interval(rd),sign_difference_CI95=interval(sd),reference_N=len(rg),core_N=len(cg),independent_bootstrap_streams=[20261011,20261012],design_limit='GEUVADIS high-PIP SNVs vs naturally occurring complete insertions;ascertainment, effect type and target selection differ;not causal OOD comparison')
 write('11_statistics','REFERENCE_OFFREF_GAP.json',gap)
 contrasts=[contrast(core[core.repeat_stratum=='SINE/Alu'],core[core.repeat_stratum=='NON_REPEAT'],'repeat_class','SINE/Alu','NON_REPEAT'),contrast(core[core.repeat_fraction>=.8],core[core.repeat_fraction<.8],'repeat_fraction','repeat-rich','low-repeat')]
 valid=[z for z in contrasts if z.get('P_bootstrap_sign') is not None]
 if valid:
  for z,q in zip(valid,multipletests([x['P_bootstrap_sign'] for x in valid],method='fdr_bh')[1]):z['BH_q']=float(q)
 write('11_statistics','PRESPECIFIED_OOD_CONTRASTS.json',contrasts)
 # Fixed gene-centered scoring: all33, unchanged selected gene and observed beta. Native/fractional scores retained separately in coverage audit.
 truth_source=pd.read_csv(P/'17_phase1r2/03_gnrs_truth/GNRS_CLASSB_PRIMARY33.tsv',sep='\t').set_index('OFFREF_ID');jobs=readjson(R/'00_protocol/PREDICTION_JOBS_FROZEN.json');orig={x['job_ID']:x for x in readjson(OLD/'00_protocol/JOBS_FROZEN.json') if x['role']=='SENTINEL' and x['panel']=='BIOLOGICAL_TRUTH'};rows=[]
 for ID,j in orig.items():
  row={'OFFREF_ID':ID,'gene':j['gene'],'beta':float(truth_source.loc[ID,'beta_presence']),'target_gene_changed':False}
  for m in ['AlphaGenome','Borzoi']:
   q=next(x for x in jobs if x['stage']=='GNRS_GENE_CENTER' and x['model']==m and x['job_ID']==ID);p=R/'08_truth_coverage'/m/(q['audit_job_ID']+'.json');z=readjson(p) if p.exists() else {'status':'MISSING'};row[m+'_status']=z['status'];row[m]=z.get('tier1_delta') if z['status']=='PASS' else None;row[m+'_mapped_bp']=z.get('REF',{}).get('mapped_bp',0);row[m+'_reason']=z.get('error')
  rows.append(row)
 g=pd.DataFrame(rows);g.to_csv(R/'08_truth_coverage/GNRS_GENE_PRESERVING_RECOVERY.tsv',sep='\t',index=False);gnrs={m:metric(g[m],g.beta,33) for m in ['AlphaGenome','Borzoi']};N=int(g[['AlphaGenome','Borzoi']].notna().all(axis=1).sum());gnrs['joint_recovered_N']=N;gnrs['coverage_verdict']='GNRS_LOW_OUTPUT_COVERAGE' if N<20 else 'GNRS_COVERAGE_RECOVERED';gnrs['original_N7_unchanged']=True;write('11_statistics','GNRS_RECOVERY_DIRECTION.json',gnrs)
 # Make temporal serving access claims from actual member predictions, not exposed enum alone.
 serving={}
 for m in ['AlphaGenome','Borzoi']:
  t=pd.read_csv(R/'11_statistics'/f'{m}_ENSEMBLE_ALL_IDENTITIES.tsv',sep='\t');serving[m]={'member_attempt_counts':t.groupby('member',dropna=False).size().to_dict(),'member_PASS_counts':t[t.status=='PASS'].groupby('member',dropna=False).size().to_dict(),'planned_members_N':4,'actual_exposed_member_predictions':bool((t.status=='PASS').any()),'biological_independence':False,'Borzoi_splits':'same test3 validation4 across four replicates' if m=='Borzoi' else 'fold-specific service;not ALL_FOLDS technical rerun variance'}
 write('03_within_model_ensembles','OFFICIAL_MEMBER_AVAILABILITY_FINAL.json',serving)
 effects=pd.read_csv(R/'09_disagreement_models/FEATURE_EFFECTS.tsv',sep='\t');supported=effects[(effects.BH_q<.05)&((effects.CI95_lower>0)|(effects.CI95_upper<0))];strong=supported.loc[supported.effect_per_standardized_covariate.abs().idxmax()].to_dict() if len(supported) else None
 ood_supported=supported[supported.feature.astype(str).str.contains('length|repeat',case=False)] if 'feature' in supported else supported.iloc[:0]
 # Gene-cluster sign-flip null respects variants in the same reference gene.
 measured_tests=[];rng=np.random.default_rng(20261011)
 for m in ['AlphaGenome','Borzoi']:
  score_column=m+'_TIER1';t=refs[[score_column,'beta','gene']].dropna();a=t[score_column].to_numpy(float);b=t.beta.to_numpy(float)
  if len(t)<3 or not np.ptp(a) or not np.ptp(b):measured_tests.append({'model':m,'P_rho':1.0,'P_sign':1.0});continue
  _,groups=np.unique(t.gene.to_numpy(),return_inverse=True);G=int(groups.max()+1);rp=[];sp=[];yr=rankdata(b);signed=(a!=0)&(b!=0)
  for step in range(0,10000,500):
   flips=rng.choice([-1.,1.],size=(500,G))[:,groups];x=flips*a;rp.extend(rowcorr(rankdata(x,axis=1),np.broadcast_to(yr,x.shape)));sp.extend(np.mean(np.sign(x[:,signed])==np.sign(b[signed]),axis=1) if signed.any() else np.full(500,np.nan))
  observed=refsum[m+'_measured_beta'];pr=(1+np.sum(np.abs(rp)>=abs(observed['rho'])))/10001;ps=(1+np.sum(np.abs(np.asarray(sp)-.5)>=abs(observed['proportion']-.5)))/10001 if observed['proportion'] is not None else 1.0;measured_tests.append({'model':m,'P_rho':float(pr),'P_sign':float(ps),'null':'10000 gene-cluster sign flips; two-sided'})
 for metric_name in ['rho','sign']:
  for z,q in zip(measured_tests,multipletests([x['P_'+metric_name] for x in measured_tests],method='fdr_bh')[1]):z['BH_q_'+metric_name]=float(q)
 write('10_nulls','REFERENCE_MEASURED_GENE_CLUSTER_NULL.json',measured_tests)
 # Scope-specific conservative operational gate, as frozen before outputs.
 rm=refsum['strict_cross_model'];adequate=rm['finite_N']>=100 and gap['reference_genes_usable_N']>=20;acceptable_agreement=adequate and rm['rho'] is not None and rm['rho']>=.2 and rm['proportion'] is not None and rm['proportion']>=.6
 measured_positive=[]
 for m in ['AlphaGenome','Borzoi']:
  z=refsum[m+'_measured_beta'];ci=z['rho_CI95'];sci=z['sign_bootstrap_CI95'];test=next(x for x in measured_tests if x['model']==m);measured_positive.append(bool((ci[0] is not None and ci[0]>0 and test['BH_q_rho']<.05) or (sci[0] is not None and sci[0]>.5 and test['BH_q_sign']<.05)))
 reference_acceptable=acceptable_agreement and any(measured_positive)
 severe_gap=bool((gap['rho_difference'] is not None and gap['rho_difference']>=.2 and gap['rho_difference_CI95'][0]>0) or (gap['sign_difference'] is not None and gap['sign_difference']>=.1 and gap['sign_difference_CI95'][0]>0))
 instability=any(z['analysis']=='PLACEMENT' and z['signed_N']>=20 and z['proportion'] is not None and 1-z['proportion']>=.2 and z['exact_CI95'][1]<.9 for z in rob)
 artifact=bool(cross['rho'] is not None and cross['rho']>=.2 and cross['proportion']>=.6 and cross['rho']-harm['original']['rho']>=.2)
 truth_supported=N>=20 and all(gnrs[m]['exact_CI95'][0] is not None and gnrs[m]['exact_CI95'][0]>.5 for m in ['AlphaGenome','Borzoi'])
 adjusted_survives=None # Not estimated: disagreement adjustment is not adjusted biological-truth validation.
 restored=artifact and truth_supported and not severe_gap and adjusted_survives
 if restored:verdict='OFFREF_ATLAS_GO_RESTORED'
 elif artifact:verdict='TECHNICAL_ARTIFACT_RESOLVED'
 elif reference_acceptable and (severe_gap or instability or len(ood_supported)>0) and cross['rho'] is not None and cross['rho']<=0:verdict='MANUSCRIPT_GO_OOD_RELIABILITY_AUDIT'
 elif adequate and not any(measured_positive) and not truth_supported:verdict='PIPELINE_OR_MODEL_CALIBRATION_FAILURE'
 else:verdict='ARCHIVE_INCONCLUSIVE'
 gate={'FINAL_VERDICT':verdict,'reference_adequate':adequate,'reference_acceptable':reference_acceptable,'reference_measured_positive_CI_models':measured_positive,'reference_agreement_acceptable':acceptable_agreement,'severe_gap':severe_gap,'placement_instability':instability,'technical_artifact_resolved':artifact,'GNRS_joint_N':N,'GNRS_supported':truth_supported,'adjusted_biological_calibration_supported':adjusted_survives,'adjustment_scope':'Full covariate models estimate model disagreement; do not infer adjusted GNRS biological validation from them','strongest_uncertainty_supported_predictor':strong,'length_repeat_supported_features_N':len(ood_supported),'Phase2':'NOT_ENTERED','provisional_manuscript_titles_permitted':verdict=='MANUSCRIPT_GO_OOD_RELIABILITY_AUDIT','source_protocol':'00_protocol/PHASE1R3_PROTOCOL.md','scope':'No special OFFREF failure claimed from missing/weak reference benchmark'};write('13_manuscript_gate','FINAL_SCIENTIFIC_ROUTE.json',gate);print('ENDPOINTS_AND_SCIENTIFIC_ROUTE_FINALIZED',verdict,flush=True)
if __name__=='__main__':main()
