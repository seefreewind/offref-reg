from pathlib import Path
import json,sys
import numpy as np,pandas as pd
from sklearn.linear_model import LogisticRegression,Ridge
from sklearn.model_selection import StratifiedKFold,KFold,cross_val_predict
from sklearn.metrics import roc_auc_score,r2_score
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parents[1];sys.path[:0]=[str(R/'00_protocol'),str(R/'11_statistics')];from audit_common import *
from statistics_core import interval

def main():
 assert_freezes();labels=readjson(R/'15_freeze/CORRECTED_LABEL_INPUT_SHA256.json');assert sha(R/labels['path'])==labels['sha256'];s=pd.read_csv(R/'06_offref_core/CORE400_DERIVED_STRATA_CORRECTED.tsv',sep='\t');cached=pd.read_csv(R/'11_statistics/CACHED_ORIGINAL_AND_TIER1_SCORES.tsv',sep='\t');s=s.merge(cached[['OFFREF_ID','AlphaGenome_tier1_cached','Borzoi_tier1_cached']],on='OFFREF_ID');a=s.AlphaGenome_tier1_cached;b=s.Borzoi_tier1_cached;finite=a.notna()&b.notna();s['disagreement_magnitude']=np.nan;s.loc[finite,'disagreement_magnitude']=abs((a[finite]-a[finite].mean())/a[finite].std(ddof=1)-(b[finite]-b[finite].mean())/b[finite].std(ddof=1));signed=finite&(a!=0)&(b!=0);s['direction_disagree']=np.nan;s.loc[signed,'direction_disagree']=(np.sign(a[signed])!=np.sign(b[signed])).astype(int)
 X=pd.DataFrame(index=s.index);X['log10_length']=np.log10(s.length);X['repeat_fraction']=s.repeat_fraction;X['GC']=s.GC;X['log10_TSS_distance']=np.log10(1+s.distance_to_TSS);X['mappability']=s.mappability;X['AF']=s.SV_frequency;X['AF_missing']=s.AF_MISSING.astype(int);X['gene_density']=s.gene_density;X['local_repeat_density']=s.local_repeat_density;X['model_output_coverage']=s.Borzoi_coverage;X['log10_mapped_exon_bp']=np.log10(1+s.Borzoi_mapped_bp);X['homology_fraction']=[sum(y-x for x,y,_ in domain(int(z.anchor),int(z.length),int(z.anchor)-262144,524288,163840,32))/196608 for z in s.itertuples()]
 for c in ['SINE/Alu','LINE','LTR','VNTR/STR','OTHER_REPEAT']:X['repeat_'+c]=(s.repeat_stratum==c).astype(int)
 for c in list(X):
  if X[c].isna().any():X[c+'_missing']=X[c].isna().astype(int);X[c]=X[c].fillna(X[c].median())
 centers=X.mean();scales=X.std(ddof=1).replace(0,1);X=(X-centers)/scales;names=list(X);rng=np.random.default_rng(20261011);rows=[];models={};B=2000
 for outcome,mask,model in [('direction_disagree',signed,LogisticRegression(C=1,max_iter=2000)),('disagreement_magnitude',finite,Ridge(alpha=1))]:
  x=X.loc[mask].to_numpy();y=s.loc[mask,outcome].to_numpy();model.fit(x,y);coef=model.coef_.ravel();boots=[]
  for _ in range(B):
   ix=rng.integers(0,len(y),len(y))
   if outcome=='direction_disagree' and len(np.unique(y[ix]))<2:continue
   model.fit(x[ix],y[ix]);boots.append(model.coef_.ravel())
  boots=np.array(boots)
  for k,n in enumerate(names):
   vals=boots[:,k];ci=interval(vals);p=min(1,2*min((np.sum(vals<=0)+1)/(len(vals)+1),(np.sum(vals>=0)+1)/(len(vals)+1)));rows.append({'outcome':outcome,'feature':n,'effect_per_standardized_covariate':float(coef[k]),'CI95_lower':ci[0],'CI95_upper':ci[1],'P_bootstrap_descriptive':p,'usable_N':len(y),'bootstrap_usable_N':len(boots)})
  if outcome=='direction_disagree':pred=cross_val_predict(LogisticRegression(C=1,max_iter=2000),x,y,cv=StratifiedKFold(5,shuffle=True,random_state=20261011),method='predict_proba')[:,1];cv={'ROC_AUC':roc_auc_score(y,pred)}
  else:pred=cross_val_predict(Ridge(alpha=1),x,y,cv=KFold(5,shuffle=True,random_state=20261011));cv={'R2':r2_score(y,pred)}
  models[outcome]={'N':len(y),'features':names,'condition_number':float(np.linalg.cond(x)),'matrix_rank':int(np.linalg.matrix_rank(x)),'prespecified_regularization':'L2 C1' if outcome=='direction_disagree' else 'ridge alpha1','cross_validation_exploratory':cv,'homology_fraction_definition':'actual shared native B output base count /196608 from frozen complete insertion and homology domains;not imputed1','constant_features':[n for n in names if X[n].std()==0]}
 q=multipletests([z['P_bootstrap_descriptive'] for z in rows],method='fdr_bh')[1]
 for z,v in zip(rows,q):z['BH_q']=float(v)
 pd.DataFrame(rows).to_csv(R/'09_disagreement_models/FEATURE_EFFECTS.tsv',sep='\t',index=False);s.to_csv(R/'09_disagreement_models/CORE400_DISAGREEMENT.tsv',sep='\t',index=False);X.to_csv(R/'09_disagreement_models/STANDARDIZED_DESIGN.tsv',sep='\t',index=False);(R/'09_disagreement_models/REGULARIZED_MODEL_QA.json').write_text(json.dumps(models,indent=2));(R/'09_disagreement_models/DESIGN_SCALING.json').write_text(json.dumps({'imputation_and_scaling_population':'all400 without outcome filtering','centers':centers.to_dict(),'scales':scales.to_dict(),'logistic_effect':'log odds of disagreement; penalized descriptive sensitivity','ridge_effect':'standardized-delta disagreement magnitude units; penalized descriptive sensitivity'},indent=2))
 # Unadjusted continuous length slope for agreement, same fixed estimator and bootstrap.
 x=X.loc[signed,['log10_length']].to_numpy();y=1-s.loc[signed,'direction_disagree'].to_numpy();m=LogisticRegression(C=1,max_iter=2000).fit(x,y);effect=float(m.coef_[0,0]);vals=[]
 for _ in range(B):
  ix=rng.integers(0,len(y),len(y));vals.append(float(LogisticRegression(C=1,max_iter=2000).fit(x[ix],y[ix]).coef_[0,0]))
 (R/'09_disagreement_models/LENGTH_CONTINUOUS_AGREEMENT.json').write_text(json.dumps({'N':len(y),'effect_log_odds_agreement_per_standardized_log10_length':effect,'CI95':interval(vals),'preset':'C1;2000 paired locus bootstrap;not monotonic claim by default'},indent=2));print('DISAGREEMENT_FEATURE_MODELS_COMPLETE',flush=True)
if __name__=='__main__':main()
