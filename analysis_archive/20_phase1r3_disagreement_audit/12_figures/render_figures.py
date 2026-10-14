from pathlib import Path
import json,sys
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','legend.frameon':False})
COL={'AlphaGenome':'#537B9D','Borzoi':'#B87654'}
def read(p):return json.loads((R/p).read_text())
def table(p):return pd.read_csv(R/p,sep='\t')
def save(fig,N,source,caption):
 stem=R/'12_figures'/f'Figure_{N:02d}';source.to_csv(str(stem)+'_source.tsv',sep='\t',index=False)
 fig.savefig(str(stem)+'.pdf')
 fig.savefig(str(stem)+'.svg')
 fig.savefig(str(stem)+'.png',dpi=300)
 plt.close(fig);return {'figure':N,'caption':caption,'source_rows':len(source),'exports':['PDF','SVG','PNG300dpi','sourceTSV']}
def scatter(ax,x,y,xlabel,ylabel,title):
 a=np.asarray(x,float);b=np.asarray(y,float);good=np.isfinite(a)&np.isfinite(b);ax.scatter(a[good],b[good],s=9,alpha=.5,c=COL['AlphaGenome'],edgecolors='none',rasterized=True);ax.axhline(0,color='.75',lw=.5);ax.axvline(0,color='.75',lw=.5);ax.set(xlabel=xlabel,ylabel=ylabel,title=f'{title}\nPaired finite n={good.sum()}')
def forest(ax,records,key,ci,title):
 for i,z in enumerate(records):
  value=z.get(key);limits=z.get(ci,[None,None]);label=z.get('label',z.get('stratum',str(i)))
  if value is not None and np.isfinite(value):
   ax.plot(value,i,'o',ms=4,color='#537B9D')
   if limits[0] is not None:ax.plot(limits,[i,i],color='#537B9D',lw=1)
  else:ax.text(.02,i,'Not estimable',fontsize=7)
 ax.set_yticks(range(len(records)),[z.get('label',z.get('stratum',str(i)))+f" (n={z.get('finite_N',0)}/{z.get('frozen_N','?')})" for i,z in enumerate(records)]);ax.axvline(0,color='.8',lw=.7);ax.set(xlabel=title);ax.invert_yaxis()
def main():
 status=read('11_statistics/ANALYSIS_STATUS.json');assert status['counts'].get('MISSING',0)==0,'Unfinished inference cannot be labeled final figures'
 core=table('11_statistics/CORE400_ALL_SCORING_DEFINITIONS.tsv');assert (core.length>0).all();harm=read('11_statistics/HARMONIZATION_AGREEMENT.json');strata=read('11_statistics/OOD_STRATA.json');refs=table('11_statistics/REFERENCE_VARIANT_SCORES.tsv');refsummary=read('11_statistics/REFERENCE_PERFORMANCE.json');gnrs=table('08_truth_coverage/GNRS_GENE_PRESERVING_RECOVERY.tsv');out=[]
 fig,ax=plt.subplots(figsize=(7.0866,3.2),layout='constrained');design=pd.DataFrame({'population':['Unbiased core','GNRS truth','Random extension','Fixed robustness subset','Reference SNVs'],'frozen_N':[400,33,67,100,500]});ax.barh(design.population,design.frozen_N,color=['#537B9D','#B87654','#A7A7A7','#88A59B','#88A59B']);ax.set(xlabel='Frozen identities',title='Prospective audit populations; robustness subset belongs to core');ax.invert_yaxis();out.append(save(fig,1,design,'Frozen populations; no outcome-based replacement.'))
 fig,axes=plt.subplots(2,2,figsize=(7.0866,5.2),layout='constrained');tiers=['original','tier1_cached','TIER2_NEAR','TIER3_GLOBAL']
 for ax,t in zip(axes.flat,tiers):scatter(ax,core['AlphaGenome_'+t],core['Borzoi_'+t],'AlphaGenome ΔRNA','Borzoi ΔRNA',t+'; rho='+('NA' if harm[t]['rho'] is None else f"{harm[t]['rho']:.3f}"))
 out.append(save(fig,2,core,'Original, exact-family, near-family and global effects. Model units differ; each effect is ALT−REF per-bp RNA. Global completeness is required.'))
 stability=[]
 for m in COL:
  for stage in ['ENSEMBLE','PLACEMENT']:
   t=table(f'11_statistics/{m}_{stage}_LOCUS_STABILITY.tsv');t['category']=stage;stability.append(t)
 st=pd.concat(stability,ignore_index=True);fig,axes=plt.subplots(1,2,figsize=(7.0866,3.4),layout='constrained')
 for ax,stage in zip(axes,['ENSEMBLE','PLACEMENT']):
  for i,m in enumerate(COL):
   t=st[(st.model==m)&(st.category==stage)];valid_col='valid_members_N' if stage=='ENSEMBLE' else 'valid_settings_N';v=t.loc[t[valid_col]>=2,'dominant_sign_fraction'].dropna();ax.scatter(np.full(len(v),i)+np.linspace(-.15,.15,len(v)),v,s=10,alpha=.5,c=COL[m]);ax.text(i,1.025,f'≥2 valid {len(v)}/100',ha='center',fontsize=7)
  ax.set(xticks=[0,1],xticklabels=list(COL),ylim=(.45,1.1),ylabel='Dominant nonzero sign fraction',title=stage+' (missing members/settings explicit)')
 out.append(save(fig,3,st,'Dominant sign fraction across available official members or fixed placements; Loci with fewer than two valid members/settings are excluded from stability marks and retained in source rows. Source includes valid member counts.'))
 fig,axes=plt.subplots(1,2,figsize=(7.0866,3.5),layout='constrained');d=table('09_disagreement_models/CORE400_DISAGREEMENT.tsv');d=core[['OFFREF_ID']].merge(d,on='OFFREF_ID',validate='one_to_one');scatter(axes[0],np.log10(core.length),d.disagreement_magnitude,'log10 insertion bp','|standardized AG − standardized B|','Joint-finite standardized discrepancy');records=[z for z in strata if z['variable']=='length_bin'];forest(axes[1],records,'rho','rho_CI95','Cross-model Spearman rho, 95% locus CI');out.append(save(fig,4,core.merge(d[['OFFREF_ID','disagreement_magnitude']],on='OFFREF_ID',validate='one_to_one'),'Length strata retain all frozen loci. Standardization uses joint-finite core, not strata.'))
 fig,axes=plt.subplots(1,2,figsize=(7.0866,4),layout='constrained');records=[z for z in strata if z['variable']=='repeat_stratum'];forest(axes[0],records,'rho','rho_CI95','Spearman rho, 95% locus CI');forest(axes[1],records,'proportion','sign_bootstrap_CI95','Signed agreement, 95% locus CI');axes[1].axvline(.5,color='.6',ls=':',lw=.8);out.append(save(fig,5,pd.json_normalize(records),'Repeat classes fixed before audit; empty strata shown.'))
 fig,axes=plt.subplots(1,2,figsize=(7.0866,3.7),layout='constrained');rec=[dict(harm['tier1_cached'],label='Complete insertions'),dict(refsummary['strict_cross_model'],label='Reference SNVs')];forest(axes[0],rec,'rho','rho_CI95','Exact-tissue Spearman rho');forest(axes[1],rec,'proportion','sign_bootstrap_CI95','Signed agreement');axes[1].axvline(.5,color='.6',ls=':',lw=.8);out.append(save(fig,6,pd.concat([core[['OFFREF_ID','gene','AlphaGenome_tier1_cached','Borzoi_tier1_cached']].rename(columns={'OFFREF_ID':'ID','AlphaGenome_tier1_cached':'AlphaGenome','Borzoi_tier1_cached':'Borzoi'}).assign(population='OFFREF'),refs[['ID','gene','AlphaGenome_TIER1','Borzoi_TIER1']].rename(columns={'AlphaGenome_TIER1':'AlphaGenome','Borzoi_TIER1':'Borzoi'}).assign(population='REFERENCE')],ignore_index=True),'Same exact RNA tier and common-exon scoring. Reference confidence intervals resample genes; core intervals resample loci. Ascertainment/effect type differ.'))
 coverage=table('08_truth_coverage/GNRS_33_COVERAGE_AUDIT.tsv');reason=read('08_truth_coverage/GNRS_COVERAGE_MECHANISMS.json')['primary_reasons'];fig,axes=plt.subplots(1,2,figsize=(7.0866,3.5),layout='constrained');axes[0].barh(list(reason),list(reason.values()),color='#A7A7A7');axes[0].set(xlabel='Original loci',title='Original 33-locus coverage causes');counts=[int(gnrs[m].notna().sum()) for m in COL];axes[1].bar(list(COL),counts,color=list(COL.values()));axes[1].axhline(20,color='.5',ls=':');axes[1].set(ylim=(0,35),ylabel='Finite gene-preserving scores /33',title='Recovery eligibility (target ≥20)');out.append(save(fig,7,coverage.merge(gnrs,on='OFFREF_ID',suffixes=('_original','_recovered')),'Frozen gene, insertion and direction unchanged. Original common N7 is retained.'))
 fig,axes=plt.subplots(1,2,figsize=(7.0866,3.5),layout='constrained');joint=int(gnrs[list(COL)].notna().all(axis=1).sum())
 for ax,m in zip(axes,COL):
  if joint>=20:scatter(ax,gnrs.beta,gnrs[m],'Measured insertion-presence beta',m+' ΔRNA',m+' gene-preserving GNRS')
  else:ax.axis('off');ax.text(.5,.6,f'GNRS_LOW_OUTPUT_COVERAGE\nJoint finite n={joint}/33\nDirectional calibration threshold:20',ha='center',va='center',transform=ax.transAxes,fontsize=9)
 out.append(save(fig,8,gnrs,'Directional scatter is displayed only if joint finite N≥20. Otherwise no directional-calibration inference is plotted.'))
 robust=read('11_statistics/PLACEMENT_AND_CONTEXT_AGREEMENT.json');fig,axes=plt.subplots(1,2,figsize=(7.0866,4.1),layout='constrained')
 for ax,analysis in zip(axes,['PLACEMENT','CONTEXT']):
  rec=[dict(z,label=z['model']+' '+str(z['alternative'])+(' masked bp' if z['model']=='Borzoi' and analysis=='CONTEXT' else ' bp')) for z in robust if z['analysis']==analysis];forest(ax,rec,'proportion','sign_bootstrap_CI95',analysis+' signed agreement vs baseline')
  ax.axvline(.5,color='.6',ls=':',lw=.8)
 out.append(save(fig,9,pd.json_normalize(robust),'All comparisons use the same fixed homologous exon intersection. Borzoi short contexts are N-mask ablations inside native524288bp.'))
 examples=d.dropna(subset=['disagreement_magnitude']).sort_values(['disagreement_magnitude','OFFREF_ID'],ascending=[False,True]).head(4);fig,axes=plt.subplots(2,2,figsize=(7.0866,4.8),layout='constrained');a=core.AlphaGenome_original;b=core.Borzoi_original;good=a.notna()&b.notna();za=(a[good]-a[good].mean())/a[good].std(ddof=1);zb=(b[good]-b[good].mean())/b[good].std(ddof=1);ex=[]
 for ax,(_,z) in zip(axes.flat,examples.iterrows()):
  idx=core.index[core.OFFREF_ID==z.OFFREF_ID][0];q=core.loc[idx];ax.bar(list(COL),[za.loc[idx],zb.loc[idx]],color=list(COL.values()));ax.axhline(0,color='.6',lw=.6);ax.set(ylabel='Joint-core standardized ΔRNA',title=f"{q.OFFREF_ID}\n{int(q.length)}bp | {q.repeat_stratum}",);ex.append({**z.to_dict(),'length':q.length,'repeat_stratum':q.repeat_stratum,'AlphaGenome_z':za.loc[idx],'Borzoi_z':zb.loc[idx],'selection':'post-hoc top absolute standardized discrepancy; illustrative'})
 out.append(save(fig,10,pd.DataFrame(ex),'Four highest-discrepancy complete insertions, selected post hoc from joint-finite core. These are illustrative and supply no independent biological validation.'))
 (R/'12_figures/FIGURE_LEGENDS_AND_QA.json').write_text(json.dumps({'status':'RENDERED_PENDING_MANUAL_FINAL_SIZE_QA','figures':out,'font_minimum_configured_pt':7,'source_integrity':'all frozen IDs retained in source tables; finite-data exclusions shown'},indent=2));print('TEN_FIGURES_RENDERED_MANUAL_QA_REQUIRED',flush=True)
if __name__=='__main__':main()
