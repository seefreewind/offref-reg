from pathlib import Path
import csv,json,hashlib,urllib.request,urllib.error,re,collections
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[1];P=R.parent;OLD=P/'19_phase1r2_execution'
def sha(s):return hashlib.sha256(s.encode()).hexdigest()
def save(name,obj):p=R/name;p.write_text(json.dumps(obj,ensure_ascii=False,indent=2))
def fetch(url,p):
 req=urllib.request.Request(url,headers={'User-Agent':'OFFREF-REG research audit','Accept':'application/json'})
 try:
  b=urllib.request.urlopen(req,timeout=40).read();p.write_bytes(b);return json.loads(b)
 except Exception as e:return {'access_error':type(e).__name__,'details':str(e)[:300],'url':url}
ag=pd.read_csv(OLD/'01_alphagenome/rna_seq_ALL_TRACKS.tsv',sep='\t').fillna('');b=pd.read_csv(P/'16_phase1r_rescue/03_borzoi/targets_human.txt',sep='\t',index_col=0).fillna('');b.index.name='track_index'
# Resolve only source experiment annotations, without any prediction/truth values.
exps=sorted(set(re.findall(r'ENCSR\w+', ' '.join(b[b.description=='RNA:GM12878'].file))));em={}
for e in exps:
 p=R/'01_model_metadata'/f'{e}.json';d=json.loads(p.read_text()) if p.exists() else fetch(f'https://www.encodeproject.org/experiments/{e}/?format=json',p);em[e]=d
 print('ENCODE_METADATA',e,'OK' if 'access_error' not in d else d['access_error'],flush=True)
rows=[]
for i,t in ag.iterrows():rows.append({'model':'AlphaGenome','track_ID':i,'track_name':t['name'],'cell_tissue':t['biosample_name'],'strand':t['strand'],'assay':'RNA_SEQ','RNA_modality':t['Assay title'],'selection_family':'polyA_positive' if 'polyA plus' in t['Assay title'] else 'total' if 'total RNA' in t['Assay title'] else 'OTHER','output_transform':'official unscaled RNA coverage','ensemble_status':'original ALL_FOLDS; FOLD_0..3 exposed; serving access to probe','training_source':t['data_source'],'ontology_curie':t['ontology_curie']})
for i,t in b.iterrows():
 desc=t.description;exp=re.findall('ENCSR\\w+',t.file);d=em.get(exp[0],{}) if exp else {};sel=d.get('library_selection_method',[]);depl=d.get('depleted_in_term_name',[]);text=' '.join(sel if isinstance(sel,list) else [sel]).lower();family='polyA_positive' if 'polya' in text and 'deplet' not in text else 'total' if ('random' in text or 'ribo' in text) and not any('polya' in str(x).lower() for x in depl) else 'UNKNOWN'
 title=d.get('assay_title','UNKNOWN');family='polyA_positive' if title=='polyA plus RNA-seq' else 'total' if title=='total RNA-seq' else family
 rows.append({'model':'Borzoi','track_ID':i,'track_name':t.identifier,'cell_tissue':desc.split(':',1)[-1],'strand':'+' if t.identifier.endswith('+') else '-' if t.identifier.endswith('-') else '.','assay':desc.split(':',1)[0],'RNA_modality':title if desc.startswith('RNA:') else 'NA','selection_family':family,'output_transform':f"legacy inverse softclip({t.clip_soft})/scale({t.scale})/power({'.75' if t.sum_stat=='sum_sqrt' else '1'}) then /32",'ensemble_status':'official replicates0–3; identical data splits','training_source':t.file,'ontology_curie':'UNKNOWN'})
full=pd.DataFrame(rows);full.to_csv(R/'02_track_harmonization/TRACK_METADATA_FULL.tsv',sep='\t',index=False)
exact={'AlphaGenome':{},'Borzoi':{}}
for m in exact:
 z=full[(full.model==m)&(full.cell_tissue=='GM12878')]
 for fam in ['polyA_positive','total']:
  exact[m][fam]=z[z.selection_family==fam].track_ID.astype(int).tolist()
shared=[fam for fam in ['polyA_positive','total'] if all(exact[m][fam] for m in exact)]
mapping={'TIER1_EXACT':{'shared_families':shared,'tracks':exact,'rule':'GM12878 identical assay/library-selection family; equal family mean; matching strand at scoring'},'TIER2_NEAR':{'rule':'All original available GM12878 RNA tracks; compatibility not same library selection'},'TIER3_GLOBAL':{'AlphaGenome':full[full.model=='AlphaGenome'].track_ID.astype(int).tolist(),'Borzoi':full[(full.model=='Borzoi')&(full.assay=='RNA')].track_ID.astype(int).tolist(),'rule':'Every official human strand-correct RNA track, equal track mean; unmatched tissue global sensitivity'},'outcomes_used_for_mapping':False}
save('02_track_harmonization/TRACK_MAPPING_FROZEN.json',mapping)
for file in ['TRACK_METADATA_FULL.tsv','TRACK_MAPPING_FROZEN.json']:
 p=R/'02_track_harmonization'/file;save('02_track_harmonization/'+file+'.sha256.json',{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'stage':'BEFORE_NEW_EFFECT_CALCULATION'})
s=pd.read_csv(OLD/'10_statistics/SENTINEL500_SCORES.tsv',sep='\t');s=s[s.panel=='UNBIASED_CLASSB_CORE'].copy();assert len(s)==400
rank=pd.concat([s[m].abs().rank(pct=True) for m in ['AlphaGenome','Borzoi']],axis=1).mean(axis=1);s['frozen_magnitude_quartile']=pd.cut(rank,[0,.25,.5,.75,1],labels=['Q1','Q2','Q3','Q4']).astype('object').fillna('MISSING');s['length_bin']=pd.cut(s.length,[49,199,499,999,1999,4999,float('inf')],labels=['50–199','200–499','500–999','1–2kb','2–5kb','>5kb']).astype(str)
def rc(t):
 t=str(t).upper()
 return 'SINE/Alu' if 'SINE' in t or 'ALU' in t else 'LINE' if 'LINE' in t else 'LTR' if 'LTR' in t else 'VNTR/STR' if any(x in t for x in ['VNTR','STR','TANDEM']) else 'NON_REPEAT' if any(x in t for x in ['NON','NONE']) else 'OTHER_REPEAT'
s['repeat_stratum']=s.repeat_class.map(rc);s['selection_hash']=s.OFFREF_ID.map(lambda x:sha('20261011|'+x));groups=[g.sort_values('selection_hash').to_dict('records') for _,g in s.groupby(['length_bin','repeat_stratum','frozen_magnitude_quartile'],observed=True)];sel=[]
while len(sel)<100:
 for g in groups:
  if g and len(sel)<100:sel.append(g.pop(0))
p=R/'04_context_robustness/SUBSET100_FROZEN.tsv';pd.DataFrame(sel).to_csv(p,sep='\t',index=False);save('04_context_robustness/SUBSET100_SHA256.json',{'N':100,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'seed':20261011,'frozen_Phase1R2_magnitude_quartiles_used':True,'new_results_used':False})
# Retain full fixed core and labels, not a new sample.
s.to_csv(R/'06_offref_core/CORE400_IDENTITIES.tsv',sep='\t',index=False);print('METADATA_MAPPING_AND_SUBSET_FROZEN',shared,flush=True)
