from pathlib import Path
import sys,json,collections,math
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'00_protocol'));from audit_common import *
assert_freezes();full=pd.read_csv(R/'02_track_harmonization/TRACK_METADATA_FULL.tsv',sep='\t').fillna('');mapping=readjson(R/'02_track_harmonization/TRACK_MAPPING_FROZEN.json');jobs=[j for j in readjson(OLD/'00_protocol/JOBS_FROZEN.json') if j['role']=='SENTINEL'];rows=[];coverage=[];truth_source=pd.read_csv(P/'17_phase1r2/03_gnrs_truth/GNRS_CLASSB_PRIMARY33.tsv',sep='\t').set_index('OFFREF_ID')
for j in jobs:
 row={'OFFREF_ID':j['job_ID'],'panel':j['panel'],'gene':j['gene'],'strand':j['strand'],'anchor0':j['anchor0'],'length':len(j['insertion'])}
 for m in ['AlphaGenome','Borzoi']:
  z=readjson(OLD/'10_statistics'/m/(j['job_ID']+'.json'));row[m+'_original']=z['common_RNA'].get('delta');row[m+'_native']=z['native_RNA'].get('delta');row[m+'_common_bp']=z['common_RNA']['mapped_bp'];row[m+'_native_bp']=z['native_RNA']['mapped_bp'];row[m+'_log1p']=z['common_RNA'].get('log1p_delta_sensitivity');score=z['common_RNA'];fams=[]
  if score['status']=='PASS':
   meta=readjson(OLD/'09_predictions'/m/(z['cache_REF']+'.json'))['metadata'];meta=[t for t in meta if t['strand']==j['strand']];deltas=score['track_deltas'];assert len(deltas)==len(meta)
   ids=[int(t['track_index']) if m=='Borzoi' else int(full[(full.model==m)&(full.track_name==t['name'])&(full.strand==t['strand'])].iloc[0].track_ID) for t in meta]
   for fam in mapping['TIER1_EXACT']['shared_families']:
    select=[k for k,x in enumerate(ids) if x in mapping['TIER1_EXACT']['tracks'][m][fam]];fams.append(float(np.mean([deltas[k] for k in select])) if select else None)
   row[m+'_tier1_cached']=float(np.mean(fams)) if fams and all(x is not None for x in fams) else None
  else:row[m+'_tier1_cached']=None
 if j['panel']=='BIOLOGICAL_TRUTH':
  exons=merge(j['exons']);ann=bool(exons and j['strand'] in ['+','-']);A=j['anchor0'];L=len(j['insertion']);bdom=domain(A,L,A-262144,524288,163840,32);rawdom=domain(A,L,A-262144,524288,163840,32,False);common=sum(max(0,min(b,y)-max(a,x)) for a,b,_ in bdom for x,y in exons);rawbp=sum(max(0,min(b,y)-max(a,x)) for a,b,_ in rawdom for x,y in exons);inside_input=sum(max(0,min(y,A+262144-L)-max(x,A-262144)) for x,y in exons)
  reason='PASS' if common else 'TRANSCRIPT_MAPPING' if not ann else 'WINDOW_CROP' if inside_input>0 else 'NO_COMMON_EXON';reason='INSUFFICIENT_HOMOLOGY' if common==0 and rawbp>0 else reason
  span=max(A+L,exons[-1][1])-min(A,exons[0][0]) if exons else None
  coverage.append({'OFFREF_ID':j['job_ID'],'frozen_gene':j['gene'],'strand':j['strand'],'source_beta':float(truth_source.loc[j['job_ID'],'beta_presence']),'primary_reason':reason,'annotation_available':ann,'model_prediction_PASS_both':True,'exon_union_bp':sum(b-a for a,b in exons),'strict_common_exon_bp':common,'fractional_common_exon_bp':rawbp,'Borzoi_input_exon_bp':inside_input,'AlphaGenome_native_exon_bp':row['AlphaGenome_native_bp'],'Borzoi_native_exon_bp':row['Borzoi_native_bp'],'insertion_and_gene_span_bp':span,'target_gene_changed':False})
 rows.append(row)
pd.DataFrame(rows).to_csv(R/'11_statistics/CACHED_ORIGINAL_AND_TIER1_SCORES.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(R/'08_truth_coverage/GNRS_33_COVERAGE_AUDIT.tsv',sep='\t',index=False)
(R/'08_truth_coverage/GNRS_COVERAGE_MECHANISMS.json').write_text(json.dumps({'frozen_N':33,'original_common_N':sum(z['strict_common_exon_bp']>0 for z in coverage),'primary_reasons':dict(collections.Counter(z['primary_reason'] for z in coverage)),'unchanged_gene_and_truth':True,'prospective_recovery_jobs':'GNRS_GENE_CENTER fixed all33, no outcome selection'},indent=2));print('CACHED_TIER1_AND_33_COVERAGE_AUDITS_COMPLETE',flush=True)
