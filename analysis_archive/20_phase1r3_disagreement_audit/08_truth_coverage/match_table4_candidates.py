from pathlib import Path
import pandas as pd,json,re
R=Path(__file__).resolve().parents[1];O=R/'08_truth_coverage'
s=pd.read_csv(O/'EXTERNAL_PUBLIC_SVEQTL_TABLE4_AB.tsv',sep='\t',dtype=str).fillna('');c=pd.read_csv(R.parent/'04_offref_catalog/OFFREF_MASTER.tsv',sep='\t',dtype=str).fillna('');c=c[c.class_A_B_C=='CLASS_B'];assert len(c)==32969;lookup={}
for z in c.to_dict('records'):lookup.setdefault((z['chrom'],int(z['anchor_left']),int(z['length'])),[]).append(z)
rows=[]
for z in s.to_dict('records'):
 vid=z['SVID'] or z['ID'];m=re.fullmatch(r'(chr[^-]+)-(\d+)-INS-(\d+)',vid)
 if not m:continue
 chrom,pos,L=m.group(1),int(m.group(2)),int(m.group(3));found={}
 for a in [pos-1,pos]:
  for x in lookup.get((chrom,a,L),[]):found[x['OFFREF_ID']]=x
 for x in found.values():rows.append({'external_variant_ID':vid,'source_sheet':z['source_sheet'],'OFFREF_ID':x['OFFREF_ID'],'external_gene':z['ensemblID'],'external_slope':z['slope'],'external_qvalue':z['qvalue'],'external_ID_POS1':pos,'catalogue_anchor_left':x['anchor_left'],'length':L,'match_category':'ANCHOR_LENGTH_CANDIDATE_SEQUENCE_UNVERIFIED','effect_allele':'UNVERIFIED','independence':'KNOWN_COHORT_OVERLAP_NO_DISJOINT_SUMMARY','use_as_direction_truth':False})
pd.DataFrame(rows,columns=list(rows[0]) if rows else ['external_variant_ID','OFFREF_ID','match_category']).to_csv(O/'EXTERNAL_ANCHOR_LENGTH_CANDIDATES.tsv',sep='\t',index=False)
p=O/'EXTERNAL_PUBLIC_TABLE_STATUS.json';a=json.loads(p.read_text());a.update(status='EXTERNAL_PUBLIC_SUMMARY_AVAILABLE_DIRECTION_UNVERIFIED_FINAL',catalogue_ClassB_N=len(c),anchor_length_candidate_rows_N=len(rows),candidate_ClassB_loci_N=len({x['OFFREF_ID'] for x in rows}),exact_sequence_match_N='UNDETERMINED',independent_direction_truth_N='UNDETERMINED',orientation_reason='No explicit effect allele in Table4. PLINK→tensorQTL dosage orientation not specified in downloaded Methods. Do not assume ALT dosage; no disjoint-donor summary.',final_audit_complete=True,search_stopped=True);p.write_text(json.dumps(a,indent=2));print('EXTERNAL_FINAL_AUDIT',len(rows),'candidate rows,',a['candidate_ClassB_loci_N'],'loci;not directional truth',flush=True)
