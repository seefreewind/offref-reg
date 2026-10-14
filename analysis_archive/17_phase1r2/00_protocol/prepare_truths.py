from pathlib import Path
import csv,json,gzip,hashlib,collections,math
R=Path(__file__).resolve().parents[1];P=R.parent;OLD=P/'16_phase1r_rescue';csv.field_size_limit(10**9)
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def write(p,rows,fields=None):
 with p.open('w') as f:w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
cov={r['OFFREF_ID']:r for r in read(OLD/'06_covariates/OFFREF_COVARIATE_COMPLETE.tsv')};seq={};v=None
for line in (OLD/'06_covariates/CLASS_B_INSERTIONS.fa').open():
 if line.startswith('>'):v=line[1:].strip();seq[v]=''
 else:seq[v]+=line.strip()
g=[r for r in read(OLD/'05_truth_gnrs/GNRS_DIRECTION_HARMONIZED.tsv') if r['class']=='CLASS_B' and r['QTL_DIRECTION_ELIGIBLE']=='TRUE'];groups=collections.defaultdict(list)
for r in g:groups[r['OFFREF_ID']].append(r)
assert len(groups)==33
fields=['OFFREF_ID','GNRS_variant_ID','GRCh38_chrom','GRCh38_anchor_0based','inserted_sequence','effect_allele','presence_allele','absence_allele','beta_presence','gene','transcript','population','distance_TSS','length','repeat_fraction','repeat_class','AF','AF_definition','orientation_evidence','eligibility_reason','beta_P','permutation_P','nominal_P','eligible_transcript_records_N']
def num(r,k):
 try:return float(r[k])
 except:return float('inf')
primary=[]
for vid,rr in sorted(groups.items()):
 a=min(rr,key=lambda r:(num(r,'beta_P'),num(r,'permutation_P'),num(r,'nominal_P'),r['gene_ID'],r['transcript_ID']));c=cov[vid]
 primary.append(dict(zip(fields,[vid,a['GNRS_ID'],a['chrom'],a['anchor_0based'],seq[vid],a['effect_allele'],a['insertion_present_allele'],a['insertion_absent_allele'],a['beta_insertion_presence'],a['gene_ID'],a['transcript_ID'],a['population'],c['distance_to_TSS'],c['length'],c['repeat_fraction'],c['repeat_class'],c['SV_frequency'],c['AF_definition'],a['orientation_proof'],'UNCHANGED_PHASE1R_FULL_DIRECTION_ELIGIBLE',a['beta_P'],a['permutation_P'],a['nominal_P'],len(rr)])))
write(R/'03_gnrs_truth/GNRS_CLASSB_PRIMARY33.tsv',primary);write(R/'03_gnrs_truth/GNRS_CLASSB_ALL_ELIGIBLE_TRANSCRIPTS.tsv',g)
(R/'03_gnrs_truth/GNRS_PREPARATION_AUDIT.json').write_text(json.dumps({'frozen_source_loci':33,'primary_loci':len(primary),'eligible_transcript_records':len(g),'aggregation':'min beta_P,permutation_P,nominal_P,gene_ID,transcript_ID;one gene association per locus','primary_truth_sha256':hashlib.sha256((R/'03_gnrs_truth/GNRS_CLASSB_PRIMARY33.tsv').read_bytes()).hexdigest(),'source_sha256':hashlib.sha256((OLD/'05_truth_gnrs/GNRS_DIRECTION_HARMONIZED.tsv').read_bytes()).hexdigest(),'prediction_informed_selection':False},indent=2))
# Cohort overlap uses exact author VCF header samples;no dosage parsing required.
with gzip.open(OLD/'05_truth_gnrs/GEUVADIS_chr1.vcf.gz','rt') as f:
 for line in f:
  if line.startswith('#CHROM'):gnrs=set(line.rstrip().split('\t')[9:]);break
mage_rows=read(R/'04_external_qtl_truth/MAGE_updated_metadata.tsv');mage={r['kgp_id'] for r in mage_rows};overlap=sorted(gnrs&mage)
write(R/'04_external_qtl_truth/MAGE_GEUVADIS_DONOR_OVERLAP.tsv',[{'donor':x,'membership':'BOTH_GNRS_GEUVADIS_AND_MAGE_SOURCE_METADATA'} for x in overlap],['donor','membership'])
(R/'04_external_qtl_truth/COHORT_INDEPENDENCE_AUDIT.json').write_text(json.dumps({'GNRS_author_analysis_VCF_sample_N':len(gnrs),'MAGE_author_updated_metadata_unique_donor_N':len(mage),'MAGE_metadata_rows':len(mage_rows),'donor_overlap_N':len(overlap),'overlap_known':True,'independent_cohort_status':'NOT_DISJOINT' if overlap else 'DISJOINT_BY_RETRIEVED_METADATA','SV_QTL_actual_analysis_subset_overlap':'UNVERIFIED_without_exact_analysis_sample_list','MAGE_sample_swap_notice':'2026-03-27 notice;2026-04-03 updated;HG00237/NA11919;Zenodo not updated per official README','independent_replication_requirement':'Need nonoverlapping donor summary or recomputed subset;different sequencing libraries alone do not establish independent people'},indent=2))
# Technical probes are committed independently of any inference.
used=set();probes=[]
conditions=[('short_nonrepeat',lambda c:50<=int(c['length'])<=199 and c['repeat_class']=='NO_DETECTED_REPEAT' and float(c['repeat_fraction'])==0),('Alu_SINE',lambda c:c['repeat_class'].startswith('SINE/Alu')),('LINE_LTR',lambda c:c['repeat_class'].split('/')[0] in ['LINE','LTR']),('insertion_1_5kb',lambda c:1000<=int(c['length'])<=5000),('difficult_repeat_rich',lambda c:int(c['length'])>=2000 and float(c['repeat_fraction'])>=.8)]
for cat,test in conditions:
 eligible=[c for c in cov.values() if c['OFFREF_ID'] not in used and test(c) and int(c['anchor'])>=524288]
 eligible.sort(key=lambda c:hashlib.sha256(('20261010|'+cat+'|'+c['OFFREF_ID']).encode()).hexdigest());assert eligible
 c=eligible[0];used.add(c['OFFREF_ID']);probes.append({'category':cat,**c,'inserted_sequence':seq[c['OFFREF_ID']],'role':'TECHNICAL_ONLY_NOT_STATISTICAL','prediction_status':'NOT_RUN_ALPHA_ACCESS_HOLD'})
write(R/'01_alphagenome/TECHNICAL_PROBE5_SELECTED.tsv',probes)
# Only existing BRFG candidate evidence is reused.
b=json.loads((OLD/'04_truth_brfg/BRFG_CARRIER_AUDIT.json').read_text());b.update(role='SECONDARY_REPRESENTATION_VALIDATION',truth='NOT_CAUSAL_EXPERIMENTAL_POSITIVES',source_file='../16_phase1r_rescue/04_truth_brfg/BRFG_GNRS_CARRIER_BRIDGE.tsv');(R/'05_brfg_secondary/SECONDARY_STATUS.json').write_text(json.dumps(b,indent=2))
print('GNRS primary',len(primary),'transcript rows',len(g),'donor overlap',len(overlap),'probes',len(probes),flush=True)
