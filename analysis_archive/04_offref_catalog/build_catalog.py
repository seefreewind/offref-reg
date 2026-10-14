#!/usr/bin/env python3
"""All HGSVC3 records retained. Provisional sequence QC is not gate qualification."""
import csv,gzip,json,re,collections,hashlib,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[1];SEED=20261009
FIELDS='OFFREF_ID source source_variant_id chrom anchor_left anchor_right GRCh38_ref inserted_sequence length SV_type sequence_resolved unique_anchor repeat_fraction repeat_class TE_family GC mappability GRCh38_presence CHM13_presence class_A_B_C allele_frequency population_AF number_haplotypes number_individuals nearest_gene distance_to_TSS QC_status exclusion_reason coordinate_system lead_sample observed_ALT_haplotypes called_haplotypes reference_SD reference_TRF annotation_support'.split()
def main():
 ann={}
 with gzip.open(ROOT/'01_sources/evidence/hgsvc3_annotation_probe.tsv.gz','rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):ann[row['ID']]=row
 writers={};handles=[]
 for label in ['MASTER','INCLUDED','EXCLUDED','PENDING']:
  f=(ROOT/f'04_offref_catalog/OFFREF_{label}.tsv').open('w');handles.append(f)
  writers[label]=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t');writers[label].writeheader()
 counts=collections.Counter();reasons=collections.Counter();samples=[];callable_by_sample=collections.Counter();missing_by_sample=collections.Counter()
 with gzip.open(ROOT/'04_offref_catalog/variants_GRCh38_sv_insdel_alt_HGSVC2024v1.0.vcf.gz','rt') as f:
  for line in f:
   if line.startswith('##'):continue
   if line.startswith('#CHROM'):samples=line.rstrip().split('\t')[9:];continue
   parts=line.rstrip('\n').split('\t');chrom,pos,vid,ref,alts,qual,filt,infos,fmt=parts[:9];info=dict((x.split('=',1)+[''])[:2] for x in infos.split(';'));gts=[x.split(':')[fmt.split(':').index('GT')] for x in parts[9:]]
   called=[]
   for sample,gt in zip(samples,gts):
    tokens=re.split(r'[|/]',gt)
    callable_by_sample[sample]+=sum(t!='.' for t in tokens);missing_by_sample[sample]+=sum(t=='.' for t in tokens)
    called.extend(t for t in tokens if t!='.')
   for ai,alt in enumerate(alts.split(','),1):
    counts['source_alleles']+=1;row={k:'UNVERIFIED' for k in FIELDS}
    row.update(OFFREF_ID='OFFREF_HGSVC3_'+hashlib.sha256(f'{chrom}:{pos}:{ref}:{alt}'.encode()).hexdigest()[:16],source='HGSVC3_2024v1.0_GRCh38',source_variant_id=vid,chrom=chrom,GRCh38_ref=ref,SV_type=info.get('SVTYPE','UNVERIFIED'),number_haplotypes=len(called),number_individuals=sum(any(t!='.' for t in re.split(r'[|/]',gt)) for gt in gts),coordinate_system='0-based_interbase_after_shared_prefix;left=right=boundary',lead_sample=info.get('SAMPLE','UNVERIFIED'),called_haplotypes=len(called),observed_ALT_haplotypes=sum(t==str(ai) for t in called),population_AF='UNVERIFIED',class_A_B_C='UNCLASSIFIED',CHM13_presence='UNVERIFIED',unique_anchor='UNVERIFIED')
    a=ann.get(vid,{})
    row.update(reference_SD=a.get('REF_SD','UNVERIFIED'),reference_TRF=a.get('REF_TRF','UNVERIFIED'),annotation_support=a.get('MERGE_SAMPLES','UNVERIFIED'))
    reasons_here=[];literal=bool(re.fullmatch('[ACGT]+',ref)) and bool(re.fullmatch('[ACGT]+',alt))
    row['sequence_resolved']='TRUE' if literal else 'FALSE'
    pref=0;suff=0
    if literal:
     while pref<min(len(ref),len(alt)) and ref[pref]==alt[pref]:pref+=1
     while suff<min(len(ref)-pref,len(alt)-pref) and ref[-1-suff]==alt[-1-suff]:suff+=1
    remaining_ref=ref[pref:len(ref)-suff if suff else len(ref)];inserted=alt[pref:len(alt)-suff if suff else len(alt)] if literal else ''
    boundary=int(pos)-1+pref;row.update(anchor_left=boundary,anchor_right=boundary,length=len(inserted),inserted_sequence=inserted or 'NA')
    pure=literal and not remaining_ref and len(inserted)>=50
    if pure:counts['candidate_pure_insertions']+=1;row.update(GC=(inserted.count('G')+inserted.count('C'))/len(inserted),GRCh38_presence='ABSENT_AT_CALLED_LOCUS_UNVERIFIED_FASTA')
    if not literal:reasons_here.append('UNRESOLVED_OR_NON_ACGT')
    elif not pure:reasons_here.append('NOT_PURE_INSERTION_GE50')
    if filt not in ['PASS','.']:reasons_here.append('SOURCE_FILTER_'+filt)
    if not re.fullmatch(r'chr([1-9]|1[0-9]|2[0-2]|X|Y)',chrom):reasons_here.append('NON_PRIMARY_ANCHOR')
    if not row['observed_ALT_haplotypes']:reasons_here.append('NO_CALLED_ALT_HAPLOTYPE')
    if literal and len(inserted)>100000:reasons_here.append('GT100KB_RECONSTRUCTION_REQUIRES_LONG_MODEL_QA')
    te=a.get('TE','')
    if te:
     row['TE_family']=te;row['repeat_class']='SINE/Alu' if 'Alu' in te else 'LINE' if re.search('L1|LINE',te) else 'LTR' if re.search('LTR|ERV',te) else 'other_repeat'
    row['allele_frequency']=sum(t==str(ai) for t in called)/len(called) if called else 'UNVERIFIED'
    if reasons_here:row.update(QC_status='EXCLUDED',exclusion_reason=';'.join(reasons_here));label='EXCLUDED';counts['excluded']+=1;reasons.update(reasons_here)
    else:row.update(QC_status='PENDING_REFERENCE_ANCHOR_CHM13_QA',exclusion_reason='NA');label='PENDING';counts['pending']+=1
    writers['MASTER'].writerow(row);writers[label].writerow(row)
 for f in handles:f.close()
 counts['high_confidence_paired']=0;counts['CLASS_A']=0;counts['CLASS_B_including_C']=0;counts['CLASS_C']=0;counts['vcf_samples']=len(samples)
 (ROOT/'04_offref_catalog/catalog_initial_counts.json').write_text(json.dumps({'counts':dict(counts),'reasons':dict(reasons),'seed':SEED},indent=2))
 with (ROOT/'04_offref_catalog/SAMPLE_CALLABILITY_AUDIT.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['sample','called_haplotype_genotypes_across_records','missing_haplotype_genotypes_across_records']);w.writerows((s,callable_by_sample[s],missing_by_sample[s]) for s in samples)
 print(json.dumps(dict(counts)))
if __name__=='__main__':main()
