"""Official donor metadata only; never infer ancestry from sample names."""
import csv,pathlib,gzip,re,json,hashlib,collections
R=pathlib.Path(__file__).resolve().parents[1];SEED=20261009
csv.field_size_limit(1000000000)
def main():
 metadata={r['sample']:r['super_population'] for r in csv.DictReader((R/'04_offref_catalog/SAMPLE_POPULATION_METADATA.tsv').open(),delimiter='\t') if r['super_population']!='UNVERIFIED'}
 d={};samples=[]
 with gzip.open(R/'04_offref_catalog/variants_GRCh38_sv_insdel_alt_HGSVC2024v1.0.vcf.gz','rt') as f:
  for line in f:
   if line.startswith('#CHROM'):samples=line.rstrip().split('\t')[9:];continue
   if line.startswith('#'):continue
   t=line.rstrip().split('\t')
   if 'SVTYPE=INS;' not in t[7] and not t[7].endswith('SVTYPE=INS'):continue
   gt_index=t[8].split(':').index('GT');gts=[re.split('[|/]',v.split(':')[gt_index]) for v in t[9:]]
   for ai,alt in enumerate(t[4].split(','),1):
    p={}
    for sample,gt in zip(samples,gts):
     if sample not in metadata:continue
     group=metadata[sample];x=p.setdefault(group,{'ALT_copies':0,'called_haplotypes':0,'called_individuals':0})
     x['ALT_copies']+=sum(g==str(ai) for g in gt);x['called_haplotypes']+=sum(g!='.' for g in gt);x['called_individuals']+=int(any(g!='.' for g in gt))
    for x in p.values():x['AF']=x['ALT_copies']/x['called_haplotypes'] if x['called_haplotypes'] else None
    vid='OFFREF_HGSVC3_'+hashlib.sha256(f'{t[0]}:{t[1]}:{t[3]}:{alt}'.encode()).hexdigest()[:16];d[vid]=json.dumps(p,separators=(',',':'),sort_keys=True)
 popcounts=collections.Counter(metadata.values());status={'seed':SEED,'verified_donors':len(metadata),'total_donors':len(samples),'unmatched_donors':len(samples)-len(metadata),'verified_superpopulation_counts':dict(popcounts),'CLASS_C_status':'NOT_EVALUABLE_NO_TWO_GROUPS_WITH_AT_LEAST20_VERIFIED_DONORS','AF_interpretation':'Discovery-cohort allele frequency;related donors possible;not population prevalence'}
 assert len([g for g,n in popcounts.items() if n>=20])<2
 for label in ['MASTER','INCLUDED','EXCLUDED','PENDING']:
  path=R/f'04_offref_catalog/OFFREF_{label}.tsv';inp=csv.DictReader(path.open(),delimiter='\t');fields=inp.fieldnames
  temp=R/f'04_offref_catalog/OFFREF_{label}.population.tmp.tsv'
  with temp.open('w') as f:
   w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
   for r in inp:
    if r['OFFREF_ID'] in d:r['population_AF']=d[r['OFFREF_ID']]
    w.writerow(r)
  temp.replace(path)
 (R/'04_offref_catalog/POPULATION_METADATA_AUDIT.json').write_text(json.dumps(status,indent=2));print(json.dumps(status),flush=True)
if __name__=='__main__':main()
