"""Prospective exact Class B matching. Run only on verified author summary+VCF inputs."""
from pathlib import Path
import csv,gzip,collections,json,argparse
R=Path(__file__).resolve().parents[1]
def normalize_ins(pos1,ref,alt):
 ref,alt=ref.upper(),alt.upper();prefix=0
 while prefix<min(len(ref),len(alt)) and ref[prefix]==alt[prefix]:prefix+=1
 rr,aa=ref[prefix:],alt[prefix:]
 while rr and aa and rr[-1]==aa[-1]:rr=rr[:-1];aa=aa[:-1]
 if rr or not aa or set(aa)-set('ACGT'):return None
 return int(pos1)-1+prefix,aa

def main():
 p=argparse.ArgumentParser();p.add_argument('--vcf',required=True);p.add_argument('--summary',required=True);p.add_argument('--orientation-proof',required=True);a=p.parse_args();proof=json.loads(Path(a.orientation_proof).read_text());assert proof['official_source_verified'] and proof['beta_dosage']=='ALT1', 'Do not guess ALT effect orientation'
 cov=list(csv.DictReader((R.parent/'16_phase1r_rescue/06_covariates/OFFREF_COVARIATE_COMPLETE.tsv').open(),delimiter='\t'));seq={};v=None
 for line in (R.parent/'16_phase1r_rescue/06_covariates/CLASS_B_INSERTIONS.fa').open():
  if line.startswith('>'):v=line[1:].strip();seq[v]=''
  else:seq[v]+=line.strip()
 exact=collections.defaultdict(list);anchors=collections.defaultdict(list)
 for c in cov:
  key=(c['chrom'],int(c['anchor']),seq[c['OFFREF_ID']]);exact[key].append(c['OFFREF_ID']);anchors[key[:2]].append((c['OFFREF_ID'],seq[c['OFFREF_ID']]))
 match={};op=gzip.open if str(a.vcf).endswith('.gz') else open
 with op(a.vcf,'rt') as f:
  for line in f:
   if line.startswith('#'):continue
   t=line.rstrip().split('\t');ch=t[0] if t[0].startswith('chr') else 'chr'+t[0];vid=t[2]
   if ',' in t[4] or '<' in t[4]:match[vid]=('AMBIGUOUS_MATCH',[],'UNRESOLVED_SYMBOLIC_OR_MULTIALLELIC');continue
   z=normalize_ins(t[1],t[3],t[4])
   if z is None:match[vid]=('NO_MATCH',[],'NOT_SEQUENCE_RESOLVED_PURE_INS');continue
   anchor,ins=z;hits=exact[(ch,anchor,ins)];key=(ch,anchor)
   if len(hits)==1:match[vid]=('EXACT_SEQUENCE_ANCHOR_MATCH',hits,'TRIMMED_VCF_INSERTION;ALT_presence;REF_absence')
   elif len(hits)>1:match[vid]=('AMBIGUOUS_MATCH',hits,'MULTIPLE_CATALOG_IDENTITIES')
   elif anchors[key]:match[vid]=('AMBIGUOUS_MATCH',[],'ANCHOR_PRESENT_BUT_SEQUENCE_DISAGREES')
   else:match[vid]=('NO_MATCH',[],'NO_EXACT_ANCHOR_SEQUENCE')
 rows=[]
 for r in csv.DictReader(Path(a.summary).open(),delimiter='\t'):
  state,hits,evidence=match.get(r['variant_ID'],('AMBIGUOUS_MATCH',[],'MISSING_VCF_RECORD'))
  for vid in hits if hits else ['NA']:
   rows.append({'OFFREF_ID':vid,'external_variant_ID':r['variant_ID'],'match_category':state,'gene':r['gene'],'effect_allele':'ALT1_verified_by_author_source','presence_allele':'trimmed_ALT_insertion','beta_presence':r['beta'],'nominal_P':r['P'],'SE':r.get('SE','NA'),'orientation_evidence':evidence+';'+proof['source'],'direction_eligible':str(state=='EXACT_SEQUENCE_ANCHOR_MATCH').upper(),'independent_cohort_eligible':str(proof.get('disjoint_donors_verified',False)).upper()})
 out=R/'04_external_qtl_truth/EXTERNAL_CLASSB_DIRECTION_TRUTH.tsv'
 with out.open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
if __name__=='__main__':main()
