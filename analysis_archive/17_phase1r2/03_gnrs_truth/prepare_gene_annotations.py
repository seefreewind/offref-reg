from pathlib import Path
import csv,gzip,collections,json
R=Path(__file__).resolve().parents[1];O=R/'03_gnrs_truth';truth=list(csv.DictReader((R/'03_gnrs_truth/GNRS_CLASSB_PRIMARY33.tsv').open(),delimiter='\t'));wanted={r['gene'].split('.')[0] for r in truth};trans={}
with gzip.open(O/'wgEncodeGencodeAttrsV20.txt.gz','rt') as f:
 for l in f:
  t=l.rstrip().split('\t')
  if len(t)>=5 and t[0].split('.')[0] in wanted:trans[t[4]]=t[0].split('.')[0]
intervals=collections.defaultdict(list)
with gzip.open(O/'wgEncodeGencodeCompV20.txt.gz','rt') as f:
 for l in f:
  t=l.rstrip().split('\t')
  if t[1] not in trans:continue
  for a,b in zip(t[9].strip(',').split(','),t[10].strip(',').split(',')):intervals[(trans[t[1]],t[2],t[3])].append((int(a),int(b)))
rows=[];audit=[]
for r in truth:
 matches=[k for k in intervals if k[0]==r['gene'].split('.')[0] and k[1]==r['GRCh38_chrom']];status='PASS' if len(matches)==1 else 'NO_NATIVE_GRCh38_GENE_MODEL' if not matches else 'AMBIGUOUS_GENE_STRAND_OR_MAPPING';n=0
 if len(matches)==1:
  k=matches[0];z=[]
  for a,b in sorted(intervals[k]):
   if z and a<=z[-1][1]:z[-1][1]=max(z[-1][1],b)
   else:z.append([a,b])
  for a,b in z:rows.append([r['OFFREF_ID'],r['gene'],k[1],k[2],a,b,'UCSC_hg38_GENCODEv20_native_exon_union;stable_ENSG_crosswalk_from_source_v19_identity']);n+=1
 audit.append([r['OFFREF_ID'],r['gene'],status,n,r['beta_P'],r['nominal_P']])
for name,data,fields in [('PRIMARY_GENE_EXONS_GRCh38.tsv',rows,['OFFREF_ID','gene','chrom','strand','start0','end0','source']),('GNRS_GENE_ANNOTATION_STATUS.tsv',audit,['OFFREF_ID','gene','annotation_status','union_exons_N','beta_P','nominal_P'])]:
 with (R/'03_gnrs_truth'/name).open('w') as f:w=csv.writer(f,delimiter='\t');w.writerow(fields);w.writerows(data)
(R/'03_gnrs_truth/GENE_ANNOTATION_STATUS.json').write_text(json.dumps({'coordinate_system':'native GRCh38','annotation':'GENCODEv20_UCSC_hg38','source_gene_ID_identity':'GENCODEv19 exact source preserved;version-independent stable ENSG crosswalk only for native gene exons','PASS_loci':sum(x[2]=='PASS' for x in audit),'UNVERIFIED_loci':sum(x[2]!='PASS' for x in audit),'gene_direction_truth_unchanged':33,'source_nominal_P_lt_005':sum(float(r['nominal_P'])<.05 for r in truth),'source_betaP_lt_005':sum(float(r['beta_P'])<.05 for r in truth),'gene_scoring_eligibility':'Distinct from truth-direction eligibility;unmapped/noncis genes remain unscored until resolved'},indent=2));print('Gene annotation',len(rows),'union exons',collections.Counter(x[2] for x in audit),'source nominal P<.05',sum(float(r['nominal_P'])<.05 for r in truth),'source betaP<.05',sum(float(r['beta_P'])<.05 for r in truth),flush=True)
