from pathlib import Path
import gzip,json,hashlib,collections
import pandas as pd
import py2bit
R=Path(__file__).resolve().parents[1];P=R.parent
# All source gene unions, exact stable IDs, native GRCh38. No model results used.
tx={};ints=collections.defaultdict(list)
with gzip.open(P/'17_phase1r2/03_gnrs_truth/wgEncodeGencodeAttrsV20.txt.gz','rt') as f:
 for l in f:
  t=l.rstrip().split('\t');tx[t[4]]=t[0].split('.')[0]
with gzip.open(P/'17_phase1r2/03_gnrs_truth/wgEncodeGencodeCompV20.txt.gz','rt') as f:
 for l in f:
  t=l.rstrip().split('\t')
  if t[1] in tx:
   for a,b in zip(t[9].strip(',').split(','),t[10].strip(',').split(',')):ints[(tx[t[1]],t[2],t[3])].append((int(a),int(b)))
union={}
for k,v in ints.items():
 z=[]
 for a,b in sorted(v):
  if z and a<=z[-1][1]:z[-1][1]=max(b,z[-1][1])
  else:z.append([a,b])
 union[k]=z
raw=pd.read_csv(R/'07_reference_benchmark/QTD000110.credible_sets.tsv.gz',sep='\t');raw=raw.drop_duplicates(['variant','gene_id']);genome=py2bit.open(str(P/'03_reference/hg38.2bit'));eligible=[];fails=collections.Counter()
for row in raw[raw.pip>=.5].to_dict('records'):
 chrom,pos,ref,alt=row['variant'].split('_');pos=int(pos);anchor=pos-1
 if len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT':fails['NOT_SNV']+=1;continue
 matches=[k for k in union if k[0]==str(row['gene_id']).split('.')[0] and k[1]==chrom]
 if len(matches)!=1:fails['NO_EXACT_GENE_ANNOTATION']+=1;continue
 if anchor<524288 or anchor+524288>genome.chroms(chrom):fails['REFERENCE_BOUNDARY']+=1;continue
 if genome.sequence(chrom,anchor,anchor+1).upper()!=ref:fails['REF_MISMATCH']+=1;continue
 k=matches[0];exons=union[k];bp=sum(max(0,min(b,anchor+98304)-max(a,anchor-98304)) for a,b in exons)
 if bp==0:fails['NO_COMMON_EXON']+=1;continue
 row.update(chrom=chrom,pos1=pos,anchor0=anchor,REF=ref,ALT=alt,gene=k[0],strand=k[2],exons=exons,common_exon_bp=bp,selection_hash=hashlib.sha256(('20261011|reference|'+row['variant']+'|'+k[0]).encode()).hexdigest());eligible.append(row)
def pick(threshold):
 selected=[];counts=collections.Counter();used=set()
 for row in sorted([x for x in eligible if x['pip']>=threshold],key=lambda x:x['selection_hash']):
  if counts[row['gene']]>=5 or row['variant'] in used:continue
  selected.append(row);counts[row['gene']]+=1;used.add(row['variant'])
  if len(selected)==500:break
 return selected
selected=pick(.9);threshold=.9
if len(selected)<500:selected=pick(.5);threshold=.5
assert len(selected)>=20 and len({z['gene'] for z in selected})>=20
jobs=[]
for row in selected:
 j={'job_ID':'REF_'+row['variant']+'_'+row['gene'],'role':'REFERENCE_EQTL','panel':'REFERENCE_EQTL','cell':'GM12878','gene':row['gene'],'strand':row['strand'],'chrom':row['chrom'],'anchor0':row['anchor0'],'insertion':'','allele_type':'SNV','pos1':row['pos1'],'REF':row['REF'],'ALT':row['ALT'],'exons':row['exons'],'beta':float(row['beta']),'PIP':float(row['pip']),'source_variant':row['variant'],'source':'eQTL Catalogue release7 GEUVADIS QTD000110','effect_allele':'ALT'};jobs.append(j)
p=R/'07_reference_benchmark/REFERENCE_PANEL_FROZEN.json';p.write_text(json.dumps(jobs,indent=2));pd.DataFrame([{k:v for k,v in z.items() if k!='exons'} for z in jobs]).to_csv(R/'07_reference_benchmark/REFERENCE_PANEL_FROZEN.tsv',sep='\t',index=False)
audit={'source_url':'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000013/QTD000110/QTD000110.credible_sets.tsv.gz','source_sha256':hashlib.sha256((R/'07_reference_benchmark/QTD000110.credible_sets.tsv.gz').read_bytes()).hexdigest(),'panel_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'N':len(jobs),'genes_N':len({z['gene'] for z in jobs}),'PIP_threshold':threshold,'source_exclusions':dict(fails),'sampling_uses_model_predictions':False,'independence':'conventional benchmark unrelated to OFFREF prediction selection;GEUVADIS donors overlap GNRS, not independent biological replication','bucket_fallback_reason':'official Borzoi bucket returned UserProjectMissing/requester-pays; public EBI fallback','allele_orientation':'ALT always effect allele per official eQTL Catalogue Data_access FAQ','before_model_inference':True}
(R/'07_reference_benchmark/REFERENCE_PANEL_FREEZE_AUDIT.json').write_text(json.dumps(audit,indent=2));print('REFERENCE_PANEL_FROZEN',len(jobs),audit['genes_N'],threshold,flush=True)
