"""50 seeded Phase0 loci; native UCSC vs generated query flank crosscheck."""
import csv,pathlib,random,subprocess,hashlib,json
R=pathlib.Path(__file__).resolve().parents[1];SEED=20261009
def fasta(p,wanted=None):
 d={};key=None
 for l in p.open():
  if l.startswith('>'):
   k=l[1:].split()[0];key=k if wanted is None or k in wanted else None
   if key:d[key]=''
  elif key:d[key]+=l.strip().upper()
 return d
def main():
 pre=dict(csv.reader((R/'04_offref_catalog/reference_precheck.tsv').open(),delimiter='\t'))
 eligible=[r for r in csv.DictReader((R/'04_offref_catalog/OFFREF_PENDING.tsv').open(),delimiter='\t') if pre.get(r['OFFREF_ID'])=='PASS'];rows=random.Random(SEED).sample(eligible,50)
 bed=R/'04_offref_catalog/NATIVE_QA50.bed';fa=R/'04_offref_catalog/NATIVE_QA50.fa';wanted=set()
 with bed.open('w') as f:
  for r in rows:
   b=int(r['anchor_left'])
   for suffix,a,z in [('L',b-2000,b),('R',b,b+2000)]:
    key=r['OFFREF_ID']+'_'+suffix;wanted.add(key);f.write(f'{r["chrom"]}\t{a}\t{z}\t{key}\n')
 with (R/'14_logs/native_reference_QA50.log').open('w') as log:subprocess.run(['/Users/zy/.codex/tools/bin/twoBitToFa',str(R/'03_reference/hg38.2bit'),str(fa),'-noMask','-bed='+str(bed)],check=True,stderr=log)
 a=fasta(fa);b=fasta(R/'04_offref_catalog/flanks_2kb.fa',wanted);out=[]
 for k in sorted(wanted):out.append([k,len(a.get(k,'')),len(b.get(k,'')),hashlib.sha256(a.get(k,'').encode()).hexdigest(),hashlib.sha256(b.get(k,'').encode()).hexdigest(),a.get(k)==b.get(k)])
 with (R/'04_offref_catalog/REFERENCE_NATIVE_CROSSCHECK_50.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['flank_ID','native_length','query_length','native_sha256','query_sha256','identical']);w.writerows(out)
 assert len(a)==len(b)==100 and all(x[-1] for x in out)
 print('PASS:50loci/100flanks identical to nativeUCSC;Phase0 only',flush=True)
if __name__=='__main__':main()
