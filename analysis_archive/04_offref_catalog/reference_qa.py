"""Exhaustive reference validation; no sampling, outcome selection, or model calls."""
import pathlib,csv,json,collections,sys,py2bit,edlib,subprocess,hashlib
csv.field_size_limit(1000000000)  # Retain excluded megabase alleles in MASTER.
ROOT=pathlib.Path(__file__).resolve().parents[1];SEED=20261009
def prepare():
 refs={k:py2bit.open(str(ROOT/f'03_reference/{k}.2bit')) for k in ['hg38']}
 for k,p in refs.items():
  with (ROOT/f'03_reference/{k}.fa').open('w') as f:
   for ch,n in p.chroms().items():
    f.write('>'+ch+'\n')
    for s in range(0,n,1000000):
     seq=p.sequence(ch,s,min(s+1000000,n));f.write('\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')
 valid=0;failed=collections.Counter()
 with (ROOT/'04_offref_catalog/flanks_2kb.fa').open('w') as out,(ROOT/'04_offref_catalog/reference_precheck.tsv').open('w') as audit:
  w=csv.writer(audit,delimiter='\t');w.writerow(['OFFREF_ID','status'])
  for r in csv.DictReader((ROOT/'04_offref_catalog/OFFREF_PENDING.tsv').open(),delimiter='\t'):
   ch=r['chrom'];b=int(r['anchor_left']);p=refs['hg38'];reason='PASS'
   if b<2000 or b+2000>p.chroms(ch):reason='FLANK_OUTSIDE_CONTIG'
   else:
    ref=r['GRCh38_ref'];start=b-len(ref)
    if p.sequence(ch,start,b).upper()!=ref:reason='VCF_REF_MISMATCH'
    else:
     left=p.sequence(ch,b-2000,b).upper();right=p.sequence(ch,b,b+2000).upper()
     if 'N' in left+right:reason='GAP_IN_FLANK'
     else:
      for suffix,seq in [('L',left),('R',right)]:out.write('>'+r['OFFREF_ID']+'_'+suffix+'\n'+seq+'\n')
      valid+=1
   w.writerow([r['OFFREF_ID'],reason]);failed[reason]+=1
 (ROOT/'04_offref_catalog/reference_precheck_counts.json').write_text(json.dumps(dict(failed),indent=2));print(valid,flush=True)
def prepare_hs1():
 p=py2bit.open(str(ROOT/'03_reference/hs1.2bit'))
 with (ROOT/'03_reference/hs1.fa').open('w') as f:
  for ch,n in p.chroms().items():
   f.write('>'+ch+'\n')
   for s in range(0,n,1000000):
    seq=p.sequence(ch,s,min(s+1000000,n));f.write('\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')
def readpaf(k):
 hits=collections.defaultdict(list)
 for line in (ROOT/f'04_offref_catalog/flanks_{k}.paf').open():
  t=line.rstrip().split('\t');tags=dict(x.split(':',2)[::2] for x in t[12:] if ':' in x)
  hit={'qs':int(t[2]),'qe':int(t[3]),'strand':t[4],'ch':t[5],'s':int(t[7]),'e':int(t[8]),'identity':int(t[9])/int(t[10]),'coverage':(int(t[3])-int(t[2]))/int(t[1]),'score':int(tags.get('AS',t[9]))}
  hits[t[0]].append(hit)
 for name,hs in hits.items():
  hs.sort(key=lambda h:h['score'],reverse=True)
 return hits
def unique(hs):
 if not hs:return None
 h=hs[0]
 if h['identity']<.95 or h['coverage']<.90:return None
 if len(hs)>1 and hs[1]['score']>=.9*h['score']:return None
 return h
def finish():
 g=readpaf('hg38');c=readpaf('hs1');gap_coords={};gapseq={}
 # Native UCSC extraction avoids the observed py2bit long-range memmove stall.
 for name in c:
  if not name.endswith('_L'):continue
  vid=name[:-2];l=unique(c.get(name,[]));rr=unique(c.get(vid+'_R',[]))
  if l and rr and l['ch']==rr['ch'] and l['strand']==rr['strand'] and l['qe']==2000 and rr['qs']==0:
   a,z=(l['e'],rr['s']) if l['strand']=='+' else (rr['e'],l['s'])
   if 0<=z-a<=200000:gap_coords[vid]=(l['ch'],a,z,l['strand'])
 bed=ROOT/'04_offref_catalog/CHM13_syntenic_gaps.bed';fa=ROOT/'04_offref_catalog/CHM13_syntenic_gaps.fa'
 with bed.open('w') as f:
  for vid,(ch,a,z,strand) in gap_coords.items():
   if z>a:f.write(f'{ch}\t{a}\t{z}\t{vid}\n')
 if bed.stat().st_size:
  with (ROOT/'14_logs/CHM13_gap_extraction.log').open('w') as log:subprocess.run(['/Users/zy/.codex/tools/bin/twoBitToFa',str(ROOT/'03_reference/hs1.2bit'),str(fa),'-noMask','-bed='+str(bed)],check=True,stderr=log)
  key=None
  for line in fa.open():
   if line.startswith('>'):key=line[1:].split()[0];gapseq[key]=''
   else:gapseq[key]+=line.strip().upper()
 pre=dict(csv.reader((ROOT/'04_offref_catalog/reference_precheck.tsv').open(),delimiter='\t'));counts=collections.Counter();reasons=collections.Counter()
 inp=csv.DictReader((ROOT/'04_offref_catalog/OFFREF_MASTER.tsv').open(),delimiter='\t');fields=inp.fieldnames+['CHM13_chrom','CHM13_gap_length','classification_evidence']
 handles={};writers={};seen=set()
 for lab in ['MASTER','INCLUDED','EXCLUDED','PENDING']:
  handles[lab]=(ROOT/f'04_offref_catalog/OFFREF_{lab}.verified.tsv').open('w');writers[lab]=csv.DictWriter(handles[lab],fieldnames=fields,delimiter='\t');writers[lab].writeheader()
 for r in inp:
  r.update(CHM13_chrom='UNVERIFIED',CHM13_gap_length='UNVERIFIED',classification_evidence='NA');lab='EXCLUDED' if r['QC_status']=='EXCLUDED' else 'PENDING';vid=r['OFFREF_ID'];reason='NA'
  if lab!='EXCLUDED':
   key=(r['chrom'],r['anchor_left'],hashlib.sha256(r['inserted_sequence'].encode()).hexdigest())
   if key in seen:reason='DUPLICATE_NORMALIZED_ALLELE';lab='EXCLUDED'
   elif pre.get(vid)!='PASS':reason=pre.get(vid,'MISSING_REFERENCE_PRECHECK');lab='EXCLUDED'
   else:
    seen.add(key)
    left=unique(g.get(vid+'_L',[]));right=unique(g.get(vid+'_R',[]));b=int(r['anchor_left'])
    expected=lambda h,s,e:h and h['ch']==r['chrom'] and h['strand']=='+' and h['s']==s and h['e']==e and h['qs']==0 and h['qe']==2000
    if not expected(left,b-2000,b) or not expected(right,b,b+2000):reason='AMBIGUOUS_OR_NONEXACT_GRCH38_ANCHOR';lab='EXCLUDED';r['unique_anchor']='FALSE'
    else:
     r.update(unique_anchor='TRUE',GRCh38_presence='ABSENT_AT_SYNTENIC_LOCUS',QC_status='HIGH_CONFIDENCE_PAIRED',exclusion_reason='NA');lab='INCLUDED';counts['high_confidence_paired']+=1
     l=unique(c.get(vid+'_L',[]));rr=unique(c.get(vid+'_R',[]))
     if l and rr and l['ch']==rr['ch'] and l['strand']==rr['strand']:
      # Inner edges must be aligned, never extrapolated from a partial query.
      boundary_ok=l['qe']==2000 and rr['qs']==0
      if l['strand']=='+':a,z=l['e'],rr['s']
      else:a,z=rr['e'],l['s']
      if boundary_ok and 0<=z-a<=200000:
       gap='' if a==z else gapseq.get(vid)
       if gap is None:raise RuntimeError('Missing native CHM13 gap sequence: '+vid)
       assert len(gap)==z-a
       gap=gap if l['strand']=='+' else gap.translate(str.maketrans('ACGT','TGCA'))[::-1]
       r.update(CHM13_chrom=l['ch'],CHM13_gap_length=len(gap))
       if 'N' not in gap and len(gap)==0:
        r['CHM13_presence']='ABSENT_CALLABLE_EMPTY_SYNTENIC_ALLELE';r['classification_evidence']='Bilateral_unique_2kb_flanks;exact_empty_gap'
        if 0<float(r['allele_frequency'])<1:r['class_A_B_C']='CLASS_B';counts['CLASS_B']+=1
       elif 'N' not in gap:
        ins=r['inserted_sequence'];al=edlib.align(ins,gap,mode='HW',task='locations',k=int(.05*len(ins)))
        if al['editDistance']>=0:
         r.update(CHM13_presence='PRESENT_SYNTENIC',class_A_B_C='CLASS_A',classification_evidence='Bilateral_unique_2kb_flanks;full_insertion_alignment_identity_ge95');counts['CLASS_A']+=1
  if lab=='EXCLUDED' and reason!='NA':r.update(QC_status='EXCLUDED',exclusion_reason=reason)
  if lab=='PENDING':reason=r['QC_status']
  writers['MASTER'].writerow(r);writers[lab].writerow(r);counts[lab]+=1;reasons[reason]+=1
 for f in handles.values():f.close()
 for lab in handles:(ROOT/f'04_offref_catalog/OFFREF_{lab}.verified.tsv').replace(ROOT/f'04_offref_catalog/OFFREF_{lab}.tsv')
 (ROOT/'04_offref_catalog/catalog_verified_counts.json').write_text(json.dumps({'counts':dict(counts),'reasons':dict(reasons),'seed':SEED,'C_status':'NOT_EVALUATED_NO_VERIFIED_POPULATION_METADATA'},indent=2));print(json.dumps(dict(counts)),flush=True)
if __name__=='__main__':globals()[sys.argv[1]]()
