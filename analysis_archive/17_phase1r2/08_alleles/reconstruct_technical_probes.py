from pathlib import Path
import csv,json,hashlib,py2bit
R=Path(__file__).resolve().parents[1];P=R.parent
probes=list(csv.DictReader((R/'01_alphagenome/TECHNICAL_PROBE5_SELECTED.tsv').open(),delimiter='\t'));assert len(probes)==5
ref=py2bit.open(str(P/'03_reference/hg38.2bit'));maps=[];audit=[]
for r in probes:
 v=r['OFFREF_ID'];ch=r['chrom'];b=int(r['anchor']);ins=r['inserted_sequence'].upper();L=len(ins);assert L==int(r['length']) and set(ins)<=set('ACGT')
 for model,W,crop in [('AlphaGenome',1048576,0),('Borzoi',524288,163840)]:
  a=b-W//2;end=a+W;assert a>=0 and end<=ref.chroms(ch) and L<W//2
  backbone=ref.sequence(ch,a,end).upper();assert len(backbone)==W
  alt=backbone[:W//2]+ins+backbone[W//2:W-L];assert len(alt)==W and alt[:W//2]==backbone[:W//2] and alt[W//2:W//2+L]==ins and alt[W//2+L:]==backbone[W//2:W-L]
  for allele,s in [('REF',backbone),('ALT',alt)]:
   p=R/'08_alleles'/f'{v}.{model}.{allele}.fa';p.write_text('>'+v+'|'+model+'|'+allele+'\n'+s+'\n')
  blocks=[('upstream',a,b,0,W//2,0,W//2,'0','HOMOLOGOUS'),('insertion',b,b,W//2,W//2,W//2,W//2+L,'NA','ALT_ONLY_EXCLUDED'),('downstream',b,end-L,W//2,W-L,W//2+L,W,'+'+str(L),'HOMOLOGOUS'),('distal_compensation',end-L,end,W-L,W,W,W,'NA','REF_ONLY_EXCLUDED')]
  for label,g0,g1,r0,r1,a0,a1,offset,state in blocks:maps.append([v,model,ch,g0,g1,r0,r1,a0,a1,label,b,L,offset,state,crop,W-crop,'ANY_BIN_INTERSECTING_ALT_INSERTION_EXCLUDED;downstream_integrate_mapped_overlap_not_equal_index'])
  audit.append({'OFFREF_ID':v,'model':model,'window_bp':W,'insertion_bp':L,'isogenic_validation':'PASS','full_insertion_preserved':True,'reference_genomic_start0':a,'reference_genomic_end0':end,'REF_sequence_sha256':hashlib.sha256(backbone.encode()).hexdigest(),'ALT_sequence_sha256':hashlib.sha256(alt.encode()).hexdigest(),'REF_N_fraction':backbone.count('N')/W,'ALT_N_fraction':alt.count('N')/W,'orientation':'exact forward-genomic frozen inserted sequence;no reverse complement','prediction_status':'NOT_RUN_ALPHA_ACCESS_HOLD','role':'TECHNICAL_ONLY'})
with (R/'08_alleles/TECHNICAL_ALLELE_HOMOLOGY_MAP.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['OFFREF_ID','model','chrom','REF_genomic_start0','REF_genomic_end0','REF_input_start0','REF_input_end0','ALT_input_start0','ALT_input_end0','homology_block','insertion_boundary0','insertion_length','downstream_offset','status','native_output_crop_start0','native_output_crop_end0','bin_rule']);w.writerows(maps)
(R/'08_alleles/TECHNICAL_ALLELE_QA.json').write_text(json.dumps({'role':'TECHNICAL_ONLY_NO_SENTINEL','loci':5,'paired_model_contexts':len(audit),'validation':audit},indent=2));print('Validated 5 technical loci,10 paired native contexts;40 homology blocks',flush=True)
