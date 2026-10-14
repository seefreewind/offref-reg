from pathlib import Path
import gzip,csv,json,math,collections,hashlib,py2bit
D=Path(__file__).resolve().parent;P=D.parents[1];SEED=20261010
rows=list(csv.DictReader(gzip.open(D/'elements.tsv.gz','rt'),delimiter='\t'));r38=[r for r in rows if r['Release']=='GRCh38'];valid=[];excluded=collections.Counter()
for r in r38:
 try:
  if r['Ref'] not in 'ACGT' or r['Alt'] not in 'ACGT' or len(r['Ref'])!=1 or len(r['Alt'])!=1 or r['Ref']==r['Alt']:excluded['NOT_SNV']+=1;continue
  if int(r['Barcodes'])<10:excluded['FEWER_THAN_10_TAGS']+=1;continue
  if int(r['DNA'])<=0 or int(r['RNA'])<=0:excluded['NONPOSITIVE_COUNTS']+=1;continue
  assert math.isfinite(float(r['Coefficient'])) and 0<=float(r['pValue'])<=1
  valid.append(r)
 except (ValueError,AssertionError):excluded['INVALID_NUMERIC']+=1
with (D/'MPRA_GRCh38_QC_VARIANTS.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(valid)
# LDLR primary + independent full-library replicate LDLR.2, explicitly HepG2 in source.
by=collections.defaultdict(dict)
for r in valid:
 if r['Element'] in ['LDLR','LDLR.2']:by[(r['Chrom'],int(r['Pos']),r['Ref'],r['Alt'])][r['Element']]=r
ref=py2bit.open(str(P/'03_reference/hg38.2bit'));paired=[];pairs=[];pools=collections.defaultdict(lambda:collections.defaultdict(list))
for key,rr in sorted(by.items()):
 if set(rr)!= {'LDLR','LDLR.2'}:continue
 ch,pos,a,b=key;ch='chr'+ch if not ch.startswith('chr') else ch;assert ref.sequence(ch,pos-1,pos).upper()==a
 x,y=rr['LDLR'],rr['LDLR.2'];fx,fy=float(x['Coefficient']),float(y['Coefficient']);px,py=float(x['pValue']),float(y['pValue']);label='STRONG_EFFECT' if px<1e-5 and py<1e-5 and fx*fy>0 else 'WEAK_UNRESOLVED_CONTROL' if px>=.1 and py>=.1 else 'INTERMEDIATE_NOT_BINARY_CONTROL'
 r={'variant_ID':f'{ch}:{pos}:{a}:{b}','chrom':ch,'pos1':pos,'REF':a,'ALT':b,'element':'LDLR','cell_type':'HepG2','measured_delta_log2':(fx+fy)/2,'replicate1_delta':fx,'replicate2_delta':fy,'replicate1_P':px,'replicate2_P':py,'label':label,'replicate_N':2,'biological_unit':'one_LDLR_promoter;variants_not_independent_elements','coordinate_REF_check':'PASS','readout':'reporter_expression_not_endogenous_gene_expression'};paired.append(r)
 if label!='INTERMEDIATE_NOT_BINARY_CONTROL':pools[a+'>'+b][label].append(r)
for mutation,z in sorted(pools.items()):
 for k in z:z[k].sort(key=lambda r:hashlib.sha256((str(SEED)+'|'+r['variant_ID']).encode()).hexdigest())
 for i,(a,b) in enumerate(zip(z['STRONG_EFFECT'],z['WEAK_UNRESOLVED_CONTROL'])):pairs.append({'pair_ID':mutation+'|'+str(i),'positive_ID':a['variant_ID'],'control_ID':b['variant_ID'],'match':'same_element_cell_type_substitution_class','control_definition':'both_source_P>=0.1;weak/unresolved_not_proven_null'})
for name,r,fields in [('MPRA_LDLR_PAIRED_VARIANTS.tsv',paired,list(paired[0])),('MPRA_MATCHED_POSITIVE_CONTROL_PAIRS.tsv',pairs,['pair_ID','positive_ID','control_ID','match','control_definition'])]:
 with (D/name).open('w') as f:w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(r)
a={'source_rows':len(rows),'GRCh38_rows':len(r38),'QC_SNV_rows':len(valid),'excluded_counts':dict(excluded),'paired_LDLR_QC_SNVs':len(paired),'binary_labels':dict(collections.Counter(r['label'] for r in paired)),'matched_pairs':len(pairs),'seed':SEED,'primary_diagnostic_subset':'LDLR plus LDLR.2 full-library replicates;HepG2 documented in main article;all QC GRCh38 source variants retained separately','independent_elements':1,'cluster_CI':'NA:one promoter;no independent-element inference','status':'MEASURED_DATA_PREPARED_MODEL_RUN_NOT_STARTED','model_development_exposure':'UNVERIFIED;not claimed held-out model validation','off_reference_truth':False}
(D/'MPRA_PREPARATION_AUDIT.json').write_text(json.dumps(a,indent=2));print(json.dumps(a),flush=True)
