from pathlib import Path
import json,hashlib,datetime
R=Path(__file__).resolve().parents[1];p=R/'00_protocol/run_predictions.py';s=p.read_text();old=hashlib.sha256(p.read_bytes()).hexdigest()
s=s.replace('from run_technical import key_from_document,Deadline','from run_technical import key_from_document,Deadline,Details')
s=s.replace("key=key_from_document();clients={}","""key=key_from_document();clients={}
  class GlobalDeadline(grpc.StreamStreamClientInterceptor):
   def intercept_stream_stream(self,continuation,d,iterator):
    return continuation(Details(d.method,600,d.metadata,d.credentials,d.wait_for_ready,d.compression),iterator)""")
s=s.replace("if job['member'] not in clients:","client_key=(job['member'],job['track_scope'])\n     if client_key not in clients:")
s=s.replace("c._channel=grpc.intercept_channel(c._channel,Deadline());clients[job['member']]=c", "c._channel=grpc.intercept_channel(c._channel,GlobalDeadline() if job['track_scope']=='GLOBAL' else Deadline());clients[client_key]=c")
s=s.replace("clients[job['member']]",'clients[client_key]')
a=s.index('      for batch in batches:');b=s.index('      assert len(observed)',a)
s=s[:a]+'''      cache_dir=outdir/'ontology_checkpoints'/job['audit_job_ID']/allele;cache_dir.mkdir(parents=True,exist_ok=True)
      queue=[list(batch) for batch in batches];completed_batches=0
      while queue:
       batch=queue.pop(0);batch_key=hashlib.sha256(json.dumps(batch).encode()).hexdigest();cp=cache_dir/(batch_key+'.json')
       if cp.exists():
        record=readjson(cp);assert record['sequence_sha256']==result['sequence_'+allele+'_sha256'] and record['intervals']==[list(x) for x in intervals];pm=record['metadata'];ps=record['summary']
       else:
        request_start=time.monotonic()
        try:
         part=clients[client_key].predict_sequence(s,organism=dna_client.Organism.HOMO_SAPIENS,requested_outputs=[dna_client.OutputType.RNA_SEQ],ontology_terms=batch,interval=aggenome.Interval(job['chrom'],job['start0'],job['start0']+job['W']))
        except Exception as e:
         if len(batch)>1:
          split=len(batch)//2;queue[:0]=[batch[:split],batch[split:]]
          (cache_dir/(batch_key+'.transport_failure.json')).write_text(json.dumps({'ontology_terms':batch,'seconds':time.monotonic()-request_start,'error_type':type(e).__name__,'error':str(e).replace(key,'[REDACTED]')[:1000],'action':'SPLIT_SAME_OFFICIAL_ONTOLOGIES_NO_TRACK_EXCLUSION'},indent=2));continue
         raise
        pv=part.rna_seq.values;pr=int(part.rna_seq.resolution);assert len(pv)*pr==job['W'];pm=[]
        for tr in part.rna_seq.metadata.to_dict('records'):
         matches=[x for x in md if x['track_name']==tr['name'] and x['strand']==tr['strand']];assert len(matches)==1;pm.append(matches[0])
        assert np.isfinite(pv).all();ps=summarize(pv,intervals,job['start0'],0,pr,allele)
        cp.write_text(json.dumps({'ontology_terms':batch,'sequence_sha256':result['sequence_'+allele+'_sha256'],'intervals':intervals,'metadata':pm,'summary':ps,'elapsed_seconds':time.monotonic()-request_start},indent=2));del pv,part
       means.extend(ps['track_means']);bps.append(ps['mapped_bp']);observed.extend(pm);completed_batches+=1
      batches_count=completed_batches
''' + s[b:]
s=s.replace("result['ontology_batch_N']=len(batches)","result['ontology_batch_N']=batches_count")
p.write_text(s)
lock=R/'15_freeze/PREDICTION_INPUT_SHA256.json';lock_before=lock.read_bytes();(R/'15_freeze/PREDICTION_INPUT_SHA256_before_global_checkpoints.json').write_bytes(lock_before);records=json.loads(lock_before)
for z in records:
 if z['path']=='00_protocol/run_predictions.py':z['sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
lock.write_text(json.dumps(records,indent=2));(R/'02_track_harmonization/GLOBAL_CHECKPOINT_TRANSPORT_REPAIR.json').write_text(json.dumps({'UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'old_runner_sha256':old,'new_runner_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'change':'600s global stream deadline; adaptive metadata-only ontology split on transport failure; per-allele per-ontology sufficient-statistic checkpoints','scoring_inputs_changed':False,'expected_RNA_tracks':667,'complete_aggregation_requires_all_official_tracks':True,'single_ontology_probe':'PASS2tracks62seconds','reason':'32-track batch previously exceeded180s; no missing-tissue imputation or partial-pool aggregate'},indent=2))
# Retain actual failed requests; retry one prospectively defined valid pilot identity.
jobs=json.loads((R/'00_protocol/PREDICTION_JOBS_FROZEN.json').read_text());jobs=sorted([j for j in jobs if j['model']=='AlphaGenome' and j['stage']=='GLOBAL_CORE'],key=lambda j:(str(j['member']),j['audit_job_ID']))
for j in jobs:
 q=R/'02_track_harmonization/AlphaGenome'/(j['audit_job_ID']+'.json')
 if q.exists():
  z=json.loads(q.read_text())
  if z['status']=='UNAVAILABLE' and z.get('error_type') not in ['ValueError']:
   hist=q.parent/'smoke_retry_history';hist.mkdir(exist_ok=True);q.rename(hist/(q.stem+'_before_checkpoints.json'))
print('GLOBAL_CHECKPOINT_REPAIR_RECORDED')
