from pathlib import Path
import sys,json,time,traceback,argparse,hashlib,re,datetime,gc
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];sys.path[:0]=[str(R/'00_protocol'),str(R/'01_model_metadata'),str(R.parent/'19_phase1r2_execution/01_alphagenome')]
from audit_common import *
import py2bit
STAGE_FOLDER={'GLOBAL_CORE':'02_track_harmonization','ENSEMBLE':'03_within_model_ensembles','PLACEMENT':'04_context_robustness','CONTEXT':'04_context_robustness','GNRS_GENE_CENTER':'08_truth_coverage','REFERENCE':'07_reference_benchmark'}
def main(model,stages=None,limit=None,shard_index=0,shards=1):
 assert_freezes();jobs=[j for j in readjson(R/'00_protocol/PREDICTION_JOBS_FROZEN.json') if j['model']==model and (not stages or j['stage'] in stages)]
 order={'GLOBAL_CORE':0,'REFERENCE':1,'GNRS_GENE_CENTER':2,'PLACEMENT':3,'CONTEXT':4,'ENSEMBLE':5};jobs.sort(key=lambda j:(str(j['member']),order[j['stage']],j['audit_job_ID']))
 jobs=[j for i,j in enumerate(jobs) if i%shards==shard_index]
 if limit:jobs=jobs[:limit]
 genome=py2bit.open(str(P/'03_reference/hg38.2bit'));full=pd.read_csv(R/'02_track_harmonization/TRACK_METADATA_FULL.tsv',sep='\t').fillna('');md=[z for z in full.to_dict('records') if z['model']==model and z['assay'] in ['RNA','RNA_SEQ']];track_table=pd.read_csv(P/'16_phase1r_rescue/03_borzoi/targets_human.txt',sep='\t',index_col=0).fillna('')
 if model=='Borzoi':
  for t in md:t.update(track_table.loc[int(t['track_ID'])].to_dict())
  from borzoi_runtime import load_model,predict,unsquash
  net=None;member=None;indices=[int(t['track_ID']) for t in md]
 else:
  from run_technical import key_from_document,Deadline,Details
  from alphagenome.models.v1 import dna_client
  from alphagenome.data import genome as aggenome
  import grpc
  key=key_from_document();clients={}
  class GlobalDeadline(grpc.StreamStreamClientInterceptor):
   def intercept_stream_stream(self,continuation,d,iterator):
    return continuation(Details(d.method,600,d.metadata,d.credentials,d.wait_for_ready,d.compression),iterator)
 oldby={j['job_ID']:j for j in readjson(OLD/'00_protocol/JOBS_FROZEN.json') if j['role']=='SENTINEL'};n=0;starttime=time.monotonic()
 for job in jobs:
  outdir=R/STAGE_FOLDER[job['stage']]/model;outdir.mkdir(exist_ok=True);path=outdir/(job['audit_job_ID']+'.json')
  if path.exists():continue
  t=time.monotonic();result={'audit_job_ID':job['audit_job_ID'],'original_ID':job['job_ID'],'model':model,'stage':job['stage'],'member':job['member'],'strand':job['strand'],'start0':job['start0'],'W':job['W'],'placement_shift':job.get('placement_shift'),'context_length':job.get('context_length'),'mask_length':job.get('mask_length'),'gene':job['gene'],'source_beta':job.get('beta'),'status':'PENDING'}
  try:
   intervals=regions(job,model,job.get('fixed_domains'))
   if not intervals:raise ValueError('NO_COMMON_EXON')
   if not job['exons'] or job['strand'] not in ['+','-']:raise ValueError('TRANSCRIPT_MAPPING')
   aseq=sequence(genome,job,model,'REF');zseq=sequence(genome,job,model,'ALT');result['sequence_REF_sha256']=hashlib.sha256(aseq.encode()).hexdigest();result['sequence_ALT_sha256']=hashlib.sha256(zseq.encode()).hexdigest();result['homologous_exon_intervals']=intervals
   # Cached native arrays reused only for exactly identical original contexts/member/available track scope.
   reuse=job['track_scope']=='LCL' and job['member']==('ALL_FOLDS' if model=='AlphaGenome' else '0') and job['job_ID'] in oldby and job['W']==(1048576 if model=='AlphaGenome' else 524288) and job['start0']==job['anchor0']-job['W']//2 and not job.get('mask_length')
   cached=readjson(OLD/'10_statistics'/model/(job['job_ID']+'.json')) if reuse else None;observed=None
   for allele,s in [('REF',aseq),('ALT',zseq)]:
    if cached:
     cid=cached['cache_'+allele];meta=readjson(OLD/'09_predictions'/model/(cid+'.json'));z=np.load(OLD/'09_predictions'/model/(cid+'.npz'));v=z['rna'];res=int(z['resolution']);z.close();observed=[]
     for tmeta in meta['metadata']:
      idx=int(tmeta['track_index']) if model=='Borzoi' else next(int(x['track_ID']) for x in md if x['track_name']==tmeta['name'] and x['strand']==tmeta['strand']);tmd=next(dict(x) for x in md if int(x['track_ID'])==idx);observed.append(tmd)
     if model=='Borzoi':v=unsquash(v,observed)
     result['cached_Phase1R2_arrays']=True
    elif model=='AlphaGenome':
     client_key=(job['member'],job['track_scope'])
     if client_key not in clients:
      c=dna_client.create(key,model_version=getattr(dna_client.ModelVersion,job['member']),timeout=30);c._channel=grpc.intercept_channel(c._channel,GlobalDeadline() if job['track_scope']=='GLOBAL' else Deadline());clients[client_key]=c
     if job['track_scope']=='GLOBAL':
      # Official ontology batches retain every RNA track without a multi-GB response.
      agmeta=pd.read_csv(OLD/'01_alphagenome/rna_seq_ALL_TRACKS.tsv',sep='\t').fillna('')
      counts=agmeta.groupby('ontology_curie',sort=True).size();batches=[];current=[];size=0
      for curie,count in counts.items():
       if current and size+count>32:batches.append(current);current=[];size=0
       current.append(curie);size+=int(count)
      if current:batches.append(current)
      observed=[];means=[];bps=[];shapes=[]
      cache_dir=outdir/'ontology_checkpoints'/job['audit_job_ID']/allele;cache_dir.mkdir(parents=True,exist_ok=True)
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
      assert len(observed)==len(md) and len({x['track_ID'] for x in observed})==len(md) and len(set(bps))==1
      result[allele]={'status':'PASS','mapped_bp':bps[0],'track_means':means};result['metadata_array_shape']=[job['W'],len(observed)];result['ontology_batch_N']=batches_count;result['global_mode']='ALL_RNA_TRACKS_BATCHED_BY_ONTOLOGY';continue
     error=None
     for attempt in range(3):
      try:
       out=clients[client_key].predict_sequence(s,organism=dna_client.Organism.HOMO_SAPIENS,requested_outputs=[dna_client.OutputType.RNA_SEQ],ontology_terms=None if job['track_scope']=='GLOBAL' else ['EFO:0002784'],interval=aggenome.Interval(job['chrom'],job['start0'],job['start0']+job['W']));error=None;break
      except Exception as e:error=e;time.sleep(2**attempt)
     if error:raise error
     v=out.rna_seq.values;res=int(out.rna_seq.resolution);assert len(v)*res==job['W'];observed=[]
     for z in out.rna_seq.metadata.to_dict('records'):
      matches=[x for x in md if x['track_name']==z['name'] and x['strand']==z['strand']];assert len(matches)==1;observed.append(matches[0])
     result['metadata_array_shape']=list(v.shape)
    else:
     if member!=job['member']:
      if net is not None:net=None;gc.collect()
      net=load_model(job['member'],indices);member=job['member']
     v=unsquash(predict(net,s),md);res=32;observed=md
    assert np.isfinite(v).all();result[allele]=summarize(v,intervals,job['start0'],163840 if model=='Borzoi' else 0,res,allele);del v
    if model=='AlphaGenome' and not cached:del out
   result['metadata']=[{k:t[k] for k in ['track_ID','track_name','strand','cell_tissue','selection_family']} for t in observed];result['status']='PASS';result['delta']=delta_score(result)
   mapping=readjson(R/'02_track_harmonization/TRACK_MAPPING_FROZEN.json')['TIER1_EXACT'];fam={k:mapping['tracks'][model][k] for k in mapping['shared_families']};result['tier1_delta']=delta_score(result,fam)
  except Exception as e:
   result['status']='UNAVAILABLE';result['error_type']=type(e).__name__;msg=str(e);msg=msg.replace(key,'[REDACTED]') if model=='AlphaGenome' else msg;result['error']=re.sub(r'AIza[\w-]+','[REDACTED]',msg)[:1500]
  result['elapsed_seconds']=round(time.monotonic()-t,2);path.write_text(json.dumps(result,indent=2));n+=1
  if n%10==0 or limit:print(model,job['stage'],n,'new jobs',result['status'],'seconds',result['elapsed_seconds'],flush=True)
 if model=='AlphaGenome':
  for c in clients.values():c._channel.close()
 print('PREDICTION_RUN_FINISHED',model,'new',n,'seconds',round(time.monotonic()-starttime),flush=True)
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--model',choices=['AlphaGenome','Borzoi'],required=True);a.add_argument('--stages',nargs='*');a.add_argument('--limit',type=int);a.add_argument('--shard-index',type=int,default=0);a.add_argument('--shards',type=int,default=1);z=a.parse_args();main(z.model,z.stages,z.limit,z.shard_index,z.shards)
