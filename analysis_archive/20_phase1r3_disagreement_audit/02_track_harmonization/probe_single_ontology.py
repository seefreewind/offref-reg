from pathlib import Path
import sys,json,time,re,hashlib
R=Path(__file__).resolve().parents[1];sys.path[:0]=[str(R/'00_protocol'),str(R.parent/'19_phase1r2_execution/01_alphagenome')]
from audit_common import *
from run_technical import key_from_document,Deadline
from alphagenome.models.v1 import dna_client
from alphagenome.data import genome as aggenome
import py2bit,grpc
jobs=[j for j in readjson(R/'00_protocol/PREDICTION_JOBS_FROZEN.json') if j['stage']=='GLOBAL_CORE' and j['model']=='AlphaGenome'];j=next(j for j in jobs if regions(j,'AlphaGenome',j['fixed_domains']));g=py2bit.open(str(R.parent/'03_reference/hg38.2bit'));s=sequence(g,j,'AlphaGenome','REF');key=key_from_document();c=dna_client.create(key,model_version=dna_client.ModelVersion.ALL_FOLDS,timeout=30);c._channel=grpc.intercept_channel(c._channel,Deadline());t=time.monotonic();z={'single_ontology':'CL:0000047','same_original_native_sequence':True,'job_ID':j['audit_job_ID']}
try:
 out=c.predict_sequence(s,organism=dna_client.Organism.HOMO_SAPIENS,requested_outputs=[dna_client.OutputType.RNA_SEQ],ontology_terms=['CL:0000047'],interval=aggenome.Interval(j['chrom'],j['start0'],j['start0']+j['W']));z.update(status='PASS',shape=list(out.rna_seq.values.shape),expected='official metadata neuronal stem cell RNA tracks;technical availability only')
except Exception as e:z.update(status='UNAVAILABLE',error_type=type(e).__name__,error=str(e).replace(key,'[REDACTED]')[:1000])
z['elapsed_seconds']=time.monotonic()-t;(R/'02_track_harmonization/SINGLE_ONTOLOGY_TRANSPORT_PROBE.json').write_text(json.dumps(z,indent=2));c._channel.close();print(z['status'],z['elapsed_seconds'],flush=True)
