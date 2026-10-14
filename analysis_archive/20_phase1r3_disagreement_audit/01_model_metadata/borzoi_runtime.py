from pathlib import Path
import collections,json
import h5py,numpy as np,torch
from borzoi_pytorch import Borzoi
from borzoi_pytorch.config_borzoi import BorzoiConfig
R=Path(__file__).resolve().parents[1];P=R.parent

def load_model(member,indices):
 torch.manual_seed(20261009);torch.set_num_threads(4);layers=collections.OrderedDict()
 weight=P/'07_model_predictions/borzoi_official_model0_best.h5' if member=='0' else R/'01_model_metadata'/f'borzoi_official_rep{member}_model0_best.h5'
 assert weight.exists(),'OFFICIAL_WEIGHT_NOT_READY'
 with h5py.File(weight,'r') as f:
  root=f['model_weights'] if 'model_weights' in f else f
  for name in root.attrs['layer_names']:
   name=name.decode() if isinstance(name,bytes) else name;g=root[name];ws=g.attrs.get('weight_names',[]);key=name.replace('sync_batch_normalization','batch_normalization')
   if len(ws):layers[key]=[torch.tensor(np.array(g[w.decode() if isinstance(w,bytes) else w])) for w in ws]
 cells=json.loads((P/'07_model_predictions/weight_conversion.ipynb').read_text())['cells'];source=next(''.join(c['source']) for c in cells if 'def transform_the_transformer' in ''.join(c.get('source',[])));source=source[:source.index('res_transformers =')];ns={'torch':torch,'np':np,'collections':collections};exec(source,ns);state={**ns['transform_the_transformer'](layers),**ns['convert_the_convs'](layers),**ns['normalize_the_norms'](layers)}
 net=Borzoi(BorzoiConfig());missing,unexpected=net.load_state_dict(state,strict=False);assert not unexpected and all(k.endswith('num_batches_tracked') for k in missing) and all(k in state for k,_ in net.named_parameters());net.eval();net.set_track_subset(indices);print('OFFICIAL_MEMBER_ALL_PARAMETERS_LOADED',member,flush=True);return net

def predict(net,sequence):
 lookup=np.zeros((256,4),np.float32)
 for i,base in enumerate('ACGT'):lookup[ord(base),i]=1
 x=torch.from_numpy(lookup[np.frombuffer(sequence.encode(),dtype=np.uint8)].T.copy())[None]
 with torch.inference_mode():out=net(x).cpu().numpy()[0].T
 assert np.isfinite(out).all() and len(out)==6144;return out

def unsquash(v,tracks):
 x=np.array(v,dtype=np.float64)
 for i,t in enumerate(tracks):
  z=x[:,i]/float(t['scale']);clip=float(t['clip_soft']);z=np.where(z>clip,(z-clip)**2+clip,z);z=z**(1/.75) if t['sum_stat']=='sum_sqrt' else z;x[:,i]=z/32
 assert np.isfinite(x).all();return x
