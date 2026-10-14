from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
from scipy.stats import spearmanr
R=Path(__file__).resolve().parents[1]
for name,h in json.loads((R/'SHA256SUMS.json').read_text()).items():
    assert hashlib.sha256((R/name).read_bytes()).hexdigest()==h,name
core=pd.read_csv(R/'source_data/core400_complete.tsv',sep='\t')
ref=pd.read_csv(R/'source_data/reference500_complete.tsv',sep='\t')
a=core.AlphaGenome_original.to_numpy(); b=core.Borzoi_original.to_numpy(); ok=np.isfinite(a)&np.isfinite(b)
assert len(core)==400 and core.gene.nunique()==387 and ok.sum()==318
assert ((np.sign(a[ok])==np.sign(b[ok]))).sum()==149
assert abs(spearmanr(a[ok],b[ok]).statistic-(-0.1113027120530779))<1e-12
assert len(ref)==500 and ref.gene.nunique()==473
assert abs(spearmanr(ref.AlphaGenome,ref.Borzoi).statistic-0.6258)<0.00005
print('PASS: file hashes, 400/387/318/149 core endpoints and 500/473 reference endpoints')
