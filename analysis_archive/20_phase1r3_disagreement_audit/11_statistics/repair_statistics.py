from pathlib import Path
R=Path(__file__).resolve().parents[1]
p=R/'11_statistics/statistics_core.py';s=p.read_text().replace('def boot(a,b,B=10000,cluster=None):','def boot(a,b,B=10000,cluster=None,seed=SEED):').replace('rng=np.random.default_rng(SEED);a=', 'rng=np.random.default_rng(seed);a=');p.write_text(s)
p=R/'11_statistics/finalize_endpoints.py';s=p.read_text().replace("B=boot(bv[:,0],bv[:,1])","B=boot(bv[:,0],bv[:,1],seed=20261012)").replace('cv=boot(cg.AlphaGenome_tier1_cached,cg.Borzoi_tier1_cached)','cv=boot(cg.AlphaGenome_tier1_cached,cg.Borzoi_tier1_cached,seed=20261012)').replace("design_limit='GEUVADIS", "independent_bootstrap_streams=[20261011,20261012],design_limit='GEUVADIS")
s=s.replace("# Gene-cluster sign-flip null respects", "ood_supported=supported[supported.feature.astype(str).str.contains('length|repeat',case=False)] if 'feature' in supported else supported.iloc[:0]\n # Gene-cluster sign-flip null respects")
s=s.replace("or strong is not None)","or len(ood_supported)>0)")
s=s.replace("'strongest_uncertainty_supported_predictor':strong,", "'strongest_uncertainty_supported_predictor':strong,'length_repeat_supported_features_N':len(ood_supported),")
p.write_text(s)
print('INDEPENDENT_BOOTSTRAP_STREAMS_AND_OOD_GATE_SCOPE_FIXED')
