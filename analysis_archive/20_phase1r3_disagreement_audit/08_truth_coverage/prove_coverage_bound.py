from pathlib import Path
import json
R=Path(__file__).resolve().parents[1];jobs=[j for j in json.loads((R/'00_protocol/PREDICTION_JOBS_FROZEN.json').read_text()) if j['stage']=='GNRS_GENE_CENTER' and j['model']=='Borzoi'];rows=[]
for j in jobs:
 A=j['anchor0'];L=len(j['insertion']);W=j['W'];crop=163840
 # Every legal window start s is in [A+L-W,A]. Union of possible reference-coordinate output support gives this relaxed bound. Ignore bin restrictions, so it can only overestimate actual recoverability.
 lo=A+L-W+crop;hi=A+W-crop-L;bp=sum(max(0,min(b,hi)-max(a,lo)) for a,b in j['exons'])
 rows.append({'OFFREF_ID':j['job_ID'],'frozen_gene':j['gene'],'insertion_bp':L,'any_legal_B_window_reachable_exon_bp_upper_bound':bp,'possible_annotation_and_output':bp>0,'complete_insertion_required':True})
N=sum(z['possible_annotation_and_output'] for z in rows);(R/'08_truth_coverage/GNRS_ARCHITECTURAL_COVERAGE_UPPER_BOUND.json').write_text(json.dumps({'frozen_N':33,'possibly_reachable_loci_upper_bound':N,'method':'union of nativeBreference-output support over every window containing whole insertion;RNA values never consulted;relaxed bin/homology upper bound','scope':'frozen target-gene native exon annotation;additional missing gene models do not count as negative biological effects','recovery_target20_feasible_under_frozen_gene_exons':N>=20,'rows':rows},indent=2));print('ARCHITECTURAL_GNRS_POSSIBLE_UPPER_BOUND',N,flush=True)
