"""Guarded prospective sampling;never produces a formal panel while any gate is HOLD."""
from pathlib import Path
import csv,json,hashlib,math,collections,random
R=Path(__file__).resolve().parents[1];SEED=20261010

def main():
 gates=list(csv.DictReader((R/'13_reports/SENTINEL_GATES.tsv').open(),delimiter='\t'))
 if len(gates)!=5 or any(g['status']!='PASS' for g in gates):raise SystemExit('NOT_RUN:formal sentinel requires all S1-S5 PASS')
 import py2bit
 ref=py2bit.open(str(R.parent/'03_reference/hg38.2bit'));cov=list(csv.DictReader((R.parent/'16_phase1r_rescue/06_covariates/OFFREF_COVARIATE_COMPLETE.tsv').open(),delimiter='\t'));by={r['OFFREF_ID']:r for r in cov};primary=list(csv.DictReader((R/'03_gnrs_truth/GNRS_CLASSB_PRIMARY33.tsv').open(),delimiter='\t'));assert len(primary)==33
 truth=[(r['OFFREF_ID'],'GNRS') for r in primary];reserved={v for v,s in truth}
 ex=list(csv.DictReader((R/'04_external_qtl_truth/EXTERNAL_CLASSB_DIRECTION_TRUTH.tsv').open(),delimiter='\t'));ex=[r for r in ex if r['direction_eligible']=='TRUE' and r['independent_cohort_eligible']=='TRUE'];ex.sort(key=lambda r:(float(r['nominal_P']),r['external_variant_ID'],r['gene']))
 for r in ex:
  if r['OFFREF_ID'] not in reserved and len(truth)<100:truth.append((r['OFFREF_ID'],'INDEPENDENT_EXTERNAL_QTL'));reserved.add(r['OFFREF_ID'])
 def compatible(r):
  b=int(r['anchor']);return r['chrom'] in ref.chroms() and b>=524288 and b+524288<=ref.chroms(r['chrom']) and int(r['length'])<262144
 eligible=[r for r in cov if r['OFFREF_ID'] not in reserved and compatible(r)];assert len(eligible)>=500
 def stratum(r):
  n=int(r['length']);l='50-199' if n<200 else '200-499' if n<500 else '500-1999' if n<2000 else '2000-5000' if n<=5000 else '>5000'
  rc=r['repeat_class'];rep='non-detected-repeat' if rc=='NO_DETECTED_REPEAT' else 'SINE' if rc.startswith('SINE') else 'LINE' if rc.startswith('LINE') else 'LTR' if rc.startswith('LTR') else 'VNTR-STR' if rc in ['VNTR','Simple_repeat'] else 'other'
  try:a=float(r['SV_frequency']);af='rare' if a<.01 else 'low' if a<.05 else 'common'
  except ValueError:af='unknown'
  t=float(r['distance_to_TSS']);ts='<10kb' if t<10000 else '10-100kb' if t<=100000 else '>100kb';return '|'.join([l,rep,af,ts])
 groups=collections.defaultdict(list)
 for r in eligible:groups[stratum(r)].append(r)
 quotas={s:math.floor(400*len(rr)/len(eligible)) for s,rr in groups.items()};order=sorted(groups,key=lambda s:(-(400*len(groups[s])/len(eligible)-quotas[s]),s))
 for s in order[:400-sum(quotas.values())]:quotas[s]+=1
 core=[];rng=random.Random(SEED)
 for s in sorted(groups):
  pool=sorted(groups[s],key=lambda r:r['OFFREF_ID']);core.extend(rng.sample(pool,quotas[s]))
 selected={r['OFFREF_ID'] for r in core}|reserved;extension=rng.sample(sorted([r for r in eligible if r['OFFREF_ID'] not in selected],key=lambda r:r['OFFREF_ID']),100-len(truth));out=[]
 for r in core:out.append({'OFFREF_ID':r['OFFREF_ID'],'panel':'UNBIASED_CLASSB_CORE','truth_source':'NONE','stratum':stratum(r),'seed':SEED})
 for v,s in truth:out.append({'OFFREF_ID':v,'panel':'BIOLOGICAL_TRUTH','truth_source':s,'stratum':stratum(by[v]),'seed':SEED})
 for r in extension:out.append({'OFFREF_ID':r['OFFREF_ID'],'panel':'RANDOM_VALIDATION_EXTENSION','truth_source':'NONE','stratum':stratum(r),'seed':SEED})
 assert len(out)==len({r['OFFREF_ID'] for r in out})==500 and len(core)==400
 p=R/'07_sentinel/SENTINEL500_FROZEN.tsv';assert not p.exists(),'Immutable panel already exists'
 with p.open('w') as f:w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
 (R/'07_sentinel/SENTINEL500_SHA256.tsv').write_text('path\tsha256\nSENTINEL500_FROZEN.tsv\t'+hashlib.sha256(p.read_bytes()).hexdigest()+'\n')
 (R/'07_sentinel/SAMPLING_AUDIT.json').write_text(json.dumps({'seed':SEED,'source_ClassB_N':len(cov),'core_eligible_N_after_truth_reservation_and_declared_window_compatibility':len(eligible),'compatibility':'canonical hg38 chromosome;full centered1Mb window;insertion<262144bp;no prediction-informed filter','truth_N':len(truth),'extension_N':len(extension),'stratum_quotas':quotas},indent=2))
if __name__=='__main__':main()
