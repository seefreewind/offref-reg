from pathlib import Path
import json,re,csv,hashlib,collections
import pandas as pd
R=Path(__file__).resolve().parents[1];P=R.parent;O=R/'08_truth_coverage'
# Convert source sheet4a/b headers exactly. Both contain measured q/slope, not predicted effects.
summary=[]
for number,label in [(2,'Table4a'),(3,'Table4b')]:
 raw=json.loads((O/f'TABLE4_SHEET{number}.json').read_text());headers={re.sub(r'\d+','',k):v for k,v in raw[0].items() if v}
 for r in raw[1:]:
  z={headers.get(re.sub(r'\d+','',k),k):v for k,v in r.items() if re.sub(r'\d+','',k) in headers};z['source_sheet']=label;summary.append(z)
s=pd.DataFrame(summary);s.to_csv(O/'EXTERNAL_PUBLIC_SVEQTL_TABLE4_AB.tsv',sep='\t',index=False)
cat=pd.read_csv(P/'04_offref_catalog/OFFREF_MASTER.tsv',sep='\t') if (P/'04_offref_catalog/OFFREF_MASTER.tsv').exists() else None
print('TABLE4_ROWS',len(s),'columns',list(s.columns),flush=True)
plain=re.sub('<[^>]+>',' ',(O/'paper.raw').read_text());segments=[]
for m in re.finditer('tensorQTL',plain):
 if m.start()>plain.find('METHODS'):segments.append(plain[max(0,m.start()-1600):m.start()+1000])
(O/'EXTERNAL_METHODS_ORIENTATION_EXCERPTS.json').write_text(json.dumps({'source':'paper.raw','targeted_methods_excerpts':segments[:4],'note':'Do not infer slope allele from VCF alone; require pipeline/author orientation evidence.'},indent=2))
(O/'EXTERNAL_PUBLIC_TABLE_STATUS.json').write_text(json.dumps({'status':'EXTERNAL_SVEQTL_PUBLIC_SUMMARY_OBTAINED_ORIENTATION_MATCHING_PENDING','source_URL':'https://www.medrxiv.org/content/medrxiv/early/2026/08/25/2026.08.21.26361050/DC6/embed/media-6.xlsx?download=true','sha256':hashlib.sha256((O/'Supplementary_Table4.xlsx').read_bytes()).hexdigest(),'summary_rows_N':len(s),'source_sheets':['Table4a','Table4b'],'beta_column':'slope','gene_column':'ensemblID','variant_column':'SVID or ID','P_column':'qvalue (adjusted); rawP andSE not provided','effect_allele_explicit_column':False,'exact_match_N':'UNDETERMINED','independent_direction_truth_N':'UNDETERMINED','known_GNRS_MAGE_shared_donors':93,'previous_Phase1R2_unchanged':True,'resource_search':'final targeted pass completed; only downloaded-file/link verification continues'},indent=2))
