from pathlib import Path
import sys,json,numpy as np
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'00_protocol'));from audit_common import *
rows=[];W=524288;start=1000000;crop=163840;res=32
for location,anchor in [('before_output',start+10000),('inside_output',start+262144),('after_output',start+500000)]:
 for L in [13,32,300]:
  j={'anchor0':anchor,'insertion':'A'*L,'W':W,'start0':start,'exons':[[start+crop-100,start+W-crop+100]],'strict':False};iv=regions(j,'Borzoi');ref=np.full((6144,1),2.);alt=np.full((6144,1),4.);a=summarize(ref,iv,start,crop,res,'REF');b=summarize(alt,iv,start,crop,res,'ALT');assert a['mapped_bp']>0 and abs(b['track_means'][0]-a['track_means'][0]-2)<1e-12;rows.append({'anchor_case':location,'insertion_bp':L,'status':'PASS','mapped_bp':a['mapped_bp']})
(R/'05_coordinate_robustness/GENE_CENTER_EDGE_QA.json').write_text(json.dumps({'status':'PASS','analytic_fixture_not_biological_data':True,'checks':rows},indent=2));print('GENE_CENTER_EDGE_COORDINATE_QA_PASS')
