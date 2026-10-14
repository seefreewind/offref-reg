from pathlib import Path
import urllib.request,urllib.error,json,re,datetime
R=Path(__file__).resolve().parents[1];O=R/'08_truth_coverage';rows=[]
urls=[('paper','https://www.medrxiv.org/content/10.64898/2026.08.21.26361050v1.full'),('supplements','https://www.medrxiv.org/content/10.64898/2026.08.21.26361050v1.supplementary-material'),('author_tree','https://api.github.com/repos/jiadong324/1KG_LongRead_SV/git/trees/main?recursive=1'),('author_readme','https://raw.githubusercontent.com/jiadong324/1KG_LongRead_SV/main/README.md'),('zenodo','https://zenodo.org/api/records/22000872'),('1KCP_portal','https://yanglab.westlake.edu.cn/1kcp/'),('1KCP_readme','https://raw.githubusercontent.com/JianYang-Lab/1KCP-analysis/main/05_eQTL/README.md')]
for name,url in urls:
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 (research reproducibility audit)'});f=urllib.request.urlopen(req,timeout=45);b=f.read();p=O/(name+'.raw');p.write_bytes(b);text=b.decode(errors='replace');links=re.findall(r'(?:href|src)=[\"\']([^\"\']+)',text);candidates=[x for x in links if re.search(r'supp|table|eqtl|\.xlsx|\.tsv|\.csv',x,re.I)];rows.append({'resource':name,'URL':url,'HTTP':f.status,'saved':str(p.relative_to(R)),'candidate_links':candidates[:100]})
  if name=='zenodo':rows[-1]['files']=[x['key'] for x in json.loads(b).get('files',[])]
  if name=='author_tree':rows[-1]['candidate_files']=[x['path'] for x in json.loads(b).get('tree',[]) if re.search('eqtl|supp|table|beta',x['path'],re.I)]
 except Exception as e:rows.append({'resource':name,'URL':url,'access_error':type(e).__name__,'reason':str(e)[:500]})
 print('FINAL_ACCESS_AUDIT',name,'OK' if 'access_error' not in rows[-1] else rows[-1]['access_error'],flush=True)
(O/'FINAL_EXTERNAL_ACCESS_EVIDENCE.json').write_text(json.dumps({'time_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'one_final_targeted_pass':True,'records':rows,'status':'REVIEW_ACTUAL_LINKS_BEFORE_FINAL_VERDICT'},indent=2))
