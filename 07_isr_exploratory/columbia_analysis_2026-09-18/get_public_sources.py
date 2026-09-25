from pathlib import Path
import urllib.request, json, re
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parent/'sources';ROOT.mkdir(exist_ok=True)
urls={
 'isr_paper.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11838536/fullTextXML',
 'isr_published.html':'https://www.nature.com/articles/s42003-026-10312-x',
 'isr_github_search.json':'https://api.github.com/search/repositories?q=GDF15+ISR',
 'hallmark_2024.gmt':'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt',
}
def fetch(item):
 name,url=item
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'Research-source-inspection'})
  data=urllib.request.urlopen(req,timeout=45).read();(ROOT/name).write_bytes(data)
  s=data.decode('utf-8','replace')
  links=sorted(set(re.findall(r'https?[^\s<>"\]]*(?:github|supplement|media)[^\s<>"\]]*',s)))
  if name.endswith('.json'):
   obj=json.loads(s);links=[x['html_url'] for x in obj.get('items',[])]
  print(name,len(data),links)
 except Exception as e:print(name,'ERROR',str(e))
with ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(fetch,urls.items()))
