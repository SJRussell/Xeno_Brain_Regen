from pathlib import Path
from lxml import html, etree
import urllib.request,json,re,sys
sys.stdout.reconfigure(encoding='utf-8')
root=Path(__file__).resolve().parent/'sources'
tree=html.fromstring((root/'isr_published.html').read_text(encoding='utf-8'))
txt='\n'.join(tree.xpath('//main//p//text()'))
(root/'isr_published_text.txt').write_text(txt,encoding='utf-8')
for key in ['Code availability','correlation matrix','scaled','regression','weights','119','Supplementary Data']:
 print('\nKEY',key)
 for m in list(re.finditer(re.escape(key),txt,re.I))[:6]:print(txt[max(0,m.start()-120):m.end()+500])
repo='mitopsychobio/2025_ISR_GDF15_JS_v2'
def get(url):return urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Research-source-inspection'}),timeout=45).read()
obj=json.loads(get('https://api.github.com/repos/'+repo));branch=obj['default_branch']
files=json.loads(get(f'https://api.github.com/repos/{repo}/git/trees/{branch}?recursive=1'))
(root/'isr_repo_tree.json').write_text(json.dumps(files,indent=2))
print('\nREPO',repo,files['sha'])
for f in files['tree']:
 if f['type']=='blob':print(f['path'],f.get('size'))
for n in [2,3,4,5,6]:
 name=f'42003_2026_10312_MOESM{n}_ESM.xlsx'
 url='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs42003-026-10312-x/MediaObjects/'+name
 try:
  data=get(url);(root/name).write_bytes(data);print(name,len(data))
 except Exception as e:print(name,str(e))
