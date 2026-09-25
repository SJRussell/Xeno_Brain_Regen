"""Insert the verified exploratory ISR result into the author's latest saved v15.

Text, legend, methods, and reference are tracked. Figure insertion is untracked.
No existing manuscript text, revision, comment, image, or citation is replaced.
"""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from copy import deepcopy
from datetime import datetime, timezone
import hashlib, json, re, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from lxml import etree as E

sys.stdout.reconfigure(encoding='utf-8')
ROOT = Path(__file__).resolve().parent
ANALYSIS = ROOT.parent / 'columbia_analysis_2026-09-18'
SRC = ROOT / 'source_v15_Sep18_snapshot.docx'
OUT = ROOT / 'Brain-Limb Communication during Limb Regeneration - v16_Sep18_ISR_tracked.docx'
ASSETS = ROOT / 'ISR_assets'
ASSETS.mkdir(exist_ok=True)
NS = {'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
      'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
      'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
      'pic':'http://schemas.openxmlformats.org/drawingml/2006/picture'}
q = lambda s: '{'+NS['w']+'}'+s
text = lambda x: ''.join(x.xpath('.//w:t[not(ancestor::w:del)]/text()', namespaces=NS))
sha = lambda b: hashlib.sha256(b).hexdigest()

# Recheck the primary contrasts and derive plotting values from the saved animal scores.
scores = pd.read_csv(ANALYSIS/'results/sample_scores.tsv', sep='\t')
tests = pd.read_csv(ANALYSIS/'results/signature_tests.tsv', sep='\t')
primary = tests[(tests.aggregation=='mean') & (tests.signature=='ISR_exclude_GDF15')]
limb = scores[(scores.tissue=='BLAST') & (scores.aggregation=='mean')].copy()
assert len(limb)==18 and limb.groupby(['treatment','time']).size().eq(3).all()
assert len(primary)==17
assert len(tests[(tests.aggregation=='mean') & (tests.family=='ISR')])==119
expected = {'CT_vs_ND_11HPA':(1.400949,.007584),
            'CT_vs_ND_24HPA':(1.955379,.001394),
            'CT_vs_ND_7DPA':(.072391,.922810),
            'INTERACTION_7DPA_vs_24HPA':(-1.882987,.010071)}
for contrast,(effect,fdr) in expected.items():
    row=primary[(primary.tissue=='BLAST') & (primary.contrast==contrast)].iloc[0]
    assert abs(row.estimate-effect)<1e-6 and abs(row.q_family-fdr)<1e-6
for tm in ['11HPA','24HPA','7DPA']:
    d=limb[limb.time==tm]
    difference=d.loc[d.treatment=='CT','ISR_exclude_GDF15'].mean()-d.loc[d.treatment=='ND','ISR_exclude_GDF15'].mean()
    assert abs(difference-expected['CT_vs_ND_'+tm][0])<1e-6

coefficients=pd.read_csv(ANALYSIS/'inputs/isr_adapted_coefficients.tsv',sep='\t')
for tissue,n in [('BLAST',101),('BRAIN',97)]:
    c=coefficients[(coefficients.tissue==tissue)&(coefficients['mode']=='exclude_GDF15')]
    assert len(c)==n and 'GDF15' not in set(c.gene)
genes=pd.read_csv(ANALYSIS/'results/target_genes_original_DESeq2.tsv',sep='\t')
g=genes[(genes.tissue=='BLAST') & (genes.contrast=='CT_vs_ND_24HPA') & genes.gene.isin(['atf3.L','atf3.S','ddit3.L'])]
assert len(g)==3 and g.log2FoldChange.gt(0).all() and g.padj.le(1.7e-6).all()

plt.rcParams.update({'font.family':'Arial','font.size':14,'axes.spines.top':False,
                     'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(6.5,3.65))
fig.subplots_adjust(left=.16,right=.98,bottom=.16,top=.91)
colors={'CT':'#BA5148','ND':'#3C6D88'}
times=['11HPA','24HPA','7DPA']
for treatment,marker,label in [('ND','o','ND'),('CT','s','MDT')]:
    means=[];xs=[]
    for i,tm in enumerate(times):
        values=limb[(limb.treatment==treatment)&(limb.time==tm)].sort_values('sample').ISR_exclude_GDF15.to_numpy()
        x=i+(-.12 if treatment=='ND' else .12)
        ax.scatter(x+np.array([-.035,0,.035]),values,s=37,marker=marker,
                   color=colors[treatment],edgecolor='white',linewidth=.5,zorder=4)
        ax.errorbar(x,values.mean(),yerr=values.std(ddof=1)/np.sqrt(3),
                    color=colors[treatment],capsize=4,linewidth=1.4,zorder=3)
        means.append(values.mean());xs.append(x)
    ax.plot(xs,means,color=colors[treatment],lw=1.5,label=label)
ax.set_xlim(-.42,2.42);ax.set_ylim(-1.85,2.5)
ax.set_xticks(range(3),['11 HPT','24 HPT','7 DPT'])
ax.set_yticks([-1,0,1,2])
ax.set_ylabel('Adapted ISR score\n(GDF15 excluded)',fontsize=14)
ax.axhline(0,color='#d2d2d2',lw=.8,zorder=0)
fig.legend(handles=[Line2D([],[],color=colors[k],marker=m,label=label,lw=1.5,markersize=5)
                    for k,m,label in [('ND','o','ND'),('CT','s','MDT')]],
           frameon=False,loc='upper center',bbox_to_anchor=(.61,1.015),ncol=2,
           handlelength=1.5,columnspacing=1.1,fontsize=13)
for i,tm in enumerate(times):
    row=primary[(primary.tissue=='BLAST') & (primary.contrast=='CT_vs_ND_'+tm)].iloc[0]
    label=f'q = {row.q_family:.4f}' if row.q_family<.01 else f'q = {row.q_family:.3f}'
    ax.plot([i-.20,i-.20,i+.20,i+.20],[1.84,1.91,1.91,1.84],color='#444444',lw=.8)
    ax.text(i,2.01,label,ha='center',va='bottom',fontsize=12.5)
for ext in ['png','pdf','svg']:
    fig.savefig(ASSETS/f'Supplementary_Figure_3_ISR.{ext}',dpi=400,facecolor='white')
plt.close(fig)
limb[['sample','treatment','time','ISR_exclude_GDF15']].to_csv(ASSETS/'Supplementary_Figure_3_source_data.tsv',sep='\t',index=False)
tests[tests.family=='ISR'].to_csv(ASSETS/'ISR_all_contrasts_and_sensitivities.tsv',sep='\t',index=False)
coefficients.to_csv(ASSETS/'ISR_adapted_coefficients_and_coverage.tsv',sep='\t',index=False)

with ZipFile(SRC) as zin:
    parts={n:zin.read(n) for n in zin.namelist()}
doc=E.fromstring(parts['word/document.xml']);original=deepcopy(doc)
body=doc.find('w:body',NS);ps=body.findall('w:p',NS)
assert len(ps)==274
assert text(ps[91]).startswith('Since brain and blastema')
assert text(ps[126]).startswith('GDF15 emerged')
assert text(ps[176]).startswith('Downstream analyses')
assert text(ps[264]).startswith('76. Dolgalev')
uid=max([int(x.get(q('id'))) for x in doc.iter() if (x.get(q('id')) or '').isdigit()]+[20000])
date=datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
new_ids=[];new_paragraphs=[];changes=[]

def revision(kind='ins'):
    global uid
    uid+=1;new_ids.append(str(uid))
    return E.Element(q(kind),{q('id'):str(uid),q('author'):'Codex',q('date'):date})

def run(s,template=None,bold=False,italic=False,superscript=False):
    r=E.Element(q('r'))
    if template is not None and template.find(q('rPr')) is not None:
        pr=deepcopy(template.find(q('rPr')));r.append(pr)
    else:pr=E.SubElement(r,q('rPr'))
    for k in ['i','iCs','b','bCs','u','vertAlign','color','highlight']:
        for e in pr.findall(q(k)):pr.remove(e)
    if bold:E.SubElement(pr,q('b'))
    if italic:E.SubElement(pr,q('i'))
    if superscript:E.SubElement(pr,q('vertAlign'),{q('val'):'superscript'})
    t=E.SubElement(r,q('t'),{'{http://www.w3.org/XML/1998/namespace}space':'preserve'});t.text=s
    return r

def add_text(p,s,template=None,bold=False):
    ins=revision()
    # Match existing nomenclature: italic Xenopus gene symbols, roman human symbol/index.
    gene_pattern=r'gdf15(?:\.[LS])?|atf3(?:\.[LS])?|ddit3(?:\.[LS])?|Xenopus'
    pattern='('+gene_pattern+r'|(?<=10)−6)'
    for segment in re.split(pattern,s):
        if segment:ins.append(run(segment,template,bold=bold,italic=bool(re.fullmatch(gene_pattern,segment)),superscript=segment=='−6'))
    p.append(ins)

def paragraph(s,template_p,template_r=None,bold=False):
    p=E.Element(q('p'));pr=deepcopy(template_p.find(q('pPr')))
    if pr is None:pr=E.Element(q('pPr'))
    p.append(pr)
    for k in ['sectPr','pPrChange']:
        for e in pr.findall(q(k)):pr.remove(e)
    rp=pr.find(q('rPr'))
    if rp is None:rp=E.SubElement(pr,q('rPr'))
    for e in list(rp):
        if e.tag in [q('ins'),q('del')]:rp.remove(e)
    rp.append(revision()) # Track the paragraph mark as well as its contents.
    if s:add_text(p,s,template_r or template_p.find('w:r',NS),bold=bold)
    new_paragraphs.append(p)
    return p

result=(
    'To test whether the limb response extended beyond gdf15, we applied an exploratory Xenopus adaptation of a human integrated stress response (ISR) index (Methods). '
    'The primary score excluded GDF15 and retained 101 of 119 ISR genes in the limb. '
    'MDT animals had higher limb scores at 11 HPT (MDT−ND difference=+1.40, BH-adjusted P=0.0076) and 24 HPT (+1.96, adjusted P=0.0014), '
    'with little treatment separation at 7 DPT (+0.07, adjusted P=0.923; Supplementary Fig. 3). '
    'A direct treatment-by-time contrast supported attenuation between 24 HPT and 7 DPT (difference in treatment effects=−1.88, adjusted P=0.010). '
    'Both atf3 homeologs and ddit3.L were also upregulated in the limb at 24 HPT (genome-wide adjusted P≤1.7×10−6). '
    'The early score elevations persisted when the higher-expressed homeolog represented each human ortholog. '
    'A corresponding brain score increase was not established at either sampled timepoint (97/119 genes; adjusted P=0.992 at 24 HPT and 0.196 at 7 DPT).'
)
p=paragraph(result,ps[91]);ps[91].addnext(p);changes.append({'section':'Results','text':result})

discussion=(
    ' The early limb ISR signature and its attenuation by day 7 place a local stress-associated transcriptional response before the later brain gdf15 increase (Supplementary Fig. 3). '
    'This pattern is consistent with an early response to regenerative induction, but the exploratory cross-species score does not establish biochemical ISR activation, a regeneration-specific mechanism, or propagation of an ISR program to the brain.'
)
add_text(ps[126],discussion,ps[126].find('w:r',NS));changes.append({'section':'Discussion','text':discussion})

heading=paragraph('Exploratory integrated stress response scoring',ps[151],bold=False)
heading.find('w:pPr/w:spacing',NS).set(q('before'),'240')
# Keep the established underlined Methods heading style.
for rp in heading.xpath('.//w:ins/w:r/w:rPr',namespaces=NS):E.SubElement(rp,q('u'),{q('val'):'single'})
ps[176].addnext(heading)
methods_prefix='The human ISR–GDF15 index of Smith et al.'
mp=paragraph(methods_prefix,ps[149]);heading.addnext(mp)

item={
    'id':'ISR_SMITH_2026','type':'article-journal',
    'title':'Mapping the GDF15 arm of the integrated stress response in human cells and tissues',
    'author':[{'given':given,'family':family} for given,family in [
        ('Janell L. M.','Smith'),('Kamaryn T.','Tanner'),('Jack','Devine'),('Anna S.','Monzel'),
        ('Taivan','Batjargal'),('Maxwell Z.','Wilson'),('Alan A.','Cohen'),('Martin','Picard')]],
    'container-title':'Communications Biology','volume':'9','issue':'1','page':'1126',
    'DOI':'10.1038/s42003-026-10312-x','URL':'https://www.nature.com/articles/s42003-026-10312-x',
    'issued':{'date-parts':[[2026,5,27]]}}
citation={'citationID':'ISRv16Sep18','properties':{'formattedCitation':r'\super 77\nosupersub{}','plainCitation':'77','noteIndex':0},
          'citationItems':[{'id':'ISR_SMITH_2026','uris':['https://doi.org/10.1038/s42003-026-10312-x'],'itemData':item}],
          'schema':'https://github.com/citation-style-language/schema/raw/master/csl-citation.json'}
ins=revision()
for kind in ['begin','code','separate','text','end']:
    if kind=='text':r=run('77',superscript=True)
    else:
        r=E.Element(q('r'))
        if kind=='code':
            e=E.SubElement(r,q('instrText'),{'{http://www.w3.org/XML/1998/namespace}space':'preserve'})
            e.text=' ADDIN ZOTERO_ITEM CSL_CITATION '+json.dumps(citation,ensure_ascii=False)+' '
        else:E.SubElement(r,q('fldChar'),{q('fldCharType'):kind})
    ins.append(r)
mp.append(ins)
methods=(
    ' was adapted using curated Xenbase human-ortholog assignments. TMM-normalized log2 counts-per-million were averaged across homeologs assigned to one human symbol and z-standardized per gene within each tissue. '
    'GDF15 was excluded from the primary score, leaving 101/119 genes in limb and 97/119 in brain. '
    'For each observed gene subset, regression scoring coefficients were reconstructed as the inverse human reference correlation matrix multiplied by the published GDF15-associated factor loadings, using the 339 fibroblast training samples and rounded loadings supplied with the published code. '
    'The adapted score was the weighted sum of Xenopus gene z-scores; its absolute scale was not compared between tissues. '
    'Ordinary linear models of animal-level scores used treatment–timepoint groups within each tissue (three independent animals per group), with two-sided t-tests for treatment, temporal and interaction contrasts. '
    'BH adjustment covered 119 tests in the exploratory ISR family: 17 contrasts across both tissues for six score variants and a canonical ISR target-gene average. '
    'Sensitivity analyses varied GDF15 inclusion, reference-SD scaling and loading-based weights; using the higher-expressed homeolog was analyzed as a separate 119-test family. '
    'The ISR analysis used R 4.6.1 and limma 3.68.5, without empirical-Bayes moderation or gene-level precision weights for the composite scores.'
)
add_text(mp,methods,ps[149].find('w:r',NS))
changes.append({'section':'Methods','text':methods_prefix+'77'+methods})

reference='77. Smith, J. L. M. et al. Mapping the GDF15 arm of the integrated stress response in human cells and tissues. Communications Biology vol. 9 1126 (2026). doi:10.1038/s42003-026-10312-x.'
ref=paragraph('',ps[264]);ins=revision();ins.append(run('77. '))
for j,segment in enumerate(re.split(r'(et al\.|Communications Biology)',reference[len('77. '):])):
    r=run(segment,italic=segment in ['et al.','Communications Biology'])
    if j==0:r.insert(1,E.Element(q('tab')))
    ins.append(r)
ref.append(ins);ps[264].addnext(ref)
changes.append({'section':'References','text':reference})

# Untracked image paragraph kept with its caption; existing images/relationships stay intact.
rels=E.fromstring(parts['word/_rels/document.xml.rels'])
rid='rIdISR20260918';assert not any(x.get('Id')==rid for x in rels)
E.SubElement(rels,'{http://schemas.openxmlformats.org/package/2006/relationships}Relationship',
             {'Id':rid,'Type':NS['r']+'/image','Target':'media/Supplementary_Figure_3_ISR.png'})
image_p=deepcopy(ps[270])
image_p.attrib.clear()
pr=image_p.find(q('pPr'));E.SubElement(pr,q('keepNext'))
for x in image_p.xpath('.//w:lastRenderedPageBreak',namespaces=NS):x.getparent().remove(x)
drawing=image_p.find('.//w:drawing',NS)
blip=drawing.find('.//a:blip',NS);blip.set('{'+NS['r']+'}embed',rid)
cx,cy=round(6.5*914400),round(3.65*914400)
for e in drawing.findall('.//wp:extent',NS)+drawing.findall('.//a:xfrm/a:ext',NS):e.set('cx',str(cx));e.set('cy',str(cy))
new_drawing_id=max(int(x.get('id')) for x in doc.findall('.//wp:docPr',NS))+1
for e in drawing.findall('.//wp:docPr',NS)+drawing.findall('.//pic:cNvPr',NS):
    e.set('id',str(new_drawing_id));e.set('name','Supplementary Figure 3 Early limb ISR')
    e.set('descr','Individual limb ISR scores for three MDT and three ND animals at 11 and 24 hours and 7 days. Early treatment differences attenuate by day 7. GDF15 is excluded from the score.')
for e in drawing.iter():
    for attr in list(e.attrib):
        if attr.rsplit('}',1)[-1] in ['anchorId','editId']:e.set(attr,'20260918')
ps[273].addnext(image_p)
legend_title='Supplementary Figure 3 | Early limb ISR-associated transcription under treatment.'
legend=(
    ' The exploratory Xenopus-adapted ISR score excludes GDF15 and uses 101 of 119 ISR genes. '
    'Points are individual animals (n=3 per treatment and timepoint); lines connect group means and error bars show s.e.m. '
    'MDT, multidrug treatment; ND, no device. Labels show BH-adjusted P values (q) for MDT versus ND from ordinary linear models of animal-level scores, corrected across the 119-test exploratory ISR family. '
    'The treatment difference decreased between 24 HPT and 7 DPT (interaction estimate=−1.88, 95% CI −2.95 to −0.82, q=0.010). '
    'Scores reflect similarity to a human-derived transcriptional signature and are not direct measurements of biochemical ISR activation.'
)
lp=paragraph(legend_title,ps[271],bold=True);add_text(lp,legend,ps[271].find('w:r',NS));image_p.addnext(lp)
changes.append({'section':'Supplementary Figure 3 legend','text':legend_title+legend})

settings=E.fromstring(parts['word/settings.xml'])
if settings.find(q('trackRevisions')) is None:settings.append(E.Element(q('trackRevisions')))
for path,root in [('word/document.xml',doc),('word/settings.xml',settings),('word/_rels/document.xml.rels',rels)]:
    parts[path]=E.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
parts['word/media/Supplementary_Figure_3_ISR.png']=(ASSETS/'Supplementary_Figure_3_ISR.png').read_bytes()
ct=E.fromstring(parts['[Content_Types].xml'])
assert any(x.get('Extension')=='png' for x in ct)

# Structural verification: all original nodes remain byte-identical except the appended Discussion text.
orig_ps=original.find('w:body',NS).findall('w:p',NS)
for i,(old,new) in enumerate(zip(orig_ps,ps)):
    if i==126:
        check=deepcopy(new)
        for child in list(check):
            if child.get(q('id')) in new_ids:check.remove(child)
        assert E.tostring(check,method='c14n',exclusive=True)==E.tostring(old,method='c14n',exclusive=True)
    else:assert E.tostring(old)==E.tostring(new),('changed original paragraph',i)
old_codes=original.xpath('//w:instrText/text()',namespaces=NS)
new_codes=doc.xpath('//w:instrText/text()',namespaces=NS)
assert [x for x in new_codes if 'ISRv16Sep18' not in x]==old_codes
original_revs={x.get(q('id')):E.tostring(x) for x in original.xpath('//w:ins|//w:del',namespaces=NS)}
current_revs={x.get(q('id')):E.tostring(x) for x in doc.xpath('//w:ins|//w:del',namespaces=NS)}
assert all(current_revs[k]==v for k,v in original_revs.items())
assert not drawing.xpath('ancestor::w:ins|ancestor::w:del',namespaces=NS)

# Reject just this edit: remove the newly tracked paragraphs and additions; image remains by request.
check=deepcopy(doc)
for p in check.findall('.//w:body/w:p',NS):
    if any(x.get(q('id')) in new_ids for x in p.findall('w:pPr/w:rPr/w:ins',NS)):
        p.getparent().remove(p)
    else:
        for x in p.xpath('.//w:ins',namespaces=NS):
            if x.get(q('id')) in new_ids:x.getparent().remove(x)
assert [text(p) for p in check.findall('.//w:body/w:p',NS) if text(p)]==[text(p) for p in orig_ps if text(p)]

with ZipFile(OUT,'w',ZIP_DEFLATED) as zout:
    for name,value in parts.items():zout.writestr(name,value)
manifest={'source':str(SRC),'source_sha256':sha(SRC.read_bytes()),'output':str(OUT),
          'output_sha256':sha(OUT.read_bytes()),'date':date,'changes':changes,
          'new_revision_ids':new_ids,'new_tracked_paragraphs':len(new_paragraphs),
          'original_revisions_preserved':len(original_revs),'original_citation_chunks_preserved':len(old_codes),
          'original_comments_part_unchanged':True,'reject_new_text_restores_source':True,
          'figure_insertion_untracked':True,'source_data_verified':True,
          'primary_statistics':primary.to_dict(orient='records')}
(ROOT/'change_manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding='utf-8')
print(OUT)
print('Verified source statistics, all original paragraphs/citations/comments/revisions, and rejection round-trip.')
