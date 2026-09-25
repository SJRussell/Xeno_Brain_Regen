from pathlib import Path
import json,hashlib,re,sys
import numpy as np
import pandas as pd
from scipy.stats import t as tdist,rankdata,false_discovery_control
from scipy.cluster.hierarchy import linkage,leaves_list
from scipy.spatial.distance import pdist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parents[1];O=Path(__file__).resolve().parent;R=O/'results';F=O/'figures';F.mkdir(exist_ok=True)
read=lambda p:pd.read_csv(p,sep='\t')
s=read(R/'signature_tests.tsv');ss=read(R/'sample_scores.tsv');fg=read(R/'hallmark_fgsea_all.tsv');checks={}
assert s[['estimate','p','ci_low','ci_high','q_family']].notna().all().all()
assert fg[['NES','pval','padj']].notna().all().all()
assert len(fg)==1300 and len(s)==2210
assert not s.duplicated(['tissue','aggregation','signature','contrast']).any()
# Independently recompute every signature estimate, standard error, p and interval from saved designs.
maxdiff={x:0. for x in ['estimate','se','p','ci_low','ci_high']}
for tissue in ['BRAIN','BLAST']:
 d=read(R/f'{tissue}_design.tsv').set_index('sample');c=read(R/f'{tissue}_contrasts.tsv').set_index('group')
 for agg in ['mean','best']:
  score=ss[(ss.tissue==tissue)&(ss.aggregation==agg)].set_index('sample').loc[d.index]
  X=d.to_numpy();df=len(X)-np.linalg.matrix_rank(X)
  for row in s[(s.tissue==tissue)&(s.aggregation==agg)].itertuples():
   y=score[row.signature].to_numpy();cv=c[row.contrast].to_numpy();b=np.linalg.lstsq(X,y,rcond=None)[0]
   se=np.sqrt(np.sum((y-X@b)**2)/df*(cv@np.linalg.inv(X.T@X)@cv));est=cv@b;p=2*tdist.sf(abs(est/se),df)
   independent=[est,se,p,est-tdist.ppf(.975,df)*se,est+tdist.ppf(.975,df)*se]
   for k,v in zip(maxdiff,independent):maxdiff[k]=max(maxdiff[k],abs(v-getattr(row,k)))
for k,v in maxdiff.items():assert v<1e-9,(k,v)
checks['independent_OLS_max_abs_errors']=maxdiff
for (agg,family),ds in s.groupby(['aggregation','family']):
 assert np.allclose(false_discovery_control(ds.p),ds.q_family,rtol=1e-10,atol=1e-12)
checks['FDR_families_per_aggregation']=s[s.aggregation=='mean'].groupby('family').size().to_dict()
checks['signature_tests']=len(s);checks['Hallmark_tests']=len(fg);checks['camera_tests']=len(read(R/'camera_tests.tsv'))
ms=ROOT/'revision_2026-09-17/Brain-Limb Communication during Limb Regeneration - v15_Sep17_tracked.docx'
checks['manuscript_sha256']=hashlib.sha256(ms.read_bytes()).hexdigest()
old=json.loads((ROOT/'revision_2026-09-17/verification.json').read_text())
assert checks['manuscript_sha256']==old['output_sha256'];checks['manuscript_unchanged']=True

# Existing comparator matrix: all 50 Hallmarks, no current transcriptomic input substitution.
H=ROOT/'supplementary_v9/supplementary_v9_staging/hallmark_50_v9redfix'
sources=[(ROOT/'tables/brain/hallmark_nes_BRAIN_24HPT_v10method.tsv','Frog24h'),(ROOT/'tables/brain/hallmark_nes_BRAIN_7DPT_v10method.tsv','Frog7d')]+[(H/f,'Mouse_'+n) for f,n in [('gsea_mPFC_full.tsv','mPFC_SNI'),('gsea_NAc_full.tsv','NAc_SNI'),('gsea_PAG_full.tsv','PAG_SNI'),('gsea_hallmark_bi_vs_nt_full.tsv','burn'),('gsea_hallmark_fa_vs_nt_full.tsv','formalin'),('gsea_hallmark_bi_vs_fa_full.tsv','burn_vs_formalin')]]
mat=pd.concat([read(p).set_index('term').NES.rename(n) for p,n in sources],axis=1).fillna(0)
assert mat.shape==(50,8)
Z=linkage(pdist(mat.to_numpy().T),method='ward');frogs=[r.tolist() for r in Z if set(r[:2])=={0,1}]
assert len(frogs)==1
checks['pain_comparator']={'shape':list(mat.shape),'frog_pair_merge_height':frogs[0][2],'frog_columns_are_sisters':True,'leaf_order':mat.columns[leaves_list(Z)].tolist()}
mat.to_csv(R/'pain_comparator_50_hallmark_matrix.tsv',sep='\t');pd.DataFrame(Z,columns=['left','right','height','members']).to_csv(R/'pain_comparator_linkage.tsv',sep='\t',index=False)

# Original leading edges, including annotations from pre-specified blood panels.
panels=json.loads((O/'inputs/panels.json').read_text());le=read(ROOT/'tables/brain/brain_LE_genes_per_hallmark_CT_vs_ND_7DPT.tsv')
heme=le[le.pathway=='HALLMARK_HEME_METABOLISM'].copy();heme['predefined_erythroid_marker']=heme.human_ortholog.isin(panels['ERYTHROID'])
heme.to_csv(R/'brain_day7_heme_leading_edge_audit.tsv',sep='\t',index=False)
checks['heme_leading_edge']={'unique_symbols':int(heme.human_ortholog.nunique()),'predefined_erythroid_symbols':sorted(heme.loc[heme.predefined_erythroid_marker,'human_ortholog'].unique().tolist()),'other_erythropoietic_regulators_present':sorted(set(heme.human_ortholog)&{'LMO2','TAL1','NFE2','GATA2'})}
leadsets={}
for tag,fn in [('24h','CT_vs_ND_24HPT'),('7d','CT_vs_ND_7DPT'),('MDT_time','CT_7v24'),('ND_time','ND_7v24')]:
 path=ROOT/f'tables/brain/brain_LE_genes_per_hallmark_{fn}.tsv'
 if path.exists():
  d=read(path);term='HALLMARK_MTORC1_SIGNALING';leadsets[tag]=set(d.loc[d.pathway==term,'human_ortholog'].dropna())
checks['mTOR_leading_edges']={k:len(v) for k,v in leadsets.items()}
if all(k in leadsets for k in ['24h','7d']):checks['mTOR_leading_edges']['24h_7d_overlap']=len(leadsets['24h']&leadsets['7d'])
(R/'mTOR_leading_edge_genes.json').write_text(json.dumps({k:sorted(v) for k,v in leadsets.items()},indent=2))

plt.rcParams.update({'font.family':'Arial','font.size':12,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
paths=['OXIDATIVE_PHOSPHORYLATION','GLYCOLYSIS','HYPOXIA','MTORC1_SIGNALING','REACTIVE_OXYGEN_SPECIES_PATHWAY','FATTY_ACID_METABOLISM','HEME_METABOLISM','UNFOLDED_PROTEIN_RESPONSE']
labels=['OXPHOS','Glycolysis','Hypoxia','mTORC1','Reactive oxygen species','Fatty-acid metabolism','Heme metabolism','Unfolded protein response']
fig,axes=plt.subplots(2,1,figsize=(15.5,10.7),gridspec_kw={'height_ratios':[1,1]},layout='constrained')
for ax,tissue in zip(axes,['BRAIN','BLAST']):
 times=['24HPA','7DPA'] if tissue=='BRAIN' else ['11HPA','24HPA','7DPA']
 contrasts=[f'CT_vs_ND_{t}' for t in times]
 for arm in ['CT','ND']:
  for a,b in ([('24HPA','7DPA')] if tissue=='BRAIN' else [('11HPA','24HPA'),('11HPA','7DPA'),('24HPA','7DPA')]):contrasts.append(f'{arm}_{b}_vs_{a}')
 ds=fg[(fg.tissue==tissue)&(fg.mapping=='curated_Xenbase')]
 m=ds.pivot(index='pathway',columns='contrast',values='NES').reindex(index=paths,columns=contrasts)
 q=ds.pivot(index='pathway',columns='contrast',values='padj').reindex(index=paths,columns=contrasts)
 im=ax.imshow(m,cmap='RdBu_r',vmin=-3.5,vmax=3.5,aspect='auto')
 for i in range(len(paths)):
  for j in range(len(contrasts)):
   v=m.iloc[i,j];ax.text(j,i,f'{v:+.2f}'+('*' if q.iloc[i,j]<.05 else ''),ha='center',va='center',fontsize=12,color='white' if abs(v)>2.1 else '#161616')
 def cn(c):return c.replace('CT_vs_ND_','MDT–ND\n').replace('CT_','MDT\n').replace('ND_','ND\n').replace('_vs_','–').replace('24HPA','24 h').replace('11HPA','11 h').replace('7DPA','7 d')
 ax.set_xticks(range(len(contrasts)),[cn(c) for c in contrasts],fontsize=11);ax.set_yticks(range(8),labels)
 ax.set_title('Brain' if tissue=='BRAIN' else 'Amputation-site tissue',loc='left',fontweight='bold',fontsize=15)
 ax.axvline(len(times)-.5,color='white',lw=3);ax.axvline(len(times)+(len(contrasts)-len(times))/2-.5,color='white',lw=3)
fig.suptitle('Exploratory metabolic pathway audit — all available pairwise contrasts',fontsize=18)
fig.colorbar(im,ax=axes,shrink=.8,label='Normalized enrichment score (NES)',pad=.025)
fig.supxlabel('Fresh GSEA from existing DESeq2 ranks; curated Xenbase mapping. * FDR < 0.05 across 50 Hallmarks per contrast.\nTemporal columns show later minus earlier. Source time labels used; time origin still requires confirmation.',fontsize=11)
for ext in ['png','pdf']:fig.savefig(F/f'metabolic_NES_audit.{ext}',dpi=180)
plt.close(fig)

fig,axes=plt.subplots(2,2,figsize=(13.5,10),layout='constrained');colors={'CT':'#BA5148','ND':'#3C6D88'}
for ax,tissue in zip(axes[0],['BRAIN','BLAST']):
 d=ss[(ss.tissue==tissue)&(ss.aggregation=='mean')];times=['24HPA','7DPA'] if tissue=='BRAIN' else ['11HPA','24HPA','7DPA']
 for arm in ['ND','CT']:
  means=[];loc=[]
  for i,tm in enumerate(times):
   v=d[(d.treatment==arm)&(d.time==tm)].ISR_exclude_GDF15.to_numpy();x=i+(-.12 if arm=='ND' else .12)
   ax.scatter(x+np.linspace(-.035,.035,len(v)),v,c=colors[arm],s=65,zorder=3);ax.errorbar(x,v.mean(),yerr=v.std(ddof=1)/np.sqrt(len(v)),color=colors[arm],capsize=4);means.append(v.mean());loc.append(x)
  ax.plot(loc,means,color=colors[arm],label='MDT' if arm=='CT' else 'ND')
 ax.set_xticks(range(len(times)),[t.replace('HPA',' h').replace('DPA',' d') for t in times]);ax.set_ylabel('Adapted ISR factor score');ax.legend(frameon=False)
 ax.set_title(('A  Brain: no clear ISR treatment effect' if tissue=='BRAIN' else 'B  Limb: early ISR signature, attenuated by day 7'),fontsize=13,loc='left',fontweight='bold')
 ax.axhline(0,color='#bbbbbb',lw=.8)
paired=read(R/'paired_gdf15_samples.tsv');x=paired.brain_L_vst.to_numpy();y=paired.blast_L_vst.to_numpy()
X=np.column_stack([np.ones(len(paired)),pd.get_dummies(paired.group,drop_first=True,dtype=float)]);rx=rankdata(x)-X@np.linalg.lstsq(X,rankdata(x),rcond=None)[0];ry=rankdata(y)-X@np.linalg.lstsq(X,rankdata(y),rcond=None)[0]
for ax,a,b,ttl in [(axes[1,0],x,y,'C  Pooled GDF15: ρ = 0.83'),(axes[1,1],rx,ry,'D  After four-group adjustment: ρ = 0.26')]:
 for arm in ['ND','CT']:
  for tm,marker in [('24HPA','o'),('7DPA','s')]:
   sel=(paired.treatment==arm)&(paired.time==tm);ax.scatter(a[sel],b[sel],c=colors[arm],marker=marker,s=65)
 xx=np.array([a.min(),a.max()]);ax.plot(xx,np.polyval(np.polyfit(a,b,1),xx),color='#666666',lw=1)
 ax.set_title(ttl,loc='left',fontsize=13,fontweight='bold')
axes[1,0].set_xlabel('Brain gdf15.L, variance-stabilized');axes[1,0].set_ylabel('Limb gdf15.L, variance-stabilized')
axes[1,1].set_xlabel('Brain ranked residual');axes[1,1].set_ylabel('Limb ranked residual')
axes[1,1].text(.04,.94,'Exact within-group permutation p = 0.488\nn = 12 paired animals',transform=axes[1,1].transAxes,va='top',fontsize=11)
fig.suptitle('Columbia follow-up: stress signatures and cross-tissue GDF15',fontsize=18)
fig.supxlabel('ISR scores exclude GDF15; 97/119 genes in brain and 101/119 in limb. Points are animals; bars are SEM.\nHuman-derived score adapted to Xenopus within each tissue; score magnitudes are not directly comparable between tissues.',fontsize=11)
for ext in ['png','pdf']:fig.savefig(F/f'ISR_and_GDF15_audit.{ext}',dpi=180)
plt.close(fig)

provenance={}
for p in [*list((O/'sources').glob('*.*')),*list((O/'inputs').glob('*.tsv')),*[ROOT/f'tables/{f}/voom_E_{t}_primary.tsv' for t,f in [('BRAIN','brain'),('BLAST','blastema')]],ROOT/'tables/sample_metadata.tsv']:
 provenance[str(p.relative_to(ROOT))]={'size':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
(O/'source_hashes.json').write_text(json.dumps(provenance,indent=2));(O/'validation.json').write_text(json.dumps(checks,indent=2))
print(json.dumps(checks,indent=2))
