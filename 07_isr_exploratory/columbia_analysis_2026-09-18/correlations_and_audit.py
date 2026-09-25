"""Independent small-sample checks and abundance audits; outputs only in this analysis folder."""
from pathlib import Path
import itertools,json,hashlib,sys,re
import numpy as np
import pandas as pd
from scipy.stats import rankdata,t as tdist,spearmanr,false_discovery_control
from scipy.cluster.hierarchy import linkage,leaves_list
from scipy.spatial.distance import pdist
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parents[1];OUT=Path(__file__).resolve().parent;RES=OUT/'results';RES.mkdir(exist_ok=True)
meta=pd.read_csv(ROOT/'tables/sample_metadata.tsv',sep='\t').set_index('sample')
def read(p):return pd.read_csv(p,sep='\t')
def design(m,mode='group'):
 if mode=='pooled':return np.ones((len(m),1))
 if mode=='additive':return np.column_stack([np.ones(len(m)),(m.treatment=='CT').astype(float),pd.get_dummies(m.time,drop_first=True,dtype=float)])
 return np.column_stack([np.ones(len(m)),pd.get_dummies(m.group,drop_first=True,dtype=float)])
def partial(x,y,m,mode='group',rank=True):
 x=np.array(x,float);y=np.array(y,float)
 if rank:x=rankdata(x);y=rankdata(y)
 X=design(m,mode);rx=x-X@np.linalg.lstsq(X,x,rcond=None)[0];ry=y-X@np.linalg.lstsq(X,y,rcond=None)[0]
 rho=float(np.corrcoef(rx,ry)[0,1]);df=len(m)-np.linalg.matrix_rank(X)-1
 p=2*tdist.sf(abs(rho)*np.sqrt(df/max(1e-15,1-rho*rho)),df)
 return rho,float(p),int(df),rx,ry
def restricted_perm(rx,ry,m):
 groups=[np.where(m.group.to_numpy()==g)[0] for g in sorted(m.group.unique())]
 poss=[list(itertools.permutations(i)) for i in groups];nperm=int(np.prod([len(x) for x in poss]))
 observed=abs(np.corrcoef(rx,ry)[0,1]);hit=0
 # Exact enumeration: n=12 -> 6^4, n=18 -> 6^6, all feasible.
 for block in itertools.product(*poss):
  idx=np.arange(len(m))
  for g,b in zip(groups,block):idx[g]=b
  r=abs(np.dot(rx,ry[idx])/np.sqrt(np.dot(rx,rx)*np.dot(ry,ry)))
  hit+=r>=observed-1e-12
 return hit/nperm,nperm
def ols_effect(y,X,c):
 b=np.linalg.lstsq(X,y,rcond=None)[0];df=len(y)-np.linalg.matrix_rank(X);rss=np.sum((y-X@b)**2)
 se=np.sqrt(rss/df*(c@np.linalg.inv(X.T@X)@c));eff=c@b;p=2*tdist.sf(abs(eff/se),df)
 return eff,se,float(p),df,eff-tdist.ppf(.975,df)*se,eff+tdist.ppf(.975,df)*se

counts={}; expr={};aud=[];sourcefiles=[]
for tissue,folder in [('BRAIN','brain'),('BLAST','blastema')]:
 ep=ROOT/f'tables/{folder}/voom_E_{tissue}_primary.tsv';cp=ROOT/f'tables/{folder}/persample_gdf15_{tissue}.tsv'
 e=read(ep).set_index('gene');c=read(cp).set_index('sample');counts[tissue]=c;expr[tissue]=e;sourcefiles.extend([ep,cp])
 percent=e.rank(pct=True)*100
 m=meta.loc[e.columns]
 for gene in ['gdf15.L','gdf15.S','gfral.L','gfral.S','ret.L','ret.S']:
  if gene not in e.index:
   aud.append({'tissue':tissue,'gene':gene,'group':'ALL','status':'absent_from_filtered_matrix; not proof of no expression'})
   continue
  for group,md in m.groupby('group'):
   ids=md.index;v=e.loc[gene,ids];cv=c.loc[ids,'normcount_'+gene] if 'normcount_'+gene in c else pd.Series(dtype=float)
   aud.append({'tissue':tissue,'gene':gene,'group':group,'n':len(ids),'status':'present_in_filtered_matrix','mean_logCPM':v.mean(),'geometric_CPM_approx':2**v.mean(),'min_logCPM':v.min(),'max_logCPM':v.max(),'n_logCPM_ge_0':int((v>=0).sum()),'mean_percentile_of_filtered_genes':percent.loc[gene,ids].mean(),'min_percentile':percent.loc[gene,ids].min(),'mean_DESeq2_normcount':cv.mean(),'min_DESeq2_normcount':cv.min(),'max_DESeq2_normcount':cv.max(),'n_normcount_ge_10':int((cv>=10).sum()) if len(cv) else np.nan})
pd.DataFrame(aud).to_csv(RES/'abundance_audit.tsv',sep='\t',index=False)

paired=meta[meta.tissue=='BRAIN'].copy();paired['key']=paired.treatment+'_'+paired.time+'_'+paired.rep
limbmeta=meta[meta.tissue=='BLAST'].copy();limbmeta['key']=limbmeta.treatment+'_'+limbmeta.time+'_'+limbmeta.rep
limbmeta=limbmeta.set_index('key').loc[paired.key]
assert len(paired)==12 and not paired.key.duplicated().any()
paired['brain_L_vst']=counts['BRAIN'].loc[paired.index,'vst_gdf15.L'].to_numpy()
paired['blast_L_vst']=counts['BLAST'].loc[limbmeta['sample'] if 'sample' in limbmeta else meta[meta.tissue=='BLAST'].assign(key=lambda x:x.treatment+'_'+x.time+'_'+x.rep).reset_index().set_index('key').loc[paired.key,'sample'],'vst_gdf15.L'].to_numpy()
limbids=meta[meta.tissue=='BLAST'].reset_index().assign(key=lambda x:x.treatment+'_'+x.time+'_'+x.rep).set_index('key').loc[paired.key,'sample'].tolist()
paired['blast_sample']=limbids
paired['brain_L_logCPM']=expr['BRAIN'].loc['gdf15.L',paired.index].to_numpy()
paired['blast_L_logCPM']=expr['BLAST'].loc['gdf15.L',limbids].to_numpy()
paired['brain_L_logcount']=np.log2(counts['BRAIN'].loc[paired.index,'normcount_gdf15.L'].to_numpy()+1)
paired['blast_L_logcount']=np.log2(counts['BLAST'].loc[limbids,'normcount_gdf15.L'].to_numpy()+1)
paired['blast_LS_logcount']=np.log2(counts['BLAST'].loc[limbids,['normcount_gdf15.L','normcount_gdf15.S']].sum(axis=1).to_numpy()+1)
paired.to_csv(RES/'paired_gdf15_samples.tsv',sep='\t')
corr=[];loo=[];cell=[]
for label,xc,yc in [('L_vst','brain_L_vst','blast_L_vst'),('L_logCPM','brain_L_logCPM','blast_L_logCPM'),('L_logcount','brain_L_logcount','blast_L_logcount'),('blast_homeolog_sum','brain_L_logcount','blast_LS_logcount')]:
 for mode in ['pooled','additive','group']:
  rho,p,df,rx,ry=partial(paired[xc],paired[yc],paired,mode)
  pp,np_=restricted_perm(rx,ry,paired) if mode=='group' else (np.nan,0)
  corr.append({'analysis':'cross_tissue','measure':label,'adjustment':mode,'n':12,'rho':rho,'p':p,'df':df,'restricted_permutation_p':pp,'n_permutations':np_})
  if mode=='group':
   for sample in paired.index:
    ds=paired.drop(sample);rr,pr,d,_,_=partial(ds[xc],ds[yc],ds,mode)
    loo.append({'measure':label,'omitted':sample,'rho':rr,'p':pr,'df':d})
 for group,ds in paired.groupby('group'):
  r,p=spearmanr(ds[xc],ds[yc]);cell.append({'measure':label,'group':group,'n':len(ds),'rho':r,'p_asymptotic_unreliable_n3':p})
pd.DataFrame(corr).to_csv(RES/'cross_tissue_correlations.tsv',sep='\t',index=False)
pd.DataFrame(loo).to_csv(RES/'cross_tissue_leave_one_out.tsv',sep='\t',index=False)
pd.DataFrame(cell).to_csv(RES/'cross_tissue_within_cell.tsv',sep='\t',index=False)

scores=read(RES/'sample_scores.tsv');associations=[];covtests=[]
markers=['ERYTHROID','LEUKOCYTE','ENDOTHELIAL','PERICYTE','MICROGLIAL_ASSOCIATED']
isrs=['ISR_exclude_GDF15','ISR_include_GDF15','ISR_referenceSD_exclude_GDF15','ISR_loading_exclude_GDF15','ISR_CANONICAL']
for tissue in ['BRAIN','BLAST']:
 for agg in ['mean','best']:
  ds=scores[(scores.tissue==tissue)&(scores.aggregation==agg)].set_index('sample')
  gdf=expr[tissue].loc['gdf15.L',ds.index].to_numpy()
  for sig in markers+isrs+['HEME_METABOLISM','UNFOLDED_PROTEIN_RESPONSE','GLYCOLYSIS','OXIDATIVE_PHOSPHORYLATION','HIF_TARGETS']:
   for mode in ['pooled','group']:
    rho,p,df,rx,ry=partial(gdf,ds[sig],ds,mode)
    pp,np_=restricted_perm(rx,ry,ds) if mode=='group' and agg=='mean' and sig=='ISR_exclude_GDF15' else (np.nan,0)
    associations.append({'tissue':tissue,'aggregation':agg,'signature':sig,'adjustment':mode,'n':len(ds),'rho':rho,'p':p,'df':df,'restricted_permutation_p':pp,'n_permutations':np_})
  if agg=='mean':
   X=design(ds);groups=sorted(ds.group.unique());descols=[groups[0]]+groups[1:]
   for sig in ['NONE']+markers:
    XX=X if sig=='NONE' else np.column_stack([X,ds[sig]])
    for tm in sorted(ds.time.unique()):
     c=np.zeros(XX.shape[1]);cg=f'CT_{tm}_{tissue}';ng=f'ND_{tm}_{tissue}'
     if cg!=groups[0]:c[groups.index(cg)]+=1
     if ng!=groups[0]:c[groups.index(ng)]-=1
     eff,se,p,df,lo,hi=ols_effect(gdf,XX,c)
     covtests.append({'tissue':tissue,'covariate':sig,'contrast':f'CT_vs_ND_{tm}','estimate_logCPM':eff,'se':se,'p':p,'df':df,'ci_low':lo,'ci_high':hi})
a=pd.DataFrame(associations);a['q_all']=false_discovery_control(a.p.to_numpy());a.to_csv(RES/'gdf15_signature_correlations.tsv',sep='\t',index=False)
b=pd.DataFrame(covtests);b['q_all']=false_discovery_control(b.p.to_numpy());b.to_csv(RES/'gdf15_covariate_sensitivity.tsv',sep='\t',index=False)

# Recalculate the current Figure 9D column clustering without importing scripts that write figures.
H=ROOT/'supplementary_v9/supplementary_v9_staging/hallmark_50_v9redfix'
print('Cross-tissue:\n',pd.DataFrame(corr).to_string(index=False))
print('ISR associations:\n',a[(a.aggregation=='mean')&a.signature.str.startswith('ISR_')].to_string(index=False))
print('GDF15 marker sensitivity:\n',b.to_string(index=False))
print('Abundance:\n',pd.DataFrame(aud).to_string(index=False))
