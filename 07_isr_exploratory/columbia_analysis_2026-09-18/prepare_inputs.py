"""Prepare frozen, auditable inputs; no manuscript writes."""
from pathlib import Path
import json, re, hashlib, sys
import numpy as np
import pandas as pd
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parents[1]; OUT=Path(__file__).resolve().parent
SRC=OUT/'sources'; INP=OUT/'inputs'; INP.mkdir(exist_ok=True)
load=pd.read_csv(SRC/'isr_repo/ISRGDF15_index/Data/gtex/FacLoads_12_factors.csv',index_col=0)
ref=pd.read_csv(SRC/'isr_training_120.csv',index_col=0)
pub=pd.read_csv(SRC/'published_scores.csv').set_index('SampleID')
ref=ref.loc[load.index,pub.index].T
z=(ref-ref.mean())/ref.std(ddof=1)
sigma=z.corr(); beta=np.linalg.solve(sigma,load.Factor1)
pred=z.to_numpy()@beta
validation={'n_genes':len(load),'n_samples':len(ref),'correlation_published':float(np.corrcoef(pred,pub.ISR_Factor1_Score)[0,1]),'rmse_published':float(np.sqrt(np.mean((pred-pub.ISR_Factor1_Score)**2))),'max_abs_difference':float(np.max(np.abs(pred-pub.ISR_Factor1_Score))),'correlation_matrix_condition':float(np.linalg.cond(sigma)),'loadings_precision':3,'training_samples':'Figure_1F published source data'}
assert validation['correlation_published']>.9999
assert validation['rmse_published']<.02
pd.DataFrame({'sample':ref.index,'published':pub.ISR_Factor1_Score,'reconstructed':pred}).to_csv(INP/'isr_training_validation.tsv',sep='\t',index=False)
sigma.to_csv(INP/'isr_reference_correlation.tsv',sep='\t')
pd.DataFrame({'loading':load.Factor1,'coefficient_full120':beta,'reference_mean':ref.mean(),'reference_sd':ref.std(ddof=1)}).to_csv(INP/'isr_reference_parameters.tsv',sep='\t')
(OUT/'isr_reference_validation.json').write_text(json.dumps(validation,indent=2))
print('REFERENCE',validation)

cols=(SRC/'xenbase_header.txt').read_text().strip().split('\t')
xmap=pd.read_csv(SRC/'GeneGenomeHumanOrtho_pivot_xl.txt',sep='\t',header=None,names=cols,dtype=str).fillna('')
# Use curated one-human-symbol mappings; record every retained/unmapped row.
xmap=xmap[['SYMBOL','ORTHOLOG_HUMAN_SYMBOL']].drop_duplicates()
valid=xmap[(xmap.ORTHOLOG_HUMAN_SYMBOL!='')&~xmap.ORTHOLOG_HUMAN_SYMBOL.str.contains('[,;| ]',regex=True)]
counts=valid.groupby('SYMBOL').ORTHOLOG_HUMAN_SYMBOL.nunique()
valid=valid[valid.SYMBOL.isin(counts[counts==1].index)]
mapping=valid.set_index('SYMBOL').ORTHOLOG_HUMAN_SYMBOL.to_dict()
panels={
 'HIF_TARGETS':'PDK1 LDHA SLC2A1 PGK1 HK2 ENO1 VEGFA EGLN3 BNIP3 CA9 SLC16A3 PFKFB3 ALDOA',
 'ISR_CANONICAL':'ATF4 ATF3 DDIT3 PPP1R15A ASNS TRIB3',
 'ERYTHROID':'ALAS2 GYPC RHAG GATA1 SPTA1 SPTB KLF1 EPB42 AHSP SLC4A1 HBA1 HBA2 HBB',
 'LEUKOCYTE':'PTPRC LCP1 CORO1A SPI1 LAPTM5 CD53 TYROBP ITGB2',
 'ENDOTHELIAL':'VWF PECAM1 CDH5 CLDN5 KDR FLT1 EMCN ESAM RAMP2',
 'PERICYTE':'PDGFRB RGS5 CSPG4 MCAM NOTCH3 ABCC9 KCNJ8 DES',
 'MICROGLIAL_ASSOCIATED':'AIF1 CSF1R CX3CR1 P2RY12 TMEM119 TREM2 TYROBP C1QA C1QB C1QC',
 # This is a descriptive receptor/channel panel, not a validated pain-state signature.
 'NOCICEPTION_RECEPTOR_CHANNEL':'TRPV1 TRPA1 ASIC1 ASIC2 ASIC3 SCN9A SCN10A SCN11A P2RX3 OPRM1 OPRD1 OPRK1 TAC1 TACR1 CALCA CALCRL NGF NTRK1 BDNF NTRK2'
}
(INP/'panels.json').write_text(json.dumps({k:v.split() for k,v in panels.items()},indent=2))
allmapping=[]; cover=[]; coefs=[]
for tissue,folder in [('BRAIN','brain'),('BLAST','blastema')]:
 e=pd.read_csv(ROOT/f'tables/{folder}/voom_E_{tissue}_primary.tsv',sep='\t',index_col=0)
 assert not e.index.duplicated().any() and np.isfinite(e.to_numpy()).all()
 m=pd.DataFrame({'gene':e.index,'human':[mapping.get(g,'') for g in e.index]})
 m['tissue']=tissue; allmapping.append(m)
 mapped=m[m.human!='']; a=e.loc[mapped.gene].copy(); a.index=mapped.human
 avg=a.groupby(level=0).mean()
 bestgenes=e.loc[mapped.gene].mean(axis=1).groupby(mapped.set_index('gene').human).idxmax()
 best=e.loc[bestgenes].copy();best.index=bestgenes.index
 avg.to_csv(INP/f'{tissue}_ortholog_mean.tsv',sep='\t',index_label='gene')
 best.to_csv(INP/f'{tissue}_ortholog_best.tsv',sep='\t',index_label='gene')
 bestgenes.to_csv(INP/f'{tissue}_best_homeolog.tsv',sep='\t',header=['gene'])
 for mode in ['exclude_GDF15','include_GDF15']:
  genes=[g for g in load.index if g in avg.index and (mode=='include_GDF15' or g!='GDF15')]
  coef=np.linalg.solve(sigma.loc[genes,genes],load.loc[genes,'Factor1'])
  for g,b in zip(genes,coef):coefs.append({'tissue':tissue,'mode':mode,'gene':g,'coefficient':b,'loading':load.loc[g,'Factor1'],'reference_sd':ref[g].std(ddof=1)})
  cover.append({'tissue':tissue,'panel':'ISR_'+mode,'available':len(genes),'total':119 if mode=='exclude_GDF15' else 120,'missing':';'.join(g for g in load.index if g not in genes and (mode=='include_GDF15' or g!='GDF15'))})
 for name,genes in panels.items():
  gs=genes.split();cover.append({'tissue':tissue,'panel':name,'available':sum(g in avg.index for g in gs),'total':len(gs),'missing':';'.join(g for g in gs if g not in avg.index)})
 print(tissue,'genes',len(e),'mapped homeologs',len(mapped),'human symbols',len(avg))
pd.concat(allmapping).to_csv(INP/'gene_mapping.tsv',sep='\t',index=False)
pd.DataFrame(cover).to_csv(INP/'coverage.tsv',sep='\t',index=False)
pd.DataFrame(coefs).to_csv(INP/'isr_adapted_coefficients.tsv',sep='\t',index=False)
print(pd.DataFrame(cover).to_string(index=False))
