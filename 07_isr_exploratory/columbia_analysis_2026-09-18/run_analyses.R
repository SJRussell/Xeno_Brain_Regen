options(stringsAsFactors=FALSE,warn=1)
out <- 'columbia_analysis_2026-09-18'; inp <- file.path(out,'inputs'); res <- file.path(out,'results')
dir.create(res,showWarnings=FALSE)
.libPaths(c(normalizePath(file.path(out,'R_library'),winslash='/'),.libPaths()))
suppressPackageStartupMessages({library(limma);library(fgsea);library(jsonlite);library(data.table)})
set.seed(20260918)
readtab <- function(p) read.delim(p,check.names=FALSE)
write <- function(x,n) write.table(x,file.path(res,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
meta <- readtab('tables/sample_metadata.tsv')
mapping <- readtab(file.path(inp,'gene_mapping.tsv'))
coef <- readtab(file.path(inp,'isr_adapted_coefficients.tsv'))
panels <- fromJSON(file.path(inp,'panels.json'))
hall <- gmtPathways(file.path(out,'sources/hallmark_2024.gmt'))
short <- sub('HALLMARK_','',names(hall)); names(hall)<-short
metabolic <- c('OXIDATIVE_PHOSPHORYLATION','GLYCOLYSIS','HYPOXIA','MTORC1_SIGNALING','REACTIVE_OXYGEN_SPECIES_PATHWAY','FATTY_ACID_METABOLISM','HEME_METABOLISM','UNFOLDED_PROTEIN_RESPONSE')
targets <- unique(c('GDF15','GFRAL','RET','HIF1A','EGLN1','EGLN2','EGLN3','VHL','HK2',unlist(panels)))
zmat <- function(e) {zz<-t(scale(t(e)));zz[!is.finite(zz)]<-0;zz}
makedesign <- function(m) {
 group<-paste(m$treatment,m$time,sep='_'); d<-model.matrix(~0+factor(group));colnames(d)<-levels(factor(group))
 times<-intersect(c('11HPA','24HPA','7DPA'),unique(m$time)); ct<-list()
 add<-function(name,terms){v<-setNames(rep(0,ncol(d)),colnames(d));for(k in names(terms))v[k]<-terms[k];ct[[name]]<<-v}
 for(tm in times) add(paste0('CT_vs_ND_',tm),setNames(c(1,-1),paste(c('CT','ND'),tm,sep='_')))
 for(pair in combn(times,2,simplify=FALSE)){
  for(tr in c('CT','ND'))add(paste(tr,pair[2],'vs',pair[1],sep='_'),setNames(c(1,-1),paste(tr,pair[c(2,1)],sep='_')))
  add(paste('INTERACTION',pair[2],'vs',pair[1],sep='_'),setNames(c(1,-1,-1,1),c(paste0('CT_',pair[2]),paste0('ND_',pair[2]),paste0('CT_',pair[1]),paste0('ND_',pair[1]))))
 }
 list(design=d,contrasts=do.call(cbind,ct))
}
scoretests<-function(scores,des,tissue,aggregation){
 # Ordinary sample-level linear models: no empirical shrinkage across arbitrary signatures.
 fit<-lmFit(t(scores),des$design); f<-contrasts.fit(fit,des$contrasts)
 rows<-list()
 for(i in seq_len(ncol(f$coefficients))){
  b<-f$coefficients[,i];se<-f$stdev.unscaled[,i]*f$sigma;df<-f$df.residual
  p<-2*pt(-abs(b/se),df);ci<-qt(.975,df)*se
  rows[[i]]<-data.frame(tissue,aggregation,signature=rownames(f$coefficients),contrast=colnames(f$coefficients)[i],estimate=b,se,df,ci_low=b-ci,ci_high=b+ci,p)
 }
 do.call(rbind,rows)
}
alltests<-list();allscores<-list();allcamera<-list();geneinter<-list();coverage<-list()
for(tissue in c('BRAIN','BLAST')) {
 e0<-readtab(file.path('tables',if(tissue=='BRAIN')'brain' else 'blastema',paste0('voom_E_',tissue,'_primary.tsv')))
 original<-as.matrix(e0[,-1]);rownames(original)<-e0[[1]]
 m<-meta[match(colnames(original),meta$sample),];stopifnot(!anyNA(m$sample),all(m$tissue==tissue))
 des<-makedesign(m);write(data.frame(sample=m$sample,des$design),paste0(tissue,'_design.tsv'));write(data.frame(group=rownames(des$contrasts),des$contrasts),paste0(tissue,'_contrasts.tsv'))
 for(aggregation in c('mean','best')){
  ed<-readtab(file.path(inp,paste0(tissue,'_ortholog_',aggregation,'.tsv')));e<-as.matrix(ed[,-1]);rownames(e)<-ed[[1]]
  stopifnot(identical(colnames(e),m$sample));z<-zmat(e)
  sets<-c(hall,panels);scores<-list()
  for(nm in names(sets)){
   gs<-intersect(sets[[nm]],rownames(e));coverage[[length(coverage)+1]]<-data.frame(tissue,aggregation,signature=nm,n=length(gs),total=length(sets[[nm]]),genes=paste(gs,collapse=';'))
   if(length(gs)>=3)scores[[nm]]<-colMeans(z[gs,,drop=FALSE])
  }
  scores[['GLYCOLYSIS_MINUS_OXPHOS']]<-scores[['GLYCOLYSIS']]-scores[['OXIDATIVE_PHOSPHORYLATION']]
  for(mode in c('exclude_GDF15','include_GDF15')){
   cf<-coef[coef$tissue==tissue&coef$mode==mode,];gs<-cf$gene
   # Missing-gene regression coefficients derive from the frozen human reference correlation.
   # Within-tissue Xenopus scaling is primary; human-SD-only scaling is a sensitivity.
   scores[[paste0('ISR_',mode)]]<-as.vector(crossprod(cf$coefficient,z[gs,,drop=FALSE]))
   scores[[paste0('ISR_referenceSD_',mode)]]<-as.vector(crossprod(cf$coefficient,sweep(e[gs,,drop=FALSE],1,cf$reference_sd,'/')))
   scores[[paste0('ISR_loading_',mode)]]<-as.vector(crossprod(cf$loading,z[gs,,drop=FALSE]))/sum(abs(cf$loading))
  }
  sm<-do.call(cbind,scores);rownames(sm)<-m$sample
  allscores[[length(allscores)+1]]<-data.frame(m,aggregation,sm)
  st<-scoretests(sm,des,tissue,aggregation)
  st$family<-ifelse(st$signature %in% c(metabolic,'GLYCOLYSIS_MINUS_OXPHOS'),'metabolic',ifelse(startsWith(st$signature,'ISR_'),'ISR','other_Hallmarks_or_panels'))
  alltests[[length(alltests)+1]]<-st
  if(aggregation=='mean'){
   ix<-lapply(sets,function(g)which(rownames(e)%in%g));ix<-ix[lengths(ix)>=3]
   for(cn in colnames(des$contrasts)){
    cm<-camera(e,index=ix,design=des$design,contrast=des$contrasts[,cn],inter.gene.cor=NA,allow.neg.cor=FALSE,trend.var=TRUE)
    allcamera[[length(allcamera)+1]]<-data.frame(tissue,contrast=cn,signature=rownames(cm),cm,check.names=FALSE)
   }
  }
 }
 # Exploratory genome-wide interaction test on available log-CPM; original precision weights unavailable.
 f<-eBayes(contrasts.fit(lmFit(original,des$design),des$contrasts),trend=TRUE)
 for(cn in grep('INTERACTION',colnames(des$contrasts),value=TRUE)){
  tt<-topTable(f,coef=cn,number=Inf,sort.by='none',confint=TRUE);tt$gene<-rownames(tt)
  mp<-mapping[mapping$tissue==tissue,];tt$human<-mp$human[match(tt$gene,mp$gene)]
  write(tt,paste0(tissue,'_',cn,'_unweighted_limma_allgenes.tsv'))
  keep<-tt$human%in%targets|toupper(sub('\\.[LS]$','',tt$gene))%in%targets
  geneinter[[length(geneinter)+1]]<-data.frame(tissue,contrast=cn,tt[keep,])
 }
 message('Completed sample signatures, camera and gene interactions: ',tissue)
}
st<-do.call(rbind,alltests)
st$q_family<-ave(st$p,interaction(st$aggregation,st$family),FUN=function(p)p.adjust(p,'BH'))
st$q_all<-ave(st$p,st$aggregation,FUN=function(p)p.adjust(p,'BH'))
write(st,'signature_tests.tsv');write(do.call(rbind,allscores),'sample_scores.tsv');write(do.call(rbind,coverage),'signature_gene_coverage.tsv')
cm<-do.call(rbind,allcamera);cm$q_all<-p.adjust(cm$PValue,'BH');write(cm,'camera_tests.tsv');write(do.call(rbind,geneinter),'target_gene_interactions.tsv')

# Re-rank existing DESeq2 tables: same genes/statistics, explicit Xenbase mapping and deterministic collapse.
depath<-'figure_rebuild_2026-07/nibi_export/fig_export'
gsea<-list();genes<-list()
for(file in list.files(depath,pattern='^de_.*\\.tsv$',full.names=TRUE)){
 stem<-sub('^de_|\\.tsv$','',basename(file));stem<-sub('\\.tsv$','',stem)
 tissue<-if(grepl('BRAIN',stem))'BRAIN' else 'BLAST';contrast<-sub(paste0('^',tissue,'_'),'',stem)
 d<-readtab(file);mp<-mapping[mapping$tissue==tissue,];d$human<-mp$human[match(d$gene,mp$gene)]
 keep<-d$human%in%targets|toupper(sub('\\.[LS]$','',d$gene))%in%targets
 genes[[length(genes)+1]]<-data.frame(tissue,contrast,d[keep,])
 for(mapmethod in c('curated_Xenbase','symbol_strip')){
  ds<-d
  if(mapmethod=='symbol_strip')ds$human<-toupper(sub('\\.[LS]$','',ds$gene))
  ds<-ds[!is.na(ds$pvalue)&!is.na(ds$log2FoldChange)&!is.na(ds$human)&ds$human!='',]
  ds$rank<-sign(ds$log2FoldChange)*-log10(pmax(ds$pvalue,1e-300))
  ds<-ds[order(-abs(ds$rank),ds$gene),];ds<-ds[!duplicated(ds$human),];r<-setNames(ds$rank,ds$human);r<-sort(r,decreasing=TRUE)
  set.seed(20260918+length(gsea));fg<-fgseaMultilevel(hall,r,minSize=10,maxSize=500,eps=1e-10,nproc=1,nPermSimple=10000)
  fg$leadingEdge<-vapply(fg$leadingEdge,paste,collapse=';',FUN.VALUE='')
  gsea[[length(gsea)+1]]<-data.frame(tissue,contrast,mapping=mapmethod,fg)
 }
 message('Completed Hallmark re-ranking: ',stem)
}
gs<-do.call(rbind,gsea);gs$q_all<-ave(gs$pval,gs$mapping,FUN=function(p)p.adjust(p,'BH'));write(gs,'hallmark_fgsea_all.tsv');write(do.call(rbind,genes),'target_genes_original_DESeq2.tsv')
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
message('ALL R ANALYSES COMPLETE')
