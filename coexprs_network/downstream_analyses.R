library(openxlsx)


expr=read.csv("../expression_metadata/VST_counts.csv",header=T, check.names = F)
rownames(expr) <- expr[,1]
expr <- expr[,2:ncol(expr)]

meta=read.xlsx("../expression_metadata/metadata.xlsx")
meta=meta[match(colnames(expr)[2:ncol(expr)],meta$subject),]

set.seed(12345)

saveto="./output/"
dir.create(saveto)

for (group in c('g1', 'g2', 'g3', 'g4')){
  setwd(group)
  saveto="./output/"
  if (!file.exists('./output')){dir.create(saveto)}
  
  sigmod=readLines("multiscale_significant.modules.txt")#significant modules#
  sigmod=lapply(sigmod,function(x) strsplit(x,"\t")[[1]])
  names(sigmod)=do.call(c,lapply(sigmod,function(x) x[1]))
  sigmod=lapply(sigmod,function(x) x[-1])
  sigmod=stack(sigmod)
  write.table(sigmod,paste(saveto,"significant_module_2column_table.txt",sep=""),sep="\t",quote=F,row.names=F)
  setwd('../../')
}

sigmod = read.table('../output/significant_module_2column_table.txt', header = T)
require(GOtest)
deg=read.csv("../deg/deg_results.csv",header=T)
deg=deg[deg$pvalue<0.05&abs(deg$log2FoldChange)>0,]
deg$Reg=ifelse(deg$log2FoldChange>0,"UP","DN")
deg_mod=GOtest(x=sigmod,deg[,c("row","Reg")],query.population = sigmod$values,
               background = "query",method="hypergeometric")
write.xlsx(deg_mod,"../output/MEGENA_mod_overlap_g1_g2_DEGs.xlsx",
            sep="\t",quote=F,row.names=F)


# ctspec_list=read.delim("E:/celltypesig/cell_type_bretigea_human_500.tsv",header=T,as.is=T)
# ct=GOtest(x=sigmod,ctspec_list,query.population = sigmod$values,background = "annotation",method="hypergeometric")
# write.table(ct,paste(saveto,"Module_overlap_celltype.xls",sep=""),sep="\t",quote=F,row.names=F)

go=msigdb.gsea(sigmod,background = "annotation",method="hypergeometric",species = "mouse")
write.table(go[go$P.adj<0.05,],paste(saveto,"MEGENA_mod_GO.xlsx",sep=""),sep="\t",quote=F,row.names=F)

tf=msigdb.gsea(sigmod,genesets = "c3.tft",background = "annotation",method = "hypergeometric",species = "mouse")
write.table(tf[tf$P.adj<0.05,],paste(saveto,"MEGENA_mod_C3.TFT.xlsx",sep=""),sep="\t",quote=F,row.names=F)

require(DGCA)
library(dplyr)

#preprocessing
meta_sub = meta[meta$virus=='g1',]
expr_sub <- expr %>% select(meta_sub$subject)

sigmod = read.table('../output/significant_module_2column_table.txt', header = T)

design_mat=model.matrix(~0+meta_sub$group)
head(design_mat)
colnames(design_mat)=c("g2","g1")
moduleDC_res=moduleDC(inputMat = expr_sub,design = design_mat,compare = c("g2","g1"),
                      genes=sigmod$values,labels=sigmod$ind,nPerm=50)
write.table(moduleDC_res,"../output/ds_MDC.xlsx",sep="\t",quote=F,row.names = F)
                      

load("MEGENA.Results.RData")

topo=MEGENA.output$module.output$module.relation
head(topo)
topo[,1]=paste0("c1_",topo[,1])
topo[,2]=paste0("c1_",topo[,2])
colnames(topo)=c("mod.parent","module.id")
topo=as.data.frame(topo)


get.topEnrich=function(enrichtable,subject,byfactor){
  s=split.data.frame(enrichtable,as.factor(enrichtable[,match(subject,colnames(enrichtable))]))
  ss=lapply(s,function(x,byfactor) {
    a=as.numeric(as.character(x[,match(byfactor,colnames(x))]))
  index=ifelse(min(a)<0.05,which.min(a),NA);
  sss=do.call(rbind,ss)
  sss=sss[!is.na(sss[,1]),]
  return(sss)
  })}
  

deg_mod1=get.topEnrich(deg_mod,subject = "Input",byfactor = "P.adj")
deg_mod1$log10FDR=ifelse(deg_mod1$Category=="UP",-log10(deg_mod1$P.adj),log10(deg_mod1$P.adj))
topo$DEG=deg_mod1$log10FDR[match(topo[,2],deg_mod1$Input)]

ct1=get.topEnrich(ct,subject = "Input",byfactor = "P.adj")
topo$celltype=ct1$Category[match(topo$module.id,ct1$Input)]

md=moduleDC_res[moduleDC_res$pVal<0.05,]
topo$MDC=md$MeDC[match(topo$module.id,md$Module)]

size=as.data.frame.vector(table(sigmod$ind))
topo$size=size[match(topo$module.id,rownames(size)),1]

go1=go[go$P.adj<0.05,]
go1=get.topEnrich(go1,subject = "Input",byfactor = "P.adj")
topo$topGO=go1$MSigDB[match(topo$module.id,go1$Input)]

tf1=tf[tf$P.adj<0.05,]
tf1=get.topEnrich(tf1,subject = "Input",byfactor = "P.adj")
topo$topTF=tf1$MSigDB[match(topo$module.id,tf1$Input)]

topo[is.na(topo)]=0
topo$Ranking=2*nrow(topo)-rank(abs(topo$DEG),ties.method = "average")-rank(abs(topo$MDC),ties.method = "average")
topo=topo[order(topo$Ranking,decreasing = F),]
topo[topo==0]=NA
write.table(topo,paste0(saveto,"module.summary.txt"),sep="\t",quote=F,row.names = F)



