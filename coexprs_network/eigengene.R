library(WGCNA) ## Eigengene calculation uses WGCNA code
library(openxlsx)
library(gplots)

setwd('')

counts=read.csv("",header=T, check.names = F)
rownames(counts) <- counts[,1]
counts <- counts[,2:ncol(counts)]

meta=read.xlsx("")
# meta=meta[match(colnames(expr)[2:ncol(expr)],meta$subject),]

set.seed(12345)

#subset metadata and counts
meta_sub = meta[meta$group=='g2',]
meta_sub = meta_sub[meta_sub$virus=='g1',]
counts_sub = t(counts[,meta_sub$subject])

#use modules with significant DEG enrichment
module_deg_enrichment = read.xlsx('')
moi = module_deg_enrichment[module_deg_enrichment$P.adj<0.05,]$Input

#sigmod is a table generated using an output file from MEGENA (multiscale_significant.modules.txt)
sigmod=readLines("multiscale_significant.modules.txt")#significant modules#
sigmod=lapply(sigmod,function(x) strsplit(x,"\t")[[1]])
names(sigmod)=do.call(c,lapply(sigmod,function(x) x[1]))
sigmod=lapply(sigmod,function(x) x[-1])
sigmod=stack(sigmod)
sigmod = read.table('significant_module_2column_table.txt', header = T)
sigmod = sigmod[sigmod$ind %in% moi,]

nSamples <- nrow(counts_sub)  ## Number of subjects
moduleColors <- sigmod$ind ## A vector of module ids corresponding to every gene in counts (note counts need to be genes in columns, subjects in rows)
## Will need counts dataframe to feature column for gene corresponding to each module it is in
expression.table <- counts_sub[,as.vector(sigmod$values)] ## counts should already be transposed before this
MEs0 = moduleEigengenes(expression.table, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

#prepare data to correlate with eigengenes
mm = meta_sub[,c("2wk_active_percentage", "day1_percent_active")]
mm$g2 = counts_sub[,'g2']

#correlate eigengenes with data, adjust p-values
moduleTraitCor = cor(MEs, mm, use = "p",method="spearman") ## Probably use spearman correlation to control for outliers, both MEs and metadata uses subjects for rows
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

adjusted.p <- data.frame(matrix(nrow=nrow(moduleTraitCor),ncol=ncol(moduleTraitCor))) ## long winded way of getting FDR-adjusted pvalues
for (i in 1:ncol(moduleTraitPvalue)){
  adjusted.p[,i] <- p.adjust(moduleTraitPvalue[,i],method="BH") ## can choose other methods like Bonferroni
  # adjusted.p <- cbind(adjusted.p,adjust)
}
#set column names and check minimum pvalues before and after correction
colnames(adjusted.p) <- colnames(moduleTraitPvalue)

#get minimum p-values by row to know which modules to remove
apply(moduleTraitPvalue,1,min)
apply(adjusted.p,2,min)

#remove modules with no significant correlations
modules_to_remove = c('MEc1_37','MEc1_204','MEc1_585', 'MEc1_1079')
moduleTraitCor = moduleTraitCor[!rownames(moduleTraitCor) %in% modules_to_remove,]
moduleTraitPvalue = moduleTraitPvalue[!rownames(moduleTraitPvalue) %in% modules_to_remove,]

#get random colors for heatmap
library(RColorBrewer)
n <- nrow(moduleTraitCor)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#get asterisks for heatmap annotation
sig_asterisks <- matrix("",nrow(moduleTraitCor),ncol(moduleTraitCor))
ast_idx = which(moduleTraitPvalue<0.05, arr.ind=T)
for (i in 1:nrow(ast_idx)){
  sig_asterisks[ast_idx[i,1], ast_idx[i,2]] <- '*'
}

ontologies = read.xlsx("../output/MEGENA_mod_GO.xlsx")
ontologies = ontologies[order(ontologies$Pvalue),]
topo = c()
for (module in rownames(moduleTraitCor)){
  c_idx = which(strsplit(module, "")[[1]] == "c")
  module = substr(module, c_idx, 100)
  ont = ontologies[ontologies$Input==module,c(1,2,10)]
  # print(ont[1:5,])
  topo = c(topo, ont$MSigDB[1])
}
write.csv(topo,'../output/eigengene_cor_heatmap_topOntologies.csv')

library(gplots)
color_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

png(paste0('../output/eigenegene_cor_heatmap.png'), unit='in', width=5.9, height=5.21, res=600)

heatmap.2(moduleTraitCor, trace = 'none', margins=c(2,15), RowSideColors = sample(col_vector, n),
          sepcolor="black", col=color_palette, density.info='none',
          colsep=1:ncol(moduleTraitCor), sepwidth=c(0.001,0.001), labRow=topo,
          rowsep=1:nrow(moduleTraitCor), cellnote = sig_asterisks, notecex=3)
dev.off()

png(paste0('../output/eigenegene_cor_heatmap2.png'), unit='in', width=5.9, height=5.21, res=600)

heatmap.2(moduleTraitCor, trace = 'none', margins=c(1,8), RowSideColors = sample(col_vector, n),sepcolor="black",
          colsep=1:ncol(moduleTraitCor), sepwidth=c(0.001,0.001), labRow=topo,
          rowsep=1:nrow(moduleTraitCor), cellnote = sig_asterisks)
dev.off()

