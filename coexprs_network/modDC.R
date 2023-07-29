require(DGCA)
library(dplyr)
library(openxlsx)

# import normalized counts
expr=read.csv("../expression_metadata/VST_counts.csv",header=T, check.names = F)

# set gene symbols as rownames
rownames(expr) <- expr[,1]
expr <- expr[,2:ncol(expr)]

# import metadata
meta=read.xlsx("../expression_metadata/metadata.xlsx")
meta=meta[match(colnames(expr)[2:ncol(expr)],meta$subject),]

set.seed(12345)

# import signiificant modules outputted from downstream_analyses.R
sigmod = read.table('significant_module_2column_table.txt', header = T)

design_mat=model.matrix(~0+meta_sub$virus)
head(design_mat)
colnames(design_mat)=c("g1","g2")
moduleDC_res=moduleDC(inputMat = expr_sub,design = design_mat,compare = c("rfp","shisa7"),
                      genes=sigmod$values,labels=sigmod$ind,nPerm=50)
write.table(moduleDC_res,"ds_MDC.xlsx",sep="\t",quote=F,row.names = F)
                      