library(data.table)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(tximport)
library(openxlsx)
library(edgeR)

meta = read.xlsx('../metadata.xlsx')
meta_sub = meta[meta$group=='1',]

files = file.path('.',list.files('.'))
colnames = c()
files_import = c()
for (file in files) {
  cutindex = unlist(gregexpr('Aligned',file))
  id = substr(file,3,(cutindex-1))
  if (id %in% meta_sub$subject){
    colnames= c(colnames,id)
    files_import = c(files_import, file)
  }
  
}

names(files_import) = colnames
txi.rsem=tximport(files_import,type='rsem', txIn=T, txOut=T)

#DESEQ
# meta$condition = paste0(meta$ID,"_", meta$Group)
meta_sub=meta_sub[match(colnames(txi.rsem$counts),meta_sub$subject),]
meta_sub$batch <- factor(meta_sub$batch)
meta_sub$virus <- factor(meta_sub$virus)

#NEXT LINE MUST OUTPUT TRUE!!!!!
print(all(meta_sub$subject == colnames(txi.rsem$counts)))
# meta$Group=factor(meta$Group)

dds = DESeqDataSetFromTximport(txi.rsem,meta_sub,~batch+virus)
#filter low-expressed genes
keep.exprs <- filterByExpr(counts(dds), group=meta_sub$virus, min.count=10)
# keep <- rowSums(counts(dds)) >= 100 #subset genes with summed counts equal to or greater than 10 
dds <- dds[keep.exprs,]
dds$virus <- relevel(dds$virus, ref = '1') #set reference group

#Run DESeq2, get results
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, tidy=T, name = 'virus_1_vs_2')

#export results
genenames = read.table('1_2Aligned.isoforms.results',header=T)
newresults = merge(res,genenames,by.x='row',by.y='transcript_id')
newresults = newresults[,c(1,8,2,3,4,5,6,7)]
write.csv(newresults, paste0('../../deg/transcript_level/deg_results.csv'))
