library(DESeq2)
library(openxlsx)
library(dplyr)
library(edgeR)
library(EnhancedVolcano)

#for volcano plotting
pval_col = 'padj' # 'pvalue' or 'padj'
threshold = 0.05

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load data
data = read.csv('../counts/raw_counts.csv', check.names = F)
#remove Undetermined sample
data = subset(data, select=-c(Undetermined))

#load metadata
meta = read.csv('../metadata/metadata_cleaned.csv')

#reorder samples in meta to be the same as data
cts = subset(data, select=-c(gene))
rownames(meta) <- meta$SampleID
meta <- meta[match(colnames(cts), rownames(meta)),]
#sanity checks for ordering; both should output TRUE
all(rownames(meta) == colnames(cts)) 
all(colnames(data[,2:ncol(data)]) == colnames(cts)) 

#remove low-expressed genes to increase statistical power and
#to help with ensembl/symbol merging (some Ensembl IDs have 
#multiple gene symbols)
keep <- filterByExpr(data[,2:ncol(data)], group=meta$Description)
data <- data[keep,]

### Annotate Ensembl IDs with gene symbols
# load the biomaRt package
library(biomaRt)
# specify the Ensembl dataset and attributes you want to retrieve
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id", "external_gene_name", "description")
# create a vector of Ensembl gene IDs you want to annotate
ensembl_ids <- data$gene
# retrieve the gene symbols using the biomaRt package
annotated <- getBM(attributes = attributes, filters = "ensembl_gene_id",
                   values = ensembl_ids, mart = ensembl)

#make rownames the gene symbols and remove the gene symbol column
exprs <- data.frame(data, check.names = F)
rownames(exprs) = exprs$gene
exprs = exprs[, 2:ncol(exprs)]

#make sure variables are categorical
meta$Description = factor(meta$Description)

#create DESeq2 object
dds <- DESeqDataSetFromMatrix(exprs, colData = meta, design = ~Description)
dds$Description <- relevel(dds$Description, ref = 'group1') #set reference group

#Run DESeq2, get results
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

#set up different comparisons
ref_column <- c("group1", "group1", "group2")
tx_column <- c("group2", "group3", "group3")
df_comp = data.frame(ref_column, tx_column)
colnames(df_comp) <- c('ref','tx')

# create a workbook object
wb <- createWorkbook()

#iterate over comparisons, save to workbook
for(row in seq(dim(df_comp)[1])){
  print(row)
  #put reference group second
  comparison = c('Description',df_comp[row,'tx'],df_comp[row,'ref'])
  res <- results(dds, comparison ,tidy=T)
  res <- res[order(res$pvalue),]
  
  res <- merge(res, annotated, all.x=T, by.x = "row", by.y = "ensembl_gene_id",
                 sort=F)
  # Move symbol column to index 2
  res <- res %>% relocate(external_gene_name, .after=row)
  
  #make worksheet and save data
  #put ref group first for clarity
  sheetname = paste0(df_comp[row,'ref'],'_vs_',df_comp[row,'tx'])
  addWorksheet(wb, sheetname)
  writeData(wb, sheetname, res)
  
  ## Volcano plot
  #create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
  #this can be achieved with nested ifelse statements
  if (pval_col=='padj'){
    res <- res[complete.cases(res[, 'padj']), ]
  }
  
  keyvals <- ifelse(
    ((res$log2FoldChange < 0) & (res[pval_col]<threshold)), 'royalblue',
    ifelse(((res$log2FoldChange > 0) & (res[pval_col]<threshold)), 'red3',
           'black'))
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'red3'] <- 'high'
  names(keyvals)[keyvals == 'black'] <- 'mid'
  names(keyvals)[keyvals == 'royalblue'] <- 'low'
  
  # some adjusted p-values may be 0. To make them able to be plotted,
  # set these one decimal place over from the smallest non-zero p-value
  # Ex. if three genes have p-values of 0, 0, 1e-100, set the first two to
  # 1e-101
  padj_col = res$padj[complete.cases(res$padj)]
  if (min(padj_col) == 0){
    smallest_nonzero_p = min(padj_col[padj_col > 0])
    overwrite = smallest_nonzero_p / 10
    # res$padj[res$padj == 0] <- overwrite
  }

  p1 = EnhancedVolcano(res,
                  lab = res$external_gene_name,
                  x = 'log2FoldChange',
                  y = pval_col,
                  selectLab = res$external_gene_name[which(names(keyvals) %in% c('high', 'low'))],
                  xlab = bquote(~Log[2]~ 'fold change'),
                  title = paste(df_comp[row,'ref'],'vs.',df_comp[row,'tx']),
                  # pCutoff = 0.0001,
                  gridlines.minor = FALSE,
                  pointSize = 1,
                  labSize = 4.5,
                  colCustom = keyvals,
                  colAlpha = 1,
                  borderWidth = 1.5,
                  cutoffLineType = 'blank',
                  drawConnectors = TRUE,
                  borderColour = 'black',
                  legendPosition = 'none',
                  caption = '',
                  max.overlaps = 3
  ) +
    #for x-axis (log2FC), add 0.1 to either side
    # xlim(min(res$log2FoldChange)-0.1,max(res$log2FoldChange)+0.1) +
    #for y-axis (-log10(pvalue)), add 0.2 on top
     ylim(-0.25,max(-log10(smallest_nonzero_p))+0.2)
  
  ggsave(paste0('../volcano_plots/',pval_col,'_',threshold,'/',df_comp[row,'ref'],'_vs_',df_comp[row,'tx'],".png"), p1, width = 8, height = 6, dpi = 300)
  
}

# save the workbook to a file
saveWorkbook(wb, "deg_results.xlsx", overwrite = T)

# VST normalization
library(tibble)

vsd <- vst(dds, blind=FALSE)
vst_cnt <- as.data.frame(assay(vsd))
vst_cnt <- rownames_to_column(vst_cnt, "ensembl")

vst_cnt <- merge(vst_cnt, annotated, all.x=T, by.x = "ensembl", by.y = "ensembl_gene_id",
             sort=F)

# Move symbol column to index 2
vst_cnt <- vst_cnt %>% relocate(external_gene_name, .after=ensembl)
write.csv(vst_cnt, file = '../counts/filtered_VST_normalized_counts.csv')
