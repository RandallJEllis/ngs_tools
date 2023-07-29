library(CARseq)
library(data.table)
require(Matrix)
library(openxlsx)

counts = read.csv('../rawcounts/rawcounts.csv', check.names=F)
rownames(counts) <- counts[,1]
counts = counts[,2:ncol(counts)]

meta = read.xlsx('../metadata.xlsx')
meta$condition = paste0(meta$ID,'_',meta$Group)

cellfrac = read.csv('../CellFractions_SVR.csv')

#Remove cells with low fractions as recommended by the developers: https://github.com/Sun-lab/CARseq/issues/3
#According to them, keeping in cells with low fractions can lead to convergence issues and many genes
#having NAs in their results
print(paste('Starting number of cell types:', dim(cellfrac)[2]-1))

#make sure you remove the ID column as this will make apply treat all values as characters instead of numerics
max_val_all_celltypes = apply(cellfrac[,2:ncol(cellfrac)], 2, max)

#remove cell types with max values below a threshold
threshold = 0.03 
print(threshold)
cells_equalTo_orAbove_threshold = cellfrac[,2:ncol(cellfrac)][, which(max_val_all_celltypes >= threshold)]
cellfrac = cbind(cellfrac$X, cells_equalTo_orAbove_threshold)
colnames(cellfrac)[1] = 'X'
print(paste('After cell type removal:', dim(cellfrac)[2]-1))
cellfrac = cellfrac[order(match(cellfrac$X, colnames(counts))),]
print(all(colnames(counts)==cellfrac$X))

# set groups
g1 = df_comp[i,1]
g2 = df_comp[i,2]
print(g1)
print(g2)

# subset metadata and expression data
print(all(colnames(counts)==cellfrac$X))
print('Confirmed samples are in same order in counts and cell fractions.')

groups=c()
for (sample in colnames(counts)){
  groups = c(groups, meta[meta$condition==sample, 'Group'])
}
print('Groups made.')

res = run_CARseq(counts, cellfrac[,2:ncol(cellfrac)], groups)
save(res, file ="carseq_rawcounts_all_genes.RData")
