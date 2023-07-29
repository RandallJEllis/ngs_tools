library(data.table)
library(DWLS)
library(stringr)
library(Seurat)
library(openxlsx)
require(Matrix)
library(biomaRt)

# import single-cell RNA-seq reference dataset
data <-fread("../single_cell_reference/Savell_et_al_2020_saline/sal_male_counts.csv")
data = data.frame(data) # convert to dataframe 
rownames(data)= data[, 1] # set gene names as row names
data = data[,2:ncol(data)] # remove gene name column

# import metadata and replace hyphens with underscores in celltype names
meta=fread('../single_cell_reference/Savell_et_al_2020_saline/sal_male_metadata.csv', header = T)
meta$celltypes = gsub('-', '_', meta$celltypes)
print('checking that order of barcodes in counts file matches order of barcodes in outputs file')

# convert cells are in correct order in counts and metadata
print(all(colnames(data)==meta$barcodes))
 
# Build signature matrix; takes several hours and requires heavy memory
# ~4-5GB per 1000 cells 
# DWLS doesn't support multi-threading so this was run on one core
print('Building signature')
Signature<-buildSignatureMatrixMAST(scdata=data,
                                   id=labels,
                                   path="./DWLS.results",
                                   diff.cutoff=0.5,
                                   pval.cutoff=0.01)

save(Signature,file="DWLS.RData")
load("DWLS.RData")

#import bulk data
counts <- read.csv("../normcounts/icer_oe_2022_normcounts.csv", check.names = F)
rownames(counts) <- make.names(counts[,1], unique = TRUE) #gene symbols
counts <- counts[,2:ncol(counts)]

#If you're using human and mouse or rat, will need to use BioMart to standardize gene symbols
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
values <- rownames(counts)
data <- getBM(

 attributes=c(
   "ensembl_gene_id",
   "mmusculus_homolog_ensembl_gene",
   "mmusculus_homolog_associated_gene_name"),

 filters = "ensembl_gene_id",
 values = values,
 mart= ensembl
)

data <- data[duplicated(data$ensembl_gene_id)==F &
              duplicated(data$mmusculus_homolog_associated_gene_name)==F ,]
rownames(data) <- data$ensembl_gene_id
data <- data[!data$mmusculus_homolog_associated_gene_name=="",]

#normalize counts
counts <- counts[rownames(data),]
CPM <-sweep(counts,2,as.numeric(colSums(counts)),FUN="/")*1000000

# 2-step process of making mouse gene symbols (avoids a formatting error)
rownames(CPM) <- rownames(counts)
counts <- CPM[rownames(Signature),] # filters  counts for genes maintained in the Signature Matrix
samples<- colnames(counts) # grabbing a vector of your subject names

# Initializing Results Data Frames for Each Algorithm that DWLS employs (DWLS being the superstar)
DWLS <- data.frame()
SVR <- data.frame()
OLS <- data.frame()

# Computing Cell Fractions for Each Sample and Adding to Results Data Frames
for(sample in samples){
  bulk <- counts[,sample]
  names(bulk) <- rownames(counts) #subsetting column deletes rownames; add gene symbols back!
  tr<-trimData(Signature,bulk)
  solDWLS<-solveDampenedWLS(tr$sig,tr$bulk)
  solSVR <- solveSVR(tr$sig,tr$bulk)
  solOLS <- solveOLS(tr$sig,tr$bulk)
  DWLS <-rbind(DWLS,solDWLS)
  SVR <- rbind(SVR,solSVR)
  OLS <- rbind(OLS,solOLS)
}

model<- list(DWLS,SVR,OLS) # List of dataframes

# Adding the rownames and columns back to each dataframe in list
L <- lapply(model, function(df)
            {
              colnames(df) <- colnames(Signature)
              rownames(df) <- samples
              df
            }
           )

# Writing each dataframe in list to file
write.csv(L[1],file="CellFractions_Day_Accumbens_ref_DWLS.csv")
write.csv(L[2],file="CellFractions_Day_Accumbens_ref_SVR.csv")
write.csv(L[3],file="CellFractions_Day_Accumbens_ref_OLS.csv")
