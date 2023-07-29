library(MEGENA)
library(Matrix)
library(openxlsx)
library(matrixStats)

options(warn=1)

# input parameters
n.cores <- 6; # number of cores/threads to call for PCP
doPar <-TRUE; # do we want to parallelize?
method = "pearson" # method for correlation. either pearson or spearman.
FDR.cutoff = 2 # FDR threshold to define significant correlations upon shuffling samples.
module.pval = 0.05 # module significance p-value. Recommended is 0.05.
hub.pval = 0.05 # connectivity significance p-value based random tetrahedral networks
cor.perm = 1; # number of permutations for calculating FDRs for all correlation pairs.
hub.perm = 100; # number of permutations for calculating connectivity significance p-value.
min.size = 50
max.size = NULL

# annotation to be done on the downstream
annot.table=NULL
id.col = 1
symbol.col= 2
###########

# import data
ds=read.csv("normcounts.csv",header=T, check.names=F)
rownames(ds) = ds[,1]
ds = ds[,2:ncol(ds)]

#optional: coefficient of variation filtering
# sd_rows <- apply(ds, 1, sd)
# mean_rows <- apply(ds, 1, mean)
# cv <- sd_rows / mean_rows * 100
# sort_cv = sort(cv, decreasing = T)
# threshold = sort_cv[as.integer(0.20 * length(sort_cv))]
# ds = ds[which(cv>=threshold),]

# #optional: retain the top 20% most variably expressed genes
# genevars = rowVars(as.matrix(ds))
# sort_genevars = sort(genevars, decreasing = T)
# threshold = sort_genevars[as.integer(0.20 * length(sort_genevars))]
# ds = ds[which(genevars>=threshold),]
# print(dim(ds))

meta=read.xlsx("metadata.xlsx")
meta$condition = paste0(meta$ID, '_', meta$Group)
rownames(meta)=meta$condition


for(group in c("g1","g2","g3")){
        
    meta_subset = meta[meta$Group==group,]
    datExpr=ds[,meta_subset$condition]
    print(dim(datExpr))
    
    #remove genes with stddev=0
    sd_rows <- apply(datExpr, 1, sd)
    datExpr = datExpr[which(sd_rows>0),]
    print(dim(datExpr))
    
    print(Sys.time())
    ####
    # calculate correlation
    ijw <- calculate.correlation(datExpr = datExpr, doPerm = cor.perm, method = method, FDR.cutoff = FDR.cutoff, saveto=".")
    saveRDS(ijw, file='ijw.RData')
    print(dim(ijw))
    print('Finished calculate.correlation')
    print(Sys.time())
    
    #### register multiple cores if needed: note that set.parallel.backend() is deprecated.
    run.par = doPar & (getDoParWorkers() == 1)
    if (run.par)
    {
      cl <- parallel::makeCluster(n.cores)
      registerDoParallel(cl)
      # check how many workers are there
      cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
    }
    
    ##### calculate PFN
    el <- calculate.PFN(ijw[,1:3],doPar = doPar,num.cores = n.cores)
    print('Finished calculate.PFN')
    write.csv(el, file="weights.for.cytoscape.csv")
    saveRDS(el, file="pfn.RData")
    rm(ijw)
    print(Sys.time())
    
    g <- graph.data.frame(el,directed = FALSE)
    saveRDS(g, "graph.RData")
    rm(el)
    
    ##### perform MCA clustering.
    MEGENA.output <- do.MEGENA(g,
                               mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = TRUE,
                               min.size = min.size,max.size = vcount(g)/2,
                               doPar = doPar,num.cores = n.cores,n.perm = hub.perm,
                               save.output = TRUE)
    print('Finished do.MEGENA')
    print(Sys.time())
    
    ###### unregister cores as these are not needed anymore.
    if (getDoParWorkers() > 1)
    {
      env <- foreach:::.foreachGlobals
      rm(list=ls(name=env), pos=env)
    }
    
    summary.output <- MEGENA.ModuleSummary(MEGENA.output,
                                           mod.pvalue = module.pval, hub.pvalue = hub.pval,
                                           min.size = min.size, max.size = vcount(g)/2,
                                           annot.table = annot.table, id.col = id.col, symbol.col = symbol.col,
                                           output.sig = TRUE)
    print('Finished MEGENA.ModuleSummary')
    print(Sys.time())

    module.output <- module_convert_to_table(MEGENA.output, mod.pval = 0.05,
                                             hub.pval = 0.05, min.size=min.size, max.size=vcount(g)/2)
    print('Finished module_convert_to_table')
    save(summary.output,MEGENA.output,module.output,g,file="MEGENA.Results.RData")
    
    rm(list = c("summary.output","module.output"))
    
    

    # prep for parallelization prior to the loop
    parallelize = TRUE
    if (foreach::getDoParWorkers() == 1 & parallelize)
    {
      cl = parallel::makeCluster(n.cores)
      doParallel::registerDoParallel(cl)
      # check how many workers are there
      cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))
    }

    # calculate PFN
    el <- calculate.PFN(ijw[,1:3],doPar = doPar,num.cores = n.cores)
    write.table(el,file = "MEGENA_Network.txt",sep = "\t",row.names = F,col.names = T,quote = F)
    rm(ijw)

    # do clustering
    g <- graph.data.frame(el,directed = F)

    MEGENA.output <- do.MEGENA(g,
                               mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = TRUE,
                               min.size = 50,max.size = vcount(g)/2,
                               doPar = TRUE,num.cores = n.cores,n.perm = 100,
                               save.output = T)

    save(MEGENA.output,file = "MEGENA_output.RData")
    setwd("../../")
    
  }


