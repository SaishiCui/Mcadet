
rm(list = ls())
gc()


library(M3Drop)
library(dplyr)
library(Seurat)
library(scry)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("M3Drop")
BiocManager::install("M3DExampleData")
BiocManager::install("scry")


# NBDrop (Negative Binomial drop)

  ## SparSim Coarse default

    for (i in 1:24) {
      load(file = paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_", i, ".RData"))
      raw_data = get(paste0("SparSim_", i)) 
      nb_raw <- NBumiConvertData(raw_data, is.counts=TRUE)
      DANB_fit <- NBumiFitModel(nb_raw)
      NBDropFS<- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thres=0.05, suppress.plot=F)
      assign(paste0("NBDrop_SparSim_",i), NBDropFS$Gene)
      save(list = paste0("NBDrop_SparSim_",i),
           file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_SparSim/NBDrop_SparSim_",
                         i , ".RData"))
      rm(list = paste0("SparSim_", i))
      rm(raw_data)
      gc()
      print(i)
    
    }


    ## SparSim Coarse 3000 genes 

for (i in 1:24) {
  load(file = paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_", i, ".RData"))
  raw_data = get(paste0("SparSim_", i)) 
  nb_raw <- NBumiConvertData(raw_data, is.counts=TRUE)
  DANB_fit <- NBumiFitModel(nb_raw)
  NBDropFS<- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thres=1, suppress.plot=F)
  assign(paste0("NBDrop_SparSim_3000_",i), NBDropFS$Gene[1:3000])
  save(list = paste0("NBDrop_SparSim_3000_",i),
       file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_SparSim/NBDrop_SparSim_3000_",
                     i , ".RData"))
  rm(list = paste0("SparSim_", i))
  rm(raw_data)
  gc()
  print(i)
  
}



    ## SparSim fine default

for (i in 1:24) {
  load(file = paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_fine_", i, ".RData"))
  raw_data = get(paste0("SparSim_fine_", i)) 
  nb_raw <- NBumiConvertData(raw_data, is.counts=TRUE)
  DANB_fit <- NBumiFitModel(nb_raw)
  NBDropFS<- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thres=0.05, suppress.plot=F)
  assign(paste0("NBDrop_SparSim_fine_",i), NBDropFS$Gene)
  save(list = paste0("NBDrop_SparSim_fine_",i),
       file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_SparSim/NBDrop_SparSim_fine_",
                     i , ".RData"))
  rm(list = paste0("SparSim_fine_", i))
  rm(raw_data)
  gc()
  print(i)
  
}




  ## SparSim fine 3000 genes

for (i in 1:24) {
  load(file = paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_fine_", i, ".RData"))
  raw_data = get(paste0("SparSim_fine_", i)) 
  nb_raw <- NBumiConvertData(raw_data, is.counts=TRUE)
  DANB_fit <- NBumiFitModel(nb_raw)
  NBDropFS<- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thres=1, suppress.plot=F)
  assign(paste0("NBDrop_SparSim_fine_3000_",i), NBDropFS$Gene[1:3000])
  save(list = paste0("NBDrop_SparSim_fine_3000_",i),
       file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_SparSim/NBDrop_SparSim_fine_3000_",
                     i , ".RData"))
  rm(list = paste0("SparSim_fine_", i))
  rm(raw_data)
  gc()
  print(i)
  
}








  ## PBMC Coarse default 

for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
  load(file = paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmcdata_",
                     donor.set[i], "_", time.set[j], ".RData"))
  raw_data = get(paste0("pbmcdata_", donor.set[i], "_", time.set[j] )) 
  nb_raw <- NBumiConvertData(raw_data, is.counts=TRUE)
  DANB_fit <- NBumiFitModel(nb_raw)
  NBDropFS<- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thres=0.05, suppress.plot=F)
  assign(paste0("NBDrop_PBMC_", donor.set[i], "_", time.set[j]), NBDropFS$Gene)
  save(list = paste0("NBDrop_PBMC_", donor.set[i], "_", time.set[j]),
       file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_PBMC/NBDrop_PBMC_",
                     donor.set[i], "_", time.set[j] , ".RData"))
  rm(list = paste0("pbmcdata_", donor.set[i], "_", time.set[j]))
  rm(raw_data)
  gc()
  print(c(i,j))
  
 }
}






 ## PBMC Coarse 3000 genes

for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file = paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmcdata_",
                       donor.set[i], "_", time.set[j], ".RData"))
    raw_data = get(paste0("pbmcdata_", donor.set[i], "_", time.set[j] )) 
    nb_raw <- NBumiConvertData(raw_data, is.counts=TRUE)
    DANB_fit <- NBumiFitModel(nb_raw)
    NBDropFS<- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thres=1, suppress.plot=F)
    assign(paste0("NBDrop_PBMC_3000_", donor.set[i], "_", time.set[j]), NBDropFS$Gene[1:3000])
    save(list = paste0("NBDrop_PBMC_3000_", donor.set[i], "_", time.set[j]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_PBMC/NBDrop_PBMC_3000_",
                       donor.set[i], "_", time.set[j] , ".RData"))
    rm(list = paste0("pbmcdata_", donor.set[i], "_", time.set[j]))
    rm(raw_data)
    gc()
    print(c(i,j))
    
  }
}



  ## PMBC fine default 

for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file = paste0("C:/Users/13541/Desktop/Thesis/pbmc/data_fine/pbmc_fine_",
                       donor.set[i], "_", time.set[j], ".RData"))
    raw_data = get(paste0("pbmc_fine_", donor.set[i], "_", time.set[j] )) 
    nb_raw <- NBumiConvertData(raw_data, is.counts=TRUE)
    DANB_fit <- NBumiFitModel(nb_raw)
    NBDropFS<- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thres=0.05, suppress.plot=F)
    assign(paste0("NBDrop_PBMC_fine_", donor.set[i], "_", time.set[j]), NBDropFS$Gene)
    save(list = paste0("NBDrop_PBMC_fine_", donor.set[i], "_", time.set[j]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_PBMC/NBDrop_PBMC_fine_",
                       donor.set[i], "_", time.set[j] , ".RData"))
    rm(list = paste0("pbmc_fine_", donor.set[i], "_", time.set[j]))
    rm(raw_data)
    gc()
    print(c(i,j))
    
  }
}




  ## PMBC fine 3000 genes

for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file = paste0("C:/Users/13541/Desktop/Thesis/pbmc/data_fine/pbmc_fine_",
                       donor.set[i], "_", time.set[j], ".RData"))
    raw_data = get(paste0("pbmc_fine_", donor.set[i], "_", time.set[j] )) 
    nb_raw <- NBumiConvertData(raw_data, is.counts=TRUE)
    DANB_fit <- NBumiFitModel(nb_raw)
    NBDropFS<- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thres= 1, suppress.plot=F)
    assign(paste0("NBDrop_PBMC_fine_3000_", donor.set[i], "_", time.set[j]), NBDropFS$Gene[1:3000])
    save(list = paste0("NBDrop_PBMC_fine_3000_", donor.set[i], "_", time.set[j]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_PBMC/NBDrop_PBMC_fine_3000_",
                       donor.set[i], "_", time.set[j] , ".RData"))
    rm(list = paste0("pbmc_fine_", donor.set[i], "_", time.set[j]))
    rm(raw_data)
    gc()
    print(c(i,j))
    
  }
}















