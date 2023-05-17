
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


# Random

## SparSim 

for (i in 1:24) {
  load(file = paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_", i, ".RData"))
  raw_data = get(paste0("SparSim_", i)) 
  set.seed( (i-1)*4+1 )
  randomFS = rownames(raw_data)[sample(c(1:nrow(raw_data)), 2000)]
  set.seed( (i-1)*4+2 )
  randomFS_3000 = rownames(raw_data)[sample(c(1:nrow(raw_data)), 3000)]
  assign(paste0("random_SparSim_",i), randomFS)
  assign(paste0("random_SparSim_3000_",i), randomFS_3000)
  save(list = paste0("random_SparSim_",i),
       file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_SparSim/random_SparSim_",
                     i , ".RData"))
  save(list = paste0("random_SparSim_3000_",i),
       file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_SparSim/random_SparSim_3000_",
                     i , ".RData"))
  rm(list = paste0("SparSim_", i))
  rm(raw_data)
  gc()
  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_fine_", i, ".RData"))
  raw_data = get(paste0("SparSim_fine_", i)) 
  set.seed((i-1)*4+3)
  randomFS = rownames(raw_data)[sample(c(1:nrow(raw_data)), 2000)]
  set.seed((i-1)*4+4)
  randomFS_3000 = rownames(raw_data)[sample(c(1:nrow(raw_data)), 3000)]
  assign(paste0("random_SparSim_fine_",i), randomFS)
  assign(paste0("random_SparSim_fine_3000_",i), randomFS_3000)
  save(list = paste0("random_SparSim_fine_",i),
       file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_SparSim/random_SparSim_fine_",
                     i , ".RData"))
  save(list = paste0("random_SparSim_fine_3000_",i),
       file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_SparSim/random_SparSim_fine_3000_",
                     i , ".RData"))
  rm(list = paste0("SparSim_fine_", i))
  rm(raw_data)
  gc()
  
  print(i)

}





## PBMC 


for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file = paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmcdata_",
                       donor.set[i], "_", time.set[j], ".RData"))
    raw_data = get(paste0("pbmcdata_", donor.set[i], "_", time.set[j])) 
    set.seed( (i-1)*12+(j-1)*4+1 )
    randomFS = rownames(raw_data)[sample(c(1:nrow(raw_data)), 2000)]
    set.seed( (i-1)*12+(j-1)*4+2 )
    randomFS_3000 = rownames(raw_data)[sample(c(1:nrow(raw_data)), 3000)]
    assign(paste0("random_PBMC_",donor.set[i], "_", time.set[j]), randomFS)
    assign(paste0("random_PBMC_3000_",donor.set[i], "_", time.set[j]), randomFS_3000)
    save(list = paste0("random_PBMC_",donor.set[i], "_", time.set[j]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_pbmc/random_PBMC_",
                       donor.set[i], "_", time.set[j] , ".RData"))
    save(list = paste0("random_PBMC_3000_",donor.set[i], "_", time.set[j]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_pbmc/random_PBMC_3000_",
                       donor.set[i], "_", time.set[j] , ".RData"))
    rm(list = paste0("pbmcdata_", donor.set[i], "_", time.set[j]))
    rm(raw_data)
    gc()
    
    

    
    load(file = paste0("C:/Users/13541/Desktop/Thesis/pbmc/data_fine/pbmc_fine_",
                       donor.set[i], "_", time.set[j], ".RData"))
    raw_data = get(paste0("pbmc_fine_", donor.set[i], "_", time.set[j])) 
    set.seed( (i-1)*12+(j-1)*4+3 )
    randomFS = rownames(raw_data)[sample(c(1:nrow(raw_data)), 2000)]
    set.seed( (i-1)*12+(j-1)*4+4 )
    randomFS_3000 = rownames(raw_data)[sample(c(1:nrow(raw_data)), 3000)]
    assign(paste0("random_PBMC_fine_", donor.set[i], "_", time.set[j]), randomFS)
    assign(paste0("random_PBMC_fine_3000_", donor.set[i], "_", time.set[j]), randomFS_3000)
    save(list = paste0("random_PBMC_fine_", donor.set[i], "_", time.set[j]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_pbmc/random_PBMC_fine_",
                       donor.set[i], "_", time.set[j] , ".RData"))
    save(list = paste0("random_PBMC_fine_3000_",donor.set[i], "_", time.set[j]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_pbmc/random_PBMC_fine_3000_",
                       donor.set[i], "_", time.set[j], ".RData"))
    rm(list = paste0("pbmc_fine_", donor.set[i], "_", time.set[j]))
    rm(raw_data)
    gc()
    
    
    print(c(i,j))
    
    
    
  }
}
















