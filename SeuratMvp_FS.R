
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


# Seurat Mvp

## SparSim Coarse default

for (i in 1:24) {
  load(file = paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_", i, ".RData"))
  raw_data = get(paste0("SparSim_", i)) 
  seurat.obj<-CreateAssayObject(counts = raw_data)
  seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj.mvp <- FindVariableFeatures(seurat.obj, selection.method = "mvp")
  SeuratMvpFS<-VariableFeatures(seurat.obj.mvp)
  assign(paste0("SeuratMvp_SparSim_",i), SeuratMvpFS)
  save(list = paste0("SeuratMvp_SparSim_",i),
       file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratMvp_SparSim_",
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
  seurat.obj<-CreateAssayObject(counts = raw_data)
  seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj.Mvp <- FindVariableFeatures(seurat.obj, selection.method = "mvp")
  SeuratMvpFS<-VariableFeatures(seurat.obj.Mvp)
  assign(paste0("SeuratMvp_SparSim_fine_",i), SeuratMvpFS)
  save(list = paste0("SeuratMvp_SparSim_fine_",i),
       file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratMvp_SparSim_fine_",
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
    seurat.obj<-CreateAssayObject(counts = raw_data)
    seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat.obj.mvp <- FindVariableFeatures(seurat.obj, selection.method = "mvp")
    SeuratMvpFS<-VariableFeatures(seurat.obj.mvp)
    assign(paste0("SeuratMvp_PBMC_", donor.set[i], "_", time.set[j]), SeuratMvpFS)
    save(list = paste0("SeuratMvp_PBMC_", donor.set[i], "_", time.set[j]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_PBMC/SeuratMvp_PBMC_",
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
    seurat.obj<-CreateAssayObject(counts = raw_data)
    seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat.obj.mvp <- FindVariableFeatures(seurat.obj, selection.method = "mvp", nfeatures = 2000)
    SeuratMvpFS<-VariableFeatures(seurat.obj.mvp)
    assign(paste0("SeuratMvp_PBMC_fine_", donor.set[i], "_", time.set[j]), SeuratMvpFS)
    save(list = paste0("SeuratMvp_PBMC_fine_", donor.set[i], "_", time.set[j]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Generating_results/OhterFS/OtherFS_results_PBMC/SeuratMvp_PBMC_fine_",
                       donor.set[i], "_", time.set[j] , ".RData"))
    rm(list = paste0("pbmc_fine_", donor.set[i], "_", time.set[j]))
    rm(raw_data)
    gc()
    print(c(i,j))
    
  }
}





