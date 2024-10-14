rm(list=ls())
gc()


library(M3Drop)
library(dplyr)
library(Seurat)
library(scry)

for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/DE_genes/DEgene_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/DE_genes/DEgene_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    
  }
}

jaccard = function(a,b){
  numerator = length(intersect(a,b))
  denominator = length(a) + length(b) - numerator
  return( numerator/denominator)
}



### SparSim Fine resolution 

for (i in 1:24) {
  
  load(file = 
         paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_fine_", 
                i, ".RData") )
  workdata = get(paste0("SparSim_fine_", i))
  
  set.seed(i)
  newdata1 = matrix(rbinom(nrow(workdata)*ncol(workdata), as.vector(as.matrix(workdata)), 0.5), nrow=nrow(workdata))
  newdata2 = workdata - newdata1
  colnames(newdata1) = colnames(workdata)
  rownames(newdata1) = rownames(workdata)
  newdata1 = as.data.frame(newdata1)
  colnames(newdata2) = colnames(workdata)
  rownames(newdata2) = rownames(workdata)
  print("Split complete")
  
  

  
  result1 = mcadet(data = newdata1,
                   n.comp= 60, run=10, n.feature=NA, nk_percent = 0.010, 
                   start_resolution = 0.5,  cell_percent =  0.005,
                   MC_iter = 50000, fdr = 0.15, seed = 1234)
  
  
  result2 = mcadet(data = newdata2,
                   n.comp= 60, run=10, n.feature=NA, nk_percent = 0.010, 
                   start_resolution = 0.5,  cell_percent =  0.005,
                   MC_iter = 50000, fdr = 0.15, seed = 1234)
  
  assign(paste0("Mcadet_SparSim_fine_stable1_", i),
         result1)
  
  assign(paste0("Mcadet_SparSim_fine_stable2_", i),
         result2)
  
  print("Mcadet Complete")
  
  
  
  result_Bre1 = BrenneckeGetVariableGenes(newdata1, suppress.plot = T, fdr=0.1, minBiolDisp=0.5)
  result_Bre2 = BrenneckeGetVariableGenes(newdata2, suppress.plot = T, fdr=0.1, minBiolDisp=0.5)
  
  assign(paste0("BreFS_SparSim_fine_stable1_", i), result_Bre1$Gene)
  assign(paste0("BreFS_SparSim_fine_stable2_", i), result_Bre2$Gene)
  
  gc()
  
  
  
  M3_data1 <- M3DropConvertData(newdata1, is.counts=TRUE)
  result_M3Drop1<- M3DropFeatureSelection(M3_data1 ,mt_method="fdr", mt_threshold=0.05)
  M3_data2 <- M3DropConvertData(newdata2, is.counts=TRUE)
  result_M3Drop2<- M3DropFeatureSelection(M3_data2 ,mt_method="fdr", mt_threshold=0.05)
  
  assign(paste0("M3Drop_SparSim_fine_stable1_", i), result_M3Drop1$Gene)
  assign(paste0("M3Drop_SparSim_fine_stable2_", i), result_M3Drop2$Gene)
  gc()
  print("M3Drop Complete")
  
  nb_raw1 <- NBumiConvertData(newdata1, is.counts=TRUE)
  DANB_fit1 <- NBumiFitModel(nb_raw1)
  result_NBDrop1<- NBumiFeatureSelectionCombinedDrop(DANB_fit1, method="fdr", qval.thres=0.05, suppress.plot=F)
  nb_raw2 <- NBumiConvertData(newdata2, is.counts=TRUE)
  DANB_fit2 <- NBumiFitModel(nb_raw2)
  result_NBDrop2<- NBumiFeatureSelectionCombinedDrop(DANB_fit2, method="fdr", qval.thres=0.05, suppress.plot=F)
  
  
  assign(paste0("NBDrop_SparSim_fine_stable1_",i), result_NBDrop1$Gene)
  assign(paste0("NBDrop_SparSim_fine_stable2_",i), result_NBDrop2$Gene)
  gc()
  
  
  
  set.seed( i )
  result_random1 = rownames(newdata1)[sample(c(1:nrow(newdata1)), 2000)]
  set.seed( i+1 )
  result_random2 = rownames(newdata2)[sample(c(1:nrow(newdata2)), 2000)]
  
  
  assign(paste0("random_SparSim_fine_stable1_", i), result_random1)
  assign(paste0("random_SparSim_fine_stable2_", i), result_random2)
  gc()
  

  
  
  scry1=devianceFeatureSelection(as.matrix(newdata1), fam = "poisson")
  scry2=devianceFeatureSelection(as.matrix(newdata2), fam = "poisson")
  
  result_scry1 = names(sort(scry1,decreasing = T)[1:2000])
  result_scry2 = names(sort(scry2,decreasing = T)[1:2000])
  
  assign(paste0("scryFS_SparSim_fine_stable1_",i), result_scry1)
  assign(paste0("scryFS_SparSim_fine_stable2_",i), result_scry2)
  gc()
  print("Scry Complete")
  
  seurat.obj1 <- CreateAssayObject(counts = newdata1)
  seurat.obj1 <- NormalizeData(seurat.obj1, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj2 <- CreateAssayObject(counts = newdata2)
  seurat.obj2 <- NormalizeData(seurat.obj2, normalization.method = "LogNormalize", scale.factor = 10000)
  
  seurat.obj.disp1 <- FindVariableFeatures(seurat.obj1, selection.method = "disp", nfeatures = 2000)
  seurat.obj.disp2 <- FindVariableFeatures(seurat.obj2, selection.method = "disp", nfeatures = 2000)
  seurat.obj.mvp1 <- FindVariableFeatures(seurat.obj1, selection.method = "mvp")
  seurat.obj.mvp2 <- FindVariableFeatures(seurat.obj2, selection.method = "mvp")
  seurat.obj.vst1 <- FindVariableFeatures(seurat.obj1, selection.method = "vst", nfeatures = 2000)
  seurat.obj.vst2 <- FindVariableFeatures(seurat.obj2, selection.method = "vst", nfeatures = 2000)
  
  
  
  
  result_SeuratDisp1 <- VariableFeatures(seurat.obj.disp1)
  result_SeuratDisp2 <- VariableFeatures(seurat.obj.disp2)
  result_SeuratMvp1  <- VariableFeatures(seurat.obj.mvp1)
  result_SeuratMvp2  <- VariableFeatures(seurat.obj.mvp2)
  result_SeuratVst1  <- VariableFeatures(seurat.obj.vst1)
  result_SeuratVst2  <- VariableFeatures(seurat.obj.vst2)
  
  
  
  
  assign(paste0("SeuratDisp_SparSim_fine_stable1_", i), result_SeuratDisp1)
  assign(paste0("SeuratDisp_SparSim_fine_stable2_", i), result_SeuratDisp2)
  assign(paste0("SeuratMvp_SparSim_fine_stable1_",  i), result_SeuratMvp1)
  assign(paste0("SeuratMvp_SparSim_fine_stable2_", i), result_SeuratMvp2)
  assign(paste0("SeuratVst_SparSim_fine_stable1_", i), result_SeuratVst1)    
  assign(paste0("SeuratVst_SparSim_fine_stable2_", i), result_SeuratVst2)    
  
  
  
  
  save(list = paste0("Mcadet_SparSim_fine_stable1_", i), file = 
paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/Mcadet_SparSim_fine_stable1_",
i, ".RData"))
  
  save(list = paste0("Mcadet_SparSim_fine_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/Mcadet_SparSim_fine_stable2_",
                i, ".RData"))
  
  
  
  
  save(list = paste0("BreFS_SparSim_fine_stable1_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/BreFS_SparSim_fine_stable1_",
                i, ".RData"))
  
  save(list = paste0("BreFS_SparSim_fine_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/BreFS_SparSim_fine_stable2_",
                i, ".RData"))
  
  
  save(list = paste0("M3Drop_SparSim_fine_stable1_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/M3Drop_SparSim_fine_stable1_",
                i, ".RData"))
  
  save(list = paste0("M3Drop_SparSim_fine_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/M3Drop_SparSim_fine_stable2_",
                i, ".RData"))
  
  
  save(list = paste0("NBDrop_SparSim_fine_stable1_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/NBDrop_SparSim_fine_stable1_",
                i, ".RData"))
  
  save(list = paste0("NBDrop_SparSim_fine_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/NBDrop_SparSim_fine_stable2_",
                i, ".RData"))
  
  
  save(list = paste0("random_SparSim_fine_stable1_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/random_SparSim_fine_stable1_",
                i, ".RData"))
  
  save(list = paste0("random_SparSim_fine_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/random_SparSim_fine_stable2_",
                i, ".RData"))
  
  
  
  save(list = paste0("scryFS_SparSim_fine_stable1_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/scryFS_SparSim_fine_stable1_",
                i, ".RData"))
  
  save(list = paste0("scryFS_SparSim_fine_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/scryFS_SparSim_fine_stable2_",
               i, ".RData"))
  
  
  save(list = paste0("SeuratDisp_SparSim_fine_stable1_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratDisp_SparSim_fine_stable1_",
                i, ".RData"))
  
  save(list = paste0("SeuratDisp_SparSim_fine_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratDisp_SparSim_fine_stable2_",
                i, ".RData"))
  
  
  
  save(list = paste0("SeuratMvp_SparSim_fine_stable1_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratMvp_SparSim_fine_stable1_",
                i, ".RData"))
  
  save(list = paste0("SeuratMvp_SparSim_fine_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratMvp_SparSim_fine_stable2_",
                i, ".RData"))
  
  
  
  save(list = paste0("SeuratVst_SparSim_fine_stable1_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratVst_SparSim_fine_stable1_",
                i, ".RData"))
  
  save(list = paste0("SeuratVst_SparSim_fine_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratVst_SparSim_fine_stable2_",
                i, ".RData"))
  
  
  
  
  rm(list =paste0("SparSim_fine_", i) )
  rm(workdata,newdata1, newdata2)
  gc()
  print(i)

}






### SparSim coarse resolution 

for (i in 1:24) {
  
  load(file = 
         paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_", 
                i, ".RData") )
  workdata = get(paste0("SparSim_", i))
  
  set.seed(i)
  newdata1 = matrix(rbinom(nrow(workdata)*ncol(workdata), as.vector(as.matrix(workdata)), 0.5), nrow=nrow(workdata))
  newdata2 = workdata - newdata1
  colnames(newdata1) = colnames(workdata)
  rownames(newdata1) = rownames(workdata)
  newdata1 = as.data.frame(newdata1)
  colnames(newdata2) = colnames(workdata)
  rownames(newdata2) = rownames(workdata)
  print("Split complete")
  

  
  
  
  
  result1 = mcadet(data = newdata1,
                   n.comp= 60, run=10, n.feature=NA, nk_percent = 0.010, 
                   start_resolution = 0.5,  cell_percent =  0.005,
                   MC_iter = 50000, fdr = 0.15, seed = 1234)
  
  
  result2 = mcadet(data = newdata2,
                   n.comp= 60, run=10, n.feature=NA, nk_percent = 0.010, 
                   start_resolution = 0.5,  cell_percent =  0.005,
                   MC_iter = 50000, fdr = 0.15, seed = 1234)
  
  assign(paste0("Mcadet_SparSim_coarse_stable1_", i),
         result1)
  
  assign(paste0("Mcadet_SparSim_coarse_stable2_", i),
         result2)
  
  print("Mcadet Complete")
  
  
  
  result_Bre1 = BrenneckeGetVariableGenes(newdata1, suppress.plot = T, fdr=0.1, minBiolDisp=0.5)
  result_Bre2 = BrenneckeGetVariableGenes(newdata2, suppress.plot = T, fdr=0.1, minBiolDisp=0.5)
  
  assign(paste0("BreFS_SparSim_coarse_stable1_", i), result_Bre1$Gene)
  assign(paste0("BreFS_SparSim_coarse_stable2_", i), result_Bre2$Gene)
  
  gc()
  
  
  
  M3_data1 <- M3DropConvertData(newdata1, is.counts=TRUE)
  result_M3Drop1<- M3DropFeatureSelection(M3_data1 ,mt_method="fdr", mt_threshold=0.05)
  M3_data2 <- M3DropConvertData(newdata2, is.counts=TRUE)
  result_M3Drop2<- M3DropFeatureSelection(M3_data2 ,mt_method="fdr", mt_threshold=0.05)
  
  assign(paste0("M3Drop_SparSim_coarse_stable1_", i), result_M3Drop1$Gene)
  assign(paste0("M3Drop_SparSim_coarse_stable2_", i), result_M3Drop2$Gene)
  gc()
  print("M3Drop Complete")
  
  nb_raw1 <- NBumiConvertData(newdata1, is.counts=TRUE)
  DANB_fit1 <- NBumiFitModel(nb_raw1)
  result_NBDrop1<- NBumiFeatureSelectionCombinedDrop(DANB_fit1, method="fdr", qval.thres=0.05, suppress.plot=F)
  nb_raw2 <- NBumiConvertData(newdata2, is.counts=TRUE)
  DANB_fit2 <- NBumiFitModel(nb_raw2)
  result_NBDrop2<- NBumiFeatureSelectionCombinedDrop(DANB_fit2, method="fdr", qval.thres=0.05, suppress.plot=F)
  
  
  assign(paste0("NBDrop_SparSim_coarse_stable1_",i), result_NBDrop1$Gene)
  assign(paste0("NBDrop_SparSim_coarse_stable2_",i), result_NBDrop2$Gene)
  gc()
  
  
  
  set.seed( i )
  result_random1 = rownames(newdata1)[sample(c(1:nrow(newdata1)), 2000)]
  set.seed(i )
  result_random2 = rownames(newdata2)[sample(c(1:nrow(newdata2)), 2000)]
  
  
  assign(paste0("random_SparSim_coarse_stable1_", i), result_random1)
  assign(paste0("random_SparSim_coarse_stable2_", i), result_random2)
  gc()
  
  
  scry1=devianceFeatureSelection(as.matrix(newdata1), fam = "poisson")
  scry2=devianceFeatureSelection(as.matrix(newdata2), fam = "poisson")
  
  result_scry1 = names(sort(scry1,decreasing = T)[1:2000])
  result_scry2 = names(sort(scry2,decreasing = T)[1:2000])
  
  assign(paste0("scryFS_SparSim_coarse_stable1_",i), result_scry1)
  assign(paste0("scryFS_SparSim_coarse_stable2_",i), result_scry2)
  gc()
  print("Scry Complete")
  
  seurat.obj1 <- CreateAssayObject(counts = newdata1)
  seurat.obj1 <- NormalizeData(seurat.obj1, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj2 <- CreateAssayObject(counts = newdata2)
  seurat.obj2 <- NormalizeData(seurat.obj2, normalization.method = "LogNormalize", scale.factor = 10000)
  
  seurat.obj.disp1 <- FindVariableFeatures(seurat.obj1, selection.method = "disp", nfeatures = 2000)
  seurat.obj.disp2 <- FindVariableFeatures(seurat.obj2, selection.method = "disp", nfeatures = 2000)
  seurat.obj.mvp1 <- FindVariableFeatures(seurat.obj1, selection.method = "mvp")
  seurat.obj.mvp2 <- FindVariableFeatures(seurat.obj2, selection.method = "mvp")
  seurat.obj.vst1 <- FindVariableFeatures(seurat.obj1, selection.method = "vst", nfeatures = 2000)
  seurat.obj.vst2 <- FindVariableFeatures(seurat.obj2, selection.method = "vst", nfeatures = 2000)
  
  
  
  
  result_SeuratDisp1 <- VariableFeatures(seurat.obj.disp1)
  result_SeuratDisp2 <- VariableFeatures(seurat.obj.disp2)
  result_SeuratMvp1  <- VariableFeatures(seurat.obj.mvp1)
  result_SeuratMvp2  <- VariableFeatures(seurat.obj.mvp2)
  result_SeuratVst1  <- VariableFeatures(seurat.obj.vst1)
  result_SeuratVst2  <- VariableFeatures(seurat.obj.vst2)
  
  
  
  
  assign(paste0("SeuratDisp_SparSim_coarse_stable1_", i), result_SeuratDisp1)
  assign(paste0("SeuratDisp_SparSim_coarse_stable2_", i), result_SeuratDisp2)
  assign(paste0("SeuratMvp_SparSim_coarse_stable1_",  i), result_SeuratMvp1)
  assign(paste0("SeuratMvp_SparSim_coarse_stable2_", i), result_SeuratMvp2)
  assign(paste0("SeuratVst_SparSim_coarse_stable1_", i), result_SeuratVst1)    
  assign(paste0("SeuratVst_SparSim_coarse_stable2_", i), result_SeuratVst2)    
  
  
  
  
  save(list = paste0("Mcadet_SparSim_coarse_stable1_", i), file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/Mcadet_SparSim_coarse_stable1_",
                i, ".RData"))
  
  save(list = paste0("Mcadet_SparSim_coarse_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/Mcadet_SparSim_coarse_stable2_",
                i, ".RData"))
  
  
  
  
  save(list = paste0("BreFS_SparSim_coarse_stable1_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/BreFS_SparSim_coarse_stable1_",
                i, ".RData"))
  
  save(list = paste0("BreFS_SparSim_coarse_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/BreFS_SparSim_coarse_stable2_",
                i, ".RData"))
  
  
  save(list = paste0("M3Drop_SparSim_coarse_stable1_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/M3Drop_SparSim_coarse_stable1_",
                i, ".RData"))
  
  save(list = paste0("M3Drop_SparSim_coarse_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/M3Drop_SparSim_coarse_stable2_",
                i, ".RData"))
  
  
  save(list = paste0("NBDrop_SparSim_coarse_stable1_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/NBDrop_SparSim_coarse_stable1_",
                i, ".RData"))
  
  save(list = paste0("NBDrop_SparSim_coarse_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/NBDrop_SparSim_coarse_stable2_",
                i, ".RData"))
  
  
  save(list = paste0("random_SparSim_coarse_stable1_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/random_SparSim_coarse_stable1_",
                i, ".RData"))
  
  save(list = paste0("random_SparSim_coarse_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/random_SparSim_coarse_stable2_",
                i, ".RData"))
  
  
  
  save(list = paste0("scryFS_SparSim_coarse_stable1_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/scryFS_SparSim_coarse_stable1_",
                i, ".RData"))
  
  save(list = paste0("scryFS_SparSim_coarse_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/scryFS_SparSim_coarse_stable2_",
                i, ".RData"))
  
  
  save(list = paste0("SeuratDisp_SparSim_coarse_stable1_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratDisp_SparSim_coarse_stable1_",
                i, ".RData"))
  
  save(list = paste0("SeuratDisp_SparSim_coarse_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratDisp_SparSim_coarse_stable2_",
                i, ".RData"))
  
  
  
  save(list = paste0("SeuratMvp_SparSim_coarse_stable1_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratMvp_SparSim_coarse_stable1_",
                i, ".RData"))
  
  save(list = paste0("SeuratMvp_SparSim_coarse_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratMvp_SparSim_coarse_stable2_",
                i, ".RData"))
  
  
  
  save(list = paste0("SeuratVst_SparSim_coarse_stable1_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratVst_SparSim_coarse_stable1_",
                i, ".RData"))
  
  save(list = paste0("SeuratVst_SparSim_coarse_stable2_", i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratVst_SparSim_coarse_stable2_",
                i, ".RData"))
  
  
  
  
  rm(list =paste0("SparSim_", i) )
  rm(workdata,newdata1, newdata2)
  gc()
  print(i)
  
}












### PBMC Fine resolution 

p = 0.1
for (i in 1:8) {
  for (j in 1:3) {
    
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
  
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/data_fine/pbmc_fine_", 
                      donor.set[i], "_", time.set[j], ".RData"))
    
    workdata = get(paste0("pbmc_fine_", donor.set[i], "_", time.set[j]))
    
    set.seed((i-1)*12+(j-1)*4+1)
    newdata1 = matrix(rbinom(nrow(workdata)*ncol(workdata), as.vector(as.matrix(workdata)), p), nrow=nrow(workdata))
    newdata2 = workdata - newdata1
    colnames(newdata1) = colnames(workdata)
    rownames(newdata1) = rownames(workdata)
    newdata1 = as.data.frame(newdata1)
    colnames(newdata2) = colnames(workdata)
    rownames(newdata2) = rownames(workdata)
    print("Split complete")
    
  
    
    result1 = mcadet(data = newdata1)
    
    
    result2 = mcadet(data = newdata2)
    
    assign(paste0("Mcadet_", p, "_stable1_", donor.set[i], "_", time.set[j]),
           result1)
    
    assign(paste0("Mcadet_", p, "_stable2_", donor.set[i], "_", time.set[j]),
           result2)
    
    print("Mcadet Complete")
    
    
    
    result_Bre1 = BrenneckeGetVariableGenes(newdata1, suppress.plot = T, fdr=0.1, minBiolDisp=0.5)
    result_Bre2 = BrenneckeGetVariableGenes(newdata2, suppress.plot = T, fdr=0.1, minBiolDisp=0.5)
    
    assign(paste0("BreFS_", p, "_stable1_", donor.set[i], "_", time.set[j]), result_Bre1$Gene)
    assign(paste0("BreFS_", p, "_stable2_", donor.set[i], "_", time.set[j]), result_Bre2$Gene)
    
    gc()
    
    
    
    M3_data1 <- M3DropConvertData(newdata1, is.counts=TRUE)
    result_M3Drop1<- M3DropFeatureSelection(M3_data1 ,mt_method="fdr", mt_threshold=0.05)
    M3_data2 <- M3DropConvertData(newdata2, is.counts=TRUE)
    result_M3Drop2<- M3DropFeatureSelection(M3_data2 ,mt_method="fdr", mt_threshold=0.05)
    
    assign(paste0("M3Drop_", p, "_stable1_",
                  donor.set[i], "_", time.set[j]), result_M3Drop1$Gene)
    assign(paste0("M3Drop_", p, "_stable2_",
                  donor.set[i], "_", time.set[j]), result_M3Drop2$Gene)
    gc()
    print("M3Drop Complete")
    
    nb_raw1 <- NBumiConvertData(newdata1, is.counts=TRUE)
    DANB_fit1 <- NBumiFitModel(nb_raw1)
    result_NBDrop1<- NBumiFeatureSelectionCombinedDrop(DANB_fit1, method="fdr", qval.thres=0.05, suppress.plot=F)
    nb_raw2 <- NBumiConvertData(newdata2, is.counts=TRUE)
    DANB_fit2 <- NBumiFitModel(nb_raw2)
    result_NBDrop2<- NBumiFeatureSelectionCombinedDrop(DANB_fit2, method="fdr", qval.thres=0.05, suppress.plot=F)
    
    
    assign(paste0("NBDrop_", p, "_stable1_",
                  donor.set[i], "_", time.set[j]), result_NBDrop1$Gene)
    assign(paste0("NBDrop_", p, "_stable2_",
                  donor.set[i], "_", time.set[j]), result_NBDrop2$Gene)
    gc()
    
    
    
    set.seed( (i-1)*12+(j-1)*4+1 )
    result_random1 = rownames(newdata1)[sample(c(1:nrow(newdata1)), 2000)]
    set.seed( (i-1)*12+(j-1)*4+2 )
    result_random2 = rownames(newdata2)[sample(c(1:nrow(newdata2)), 2000)]
    
    
    assign(paste0("random_", p, "_stable1_",
                  donor.set[i], "_", time.set[j]), result_random1)
    assign(paste0("random_", p, "_stable2_",
                  donor.set[i], "_", time.set[j]), result_random2)
    gc()
    
    
    scry1=devianceFeatureSelection(as.matrix(newdata1), fam = "poisson")
    scry2=devianceFeatureSelection(as.matrix(newdata2), fam = "poisson")
    
    result_scry1 = names(sort(scry1,decreasing = T)[1:2000])
    result_scry2 = names(sort(scry2,decreasing = T)[1:2000])
    
    assign(paste0("scryFS_", p, "_stable1_", 
                  donor.set[i], "_", time.set[j]), result_scry1)
    assign(paste0("scryFS_", p, "_stable2_", 
                  donor.set[i], "_", time.set[j]), result_scry2)
    gc()
    print("Scry Complete")
    
    seurat.obj1 <- CreateAssayObject(counts = newdata1)
    seurat.obj1 <- NormalizeData(seurat.obj1, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat.obj2 <- CreateAssayObject(counts = newdata2)
    seurat.obj2 <- NormalizeData(seurat.obj2, normalization.method = "LogNormalize", scale.factor = 10000)
    
    seurat.obj.disp1 <- FindVariableFeatures(seurat.obj1, selection.method = "disp", nfeatures = 2000)
    seurat.obj.disp2 <- FindVariableFeatures(seurat.obj2, selection.method = "disp", nfeatures = 2000)
    seurat.obj.mvp1 <- FindVariableFeatures(seurat.obj1, selection.method = "mvp")
    seurat.obj.mvp2 <- FindVariableFeatures(seurat.obj2, selection.method = "mvp")
    seurat.obj.vst1 <- FindVariableFeatures(seurat.obj1, selection.method = "vst", nfeatures = 2000)
    seurat.obj.vst2 <- FindVariableFeatures(seurat.obj2, selection.method = "vst", nfeatures = 2000)
    
    
    
    
    result_SeuratDisp1 <- VariableFeatures(seurat.obj.disp1)
    result_SeuratDisp2 <- VariableFeatures(seurat.obj.disp2)
    result_SeuratMvp1  <- VariableFeatures(seurat.obj.mvp1)
    result_SeuratMvp2  <- VariableFeatures(seurat.obj.mvp2)
    result_SeuratVst1  <- VariableFeatures(seurat.obj.vst1)
    result_SeuratVst2  <- VariableFeatures(seurat.obj.vst2)
    
    
    
    
    assign(paste0("SeuratDisp_", p, "_stable1_", donor.set[i], "_", time.set[j]), result_SeuratDisp1)
    assign(paste0("SeuratDisp_", p, "_stable2_", donor.set[i], "_", time.set[j]), result_SeuratDisp2)
    assign(paste0("SeuratMvp_", p, "_stable1_", donor.set[i], "_", time.set[j]), result_SeuratMvp1)
    assign(paste0("SeuratMvp_", p, "_stable2_", donor.set[i], "_", time.set[j]), result_SeuratMvp2)
    assign(paste0("SeuratVst_", p, "_stable1_", donor.set[i], "_", time.set[j]), result_SeuratVst1)    
    assign(paste0("SeuratVst_", p, "_stable2_", donor.set[i], "_", time.set[j]), result_SeuratVst2)    
    
    
    
    
    
    rm(list =paste0("pbmc_fine_", donor.set[i], "_", time.set[j]) )
    rm(workdata,newdata1, newdata2)
    gc()
    print(c(i,j))
    
    
  }
}




### PBMC coarse resolution 


for (i in 1:8) {
  for (j in 1:3) {
    
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    
    load(file =paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmcdata_", 
                      donor.set[i], "_", time.set[j], ".RData"))
    
    workdata = get(paste0("pbmcdata_", donor.set[i], "_", time.set[j]))
    set.seed((i-1)*12+(j-1)*4+1)
    newdata1 = matrix(rbinom(nrow(workdata)*ncol(workdata), as.vector(as.matrix(workdata)), 0.5), nrow=nrow(workdata))
    newdata2 = workdata - newdata1
    colnames(newdata1) = colnames(workdata)
    rownames(newdata1) = rownames(workdata)
    newdata1 = as.data.frame(newdata1)
    colnames(newdata2) = colnames(workdata)
    rownames(newdata2) = rownames(workdata)
    print("Split complete")
    
    
    result1 = mcadet(data = newdata1,
                     n.comp= 60, run=10, n.feature=NA, nk_percent = 0.010, 
                     start_resolution = 0.5,  cell_percent =  0.005,
                     MC_iter = 50000, fdr = 0.15, seed = 1234)
    
    
    result2 = mcadet(data = newdata2,
                     n.comp= 60, run=10, n.feature=NA, nk_percent = 0.010, 
                     start_resolution = 0.5,  cell_percent =  0.005,
                     MC_iter = 50000, fdr = 0.15, seed = 1234)
    
    assign(paste0("Mcadet_PBMC_coarse_stable1_", donor.set[i], "_", time.set[j]),
           result1)
    
    assign(paste0("Mcadet_PBMC_coarse_stable2_",donor.set[i], "_", time.set[j]),
           result2)
    
    print("Mcadet Complete")
    
    
    
    result_Bre1 = BrenneckeGetVariableGenes(newdata1, suppress.plot = T, fdr=0.1, minBiolDisp=0.5)
    result_Bre2 = BrenneckeGetVariableGenes(newdata2, suppress.plot = T, fdr=0.1, minBiolDisp=0.5)
    
    assign(paste0("BreFS_PBMC_coarse_stable1_", donor.set[i], "_", time.set[j]), result_Bre1$Gene)
    assign(paste0("BreFS_PBMC_coarse_stable2_", donor.set[i], "_", time.set[j]), result_Bre2$Gene)
    
    gc()
    
    
    
    M3_data1 <- M3DropConvertData(newdata1, is.counts=TRUE)
    result_M3Drop1<- M3DropFeatureSelection(M3_data1 ,mt_method="fdr", mt_threshold=0.05)
    M3_data2 <- M3DropConvertData(newdata2, is.counts=TRUE)
    result_M3Drop2<- M3DropFeatureSelection(M3_data2 ,mt_method="fdr", mt_threshold=0.05)
    
    assign(paste0("M3Drop_PBMC_coarse_stable1_",
                  donor.set[i], "_", time.set[j]), result_M3Drop1$Gene)
    assign(paste0("M3Drop_PBMC_coarse_stable2_",
                  donor.set[i], "_", time.set[j]), result_M3Drop2$Gene)
    gc()
    print("M3Drop Complete")
    
    nb_raw1 <- NBumiConvertData(newdata1, is.counts=TRUE)
    DANB_fit1 <- NBumiFitModel(nb_raw1)
    result_NBDrop1<- NBumiFeatureSelectionCombinedDrop(DANB_fit1, method="fdr", qval.thres=0.05, suppress.plot=F)
    nb_raw2 <- NBumiConvertData(newdata2, is.counts=TRUE)
    DANB_fit2 <- NBumiFitModel(nb_raw2)
    result_NBDrop2<- NBumiFeatureSelectionCombinedDrop(DANB_fit2, method="fdr", qval.thres=0.05, suppress.plot=F)

    
    assign(paste0("NBDrop_PBMC_coarse_stable1_",
                  donor.set[i], "_", time.set[j]), result_NBDrop1$Gene)
    assign(paste0("NBDrop_PBMC_coarse_stable2_",
                  donor.set[i], "_", time.set[j]), result_NBDrop2$Gene)
    gc()
    
    
    
    set.seed( (i-1)*12+(j-1)*4+1 )
    result_random1 = rownames(newdata1)[sample(c(1:nrow(newdata1)), 2000)]
    set.seed( (i-1)*12+(j-1)*4+2 )
    result_random2 = rownames(newdata2)[sample(c(1:nrow(newdata2)), 2000)]
    
    
    assign(paste0("random_PBMC_coarse_stable1_",
                  donor.set[i], "_", time.set[j]), result_random1)
    assign(paste0("random_PBMC_coarse_stable2_",
                  donor.set[i], "_", time.set[j]), result_random2)
    gc()
    
    
    scry1=devianceFeatureSelection(as.matrix(newdata1), fam = "poisson")
    scry2=devianceFeatureSelection(as.matrix(newdata2), fam = "poisson")
    
    result_scry1 = names(sort(scry1,decreasing = T)[1:2000])
    result_scry2 = names(sort(scry2,decreasing = T)[1:2000])
    
    assign(paste0("scryFS_PBMC_coarse_stable1_", 
                  donor.set[i], "_", time.set[j]), result_scry1)
    assign(paste0("scryFS_PBMC_coarse_stable2_", 
                  donor.set[i], "_", time.set[j]), result_scry2)
    gc()
    print("Scry Complete")
    
    seurat.obj1 <- CreateAssayObject(counts = newdata1)
    seurat.obj1 <- NormalizeData(seurat.obj1, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat.obj2 <- CreateAssayObject(counts = newdata2)
    seurat.obj2 <- NormalizeData(seurat.obj2, normalization.method = "LogNormalize", scale.factor = 10000)

    seurat.obj.disp1 <- FindVariableFeatures(seurat.obj1, selection.method = "disp", nfeatures = 2000)
    seurat.obj.disp2 <- FindVariableFeatures(seurat.obj2, selection.method = "disp", nfeatures = 2000)
    seurat.obj.mvp1 <- FindVariableFeatures(seurat.obj1, selection.method = "mvp")
    seurat.obj.mvp2 <- FindVariableFeatures(seurat.obj2, selection.method = "mvp")
    seurat.obj.vst1 <- FindVariableFeatures(seurat.obj1, selection.method = "vst", nfeatures = 2000)
    seurat.obj.vst2 <- FindVariableFeatures(seurat.obj2, selection.method = "vst", nfeatures = 2000)
    
    
  
  
    result_SeuratDisp1 <- VariableFeatures(seurat.obj.disp1)
    result_SeuratDisp2 <- VariableFeatures(seurat.obj.disp2)
    result_SeuratMvp1  <- VariableFeatures(seurat.obj.mvp1)
    result_SeuratMvp2  <- VariableFeatures(seurat.obj.mvp2)
    result_SeuratVst1  <- VariableFeatures(seurat.obj.vst1)
    result_SeuratVst2  <- VariableFeatures(seurat.obj.vst2)
    
    
    
    
    assign(paste0("SeuratDisp_PBMC_coarse_stable1_", donor.set[i], "_", time.set[j]), result_SeuratDisp1)
    assign(paste0("SeuratDisp_PBMC_coarse_stable2_", donor.set[i], "_", time.set[j]), result_SeuratDisp2)
    assign(paste0("SeuratMvp_PBMC_coarse_stable1_", donor.set[i], "_", time.set[j]), result_SeuratMvp1)
    assign(paste0("SeuratMvp_PBMC_coarse_stable2_", donor.set[i], "_", time.set[j]), result_SeuratMvp2)
    assign(paste0("SeuratVst_PBMC_coarse_stable1_", donor.set[i], "_", time.set[j]), result_SeuratVst1)    
    assign(paste0("SeuratVst_PBMC_coarse_stable2_", donor.set[i], "_", time.set[j]), result_SeuratVst2)    
    
    
    
    
    
    

    
    
    save(list = paste0("Mcadet_PBMC_coarse_stable1_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/Mcadet_PBMC_coarse_stable1_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    save(list = paste0("Mcadet_PBMC_coarse_stable2_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/Mcadet_PBMC_coarse_stable2_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    
    
    
    save(list = paste0("BreFS_PBMC_coarse_stable1_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/BreFS_PBMC_coarse_stable1_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    save(list = paste0("BreFS_PBMC_coarse_stable2_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/BreFS_PBMC_coarse_stable2_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    
    save(list = paste0("M3Drop_PBMC_coarse_stable1_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/M3Drop_PBMC_coarse_stable1_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    save(list = paste0("M3Drop_PBMC_coarse_stable2_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/M3Drop_PBMC_coarse_stable2_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    
    save(list = paste0("NBDrop_PBMC_coarse_stable1_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/NBDrop_PBMC_coarse_stable1_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    save(list = paste0("NBDrop_PBMC_coarse_stable2_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/NBDrop_PBMC_coarse_stable2_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    
    save(list = paste0("random_PBMC_coarse_stable1_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/random_PBMC_coarse_stable1_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    save(list = paste0("random_PBMC_coarse_stable2_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/random_PBMC_coarse_stable2_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    
    
    save(list = paste0("scryFS_PBMC_coarse_stable1_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/scryFS_PBMC_coarse_stable1_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    save(list = paste0("scryFS_PBMC_coarse_stable2_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/scryFS_PBMC_coarse_stable2_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    
    save(list = paste0("SeuratDisp_PBMC_coarse_stable1_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratDisp_PBMC_coarse_stable1_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    save(list = paste0("SeuratDisp_PBMC_coarse_stable2_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratDisp_PBMC_coarse_stable2_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    
    
    save(list = paste0("SeuratMvp_PBMC_coarse_stable1_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratMvp_PBMC_coarse_stable1_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    save(list = paste0("SeuratMvp_PBMC_coarse_stable2_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratMvp_PBMC_coarse_stable2_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    
    
    save(list = paste0("SeuratVst_PBMC_coarse_stable1_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratVst_PBMC_coarse_stable1_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    save(list = paste0("SeuratVst_PBMC_coarse_stable2_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratVst_PBMC_coarse_stable2_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    
    
    
    rm(list =paste0("pbmcdata_", donor.set[i], "_", time.set[j]) )
    rm(workdata,newdata1, newdata2)
    gc()
    print(c(i,j))
    
    
  }
}








## generate stability data frame


true_SparSim_vector = c()
for(i in 1:2000){
  true_SparSim_vector = append(true_SparSim_vector, paste0("Gene", i))
  
}




stability_output <- function(data, resolution, p){
Jaccard = c()
Method  = c()
group = c()

  if(data == "SparSim"){
    
  for (i in 1:24) {
    
  mcadet1 <- get(paste0("Mcadet_SparSim_" , resolution , "_stable1_", i))$gene
  mcadet2 <- get(paste0("Mcadet_SparSim_" , resolution , "_stable2_", i))$gene
  mcadet_jac = jaccard(mcadet1, mcadet2)
  Jaccard = append(Jaccard, mcadet_jac)
  Method = append(Method, "Mcadet")
  group = append(group, "Between split datasets")
  
  Jaccard = append(Jaccard, jaccard(mcadet1, true_SparSim_vector) )
  Method = append(Method, "Mcadet")
  group = append(group, "Between Split 1 and HVGs")
  
  Jaccard = append(Jaccard, jaccard(mcadet2, true_SparSim_vector) )
  Method = append(Method, "Mcadet")
  group = append(group, "Between Split 2 and HVGs")
  
  

  
  vst1 <- get(paste0("SeuratVst_SparSim_" , resolution , "_stable1_", i))
  vst2 <- get(paste0("SeuratVst_SparSim_" , resolution , "_stable2_", i))
  vst_jac = jaccard(vst1, vst2)
  Jaccard = append(Jaccard, vst_jac)
  Method = append(Method, "Seurat Vst")
  group = append(group, "Between split datasets")
  
  Jaccard = append(Jaccard, jaccard(vst1, true_SparSim_vector) )
  Method = append(Method, "Seurat Vst")
  group = append(group, "Between Split 1 and HVGs")
  
  Jaccard = append(Jaccard, jaccard(vst2, true_SparSim_vector) )
  Method = append(Method, "Seurat Vst")
  group = append(group, "Between Split 2 and HVGs")
  
  
  
  
  disp1 <- get(paste0("SeuratDisp_SparSim_" , resolution , "_stable1_", i))
  disp2 <- get(paste0("SeuratDisp_SparSim_" , resolution , "_stable2_", i))
  disp_jac = jaccard(disp1, disp2)
  Jaccard = append(Jaccard, disp_jac)
  Method = append(Method, "Seurat Disp")
  group = append(group, "Between split datasets")
  
  Jaccard = append(Jaccard, jaccard(disp1, true_SparSim_vector) )
  Method = append(Method, "Seurat Disp")
  group = append(group, "Between Split 1 and HVGs")
  
  Jaccard = append(Jaccard, jaccard(disp2, true_SparSim_vector) )
  Method = append(Method, "Seurat Disp")
  group = append(group, "Between Split 2 and HVGs")
  
  
  
  mvp1 <- get(paste0("SeuratMvp_SparSim_" , resolution , "_stable1_", i))
  mvp2 <- get(paste0("SeuratMvp_SparSim_" , resolution , "_stable2_", i))
  mvp_jac = jaccard(mvp1, mvp2)
  Jaccard = append(Jaccard, mvp_jac)
  Method = append(Method, "Seurat Mvp")
  group = append(group, "Between split datasets")
  
  Jaccard = append(Jaccard, jaccard(mvp1, true_SparSim_vector) )
  Method = append(Method, "Seurat Mvp")
  group = append(group, "Between Split 1 and HVGs")
  
  Jaccard = append(Jaccard, jaccard(mvp2, true_SparSim_vector) )
  Method = append(Method, "Seurat Mvp")
  group = append(group, "Between Split 2 and HVGs")
  
  
  NBdrop1 <- get(paste0("NBDrop_SparSim_" , resolution , "_stable1_", i))
  NBdrop2 <- get(paste0("NBDrop_SparSim_" , resolution , "_stable2_", i))
  NBdrop_jac = jaccard(NBdrop1, NBdrop2)
  Jaccard = append(Jaccard, NBdrop_jac)
  Method = append(Method, "NBDrop")
  group = append(group, "Between split datasets")
  
  Jaccard = append(Jaccard, jaccard(NBdrop1, true_SparSim_vector) )
  Method = append(Method, "NBDrop")
  group = append(group, "Between Split 1 and HVGs")
  
  Jaccard = append(Jaccard, jaccard(NBdrop2, true_SparSim_vector) )
  Method = append(Method, "NBDrop")
  group = append(group, "Between Split 2 and HVGs")
  
  
  M3Drop1 <- get(paste0("M3Drop_SparSim_" , resolution , "_stable1_", i))
  M3Drop2 <- get(paste0("M3Drop_SparSim_" , resolution , "_stable2_", i))
  M3Drop_jac = jaccard(M3Drop1, M3Drop2)
  Jaccard = append(Jaccard, M3Drop_jac)
  Method = append(Method, "M3Drop")
  group = append(group, "Between split datasets")
  
  Jaccard = append(Jaccard, jaccard(M3Drop1, true_SparSim_vector) )
  Method = append(Method, "M3Drop")
  group = append(group, "Between Split 1 and HVGs")
  
  Jaccard = append(Jaccard, jaccard(M3Drop2, true_SparSim_vector) )
  Method = append(Method, "M3Drop")
  group = append(group, "Between Split 2 and HVGs")
  
  
  Brennecke1 <- get(paste0("BreFS_SparSim_" , resolution , "_stable1_", i))
  Brennecke2 <- get(paste0("BreFS_SparSim_" , resolution , "_stable2_", i))
  Brennecke_jac = jaccard(Brennecke1, Brennecke2)
  Jaccard = append(Jaccard, Brennecke_jac)
  Method = append(Method, "Brennecke")
  group = append(group, "Between split datasets")
  
  Jaccard = append(Jaccard, jaccard(Brennecke1, true_SparSim_vector) )
  Method = append(Method, "Brennecke")
  group = append(group, "Between Split 1 and HVGs")
  
  Jaccard = append(Jaccard, jaccard(Brennecke2, true_SparSim_vector) )
  Method = append(Method, "Brennecke")
  group = append(group, "Between Split 2 and HVGs")
  
  
  
  scry1 <- get(paste0("scryFS_SparSim_" , resolution , "_stable1_", i))
  scry2 <- get(paste0("scryFS_SparSim_" , resolution , "_stable2_", i))
  scry_jac = jaccard(scry1, scry2)
  Jaccard = append(Jaccard, scry_jac)
  Method = append(Method, "Scry")
  group = append(group, "Between split datasets")
  
  Jaccard = append(Jaccard, jaccard(scry1, true_SparSim_vector) )
  Method = append(Method, "Scry")
  group = append(group, "Between Split 1 and HVGs")
  
  Jaccard = append(Jaccard, jaccard(scry2, true_SparSim_vector) )
  Method = append(Method, "Scry")
  group = append(group, "Between Split 2 and HVGs")
  
  
  random1 <- get(paste0("random_SparSim_" , resolution , "_stable1_", i))
  random2 <- get(paste0("random_SparSim_" , resolution , "_stable2_", i))
  random_jac = jaccard(random1, random2)
  Jaccard = append(Jaccard, random_jac)
  Method = append(Method, "Random")
  group = append(group, "Between split datasets")
  
  Jaccard = append(Jaccard, jaccard(random1, true_SparSim_vector) )
  Method = append(Method, "Random")
  group = append(group, "Between Split 1 and HVGs")
  
  Jaccard = append(Jaccard, jaccard(random2, true_SparSim_vector) )
  Method = append(Method, "Random")
  group = append(group, "Between Split 2 and HVGs")
  
  
  }
    
    

  }else{
    
    
    for (i in 1:8) {
      for (j in 1:3) {
        
      donor.set<-c(1,2,3,4,5,6,7,8)
      time.set<-c(0,3,7)
      if(resolution == "coarse"){
        pbmc_true = get(paste0("DEgene_", donor.set[i], "_", time.set[j] )) 
      }else{
        pbmc_true = get(paste0("DEgene_", resolution, "_", donor.set[i], "_", time.set[j] )) 
      }
      

      mcadet1 <- get(paste0("Mcadet_" ,  p, "_stable1_", donor.set[i], "_", time.set[j]))$gene
      mcadet2 <- get(paste0("Mcadet_" ,  p, "_stable2_", donor.set[i], "_", time.set[j]))$gene
      mcadet_jac = jaccard(mcadet1, mcadet2)
      Jaccard = append(Jaccard, mcadet_jac)
      Method = append(Method, "Mcadet")
      group = append(group, "Between split datasets")
      
      Jaccard = append(Jaccard, jaccard(mcadet1, pbmc_true))
      Method = append(Method, "Mcadet")
      group = append(group, "Between Split 1 and HVGs")
      
      Jaccard = append(Jaccard, jaccard(mcadet2, pbmc_true))
      Method = append(Method, "Mcadet")
      group = append(group, "Between Split 2 and HVGs")
      
      
      
      vst1 <- get(paste0("SeuratVst_" , p, "_stable1_", donor.set[i], "_", time.set[j]))
      vst2 <- get(paste0("SeuratVst_" , p, "_stable2_", donor.set[i], "_", time.set[j]))
      vst_jac = jaccard(vst1, vst2)
      Jaccard = append(Jaccard, vst_jac)
      Method = append(Method, "Seurat Vst")
      group = append(group, "Between split datasets")
      
      Jaccard = append(Jaccard, jaccard(vst1, pbmc_true))
      Method = append(Method, "Seurat Vst")
      group = append(group, "Between Split 1 and HVGs")
      
      Jaccard = append(Jaccard, jaccard(vst2, pbmc_true))
      Method = append(Method, "Seurat Vst")
      group = append(group, "Between Split 2 and HVGs")
      
      
      disp1 <- get(paste0("SeuratDisp_" , p, "_stable1_", donor.set[i], "_", time.set[j]))
      disp2 <- get(paste0("SeuratDisp_" , p, "_stable2_", donor.set[i], "_", time.set[j]))
      disp_jac = jaccard(disp1, disp2)
      Jaccard = append(Jaccard, disp_jac)
      Method = append(Method, "Seurat Disp")
      group = append(group, "Between split datasets")
      
      Jaccard = append(Jaccard, jaccard(disp1, pbmc_true))
      Method = append(Method, "Seurat Disp")
      group = append(group, "Between Split 1 and HVGs")
      
      Jaccard = append(Jaccard, jaccard(disp2, pbmc_true))
      Method = append(Method, "Seurat Disp")
      group = append(group, "Between Split 2 and HVGs")
      
      
      
      mvp1 <- get(paste0("SeuratMvp_" , p , "_stable1_", donor.set[i], "_", time.set[j]))
      mvp2 <- get(paste0("SeuratMvp_" , p , "_stable2_", donor.set[i], "_", time.set[j]))
      mvp_jac = jaccard(mvp1, mvp2)
      Jaccard = append(Jaccard, mvp_jac)
      Method = append(Method, "Seurat Mvp")
      group = append(group, "Between split datasets")
      
      Jaccard = append(Jaccard, jaccard(mvp1, pbmc_true))
      Method = append(Method, "Seurat Mvp")
      group = append(group, "Between Split 1 and HVGs")
      
      Jaccard = append(Jaccard, jaccard(mvp2, pbmc_true))
      Method = append(Method, "Seurat Mvp")
      group = append(group, "Between Split 2 and HVGs")
      
      NBdrop1 <- get(paste0("NBDrop_" , p , "_stable1_", donor.set[i], "_", time.set[j]))
      NBdrop2 <- get(paste0("NBDrop_" , p , "_stable2_", donor.set[i], "_", time.set[j]))
      NBdrop_jac = jaccard(NBdrop1, NBdrop2)
      Jaccard = append(Jaccard, NBdrop_jac)
      Method = append(Method, "NBDrop")
      group = append(group, "Between split datasets")
      
      Jaccard = append(Jaccard, jaccard(NBdrop1, pbmc_true))
      Method = append(Method, "NBDrop")
      group = append(group, "Between Split 1 and HVGs")
      
      Jaccard = append(Jaccard, jaccard(NBdrop2, pbmc_true))
      Method = append(Method, "NBDrop")
      group = append(group, "Between Split 2 and HVGs")
      
      M3Drop1 <- get(paste0("M3Drop_" , p , "_stable1_", donor.set[i], "_", time.set[j]))
      M3Drop2 <- get(paste0("M3Drop_" , p , "_stable2_", donor.set[i], "_", time.set[j]))
      M3Drop_jac = jaccard(M3Drop1, M3Drop2)
      Jaccard = append(Jaccard, M3Drop_jac)
      Method = append(Method, "M3Drop")
      group = append(group, "Between split datasets")
      
      Jaccard = append(Jaccard, jaccard(M3Drop1, pbmc_true))
      Method = append(Method, "M3Drop")
      group = append(group, "Between Split 1 and HVGs")
      
      Jaccard = append(Jaccard, jaccard(M3Drop2, pbmc_true))
      Method = append(Method, "M3Drop")
      group = append(group, "Between Split 2 and HVGs")
      
      
      Brennecke1 <- get(paste0("BreFS_" , p , "_stable1_", donor.set[i], "_", time.set[j]))
      Brennecke2 <- get(paste0("BreFS_" , p , "_stable2_", donor.set[i], "_", time.set[j]))
      Brennecke_jac = jaccard(Brennecke1, Brennecke2)
      Jaccard = append(Jaccard, Brennecke_jac)
      Method = append(Method, "Brennecke")
      group = append(group, "Between split datasets")
      
      Jaccard = append(Jaccard, jaccard(Brennecke1, pbmc_true))
      Method = append(Method, "Brennecke")
      group = append(group, "Between Split 1 and HVGs")
      
      Jaccard = append(Jaccard, jaccard(Brennecke2, pbmc_true))
      Method = append(Method, "Brennecke")
      group = append(group, "Between Split 2 and HVGs")
      
      
      
      scry1 <- get(paste0("scryFS_" , p , "_stable1_", donor.set[i], "_", time.set[j]))
      scry2 <- get(paste0("scryFS_" , p , "_stable2_", donor.set[i], "_", time.set[j]))
      scry_jac = jaccard(scry1, scry2)
      Jaccard = append(Jaccard, scry_jac)
      Method = append(Method, "Scry")
      group = append(group, "Between split datasets")
      
      Jaccard = append(Jaccard, jaccard(scry1, pbmc_true))
      Method = append(Method, "Scry")
      group = append(group, "Between Split 1 and HVGs")
      
      Jaccard = append(Jaccard, jaccard(scry2, pbmc_true))
      Method = append(Method, "Scry")
      group = append(group, "Between Split 2 and HVGs")
      
      
      
      random1 <- get(paste0("random_" , p , "_stable1_", donor.set[i], "_", time.set[j]))
      random2 <- get(paste0("random_" , p , "_stable2_", donor.set[i], "_", time.set[j]))
      random_jac = jaccard(random1, random2)
      Jaccard = append(Jaccard, random_jac)
      Method = append(Method, "Random")
      group = append(group, "Between split datasets")
      
      Jaccard = append(Jaccard, jaccard(random1, pbmc_true))
      Method = append(Method, "Random")
      group = append(group, "Between Split 1 and HVGs")
      
      Jaccard = append(Jaccard, jaccard(random2, pbmc_true))
      Method = append(Method, "Random")
      group = append(group, "Between Split 2 and HVGs")
      
      
      
     }
    }
  }
  
   


 out <- data.frame(Jaccard, Method, group)
 return(out)
}



stability_df_0.5 <- stability_output("PBMC", "fine", p = 0.5)
stability_df_0.4 <- stability_output("PBMC", "fine", p = 0.4)
stability_df_0.3 <- stability_output("PBMC", "fine", p = 0.3)
stability_df_0.2 <- stability_output("PBMC", "fine", p = 0.2)
stability_df_0.1 <- stability_output("PBMC", "fine", p = 0.1)


### load data

for (i in 1:24) {
  
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/Mcadet_SparSim_fine_stable1_",
                i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/Mcadet_SparSim_fine_stable2_",
                     i, ".RData"))  
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/Mcadet_SparSim_coarse_stable1_",
                   i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/Mcadet_SparSim_coarse_stable2_",
                   i, ".RData"))   
  
  
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/NBDrop_SparSim_fine_stable1_",
                   i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/NBDrop_SparSim_fine_stable2_",
                   i, ".RData"))  
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/NBDrop_SparSim_coarse_stable1_",
                   i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/NBDrop_SparSim_coarse_stable2_",
                   i, ".RData"))  
  

load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/M3Drop_SparSim_fine_stable1_",
                   i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/M3Drop_SparSim_fine_stable2_",
                   i, ".RData"))  
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/M3Drop_SparSim_coarse_stable1_",
                   i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/M3Drop_SparSim_coarse_stable2_",
                   i, ".RData")) 



load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/BreFS_SparSim_fine_stable1_",
                   i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/BreFS_SparSim_fine_stable2_",
                   i, ".RData"))  
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/BreFS_SparSim_coarse_stable1_",
                   i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/BreFS_SparSim_coarse_stable2_",
                   i, ".RData")) 

  



load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/scryFS_SparSim_fine_stable1_",
                   i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/scryFS_SparSim_fine_stable2_",
                   i, ".RData"))  
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/scryFS_SparSim_coarse_stable1_",
                   i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/scryFS_SparSim_coarse_stable2_",
                   i, ".RData")) 




load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratDisp_SparSim_fine_stable1_",
                   i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratDisp_SparSim_fine_stable2_",
                   i, ".RData"))  
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratDisp_SparSim_coarse_stable1_",
                   i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratDisp_SparSim_coarse_stable2_",
                   i, ".RData")) 




load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratMvp_SparSim_fine_stable1_",
                   i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratMvp_SparSim_fine_stable2_",
                   i, ".RData"))  
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratMvp_SparSim_coarse_stable1_",
                   i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratMvp_SparSim_coarse_stable2_",
                   i, ".RData"))



load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratVst_SparSim_fine_stable1_",
                   i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratVst_SparSim_fine_stable2_",
                   i, ".RData"))  
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratVst_SparSim_coarse_stable1_",
                   i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratVst_SparSim_coarse_stable2_",
                   i, ".RData"))




load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/random_SparSim_fine_stable1_",
                   i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/random_SparSim_fine_stable2_",
                   i, ".RData"))  
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/random_SparSim_coarse_stable1_",
                   i, ".RData"))
load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/random_SparSim_coarse_stable2_",
                   i, ".RData"))

}







for (i in 1:8) {
  for (j in 1:3) {
    
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)

  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/Mcadet_PBMC_fine_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/Mcadet_PBMC_fine_stable2_",
                     donor.set[i], "_", time.set[j], ".RData"))  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/Mcadet_PBMC_coarse_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/Mcadet_PBMC_coarse_stable2_",
                     donor.set[i], "_", time.set[j], ".RData"))   
  
  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/NBDrop_PBMC_fine_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/NBDrop_PBMC_fine_stable2_",
                     donor.set[i], "_", time.set[j], ".RData"))  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/NBDrop_PBMC_coarse_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/NBDrop_PBMC_coarse_stable2_",
                     donor.set[i], "_", time.set[j], ".RData"))  
  
  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/M3Drop_PBMC_fine_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/M3Drop_PBMC_fine_stable2_",
                     donor.set[i], "_", time.set[j], ".RData"))  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/M3Drop_PBMC_coarse_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/M3Drop_PBMC_coarse_stable2_",
                     donor.set[i], "_", time.set[j], ".RData")) 
  
  
  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/BreFS_PBMC_fine_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/BreFS_PBMC_fine_stable2_",
                     donor.set[i], "_", time.set[j], ".RData"))  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/BreFS_PBMC_coarse_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/BreFS_PBMC_coarse_stable2_",
                     donor.set[i], "_", time.set[j], ".RData")) 
  
  
  
  
  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/scryFS_PBMC_fine_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/scryFS_PBMC_fine_stable2_",
                     donor.set[i], "_", time.set[j], ".RData"))  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/scryFS_PBMC_coarse_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/scryFS_PBMC_coarse_stable2_",
                     donor.set[i], "_", time.set[j], ".RData")) 
  
  
  
  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratDisp_PBMC_fine_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratDisp_PBMC_fine_stable2_",
                     donor.set[i], "_", time.set[j], ".RData"))  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratDisp_PBMC_coarse_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratDisp_PBMC_coarse_stable2_",
                     donor.set[i], "_", time.set[j], ".RData")) 
  
  
  
  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratMvp_PBMC_fine_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratMvp_PBMC_fine_stable2_",
                     donor.set[i], "_", time.set[j], ".RData"))  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratMvp_PBMC_coarse_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratMvp_PBMC_coarse_stable2_",
                     donor.set[i], "_", time.set[j], ".RData"))
  
  
  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratVst_PBMC_fine_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratVst_PBMC_fine_stable2_",
                     donor.set[i], "_", time.set[j], ".RData"))  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratVst_PBMC_coarse_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/SeuratVst_PBMC_coarse_stable2_",
                     donor.set[i], "_", time.set[j], ".RData"))
  
  
  
  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/random_PBMC_fine_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/random_PBMC_fine_stable2_",
                     donor.set[i], "_", time.set[j], ".RData"))  
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/random_PBMC_coarse_stable1_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/stability_results/random_PBMC_coarse_stable2_",
                     donor.set[i], "_", time.set[j], ".RData"))
  
 }
}





## generate output boxplot figure

generate_figure_boxplot <- function(DF, plot_title){


DF$Method = factor(DF$Method,
                           levels = c("Mcadet",
                                      "Seurat Vst",
                                      "Seurat Disp",
                                      "Seurat Mvp",
                                      "NBDrop",
                                      "M3Drop",
                                      "Brennecke",
                                      "Scry",
                                      "Random"))

library(ggplot2)
library(ggpubr)

cbPalette <- c("#cc340c", "#CC79A7", "#56B4E9",
               "#009E73", "#E69F00", "#F0E442", 
               "#999999","#0072B2","#D55E00")


out_figure = ggplot(DF, 
                    aes(Method, Jaccard, color = Method )) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  guides(fill = "none") +
  labs(x = "", y = "Normalized Jaccard Similarity", 
       title = plot_title )+theme_bw()+ 
  scale_color_manual(values = cbPalette)+
  stat_compare_means(label = "p.signif", method = "t.test",label.y = 3.5,
                     ref.group = ".all.")+
  geom_hline(yintercept = mean(DF$Jaccard),
             linetype=2)+
  theme( panel.grid = element_blank(),  
         axis.title = element_text(face ="bold", size = 16),
         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
         legend.position = "none",
         axis.text.x = element_text(face = "bold", size = 14, color = "black"), 
         axis.text.y = element_text(face = "bold", size = 12))

return(out_figure)

}



generate_figure_boxplot(DF =stability_SparSim_coarse_df, plot_title = "Simulation coarse datasets" )
generate_figure_boxplot(DF =stability_SparSim_fine_df, plot_title = "Simulation fine datasets" )
generate_figure_boxplot(DF =stability_PBMC_coarse_df, plot_title = "PBMC coarse datasets" )
generate_figure_boxplot(DF =stability_PBMC_fine_df, plot_title = "PBMC fine datasets" )







### generate barplot


generate_barplot_df <- function(DF){
  
  
  mean_jaccard = c()
  sd_jaccard = c()
  
  
  
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Mcadet" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Mcadet" & DF$group == "Between Split 1 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "NBDrop" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "NBDrop" & DF$group == "Between Split 1 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "M3Drop" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "M3Drop" & DF$group == "Between Split 1 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Brennecke" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Brennecke" & DF$group == "Between Split 1 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Scry" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Scry" & DF$group == "Between Split 1 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Seurat Disp" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Seurat Disp" & DF$group == "Between Split 1 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Seurat Mvp" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Seurat Mvp" & DF$group == "Between Split 1 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Seurat Vst" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Seurat Vst" & DF$group == "Between Split 1 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Random" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Random" & DF$group == "Between Split 1 and HVGs"]))
  
  
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Mcadet" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Mcadet" & DF$group == "Between Split 2 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "NBDrop" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "NBDrop" & DF$group == "Between Split 2 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "M3Drop" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "M3Drop" & DF$group == "Between Split 2 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Brennecke" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Brennecke" & DF$group == "Between Split 2 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Scry" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Scry" & DF$group == "Between Split 2 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Seurat Disp" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Seurat Disp" & DF$group == "Between Split 2 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Seurat Mvp" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Seurat Mvp" & DF$group == "Between Split 2 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Seurat Vst" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Seurat Vst" & DF$group == "Between Split 2 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Random" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Random" & DF$group == "Between Split 2 and HVGs"]))  
  
  
  DF_bar = data.frame("mean" = mean_jaccard, "sd" = sd_jaccard, "Method" = 
                        c("Mcadet", "NBDrop", "M3Drop", "Brennecke", "Scry", "Seurat Disp",
                          "Seurat Mvp", "Seurat Vst", "Random"),
                      "group" = rep(c("Between Split 1 and semi-true HVGs",
                                      "Between Split 2 and semi-true HVGs"), each = 9 ) )

 return(DF_bar)
}


barplot_df_0.5 <- generate_barplot_df(DF = stability_df_0.5)
barplot_df_0.4 <- generate_barplot_df(DF = stability_df_0.4)
barplot_df_0.3 <- generate_barplot_df(DF = stability_df_0.3)
barplot_df_0.2 <- generate_barplot_df(DF = stability_df_0.2)
barplot_df_0.1 <- generate_barplot_df(DF = stability_df_0.1)





save(barplot_df_0.5, file = "C:/Users/css22/Desktop/Thesis1/Consistency_improve/barplot_df_0.5.RData")
save(barplot_df_0.4, file = "C:/Users/css22/Desktop/Thesis1/Consistency_improve/barplot_df_0.4.RData")
save(barplot_df_0.3, file = "C:/Users/css22/Desktop/Thesis1/Consistency_improve/barplot_df_0.3.RData")
save(barplot_df_0.2, file = "C:/Users/css22/Desktop/Thesis1/Consistency_improve/barplot_df_0.2.RData")
save(barplot_df_0.1, file = "C:/Users/css22/Desktop/Thesis1/Consistency_improve/barplot_df_0.1.RData")


barplot_df = rbind(barplot_df_0.1, barplot_df_0.2, barplot_df_0.3, barplot_df_0.4, barplot_df_0.5)
barplot_df$"probability" = rep( c(0.1,0.2,0.3,0.4,0.5), each= 18 )
barplot_df$Method = factor(barplot_df$Method, levels = c("Mcadet", "Scry", "Seurat Disp", "Seurat Vst", "Seurat Mvp",
                                                         "Brennecke", "M3Drop", "NBDrop", "Random"))





library(ggplot2)
library(ggpubr)

cbPalette <- c("#cc340c", "#CC79A7", "#56B4E9",
               "#009E73", "#E69F00", "#F0E442", 
               "#999999","#0072B2","#D55E00")

ggplot(barplot_df[barplot_df$probability == 0.5,], 
                      aes(x= Method, y = mean, fill = group )) +
    geom_bar(width = 0.5, stat="identity", position = position_dodge(),
             colour = "black", linewidth = 1) +
     geom_errorbar(aes(ymax = mean+sd, ymin= mean-sd),
                                          width = 0.4, position = position_dodge(0.5))+
    labs(x = "", y = "Mean Jaccard Similarity", 
         title = "" )+theme_bw() + 
  scale_fill_manual(values = c("#7FC97F", "#BEAED4"), name = "Group",
                    labels = c(
                               "Between data 1 and semi-ground truth", 
                               "Between data 2 and semi-ground truth"))+
    geom_hline(yintercept = 0, linewidth = 0.8, color = "blue", linetype=1)+
    theme( panel.grid = element_blank(),  
           axis.title = element_text(face ="bold", size = 16),
           plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
           legend.position = "bottom",
           legend.text = element_text(size = 10, face = "bold"),
           legend.title = element_blank(),  
           legend.background = element_rect(color = "black", 
                                            linewidth = 1.5,
                                            fill = "lightblue"),
           axis.text.x = element_text(face = "bold", size = 14, color = "black"), 
           axis.text.y = element_text(face = "bold", size = 12))





probability_labels <- function(variable, value) {
  paste("Probability =", value)
}

ggplot(barplot_df[barplot_df$probability != 0.5,], 
       aes(x= Method, y = mean, fill = group )) +
  geom_bar(width = 0.3, stat="identity", position = position_dodge(), 
           colour = "black", linewidth = 1) +
  geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd),
                width = 0.2, position = position_dodge(0.3)) +
  labs(x = "", y = "Mean Jaccard Similarity", 
       title = "" ) +
  theme_bw() + 
  scale_fill_manual(values = c("#7FC97F", "#BEAED4"), name = "Group",
                    labels = c(
                      "Between data 1 and semi-ground truth", 
                      "Between data 2 and semi-ground truth")) +
  geom_hline(yintercept = 0, linewidth = 0.8, color = "blue", linetype = 1) +
  theme(panel.grid = element_blank(),  
        axis.title = element_text(face = "bold", size = 16),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "bottom",
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(),  
        legend.background = element_rect(color = "black", 
                                         linewidth = 1.5,
                                         fill = "lightblue"),
        axis.text.x = element_text(face = "bold", size = 14, color = "black"), 
        axis.text.y = element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold", size = 14)) +  
  facet_wrap(~ probability, nrow = 5, labeller = probability_labels)  




DF = stability_PBMC_fine_df
plot_title = "1"
show_legend = TRUE
library(RColorBrewer)
brewer.pal(8,"Accent")


generate_figure_barplot <- function(DF, plot_title, show_legend){
  
  
  mean_jaccard = c()
  sd_jaccard = c()
  

  
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Mcadet" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Mcadet" & DF$group == "Between Split 1 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "NBDrop" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "NBDrop" & DF$group == "Between Split 1 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "M3Drop" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "M3Drop" & DF$group == "Between Split 1 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Brennecke" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Brennecke" & DF$group == "Between Split 1 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Scry" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Scry" & DF$group == "Between Split 1 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Seurat Disp" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Seurat Disp" & DF$group == "Between Split 1 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Seurat Mvp" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Seurat Mvp" & DF$group == "Between Split 1 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Seurat Vst" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Seurat Vst" & DF$group == "Between Split 1 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Random" & DF$group == "Between Split 1 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Random" & DF$group == "Between Split 1 and HVGs"]))
  
  
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Mcadet" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Mcadet" & DF$group == "Between Split 2 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "NBDrop" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "NBDrop" & DF$group == "Between Split 2 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "M3Drop" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "M3Drop" & DF$group == "Between Split 2 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Brennecke" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Brennecke" & DF$group == "Between Split 2 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Scry" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Scry" & DF$group == "Between Split 2 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Seurat Disp" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Seurat Disp" & DF$group == "Between Split 2 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Seurat Mvp" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Seurat Mvp" & DF$group == "Between Split 2 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Seurat Vst" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Seurat Vst" & DF$group == "Between Split 2 and HVGs"]))
  mean_jaccard = append(mean_jaccard, mean(DF$Jaccard[DF$Method == "Random" & DF$group == "Between Split 2 and HVGs"]))
  sd_jaccard = append(sd_jaccard, sd(DF$Jaccard[DF$Method == "Random" & DF$group == "Between Split 2 and HVGs"]))  
  
  
  DF_bar = data.frame("mean" = mean_jaccard, "sd" = sd_jaccard, "Method" = 
               c("Mcadet", "NBDrop", "M3Drop", "Brennecke", "Scry", "Seurat Disp",
                 "Seurat Mvp", "Seurat Vst", "Random"),
             "group" = rep(c("Between split datasets", "Between Split 1 and HVGs",
                             "Between Split 2 and HVGs"), each = 9 ) )
  
  DF_bar$Method = factor(DF_bar$Method,
                     levels = c("Mcadet",
                                "Seurat Vst",
                                "Seurat Disp",
                                "Seurat Mvp",
                                "NBDrop",
                                "M3Drop",
                                "Brennecke",
                                "Scry",
                                "Random"))
  
  DF_bar$group = factor(DF_bar$group,
                    levels = c("Between split datasets","Between Split 1 and HVGs",
                               "Between Split 2 and HVGs"))
  
  library(ggplot2)
  library(ggpubr)
  
  cbPalette <- c("#cc340c", "#CC79A7", "#56B4E9",
                 "#009E73", "#E69F00", "#F0E442", 
                 "#999999","#0072B2","#D55E00")
  if(!show_legend){
  out_figure = ggplot(DF_bar, 
                      aes(x= Method, y = mean, fill = group )) +
    geom_bar(width = 0.5, stat="identity", position = position_dodge() ) +
    guides(fill = "none") + geom_errorbar(aes(ymax = mean+sd, ymin= mean-sd),
                            width = 0.3, position = position_dodge(0.5))+
    labs(x = "", y = "Mean Jaccard Similarity (Normalized)", 
         title = plot_title )+theme_bw() + 
    scale_fill_brewer(palette = "Accent") +
    geom_hline(yintercept = 0, linewidth = 0.8, color = "blue", linetype=1)+
    theme( panel.grid = element_blank(),  
           axis.title = element_text(face ="bold", size = 16),
           plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
           legend.position = "none",
           axis.text.x = element_text(face = "bold", size = 14, color = "black"), 
           axis.text.y = element_text(face = "bold", size = 12))
  }else{
    
    out_figure = ggplot(DF_bar, 
                        aes(x= Method, y = mean, fill = group )) +
      geom_bar(width = 0.5, stat="identity", position = position_dodge() ) +
      geom_errorbar(aes(ymax = mean+sd, ymin= mean-sd),
                                            width = 0.3, position = position_dodge(0.5))+
      labs(x = "", y = "Mean Jaccard Similarity (Normalized)", 
           title = plot_title )+theme_bw() + 
      scale_fill_manual(values = c("#7FC97F", "#BEAED4", "#FDC086"), name = "Group",
                        labels = c("Between data 1 and 2",
                                   "Between data 1 and ground truth", 
                                   "Between data 2 and ground truth"))+
      geom_hline(yintercept = 0, linewidth = 0.8, color = "blue", linetype=1)+
      theme( panel.grid = element_blank(),  
             axis.title = element_text(face ="bold", size = 16),
             plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
             legend.position = "right",
             legend.text = element_text(size = 10, face = "bold"),
             legend.title = element_text(size = 12, face = "bold"),
             legend.background = element_rect(color = "black", 
                                              linewidth = 1.5,
                                              fill = "lightblue"),
             axis.text.x = element_text(face = "bold", size = 14, color = "black"), 
             axis.text.y = element_text(face = "bold", size = 12))
    
  }
  
  return(out_figure)
  
}


generate_figure_barplot(stability_SparSim_coarse_df, "Simulation coarse datasets", FALSE)
generate_figure_barplot(stability_SparSim_fine_df, "Simulation fine datasets", FALSE)
generate_figure_barplot(stability_PBMC_coarse_df, "PBMC coarse datasets",FALSE)
generate_figure_barplot(stability_PBMC_fine_df, "PBMC fine datasets", FALSE)
generate_figure_barplot(stability_PBMC_fine_df, "PBMC fine datasets", TRUE)





