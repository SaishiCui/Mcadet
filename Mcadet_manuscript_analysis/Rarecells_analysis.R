rm(list=ls())
gc()
################ load data #########################################

library(M3Drop)
library(dplyr)
library(Seurat)
library(scry)

for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =
           paste0("C:/Users/css22/Desktop/Thesis1/pbmc/DE_genes/rank_df_", donor.set[i], "_", time.set[j], ".RData"))
    load(file =
           paste0("C:/Users/css22/Desktop/Thesis1/pbmc/DE_genes/rank_df_fine_", donor.set[i], "_", time.set[j], ".RData"))
    
  }
}




## generate HVGs for rare PBMC coarse (B cells) datasets ##


for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    rank_df <- get(paste0( "rank_df_", donor.set[i], "_", time.set[j] ))
    cluster_label <- rank_df$cluster=="B"
    extract_rank_df<- rank_df[cluster_label,]
    order_by_padjust<-order(extract_rank_df[,"p_val_adj"],decreasing = F)
    
    sig_label <- extract_rank_df[order_by_padjust,"p_val_adj"]<0.05
    extract_DEgene <- extract_rank_df[order_by_padjust[sig_label],"gene"]
    assign(paste0("rare_PBMC_coarse_HVG_", donor.set[i], "_", time.set[j]),
           extract_DEgene)
    
    save(list = paste0("rare_PBMC_coarse_HVG_", donor.set[i], "_", time.set[j]),
         file = paste0("C:/Users/css22/Desktop/Thesis1/Results/rare_results/rare_PBMC_coarse_HVG_", 
                       donor.set[i], "_", time.set[j], ".RData"))  
    print(c(i,j))
    
  }
}





## generate HVGs for rare PBMC fine (dnT cells) datasets ##


for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    rank_df <- get(paste0( "rank_df_fine_", donor.set[i], "_", time.set[j] ))
    cluster_label <- rank_df$cluster=="dnT"
    extract_rank_df<- rank_df[cluster_label,]
    order_by_padjust<-order(extract_rank_df[,"p_val_adj"],decreasing = F)
    
    sig_label <- extract_rank_df[order_by_padjust,"p_val_adj"]<0.05
    extract_DEgene <- extract_rank_df[order_by_padjust[sig_label],"gene"]
    assign(paste0("rare_PBMC_fine_HVG_", donor.set[i], "_", time.set[j]),
           extract_DEgene)
    
    save(list = paste0("rare_PBMC_fine_HVG_", donor.set[i], "_", time.set[j]),
         file = paste0("C:/Users/css22/Desktop/Thesis1/Results/rare_results/rare_PBMC_fine_HVG_", 
                       donor.set[i], "_", time.set[j], ".RData"))  
    print(c(i,j))
    
  }
}


## HVGs for dnT cells for sample donor 5, t = 7 and donor 6, t = 1, should be re-designed.

rank_df <- get(paste0( "rank_df_fine_", 5, "_", 0 ))
cluster_label <- rank_df$cluster=="dnT"
extract_rank_df<- rank_df[cluster_label,]
order_by_padjust<-order(extract_rank_df[,"p_val_adj"],decreasing = F)
sig_label <- extract_rank_df[order_by_padjust,"p_val_adj"]<0.05
extract_DEgene1 <- extract_rank_df[order_by_padjust[sig_label],"gene"]
rank_df <- get(paste0( "rank_df_fine_", 5, "_", 3 ))
cluster_label <- rank_df$cluster=="dnT"
extract_rank_df<- rank_df[cluster_label,]
order_by_padjust<-order(extract_rank_df[,"p_val_adj"],decreasing = F)
sig_label <- extract_rank_df[order_by_padjust,"p_val_adj"]<0.05
extract_DEgene2 <- extract_rank_df[order_by_padjust[sig_label],"gene"]


assign(paste0("rare_PBMC_fine_HVG_", 5, "_", 7),
       union(extract_DEgene1, extract_DEgene2))

save(list = paste0("rare_PBMC_fine_HVG_",5, "_", 7),
     file = paste0("C:/Users/css22/Desktop/Thesis1/Results/rare_results/rare_PBMC_fine_HVG_", 
                   5, "_", 7, ".RData")) 



rank_df <- get(paste0( "rank_df_fine_", 6, "_", 3 ))
cluster_label <- rank_df$cluster=="dnT"
extract_rank_df<- rank_df[cluster_label,]
order_by_padjust<-order(extract_rank_df[,"p_val_adj"],decreasing = F)
sig_label <- extract_rank_df[order_by_padjust,"p_val_adj"]<0.05
extract_DEgene1 <- extract_rank_df[order_by_padjust[sig_label],"gene"]
rank_df <- get(paste0( "rank_df_fine_", 6, "_", 7 ))
cluster_label <- rank_df$cluster=="dnT"
extract_rank_df<- rank_df[cluster_label,]
order_by_padjust<-order(extract_rank_df[,"p_val_adj"],decreasing = F)
sig_label <- extract_rank_df[order_by_padjust,"p_val_adj"]<0.05
extract_DEgene2 <- extract_rank_df[order_by_padjust[sig_label],"gene"]


assign(paste0("rare_PBMC_fine_HVG_", 6, "_", 0),
       union(extract_DEgene1, extract_DEgene2))

save(list = paste0("rare_PBMC_fine_HVG_",6, "_", 0),
     file = paste0("C:/Users/css22/Desktop/Thesis1/Results/rare_results/rare_PBMC_fine_HVG_", 
                   6, "_",  0, ".RData")) 


########## load data


for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)

    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/pbmc_results/Mcadet_pbmc_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/BreFS_PBMC_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/M3Drop_PBMC_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/NBDrop_PBMC_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/scryFS_PBMC_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/SeuratDisp_PBMC_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/SeuratMvp_PBMC_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/SeuratVst_PBMC_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/random_PBMC_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    
    
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/pbmc_results/Mcadet_pbmc_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/BreFS_PBMC_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/M3Drop_PBMC_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/NBDrop_PBMC_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/scryFS_PBMC_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/SeuratDisp_PBMC_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/SeuratMvp_PBMC_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/SeuratVst_PBMC_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/random_PBMC_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    
    

  }
}


for (i in 1:24) {
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/SparSim_results/Mcadet_SparSim_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/SparSim_results/Mcadet_SparSim_fine_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/BreFS_SparSim_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/BreFS_SparSim_fine_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/M3Drop_SparSim_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/M3Drop_SparSim_fine_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/NBDrop_SparSim_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/NBDrop_SparSim_fine_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/scryFS_SparSim_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/scryFS_SparSim_fine_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratDisp_SparSim_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratDisp_SparSim_fine_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratMvp_SparSim_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratMvp_SparSim_fine_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratVst_SparSim_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratVst_SparSim_fine_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/random_SparSim_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/random_SparSim_fine_", i, ".RData"))
  
}



## generate HVGs for rare both SparSim fine and coarse (type 1 cells) datasets ##

rare_SparSim_HVG = c()
for (i in 1:1600) {
  rare_SparSim_HVG<-append(rare_SparSim_HVG, paste0("Gene", i ))
}

rare_SparSim_HVG = rare_SparSim_HVG[c(c(1:20), c(201:400), c(801:1000), c(1401:1600))]




F1_score = function(a,b){
  recall = length(intersect(a,b))/length(a)
  precision = length(intersect(a,b))/length(a)
  return(2/(1/recall + 1/precision))
}


  ## fine 
F1_rare_SparSim_fine_list = c()
for(i in 1:24){
  
a= c(F1_score(rare_SparSim_HVG, get(paste0("Mcadet_SparSim_fine_", i))$gene),
     F1_score(rare_SparSim_HVG, get(paste0("NBDrop_SparSim_fine_", i))),
     F1_score(rare_SparSim_HVG, get(paste0("M3Drop_SparSim_fine_", i))),
     F1_score(rare_SparSim_HVG, get(paste0("BreFS_SparSim_fine_", i))),
     F1_score(rare_SparSim_HVG, get(paste0("scryFS_SparSim_fine_", i))),
     F1_score(rare_SparSim_HVG, get(paste0("SeuratDisp_SparSim_fine_", i))),
     F1_score(rare_SparSim_HVG, get(paste0("SeuratMvp_SparSim_fine_", i))),
     F1_score(rare_SparSim_HVG, get(paste0("SeuratVst_SparSim_fine_", i))),
     F1_score(rare_SparSim_HVG, get(paste0("random_SparSim_fine_", i))))

F1_rare_SparSim_fine_list <- append(F1_rare_SparSim_fine_list, a)

}

Method = rep(c("Mcadet", "NBDrop", "M3Drop", "Brennecke", "Scry", "Seurat Disp",
           "Seurat Mvp", "Seurat Vst", "Random"), 24)

F1_rare_SparSim_fine_df <- 
  data.frame("F1" = F1_rare_SparSim_fine_list, Method)



  ## coarse

F1_rare_SparSim_coarse_list = c()
for(i in 1:24){
  
  a= c(F1_score(rare_SparSim_HVG, get(paste0("Mcadet_SparSim_", i))$gene),
       F1_score(rare_SparSim_HVG, get(paste0("NBDrop_SparSim_", i))),
       F1_score(rare_SparSim_HVG, get(paste0("M3Drop_SparSim_", i))),
       F1_score(rare_SparSim_HVG, get(paste0("BreFS_SparSim_", i))),
       F1_score(rare_SparSim_HVG, get(paste0("scryFS_SparSim_", i))),
       F1_score(rare_SparSim_HVG, get(paste0("SeuratDisp_SparSim_", i))),
       F1_score(rare_SparSim_HVG, get(paste0("SeuratMvp_SparSim_", i))),
       F1_score(rare_SparSim_HVG, get(paste0("SeuratVst_SparSim_", i))),
       F1_score(rare_SparSim_HVG, get(paste0("random_SparSim_", i))))
  
  F1_rare_SparSim_coarse_list <- append(F1_rare_SparSim_coarse_list, a)
  
}

Method = rep(c("Mcadet", "NBDrop", "M3Drop", "Brennecke", "Scry", "Seurat Disp",
               "Seurat Mvp", "Seurat Vst", "Random"), 24)

F1_rare_SparSim_coarse_df <-
  data.frame("F1" = F1_rare_SparSim_coarse_list, Method)






## generate HVGs for rare PBMC fine (dnT cells) datasets ##

F1_rare_PBMC_fine_list = c()
for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    
    load(file = paste0("C:/Users/css22/Desktop/Thesis1/Results/rare_results/rare_PBMC_fine_HVG_", 
                  donor.set[i], "_", time.set[j], ".RData"))
    
    rare_PBMC_HVG = get(paste0("rare_PBMC_fine_HVG_", donor.set[i], "_", time.set[j]))
    
    
    a= c(F1_score(rare_PBMC_HVG, get(paste0("Mcadet_pbmc_", donor.set[i], "_", time.set[j]))$gene),
         F1_score(rare_PBMC_HVG, get(paste0("NBDrop_PBMC_", donor.set[i], "_", time.set[j]))),
         F1_score(rare_PBMC_HVG, get(paste0("M3Drop_PBMC_", donor.set[i], "_", time.set[j]))),
         F1_score(rare_PBMC_HVG, get(paste0("BreFS_PBMC_",  donor.set[i], "_", time.set[j]))),
         F1_score(rare_PBMC_HVG, get(paste0("scryFS_PBMC_", donor.set[i], "_", time.set[j]))),
         F1_score(rare_PBMC_HVG, get(paste0("SeuratDisp_PBMC_", donor.set[i], "_", time.set[j]))),
         F1_score(rare_PBMC_HVG, get(paste0("SeuratMvp_PBMC_", donor.set[i], "_", time.set[j]))),
         F1_score(rare_PBMC_HVG, get(paste0("SeuratVst_PBMC_", donor.set[i], "_", time.set[j]))),
         F1_score(rare_PBMC_HVG, get(paste0("random_PBMC_", donor.set[i], "_", time.set[j]))))
    
    F1_rare_PBMC_fine_list <- append(F1_rare_PBMC_fine_list, a)
    
    
    }
}


F1_rare_PBMC_fine_df <-
  data.frame("F1" = F1_rare_PBMC_fine_list, Method)



#### generate new results for PBMC coarse rare datasets (B cells) #####


for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file = paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmcdata_", donor.set[i], "_", time.set[j], ".RData" ) )   
    load(file = paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmccell1_", donor.set[i], "_", time.set[j], ".RData" ) )
    pbmcdata = get(paste0("pbmcdata_", donor.set[i], "_", time.set[j]))
    pbmclabel = get(paste0("pbmccell1_", donor.set[i], "_", time.set[j]))
    
    which_B <- which(pbmclabel == "B")
    total_B <- length(which(pbmclabel== "B"))

    discard_B <- sample(which_B,size =  total_B - 60, replace = F)

    assign(paste0("pbmc_rareB_", donor.set[i], "_", time.set[j]), 
           pbmcdata[, -discard_B])
    
    save(list = paste0("pbmc_rareB_", donor.set[i], "_", time.set[j]),
    file = paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/pbmc_rareB_",
           donor.set[i], "_", time.set[j]))
    
    rm(pbmcdata,pbmclabel)
    rm(list = paste0("pbmcdata_", donor.set[i], "_", time.set[j]))
    rm(list = paste0("pbmc_rareB_", donor.set[i], "_", time.set[j]) )
    print(c(i,j))
    gc()
    
   }
  }









for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    
    load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/pbmc_rareB_",
                       donor.set[i], "_", time.set[j]))
    
    rare_data = get(paste0("pbmc_rareB_", donor.set[i], "_", time.set[j]))
    result = mcadet(data = rare_data,
                    n.comp= 60, run=10, n.feature=NA, nk_percent = 0.010, 
                    start_resolution = 0.5,  cell_percent =  0.005,
                    MC_iter = 50000, fdr = 0.15, seed = 1234)
    
    assign(paste0("Mcadet_pbmc_rareB_", donor.set[i], "_", time.set[j]),
           result)
    
    
    result_Bre = BrenneckeGetVariableGenes(rare_data, suppress.plot = T, fdr=0.1, minBiolDisp=0.5)

    assign(paste0("BreFS_pbmc_rareB_", donor.set[i], "_", time.set[j]), 
           result_Bre$Gene)
    gc()
    
   
    M3_data <- M3DropConvertData(rare_data, is.counts=TRUE)
    result_M3Drop<- M3DropFeatureSelection(M3_data ,mt_method="fdr", mt_threshold=0.05)
    assign(paste0("M3Drop_pbmc_rareB_", donor.set[i], "_", time.set[j]), result_M3Drop$Gene)
    gc()
    
    
    nb_raw <- NBumiConvertData(rare_data, is.counts=TRUE)
    DANB_fit <- NBumiFitModel(nb_raw)
    result_NBDrop<- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thres=0.05, suppress.plot=F)
    
    assign(paste0("NBDrop_pbmc_rareB_",donor.set[i], "_", time.set[j]), result_NBDrop$Gene)
 
    
    
    set.seed((i-1)*12+(j-1)*4+1)
    result_random = rownames(rare_data)[sample(c(1:nrow(rare_data)), 2000)]
    assign(paste0("random_pbmc_rareB_",donor.set[i], "_", time.set[j]), 
           result_random)
    gc()
    
 
    scry=devianceFeatureSelection(as.matrix(rare_data), fam = "poisson")
    result_scry = names(sort(scry, decreasing = T)[1:2000])
    assign(paste0("scryFS_pbmc_rareB_", donor.set[i], "_", time.set[j]), result_scry)
    gc()

    
    seurat.obj <- CreateAssayObject(counts = rare_data)
    seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat.obj.disp <- FindVariableFeatures(seurat.obj, selection.method = "disp", nfeatures = 2000)
    seurat.obj.mvp <- FindVariableFeatures(seurat.obj, selection.method = "mvp")
    seurat.obj.vst <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
   

    result_SeuratDisp <- VariableFeatures(seurat.obj.disp)
    result_SeuratMvp  <- VariableFeatures(seurat.obj.mvp)
    result_SeuratVst  <- VariableFeatures(seurat.obj.vst)


    assign(paste0("SeuratDisp_pbmc_rareB_", donor.set[i], "_", time.set[j]), result_SeuratDisp)
    assign(paste0("SeuratMvp_pbmc_rareB_", donor.set[i], "_", time.set[j]), result_SeuratMvp)
    assign(paste0("SeuratVst_pbmc_rareB_", donor.set[i], "_", time.set[j]), result_SeuratVst)    
    
    
    save(list = paste0("Mcadet_pbmc_rareB_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/Mcadet_pbmc_rareB_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    
    save(list = paste0("BreFS_pbmc_rareB_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/BreFS_pbmc_rareB_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    
    save(list = paste0("M3Drop_pbmc_rareB_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/M3Drop_pbmc_rareB_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    save(list = paste0("NBDrop_pbmc_rareB_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/NBDrop_pbmc_rareB_",
                  donor.set[i], "_", time.set[j], ".RData"))
   
    save(list = paste0("scryFS_pbmc_rareB_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/scryFS_pbmc_rareB_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    
    save(list = paste0("SeuratDisp_pbmc_rareB_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/SeuratDisp_pbmc_rareB_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    
    save(list = paste0("SeuratMvp_pbmc_rareB_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/SeuratMvp_pbmc_rareB_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    
    save(list = paste0("SeuratVst_pbmc_rareB_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/SeuratVst_pbmc_rareB_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    save(list = paste0("random_pbmc_rareB_", donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/random_pbmc_rareB_",
                  donor.set[i], "_", time.set[j], ".RData"))
    
    
    print(c(i,j))
  }
  
  }  




### load data


for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)

load(file = 
    paste0("C:/Users/css22/Desktop/Thesis1/Results/rare_results/Mcadet_pbmc_rareB_",
    donor.set[i], "_", time.set[j], ".RData"))
load(file = 
       paste0("C:/Users/css22/Desktop/Thesis1/Results/rare_results/Mcadet_pbmc_rareB_",
              donor.set[i], "_", time.set[j], ".RData"))



load(file = 
    paste0("C:/Users/css22/Desktop/Thesis1/Results/rare_results/BreFS_pbmc_rareB_",
              donor.set[i], "_", time.set[j], ".RData"))


load(file = 
    paste0("C:/Users/css22/Desktop/Thesis1/Results/rare_results/M3Drop_pbmc_rareB_",
    donor.set[i], "_", time.set[j], ".RData"))


load(file = 
    paste0("C:/Users/css22/Desktop/Thesis1/Results/rare_results/NBDrop_pbmc_rareB_",
    donor.set[i], "_", time.set[j], ".RData"))

load(file = 
     paste0("C:/Users/css22/Desktop/Thesis1/Results/rare_results/scryFS_pbmc_rareB_",
     donor.set[i], "_", time.set[j], ".RData"))


load(file = 
     paste0("C:/Users/css22/Desktop/Thesis1/Results/rare_results/SeuratDisp_pbmc_rareB_",
     donor.set[i], "_", time.set[j], ".RData"))


load(file = 
     paste0("C:/Users/css22/Desktop/Thesis1/Results/rare_results/SeuratMvp_pbmc_rareB_",
     donor.set[i], "_", time.set[j], ".RData"))


load(file = 
     paste0("C:/Users/css22/Desktop/Thesis1/Results/rare_results/SeuratVst_pbmc_rareB_",
     donor.set[i], "_", time.set[j], ".RData"))

load(file = 
    paste0("C:/Users/css22/Desktop/Thesis1/Results/rare_results/random_pbmc_rareB_",
    donor.set[i], "_", time.set[j], ".RData"))

}}







F1_rare_PBMC_coarse_list = c()
for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    
    load(file = paste0("C:/Users/css22/Desktop/Thesis1/Results/rare_results/rare_PBMC_coarse_HVG_", 
                       donor.set[i], "_", time.set[j], ".RData"))
    
    rare_PBMC_HVG = get(paste0("rare_PBMC_coarse_HVG_", donor.set[i], "_", time.set[j]))
    
    
    a= c(F1_score(rare_PBMC_HVG, get(paste0("Mcadet_pbmc_rareB_", donor.set[i], "_", time.set[j]))$gene),
         F1_score(rare_PBMC_HVG, get(paste0("NBDrop_pbmc_rareB_", donor.set[i], "_", time.set[j]))),
         F1_score(rare_PBMC_HVG, get(paste0("M3Drop_pbmc_rareB_", donor.set[i], "_", time.set[j]))),
         F1_score(rare_PBMC_HVG, get(paste0("BreFS_pbmc_rareB_",  donor.set[i], "_", time.set[j]))),
         F1_score(rare_PBMC_HVG, get(paste0("scryFS_pbmc_rareB_", donor.set[i], "_", time.set[j]))),
         F1_score(rare_PBMC_HVG, get(paste0("SeuratDisp_pbmc_rareB_", donor.set[i], "_", time.set[j]))),
         F1_score(rare_PBMC_HVG, get(paste0("SeuratMvp_pbmc_rareB_", donor.set[i], "_", time.set[j]))),
         F1_score(rare_PBMC_HVG, get(paste0("SeuratVst_pbmc_rareB_", donor.set[i], "_", time.set[j]))),
         F1_score(rare_PBMC_HVG, get(paste0("random_pbmc_rareB_", donor.set[i], "_", time.set[j]))))
    
    F1_rare_PBMC_coarse_list <- append(F1_rare_PBMC_coarse_list, a)
    
    
  }
}

F1_rare_PBMC_coarse_df <-
  data.frame("F1" = F1_rare_PBMC_coarse_list, Method)



### generate barplot for the results




generate_rarefigure_barplot <- function(DF, label, title, y_title){
  


  
  
  library(ggplot2)
  library(ggpubr)
  

  

  
  for (method in c("Seurat Vst",
                   "Seurat Disp",
                   "Seurat Mvp",
                   "NBDrop",
                   "M3Drop",
                   "Brennecke",
                   "Scry",
                   "Random")) {
    group1 <-  DF[DF$Method == "Mcadet",]$F1
    group2 <-  DF[DF$Method == method,]$F1
    
    test_result <- t.test(group1, group2, alternative = "greater")
    
    if (test_result$p.value >= 0.05) {
      DF$Significance[DF$Method == method] <- "NS"
    } else if( test_result$p.value < 0.05 & test_result$p.value >= 0.01 ){
      DF$Significance[DF$Method == method] <- "*"
    }else if( test_result$p.value < 0.01 & test_result$p.value >= 0.001 ){
      DF$Significance[DF$Method == method] <- "**"
    }else if( test_result$p.value < 0.001){
      DF$Significance[DF$Method == method] <- "***"
    }
    
    
  }
  
  
  significance_data <- unique(DF[, c("Method", "Significance")])
  significance_data$y <- 0.65
  significance_data$Significance[1] = ""
  significance_data$mean = c(mean(DF$F1[DF$Method == "Mcadet"]),
                             mean(DF$F1[DF$Method == "NBDrop"]),
                             mean(DF$F1[DF$Method == "M3Drop"]),
                             mean(DF$F1[DF$Method == "Brennecke"]),
                             mean(DF$F1[DF$Method == "Scry"]),
                             mean(DF$F1[DF$Method == "Seurat Disp"]),
                             mean(DF$F1[DF$Method == "Seurat Mvp"]),
                             mean(DF$F1[DF$Method == "Seurat Vst"]),
                             mean(DF$F1[DF$Method == "Random"]))
  
  significance_data$sd = c(sd(DF$F1[DF$Method == "Mcadet"]),
                           sd(DF$F1[DF$Method == "NBDrop"]),
                           sd(DF$F1[DF$Method == "M3Drop"]),
                           sd(DF$F1[DF$Method == "Brennecke"]),
                           sd(DF$F1[DF$Method == "Scry"]),
                           sd(DF$F1[DF$Method == "Seurat Disp"]),
                           sd(DF$F1[DF$Method == "Seurat Mvp"]),
                           sd(DF$F1[DF$Method == "Seurat Vst"]),
                           sd(DF$F1[DF$Method == "Random"]))
  
  
  
  significance_data$Method = factor(c("Mcadet",
                                      "Seurat Vst",
                                      "Seurat Disp",
                                      "Seurat Mvp",
                                      "NBDrop",
                                      "M3Drop",
                                      "Brennecke",
                                      "Scry",
                                      "Random"),
                     levels = c("Mcadet",
                                "Seurat Vst",
                                "Seurat Disp",
                                "Seurat Mvp",
                                "NBDrop",
                                "M3Drop",
                                "Brennecke",
                                "Scry",
                                "Random"))

  
  method_colors <- c("Mcadet" = "#D73027",
                     "Brennecke" = "#E0F3F8", "M3Drop" = "#FFFFBF", "NBDrop" = "#FDAE61", 
                     "Scry" = "#F46D43", "Seurat Disp" = "#A6D96A", "Seurat Mvp" = "#66BD63",
                     "Seurat Vst" = "#1A9850")
  
  
  
  out_figure<- ggplot(significance_data, aes(x = Method, y = mean, fill = Method)) +
    geom_bar(stat = "identity", width = 0.5, color = "black", alpha = 0.8, linewidth = 1) +
    geom_errorbar(aes(ymin = mean - 0.5 * sd, ymax = mean + 0.5 * sd), width = 0.2, color = "black", linewidth = 1) +
    labs(title = title,
         x = "",
         y = y_title) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
      axis.title = element_text(face = "bold", size = 13), 
      axis.text = element_text(face = "bold", size = 12), 
      axis.line = element_line(color = "black"), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.background = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),  
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      legend.title = element_blank(),
      legend.position = "none"
    ) + 
    geom_hline(yintercept = seq(0, 0.5, by = 0.1), linetype = "dashed", color = "black", size = 0.3) + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
    scale_fill_manual(values = method_colors) +
    geom_text(aes(x = Method, y = y, label = Significance), 
              size = 5, color = "black", fontface = "bold") +
    annotate("text", x = -Inf, y = Inf, label = label, hjust = -1, vjust = 1.5, size = 8, fontface = "bold")
  
  
  
  
  
  return(out_figure)
  
}



bar_rare_PBMC_coarse <- generate_rarefigure_barplot(DF = F1_rare_PBMC_coarse_df,
                                                    label = "A",
                                                    title = "PBMC coarse-resolution datasets",
                                                    y_title = "Mean F1 Score")

bar_rare_PBMC_fine <- generate_rarefigure_barplot(DF = F1_rare_PBMC_fine_df,
                                                    label = "B",
                                                    title = "PBMC fine-resolution datasets",
                                                    y_title = "Mean F1 Score")

bar_rare_simulated_coarse <- generate_rarefigure_barplot(DF = F1_rare_SparSim_coarse_df,
                                                    label = "C",
                                                    title = "Simulated coarse-resolution datasets",
                                                    y_title = "Mean F1 Score")

bar_rare_simulated_fine <- generate_rarefigure_barplot(DF = F1_rare_SparSim_fine_df,
                                                  label = "D",
                                                  title = "Simulated fine-resolution datasets",
                                                  y_title = "Mean F1 Score")


library(patchwork)

(bar_rare_PBMC_coarse + bar_rare_PBMC_fine) / (bar_rare_simulated_coarse + bar_rare_simulated_fine)

