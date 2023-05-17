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
           paste0("C:/Users/13541/Desktop/Thesis/pbmc/DE_genes/rank_df_", donor.set[i], "_", time.set[j], ".RData"))
    load(file =
           paste0("C:/Users/13541/Desktop/Thesis/pbmc/DE_genes/rank_df_fine_", donor.set[i], "_", time.set[j], ".RData"))
    
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
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/rare_PBMC_coarse_HVG_", 
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
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/rare_PBMC_fine_HVG_", 
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
     file = paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/rare_PBMC_fine_HVG_", 
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
     file = paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/rare_PBMC_fine_HVG_", 
                   6, "_",  0, ".RData")) 




## generate HVGs for rare both SparSim fine and coarse (type 1 cells) datasets ##

rare_SparSim_HVG = c()
for (i in 1:1600) {
  rare_SparSim_HVG<-append(rare_SparSim_HVG, paste0("Gene", i ))
}

rare_SparSim_HVG = rare_SparSim_HVG[c(c(1:20), c(201:400), c(801:1000), c(1401:1600))]




propotion = function(a,b){
  return(length(intersect(a,b))/length(a))

}

  ## fine 
propotion_rare_SparSim_fine_list = c()
for(i in 1:24){
  
a= c(propotion(rare_SparSim_HVG, get(paste0("Mcadet_SparSim_fine_", i))$gene),
     propotion(rare_SparSim_HVG, get(paste0("NBDrop_SparSim_fine_", i))),
     propotion(rare_SparSim_HVG, get(paste0("M3Drop_SparSim_fine_", i))),
     propotion(rare_SparSim_HVG, get(paste0("BreFS_SparSim_fine_", i))),
     propotion(rare_SparSim_HVG, get(paste0("scryFS_SparSim_fine_", i))),
     propotion(rare_SparSim_HVG, get(paste0("SeuratDisp_SparSim_fine_", i))),
     propotion(rare_SparSim_HVG, get(paste0("SeuratMvp_SparSim_fine_", i))),
     propotion(rare_SparSim_HVG, get(paste0("SeuratVst_SparSim_fine_", i))),
     propotion(rare_SparSim_HVG, get(paste0("random_SparSim_fine_", i))))

propotion_rare_SparSim_fine_list <- append(propotion_rare_SparSim_fine_list, a)

}

Method = rep(c("Mcadet", "NBDrop", "M3Drop", "Brennecke", "Scry", "Seurat Disp",
           "Seurat Mvp", "Seurat Vst", "Random"), 24)

propotion_rare_SparSim_fine_df <- 
  data.frame("proportion" = propotion_rare_SparSim_fine_list, Method)



  ## coarse

propotion_rare_SparSim_coarse_list = c()
for(i in 1:24){
  
  a= c(propotion(rare_SparSim_HVG, get(paste0("Mcadet_SparSim_", i))$gene),
       propotion(rare_SparSim_HVG, get(paste0("NBDrop_SparSim_", i))),
       propotion(rare_SparSim_HVG, get(paste0("M3Drop_SparSim_", i))),
       propotion(rare_SparSim_HVG, get(paste0("BreFS_SparSim_", i))),
       propotion(rare_SparSim_HVG, get(paste0("scryFS_SparSim_", i))),
       propotion(rare_SparSim_HVG, get(paste0("SeuratDisp_SparSim_", i))),
       propotion(rare_SparSim_HVG, get(paste0("SeuratMvp_SparSim_", i))),
       propotion(rare_SparSim_HVG, get(paste0("SeuratVst_SparSim_", i))),
       propotion(rare_SparSim_HVG, get(paste0("random_SparSim_", i))))
  
  propotion_rare_SparSim_coarse_list <- append(propotion_rare_SparSim_coarse_list, a)
  
}

Method = rep(c("Mcadet", "NBDrop", "M3Drop", "Brennecke", "Scry", "Seurat Disp",
               "Seurat Mvp", "Seurat Vst", "Random"), 24)

propotion_rare_SparSim_coarse_df <-
  data.frame("proportion" = propotion_rare_SparSim_coarse_list, Method)






## generate HVGs for rare PBMC fine (dnT cells) datasets ##

propotion_rare_PBMC_fine_list = c()
for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    
    load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/rare_PBMC_fine_HVG_", 
                  donor.set[i], "_", time.set[j], ".RData"))
    
    rare_PBMC_HVG = get(paste0("rare_PBMC_fine_HVG_", donor.set[i], "_", time.set[j]))
    
    
    a= c(propotion(rare_PBMC_HVG, get(paste0("Mcadet_pbmc_", donor.set[i], "_", time.set[j]))$gene),
         propotion(rare_PBMC_HVG, get(paste0("NBDrop_PBMC_", donor.set[i], "_", time.set[j]))),
         propotion(rare_PBMC_HVG, get(paste0("M3Drop_PBMC_", donor.set[i], "_", time.set[j]))),
         propotion(rare_PBMC_HVG, get(paste0("BreFS_PBMC_",  donor.set[i], "_", time.set[j]))),
         propotion(rare_PBMC_HVG, get(paste0("scryFS_PBMC_", donor.set[i], "_", time.set[j]))),
         propotion(rare_PBMC_HVG, get(paste0("SeuratDisp_PBMC_", donor.set[i], "_", time.set[j]))),
         propotion(rare_PBMC_HVG, get(paste0("SeuratMvp_PBMC_", donor.set[i], "_", time.set[j]))),
         propotion(rare_PBMC_HVG, get(paste0("SeuratVst_PBMC_", donor.set[i], "_", time.set[j]))),
         propotion(rare_PBMC_HVG, get(paste0("random_PBMC_", donor.set[i], "_", time.set[j]))))
    
    propotion_rare_PBMC_fine_list <- append(propotion_rare_PBMC_fine_list, a)
    
    
    }
}


propotion_rare_PBMC_fine_df <-
  data.frame("proportion" = propotion_rare_PBMC_fine_list, Method)



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
    paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/Mcadet_pbmc_rareB_",
    donor.set[i], "_", time.set[j], ".RData"))


load(file = 
    paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/BreFS_pbmc_rareB_",
              donor.set[i], "_", time.set[j], ".RData"))


load(file = 
    paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/M3Drop_pbmc_rareB_",
    donor.set[i], "_", time.set[j], ".RData"))


load(file = 
    paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/NBDrop_pbmc_rareB_",
    donor.set[i], "_", time.set[j], ".RData"))

load(file = 
     paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/scryFS_pbmc_rareB_",
     donor.set[i], "_", time.set[j], ".RData"))


load(file = 
     paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/SeuratDisp_pbmc_rareB_",
     donor.set[i], "_", time.set[j], ".RData"))


load(file = 
     paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/SeuratMvp_pbmc_rareB_",
     donor.set[i], "_", time.set[j], ".RData"))


load(file = 
     paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/SeuratVst_pbmc_rareB_",
     donor.set[i], "_", time.set[j], ".RData"))

load(file = 
    paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/random_pbmc_rareB_",
    donor.set[i], "_", time.set[j], ".RData"))

}}







propotion_rare_PBMC_coarse_list = c()
for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    
    load(file = paste0("C:/Users/13541/Desktop/Thesis/Results/rare_results/rare_PBMC_coarse_HVG_", 
                       donor.set[i], "_", time.set[j], ".RData"))
    
    rare_PBMC_HVG = get(paste0("rare_PBMC_coarse_HVG_", donor.set[i], "_", time.set[j]))
    
    
    a= c(propotion(rare_PBMC_HVG, get(paste0("Mcadet_pbmc_rareB_", donor.set[i], "_", time.set[j]))$gene),
         propotion(rare_PBMC_HVG, get(paste0("NBDrop_pbmc_rareB_", donor.set[i], "_", time.set[j]))),
         propotion(rare_PBMC_HVG, get(paste0("M3Drop_pbmc_rareB_", donor.set[i], "_", time.set[j]))),
         propotion(rare_PBMC_HVG, get(paste0("BreFS_pbmc_rareB_",  donor.set[i], "_", time.set[j]))),
         propotion(rare_PBMC_HVG, get(paste0("scryFS_pbmc_rareB_", donor.set[i], "_", time.set[j]))),
         propotion(rare_PBMC_HVG, get(paste0("SeuratDisp_pbmc_rareB_", donor.set[i], "_", time.set[j]))),
         propotion(rare_PBMC_HVG, get(paste0("SeuratMvp_pbmc_rareB_", donor.set[i], "_", time.set[j]))),
         propotion(rare_PBMC_HVG, get(paste0("SeuratVst_pbmc_rareB_", donor.set[i], "_", time.set[j]))),
         propotion(rare_PBMC_HVG, get(paste0("random_pbmc_rareB_", donor.set[i], "_", time.set[j]))))
    
    propotion_rare_PBMC_coarse_list <- append(propotion_rare_PBMC_coarse_list, a)
    
    
  }
}

propotion_rare_PBMC_coarse_df <-
  data.frame("proportion" = propotion_rare_PBMC_coarse_list, Method)



### generate barplot for the results




generate_rarefigure_barplot <- function(DF, plot_title){
  
  
  mean_proportion = c()
  sd_proportion = c()
  
  mean_proportion = append(mean_proportion, mean(DF$proportion[DF$Method == "Mcadet"]))
  sd_proportion = append(sd_proportion, sd(DF$proportion[DF$Method == "Mcadet"]))
  mean_proportion = append(mean_proportion, mean(DF$proportion[DF$Method == "NBDrop"]))
  sd_proportion = append(sd_proportion, sd(DF$proportion[DF$Method == "NBDrop"]))
  mean_proportion = append(mean_proportion, mean(DF$proportion[DF$Method == "M3Drop"]))
  sd_proportion = append(sd_proportion, sd(DF$proportion[DF$Method == "M3Drop"]))
  mean_proportion = append(mean_proportion, mean(DF$proportion[DF$Method == "Brennecke"]))
  sd_proportion = append(sd_proportion, sd(DF$proportion[DF$Method == "Brennecke"]))
  mean_proportion = append(mean_proportion, mean(DF$proportion[DF$Method == "Scry"]))
  sd_proportion = append(sd_proportion, sd(DF$proportion[DF$Method == "Scry"]))
  mean_proportion = append(mean_proportion, mean(DF$proportion[DF$Method == "Seurat Disp"]))
  sd_proportion = append(sd_proportion, sd(DF$proportion[DF$Method == "Seurat Disp"]))
  mean_proportion = append(mean_proportion, mean(DF$proportion[DF$Method == "Seurat Mvp"]))
  sd_proportion = append(sd_proportion, sd(DF$proportion[DF$Method == "Seurat Mvp"]))
  mean_proportion = append(mean_proportion, mean(DF$proportion[DF$Method == "Seurat Vst"]))
  sd_proportion = append(sd_proportion, sd(DF$proportion[DF$Method == "Seurat Vst"]))
  mean_proportion = append(mean_proportion, mean(DF$proportion[DF$Method == "Random"]))
  sd_proportion = append(sd_proportion, sd(DF$proportion[DF$Method == "Random"]))
  
  
  DF_bar = data.frame("mean" = mean_proportion, "sd" = sd_proportion, "Method" = 
                        c("Mcadet", "NBDrop", "M3Drop", "Brennecke", "Scry", "Seurat Disp",
                          "Seurat Mvp", "Seurat Vst", "Random") )
  
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
  
  
  library(ggplot2)
  library(ggpubr)
  
  cbPalette <- c("#cc340c", "#CC79A7", "#56B4E9",
                 "#009E73", "#E69F00", "#F0E442", 
                 "#999999","#0072B2","#D55E00")
  
  out_figure = ggplot(DF_bar, 
                      aes(x= Method, y = mean, fill = Method )) +
    geom_bar(width = 0.5, stat="identity") +
    guides(fill = "none") + geom_errorbar(aes(ymax = mean+sd, ymin= mean-sd),
                                          width = 0.3)+
    labs(x = "", y = "Mean proportion", 
         title = plot_title )+theme_bw() + 
    scale_color_manual(values = cbPalette)+
    geom_hline(yintercept = mean(DF_bar$mean), linewidth = 0.6, color = "blue", linetype=2)+
    theme( panel.grid = element_blank(),  
           axis.title = element_text(face ="bold", size = 16),
           plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
           legend.position = "none",
           axis.text.x = element_text(face = "bold", size = 14, color = "black"), 
           axis.text.y = element_text(face = "bold", size = 12))
  
  return(out_figure)
  
}


generate_rarefigure_barplot(propotion_rare_SparSim_coarse_df, "Simulation coarse datasets")
generate_rarefigure_barplot(propotion_rare_SparSim_fine_df, "Simulation fine datasets")
generate_rarefigure_barplot(propotion_rare_PBMC_coarse_df, "PBMC coarse datasets")
generate_rarefigure_barplot(propotion_rare_PBMC_fine_df, "PBMC fine datasets")


