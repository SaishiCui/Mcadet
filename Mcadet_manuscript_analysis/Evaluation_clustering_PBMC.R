rm(list = ls())
gc()

library(Seurat)
library(irlba)
library(cccd)

## evaluation function function  ###


eva<-function(genelist, cell.info, gene.name, workdata, n_components, 
              kmeans.centers, nstart, k){
  
  library("umap")
  library("uwot")
  library("cluster")
  library("factoextra")
  library("funtimes")
  library("mclust")
  library("aricode")
  
  label<-NULL
  for (i in 1:length(genelist)) {
    label[i]=which(gene.name==genelist[i]) 
  }
  
  new_data<-workdata[label,]
  
  set.seed(1234)
  seurat.obj<-CreateAssayObject(counts = new_data)
  seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  
  fastpca_data<- as.matrix(seurat.obj@data)
  
  
  pc_df <- as.data.frame(
    prcomp_irlba(fastpca_data, n = n_components,
                 retx = TRUE, center = TRUE, scale. = TRUE)$rotation)

  
  
  partition.pc<-kmeans(pc_df, centers = kmeans.centers, nstart = nstart)
  
  
  sil= silhouette(partition.pc$cluster, dist(pc_df))
  
  sil.mean.cluster<-NULL
  for (m in 1:kmeans.centers) {
    sil.mean.cluster[m]<-mean(sil[sil[,1]==m,3])
  }
  
  
  knn.graph.cells=nng(pc_df, mutual = F, k= k)   ## NNgraph algorithm based on Euclidean distance
  adjacency.matrix <- igraph::as_adjacency_matrix(knn.graph.cells)
  
  knn_ret_list<-NULL
  for (n in 1:length(cell.info)) {
    knn_ret<-mean(cell.info[n]==cell.info[which(adjacency.matrix[n,]==1)])
    knn_ret_list<-append(knn_ret_list, knn_ret)
  }
  
  sil.value<-mean(sil.mean.cluster)
  purity.value<-purity(partition.pc$cluster, as.factor(cell.info))$pur
  ari.value<-adjustedRandIndex(partition.pc$cluster,as.factor(cell.info))
  NMI.value<- NMI(partition.pc$cluster,as.factor(cell.info))
  knnretention.value <- mean(knn_ret_list)
  
  
  return(c(sil.value, purity.value, ari.value, NMI.value, knnretention.value))
}




### PBMC fine 3456 runs (24 datasets * 9 methods * 4 pcs * 4 k) 



for (k in 1:4) {
 for (p in 1:4) {
   
  Mcadetlist <- c()
  NBDroplist <- c()
  M3Droplist <- c()
  SeuratDisplist <- c()
  SeuratMvplist <- c()
  SeuratVstlist <- c()
  scryFSlist <- c()
  BreFSlist <- c()
  randomlist <- c()
  
  for (i in 1:8) {
    for (j in 1:3) {
      
      donor.set<-c(1,2,3,4,5,6,7,8)
      time.set<-c(0,3,7)
      k.set = c(15,20,25,30)
      pc.set = c(10,15,20,25)
      
      

      load(file =paste0("C:/Users/13541/Desktop/Thesis/pbmc/data_fine/pbmc_fine_", 
                        donor.set[i], "_", time.set[j], ".RData"))
      load(file =paste0("C:/Users/13541/Desktop/Thesis/pbmc/data_fine/pbmc_fine_clabel_", 
                        donor.set[i], "_", time.set[j], ".RData"))
      
      cell.info<-get(paste0("pbmc_fine_clabel_", donor.set[i], "_", time.set[j]))
      workdata=get(paste0("pbmc_fine_", donor.set[i], "_", time.set[j]))
      gene.name<-rownames(workdata)  
      
      genelist.Mcadet<-get(paste0("Mcadet_pbmc_fine_", donor.set[i], "_", time.set[j]))$gene
      genelist.NBDrop<-get(paste0("NBDrop_PBMC_fine_", donor.set[i], "_", time.set[j]))
      genelist.M3Drop<-get(paste0("M3Drop_PBMC_fine_", donor.set[i], "_", time.set[j]))
      genelist.SeuratDisp<-get(paste0("SeuratDisp_PBMC_fine_", donor.set[i], "_", time.set[j]))
      genelist.SeuratMvp<-get(paste0("SeuratMvp_PBMC_fine_", donor.set[i], "_", time.set[j]))
      genelist.SeuratVst<-get(paste0("SeuratVst_PBMC_fine_", donor.set[i], "_", time.set[j]))
      genelist.scryFS<-get(paste0("scryFS_PBMC_fine_", donor.set[i], "_", time.set[j]))
      genelist.BreFS<-get(paste0("BreFS_PBMC_fine_", donor.set[i], "_", time.set[j]))
      genelist.random<-get(paste0("random_PBMC_fine_", donor.set[i], "_", time.set[j]))
      


      Mcadetlist <- append(Mcadetlist,
             eva(genelist.Mcadet, cell.info, gene.name, workdata,
                 n_components=pc.set[p], 
                 kmeans.centers= 11, 
                 nstart = 10, 
                 k = k.set[k]))
      
      NBDroplist <- append(NBDroplist,
             eva(genelist.NBDrop, cell.info, gene.name, workdata,
                 n_components=pc.set[p], 
                 kmeans.centers= 11, 
                 nstart = 10, 
                 k = k.set[k]))
      
      M3Droplist <- append(M3Droplist,
                    eva(genelist.M3Drop, cell.info, gene.name, workdata,
                    n_components=pc.set[p], 
                    kmeans.centers= 11, 
                    nstart = 10, 
                    k = k.set[k]))
      

      SeuratDisplist <- append(SeuratDisplist,
             eva(genelist.SeuratDisp, cell.info, gene.name, workdata,
                 n_components=pc.set[p], 
                 kmeans.centers= 11, 
                 nstart = 10, 
                 k = k.set[k]))
      
      
      SeuratVstlist <- append(SeuratVstlist,
             eva(genelist.SeuratVst, cell.info, gene.name, workdata,
                 n_components=pc.set[p], 
                 kmeans.centers= 11, 
                 nstart = 10, 
                 k = k.set[k]))
      
      SeuratMvplist <- append(SeuratMvplist,
             eva(genelist.SeuratMvp, cell.info, gene.name, workdata,
                 n_components=pc.set[p], 
                 kmeans.centers= 11, 
                 nstart = 10, 
                 k = k.set[k]))
      
      scryFSlist <- append(scryFSlist,
             eva(genelist.scryFS, cell.info, gene.name, workdata,
                 n_components=pc.set[p], 
                 kmeans.centers= 11, 
                 nstart = 10, 
                 k = k.set[k]))      

      BreFSlist <- append(BreFSlist,
             eva(genelist.BreFS, cell.info, gene.name, workdata,
                 n_components=pc.set[p], 
                 kmeans.centers= 11, 
                 nstart = 10, 
                 k = k.set[k]))    
      
      randomlist <- append(randomlist,
             eva(genelist.random, cell.info, gene.name, workdata,
                 n_components=pc.set[p], 
                 kmeans.centers= 11, 
                 nstart = 10, 
                 k = k.set[k]))      
      
      
      rm(list = paste0("pbmc_fine_", donor.set[i], "_", time.set[j]))
      gc()
      
      

      print(c(k,p,i,j))

      
   }
  }
  
  assign(paste0("Mcadet_pbmcfineeva_",
                "PC_", pc.set[p], "_", "k_", k.set[k] ),
         t(matrix(Mcadetlist, nrow = 5, ncol = 24))) 
  
  assign(paste0("NBDrop_pbmcfineeva_",
                "PC_", pc.set[p], "_", "k_", k.set[k] ),
         t(matrix(NBDroplist, nrow = 5, ncol = 24))) 
  
  assign(paste0("M3Drop_pbmcfineeva_",
                "PC_", pc.set[p], "_", "k_", k.set[k] ),
         t(matrix(M3Droplist, nrow = 5, ncol = 24))) 
  
  
  assign(paste0("SeuratDisp_pbmcfineeva_",
                "PC_", pc.set[p], "_", "k_", k.set[k] ),
         t(matrix(SeuratDisplist, nrow = 5, ncol = 24))) 
  
  
  assign(paste0("SeuratMvp_pbmcfineeva_",
                "PC_", pc.set[p], "_", "k_", k.set[k] ),
         t(matrix(SeuratMvplist, nrow = 5, ncol = 24))) 
  
  assign(paste0("SeuratVst_pbmcfineeva_",
                "PC_", pc.set[p], "_", "k_", k.set[k] ),
         t(matrix(SeuratVstlist, nrow = 5, ncol = 24))) 
  
  assign(paste0("scryFS_pbmcfineeva_",
                "PC_", pc.set[p], "_", "k_", k.set[k] ),
         t(matrix(scryFSlist, nrow = 5, ncol = 24))) 
  
  assign(paste0("BreFS_pbmcfineeva_",
                "PC_", pc.set[p], "_", "k_", k.set[k] ),
         t(matrix(BreFSlist, nrow = 5, ncol = 24))) 
  
  assign(paste0("random_pbmcfineeva_",
                "PC_", pc.set[p], "_", "k_", k.set[k] ),
         t(matrix(randomlist, nrow = 5, ncol = 24))) 
  
  
  
  save(list=paste0("Mcadet_pbmcfineeva_",
                   "PC_", pc.set[p], "_", "k_", k.set[k] ),
       file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/Mcadet_pbmcfineeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
  
  
  save(list=paste0("NBDrop_pbmcfineeva_",
                   "PC_", pc.set[p], "_", "k_", k.set[k] ),
       file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/NBDrop_pbmcfineeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
  
  save(list=paste0("M3Drop_pbmcfineeva_",
                   "PC_", pc.set[p], "_", "k_", k.set[k] ),
       file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/M3Drop_pbmcfineeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
  
  
  save(list=paste0("SeuratDisp_pbmcfineeva_",
                   "PC_", pc.set[p], "_", "k_", k.set[k] ),
       file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/SeuratDisp_pbmcfineeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
  
  
  save(list=paste0("SeuratMvp_pbmcfineeva_",
                   "PC_", pc.set[p], "_", "k_", k.set[k] ),
       file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/SeuratMvp_pbmcfineeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
  
  
  save(list=paste0("SeuratVst_pbmcfineeva_",
                   "PC_", pc.set[p], "_", "k_", k.set[k] ),
       file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/SeuratVst_pbmcfineeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
  
  
  save(list=paste0("scryFS_pbmcfineeva_",
                   "PC_", pc.set[p], "_", "k_", k.set[k] ),
       file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/scryFS_pbmcfineeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
  
  
  save(list=paste0("BreFS_pbmcfineeva_",
                   "PC_", pc.set[p], "_", "k_", k.set[k] ),
       file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/BreFS_pbmcfineeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
  
  
  save(list=paste0("random_pbmcfineeva_",
                   "PC_", pc.set[p], "_", "k_", k.set[k] ),
       file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/random_pbmcfineeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
  
 }
}







### PBMC coarse 3456 runs (24 datasets * 9 methods * 4 pcs * 4 k) 



for (k in 1:4) {
  for (p in 1:4) {
    
    Mcadetlist <- c()
    NBDroplist <- c()
    M3Droplist <- c()
    SeuratDisplist <- c()
    SeuratMvplist <- c()
    SeuratVstlist <- c()
    scryFSlist <- c()
    BreFSlist <- c()
    randomlist <- c()
    
    for (i in 1:8) {
      for (j in 1:3) {
        
        donor.set<-c(1,2,3,4,5,6,7,8)
        time.set<-c(0,3,7)
        k.set = c(15,20,25,30)
        pc.set = c(10,15,20,25)
        
        
        
        load(file =paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmcdata_", 
                          donor.set[i], "_", time.set[j], ".RData"))
        load(file =paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmccell1_", 
                          donor.set[i], "_", time.set[j], ".RData"))
        
        cell.info<-get(paste0("pbmccell1_", donor.set[i], "_", time.set[j]))
        workdata=get(paste0("pbmcdata_", donor.set[i], "_", time.set[j]))
        gene.name<-rownames(workdata)  
        
        genelist.Mcadet<-get(paste0("Mcadet_pbmc_", donor.set[i], "_", time.set[j]))$gene
        genelist.NBDrop<-get(paste0("NBDrop_PBMC_", donor.set[i], "_", time.set[j]))
        genelist.M3Drop<-get(paste0("M3Drop_PBMC_", donor.set[i], "_", time.set[j]))
        genelist.SeuratDisp<-get(paste0("SeuratDisp_PBMC_", donor.set[i], "_", time.set[j]))
        genelist.SeuratMvp<-get(paste0("SeuratMvp_PBMC_", donor.set[i], "_", time.set[j]))
        genelist.SeuratVst<-get(paste0("SeuratVst_PBMC_", donor.set[i], "_", time.set[j]))
        genelist.scryFS<-get(paste0("scryFS_PBMC_", donor.set[i], "_", time.set[j]))
        genelist.BreFS<-get(paste0("BreFS_PBMC_", donor.set[i], "_", time.set[j]))
        genelist.random<-get(paste0("random_PBMC_", donor.set[i], "_", time.set[j]))
        
        
        
        Mcadetlist <- append(Mcadetlist,
                             eva(genelist.Mcadet, cell.info, gene.name, workdata,
                                 n_components=pc.set[p], 
                                 kmeans.centers= 8, 
                                 nstart = 10, 
                                 k = k.set[k]))
        
        NBDroplist <- append(NBDroplist,
                             eva(genelist.NBDrop, cell.info, gene.name, workdata,
                                 n_components=pc.set[p], 
                                 kmeans.centers= 8, 
                                 nstart = 10, 
                                 k = k.set[k]))
        
        M3Droplist <- append(M3Droplist,
                             eva(genelist.M3Drop, cell.info, gene.name, workdata,
                                 n_components=pc.set[p], 
                                 kmeans.centers= 8, 
                                 nstart = 10, 
                                 k = k.set[k]))
        
        
        SeuratDisplist <- append(SeuratDisplist,
                                 eva(genelist.SeuratDisp, cell.info, gene.name, workdata,
                                     n_components=pc.set[p], 
                                     kmeans.centers= 8, 
                                     nstart = 10, 
                                     k = k.set[k]))
        
        
        SeuratVstlist <- append(SeuratVstlist,
                                eva(genelist.SeuratVst, cell.info, gene.name, workdata,
                                    n_components=pc.set[p], 
                                    kmeans.centers= 8, 
                                    nstart = 10, 
                                    k = k.set[k]))
        
        SeuratMvplist <- append(SeuratMvplist,
                                eva(genelist.SeuratMvp, cell.info, gene.name, workdata,
                                    n_components=pc.set[p], 
                                    kmeans.centers= 8, 
                                    nstart = 10, 
                                    k = k.set[k]))
        
        scryFSlist <- append(scryFSlist,
                             eva(genelist.scryFS, cell.info, gene.name, workdata,
                                 n_components=pc.set[p], 
                                 kmeans.centers= 8, 
                                 nstart = 10, 
                                 k = k.set[k]))      
        
        BreFSlist <- append(BreFSlist,
                            eva(genelist.BreFS, cell.info, gene.name, workdata,
                                n_components=pc.set[p], 
                                kmeans.centers= 8, 
                                nstart = 10, 
                                k = k.set[k]))    
        
        randomlist <- append(randomlist,
                             eva(genelist.random, cell.info, gene.name, workdata,
                                 n_components=pc.set[p], 
                                 kmeans.centers= 8, 
                                 nstart = 10, 
                                 k = k.set[k]))      
        
        
        rm(list = paste0("pbmcdata_", donor.set[i], "_", time.set[j]))
        gc()
        
        
        
        print(c(k,p,i,j))
        
        
      }
    }
    
    assign(paste0("Mcadet_pbmccoarseeva_",
                  "PC_", pc.set[p], "_", "k_", k.set[k] ),
           t(matrix(Mcadetlist, nrow = 5, ncol = 24))) 
    
    assign(paste0("NBDrop_pbmccoarseeva_",
                  "PC_", pc.set[p], "_", "k_", k.set[k] ),
           t(matrix(NBDroplist, nrow = 5, ncol = 24))) 
    
    assign(paste0("M3Drop_pbmccoarseeva_",
                  "PC_", pc.set[p], "_", "k_", k.set[k] ),
           t(matrix(M3Droplist, nrow = 5, ncol = 24))) 
    
    
    assign(paste0("SeuratDisp_pbmccoarseeva_",
                  "PC_", pc.set[p], "_", "k_", k.set[k] ),
           t(matrix(SeuratDisplist, nrow = 5, ncol = 24))) 
    
    
    assign(paste0("SeuratMvp_pbmccoarseeva_",
                  "PC_", pc.set[p], "_", "k_", k.set[k] ),
           t(matrix(SeuratMvplist, nrow = 5, ncol = 24))) 
    
    assign(paste0("SeuratVst_pbmccoarseeva_",
                  "PC_", pc.set[p], "_", "k_", k.set[k] ),
           t(matrix(SeuratVstlist, nrow = 5, ncol = 24))) 
    
    assign(paste0("scryFS_pbmccoarseeva_",
                  "PC_", pc.set[p], "_", "k_", k.set[k] ),
           t(matrix(scryFSlist, nrow = 5, ncol = 24))) 
    
    assign(paste0("BreFS_pbmccoarseeva_",
                  "PC_", pc.set[p], "_", "k_", k.set[k] ),
           t(matrix(BreFSlist, nrow = 5, ncol = 24))) 
    
    assign(paste0("random_pbmccoarseeva_",
                  "PC_", pc.set[p], "_", "k_", k.set[k] ),
           t(matrix(randomlist, nrow = 5, ncol = 24))) 
    
    
    
    save(list=paste0("Mcadet_pbmccoarseeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k] ),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/Mcadet_pbmccoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("NBDrop_pbmccoarseeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k] ),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/NBDrop_pbmccoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    save(list=paste0("M3Drop_pbmccoarseeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k] ),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/M3Drop_pbmccoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("SeuratDisp_pbmccoarseeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k] ),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/SeuratDisp_pbmccoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("SeuratMvp_pbmccoarseeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k] ),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/SeuratMvp_pbmccoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("SeuratVst_pbmccoarseeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k] ),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/SeuratVst_pbmccoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("scryFS_pbmccoarseeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k] ),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/scryFS_pbmccoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("BreFS_pbmccoarseeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k] ),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/BreFS_pbmccoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    save(list=paste0("random_pbmccoarseeva_",
                     "PC_", pc.set[p], "_", "k_", k.set[k] ),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/random_pbmccoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
  }
}




## loading data
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/Mcadet_pbmcfineeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/NBDrop_pbmcfineeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/M3Drop_pbmcfineeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/SeuratDisp_pbmcfineeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/SeuratMvp_pbmcfineeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/SeuratVst_pbmcfineeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/scryFS_pbmcfineeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/BreFS_pbmcfineeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/random_pbmcfineeva_PC_15_k_30.RData")

load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/Mcadet_pbmccoarseeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/NBDrop_pbmccoarseeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/M3Drop_pbmccoarseeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/SeuratDisp_pbmccoarseeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/SeuratMvp_pbmccoarseeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/SeuratVst_pbmccoarseeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/scryFS_pbmccoarseeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/BreFS_pbmccoarseeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/random_pbmccoarseeva_PC_15_k_30.RData")





## boxplot

generate_figure <- function(number, input_y, input_title, label, fine,  y_position){
  if(fine){
    Df_pbmc = data.frame("Metric" = c(Mcadet_pbmcfineeva_PC_15_k_30[,number],
                                      NBDrop_pbmcfineeva_PC_15_k_30[,number],
                                      M3Drop_pbmcfineeva_PC_15_k_30[,number],
                                      SeuratDisp_pbmcfineeva_PC_15_k_30[,number],
                                      SeuratMvp_pbmcfineeva_PC_15_k_30[,number],
                                      SeuratVst_pbmcfineeva_PC_15_k_30[,number],
                                      scryFS_pbmcfineeva_PC_15_k_30[,number],
                                      BreFS_pbmcfineeva_PC_15_k_30[,number],
                                      random_pbmcfineeva_PC_15_k_30[,number]),
                         "Method" = c(rep("Mcadet", 24),
                                      rep("NBDrop", 24),
                                      rep("M3Drop", 24),
                                      rep("Brennecke", 24),
                                      rep("Seurat Disp", 24),
                                      rep("Seurat Mvp", 24),
                                      rep("Seurat Vst", 24),
                                      rep("Scry", 24),
                                      rep("Random", 24)))
    
  }else{
    Df_pbmc = data.frame("Metric" = c(Mcadet_pbmccoarseeva_PC_15_k_30[,number],
                                      NBDrop_pbmccoarseeva_PC_15_k_30[,number],
                                      M3Drop_pbmccoarseeva_PC_15_k_30[,number],
                                      SeuratDisp_pbmccoarseeva_PC_15_k_30[,number],
                                      SeuratMvp_pbmccoarseeva_PC_15_k_30[,number],
                                      SeuratVst_pbmccoarseeva_PC_15_k_30[,number],
                                      scryFS_pbmccoarseeva_PC_15_k_30[,number],
                                      BreFS_pbmccoarseeva_PC_15_k_30[,number],
                                      random_pbmccoarseeva_PC_15_k_30[,number]),
                         "Method" = c(rep("Mcadet", 24),
                                      rep("NBDrop", 24),
                                      rep("M3Drop", 24),
                                      rep("Brennecke", 24),
                                      rep("Seurat Disp", 24),
                                      rep("Seurat Mvp", 24),
                                      rep("Seurat Vst", 24),
                                      rep("Scry", 24),
                                      rep("Random", 24)))
    
  }
  
  
  
  
  Df_pbmc$Method = factor(Df_pbmc$Method,
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
  
  
  df = Df_pbmc
  for (method in c("Seurat Vst",
                   "Seurat Disp",
                   "Seurat Mvp",
                   "NBDrop",
                   "M3Drop",
                   "Brennecke",
                   "Scry",
                   "Random")) {
    group1 <- df[df$Method == "Mcadet",]$Metric
    group2 <- df[df$Method == method,]$Metric
    
    test_result <- t.test(group1, group2, alternative = "greater")
    
    if (test_result$p.value >= 0.05) {
      df$Significance[df$Method == method] <- "NS"
    } else if( test_result$p.value < 0.05 & test_result$p.value >= 0.01 ){
      df$Significance[df$Method == method] <- "*"
    }else if( test_result$p.value < 0.01 & test_result$p.value >= 0.001 ){
      df$Significance[df$Method == method] <- "**"
    }else if( test_result$p.value < 0.001){
      df$Significance[df$Method == method] <- "***"
    }
    
    
  }
  
  
  
  significance_data <- unique(df[, c("Method", "Significance")])
  significance_data$y <-  y_position
  significance_data$Significance[1] = ""
  
  
  
  out_figure <- ggplot(df,aes(Method, Metric, color = Method, fill = Method)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.6, linewidth = 1.1) +
    geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 1.3, alpha = 0.9) +
    guides(fill = "none") +
    labs(x = "", y = input_y , title = input_title ) +
    theme_bw() + 
    scale_color_manual(values = cbPalette) +
    scale_fill_manual(values = cbPalette)+
    geom_hline(yintercept = mean(df$Metric), linetype = 2, linewidth = 1) +
    theme(
      panel.grid = element_blank(),  
      axis.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      legend.position = "none",
      axis.text.x = element_text(face = "bold", size = 13, color = "black"), 
      axis.text.y = element_text(face = "bold", size = 13),
      axis.title.y = element_text(face = "bold", color = "black", size = 15),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.9) 
    )+
    annotate("text", x = -Inf, y = Inf, label = label, hjust = -1, vjust = 1.5, size = 8, fontface = "bold")+
    geom_text(data = significance_data, aes(x = Method, y = y, label = Significance), 
              size = 5, color = "black", fontface = "bold") 
  
  
  
  
  return(out_figure)
}




c(sil.value, purity.value, ari.value, NMI.value, knnretention.value)




Box_nhit_PBMC_coarse = generate_figure(fine = F, number = 5,
                             input_y = "Neighborhood hit",  y_position = 1,
                             input_title = "PBMC coarse-resolution datasets", label = "A")



Box_nhit_PBMC_fine = generate_figure(fine = T, number = 5,
                                      input_y = "",   y_position = 1,
                                      input_title = "PBMC fine-resolution datasets", label = "B")







######## Average all metrics





Df_PBMC = data.frame("Metric" = c(apply(Mcadet_pbmcfineeva_PC_15_k_30, 1, mean),
                                       apply(NBDrop_pbmcfineeva_PC_15_k_30, 1, mean),
                                       apply(M3Drop_pbmcfineeva_PC_15_k_30, 1, mean),
                                       apply(SeuratDisp_pbmcfineeva_PC_15_k_30, 1, mean),
                                       apply(SeuratMvp_pbmcfineeva_PC_15_k_30, 1, mean),
                                       apply(SeuratVst_pbmcfineeva_PC_15_k_30, 1, mean),
                                       apply(scryFS_pbmcfineeva_PC_15_k_30, 1, mean),
                                       apply(BreFS_pbmcfineeva_PC_15_k_30, 1, mean),
                                       apply(random_pbmcfineeva_PC_15_k_30, 1, mean)),
                          "Method" = c(rep("Mcadet", 24),
                                       rep("NBDrop", 24),
                                       rep("M3Drop", 24),
                                       rep("Brennecke", 24),
                                       rep("Seurat Disp", 24),
                                       rep("Seurat Mvp", 24),
                                       rep("Seurat Vst", 24),
                                       rep("Scry", 24),
                                       rep("Random", 24)))







Df_PBMC$Method =      factor(Df_PBMC$Method,
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


df = Df_PBMC 
for (method in c("Seurat Vst",
                 "Seurat Disp",
                 "Seurat Mvp",
                 "NBDrop",
                 "M3Drop",
                 "Brennecke",
                 "Scry",
                 "Random")) {
  group1 <- df[df$Method == "Mcadet",]$Metric
  group2 <- df[df$Method == method,]$Metric
  
  test_result <- t.test(group1, group2, alternative = "greater")
  
  if (test_result$p.value >= 0.05) {
    df$Significance[df$Method == method] <- "NS"
  } else if( test_result$p.value < 0.05 & test_result$p.value >= 0.01 ){
    df$Significance[df$Method == method] <- "*"
  }else if( test_result$p.value < 0.01 & test_result$p.value >= 0.001 ){
    df$Significance[df$Method == method] <- "**"
  }else if( test_result$p.value < 0.001){
    df$Significance[df$Method == method] <- "***"
  }
  
  
}



significance_data <- unique(df[, c("Method", "Significance")])
significance_data$y <-  y_position
significance_data$Significance[1] = ""



Box_ave_PBMC_fine <- ggplot(df,aes(Method, Metric, color = Method, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.6, linewidth = 1.1) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 1.3, alpha = 0.9) +
  guides(fill = "none") +
  labs(x = "", y = "" , title = "PBMC fine-resolution datasets" ) +
  theme_bw() + 
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette)+
  geom_hline(yintercept = mean(df$Metric), linetype = 2, linewidth = 1) +
  theme(
    panel.grid = element_blank(),  
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.position = "none",
    axis.text.x = element_text(face = "bold", size = 13, color = "black"), 
    axis.text.y = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", color = "black", size = 15),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.9) 
  )+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -1, vjust = 1.5, size = 8, fontface = "bold")+
  geom_text(data = significance_data, aes(x = Method, y = 1, label = Significance), 
            size = 5, color = "black", fontface = "bold") 




##################### Coarse 



Df_PBMC = data.frame("Metric" = c(apply(Mcadet_pbmccoarseeva_PC_15_k_30, 1, mean),
                                       apply(NBDrop_pbmccoarseeva_PC_15_k_30, 1, mean),
                                       apply(M3Drop_pbmccoarseeva_PC_15_k_30, 1, mean),
                                       apply(SeuratDisp_pbmccoarseeva_PC_15_k_30, 1, mean),
                                       apply(SeuratMvp_pbmccoarseeva_PC_15_k_30, 1, mean),
                                       apply(SeuratVst_pbmccoarseeva_PC_15_k_30, 1, mean),
                                       apply(scryFS_pbmccoarseeva_PC_15_k_30, 1, mean),
                                       apply(BreFS_pbmccoarseeva_PC_15_k_30, 1, mean),
                                       apply(random_pbmccoarseeva_PC_15_k_30, 1, mean)),
                          "Method" = c(rep("Mcadet", 24),
                                       rep("NBDrop", 24),
                                       rep("M3Drop", 24),
                                       rep("Brennecke", 24),
                                       rep("Seurat Disp", 24),
                                       rep("Seurat Mvp", 24),
                                       rep("Seurat Vst", 24),
                                       rep("Scry", 24),
                                       rep("Random", 24)))



Df_PBMC$Method = factor(Df_PBMC$Method,
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


df = Df_PBMC
for (method in c("Seurat Vst",
                 "Seurat Disp",
                 "Seurat Mvp",
                 "NBDrop",
                 "M3Drop",
                 "Brennecke",
                 "Scry",
                 "Random")) {
  group1 <- df[df$Method == "Mcadet",]$Metric
  group2 <- df[df$Method == method,]$Metric
  
  test_result <- t.test(group1, group2, alternative = "greater")
  
  if (test_result$p.value >= 0.05) {
    df$Significance[df$Method == method] <- "NS"
  } else if( test_result$p.value < 0.05 & test_result$p.value >= 0.01 ){
    df$Significance[df$Method == method] <- "*"
  }else if( test_result$p.value < 0.01 & test_result$p.value >= 0.001 ){
    df$Significance[df$Method == method] <- "**"
  }else if( test_result$p.value < 0.001){
    df$Significance[df$Method == method] <- "***"
  }
  
  
}



significance_data <- unique(df[, c("Method", "Significance")])
significance_data$y <-  y_position
significance_data$Significance[1] = ""



Box_ave_PBMC_coarse <- ggplot(df,aes(Method, Metric, color = Method, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.6, linewidth = 1.1) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 1.3, alpha = 0.9) +
  guides(fill = "none") +
  labs(x = "", y = "Averaged clustering metrics" , title = "PBMC coarse-resolution datasets" ) +
  theme_bw() + 
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette)+
  geom_hline(yintercept = mean(df$Metric), linetype = 2, linewidth = 1) +
  theme(
    panel.grid = element_blank(),  
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.position = "none",
    axis.text.x = element_text(face = "bold", size = 13, color = "black"), 
    axis.text.y = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", color = "black", size = 15),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.9) 
  )+
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -1, vjust = 1.5, size = 8, fontface = "bold")+
  geom_text(data = significance_data, aes(x = Method, y = 1, label = Significance), 
            size = 5, color = "black", fontface = "bold") 


(Box_ave_PBMC_coarse + Box_ave_PBMC_fine) /
(Box_ave_simulated_coarse + Box_ave_simulated_fine)






















### PBMC Linechart of ARI NMI Purity Silhoutte ###


  ## Silhoutte

silmean<-NULL
silsd<-NULL
for (i in 1:6) {
  number.set<-c(100,300,600,1000,1500,2000)
  df<-get(paste0("Silhoutte", "_pbmc_", number.set[i]  , ".df"))
  silmean<-append(silmean,tapply(df$Y,  factor(df$method), mean))
  silsd<-append(silsd,tapply(df$Y,  factor(df$method), sd))
}

line.chart.sil.df<-data.frame(silmean, silsd, "method"=rep(c("Mcadet_fine", "Mcadet", "NBdrop",
                                                             "Brennecke", "Seurat", "Disp", "Random"),6),
                              "n_genes"=rep(c(100,300,600,1000,1500,2000), each=7))


line.chart.sil.df$method<-factor(line.chart.sil.df$method, 
                                 levels = c("Mcadet_fine","Mcadet", "Seurat", "NBdrop", "Disp", "Brennecke", "Random"))


library(ggplot2)

  
a=ggplot(data = line.chart.sil.df, mapping = aes(x = n_genes, y =silmean , 
                                                 colour = method, shape=method,linetype = method))+
  geom_line(size=0.8)+geom_point()+theme_bw()+theme(legend.position="none")+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Mean silhoutte score")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(panel.grid.major = element_blank(),
        panel.background = element_blank(),axis.line = element_blank())+
  scale_x_continuous(limits=c(100,2000), breaks = c(100,300,600,1000,1500,2000))











  ## Purity

puritymean<-NULL
puritysd<-NULL
for (i in 1:6) {
  number.set<-c(100,300,600,1000,1500,2000)
  df<-get(paste0("Purity", "_pbmc_", number.set[i]  , ".df"))
  puritymean<-append(puritymean,tapply(df$Y,  factor(df$method), mean))
  puritysd<-append(puritysd,tapply(df$Y,  factor(df$method), sd))
}

line.chart.purity.df<-data.frame(puritymean, puritysd, "method"=rep(c("Mcadet_fine","Mcadet", "NBdrop",
                                                                      "Brennecke", "Seurat", "Disp", "Random"),6),
                                 "n_genes"=rep(c(100,300,600,1000,1500,2000), each=7))


line.chart.purity.df$method<-factor(line.chart.purity.df$method, 
                                    levels = c("Mcadet_fine","Mcadet","Seurat", "NBdrop", "Disp","Brennecke","Random"))



library(ggplot2)

b=ggplot(data = line.chart.purity.df, mapping = aes(x = n_genes, y =puritymean , 
                                                    colour = method, shape=method,linetype = method))+
  geom_line(size=0.8)+geom_point()+theme_bw()+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Mean purity")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.line = element_blank(),panel.background = element_blank(),
        panel.grid.major = element_blank())+
  scale_x_continuous(limits=c(100,2000), breaks = c(100,300,600,1000,1500,2000))






  ## ARI

arimean<-NULL
arisd<-NULL
for (i in 1:6) {
  number.set<-c(100,300,600,1000,1500,2000)
  df<-get(paste0("ARI", "_pbmc_", number.set[i]  , ".df"))
  arimean<-append(arimean,tapply(df$Y,  factor(df$method), mean))
  arisd<-append(arisd,tapply(df$Y,  factor(df$method), sd))
}

line.chart.ari.df<-data.frame(arimean, arisd, "method"=rep(c("Mcadet_fine","Mcadet", "NBdrop",
                                                             "Brennecke", "Seurat", "Disp", "Random"),6),
                              "n_genes"=rep(c(100,300,600,1000,1500,2000), each=7))


line.chart.ari.df$method<-factor(line.chart.ari.df$method, 
                                 levels = c("Mcadet_fine","Mcadet","Seurat", "NBdrop",  "Disp","Brennecke","Random"))




library(ggplot2)

c=ggplot(data = line.chart.ari.df, mapping = aes(x = n_genes, y =arimean , 
                                                 colour = method, shape=method,linetype = method))+
  geom_line(size=0.8)+geom_point()+theme_bw()+theme(legend.position="none")+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Mean ARI")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(panel.background = element_blank(),axis.line = element_blank(),
        panel.grid.major = element_blank())+
  scale_x_continuous(limits=c(100,2000), breaks = c(100,300,600,1000,1500,2000))







  ## NMI

nmimean<-NULL
nmisd<-NULL
for (i in 1:6) {
  number.set<-c(100,300,600,1000,1500,2000)
  df<-get(paste0("NMI", "_pbmc_", number.set[i]  , ".df"))
  nmimean<-append(nmimean,tapply(df$Y,  factor(df$method), mean))
  nmisd<-append(nmisd,tapply(df$Y,  factor(df$method), sd))
}

line.chart.nmi.df<-data.frame(nmimean, nmisd, "method"=rep(c( "Mcadet_fine","Mcadet", "NBdrop",
                                                              "Brennecke", "Seurat", "Disp", "Random"),6),
                              "n_genes"=rep(c(100,300,600,1000,1500,2000), each=7))


line.chart.nmi.df$method<-factor(line.chart.nmi.df$method, 
                                 levels = c("Mcadet_fine","Mcadet","Seurat", "NBdrop", "Disp", "Brennecke","Random"))



library(ggplot2)

d=ggplot(data = line.chart.nmi.df, mapping = aes(x = n_genes, y =nmimean , 
                                                 colour = method, shape=method,linetype = method))+
  geom_line(size=0.8)+geom_point()+theme_bw()+theme(legend.position="none")+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Mean NMI")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(panel.background = element_blank(), axis.line = element_blank(),
        panel.grid.major = element_blank())+
  scale_x_continuous(limits=c(100,2000), breaks = c(100,300,600,1000,1500,2000))




library('patchwork')

(jac_pbmc_linechart_p+knn.only)/(a+b)/(c+d)








