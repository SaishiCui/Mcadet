rm(list = ls())
gc()

library(Seurat)
library(irlba)
library(cccd)

cell.info<-c(rep("Type 1",  60),
             rep("Type 2",  180),
             rep("Type 3",  240),
             rep("Type 4",  360),
             rep("Type 5",  420),
             rep("Type 6",  480),
             rep("Type 7",  600),
             rep("Type 8",  960),
             rep("Type 9",  1200),
             rep("Type 10", 1500))


## evaluation function function  ###


genelist = gene
cell.info = cell.info
gene.name = gene.name
workdata = workdata
n_components= 15 
kmeans.centers= 11
nstart = 10
k = 30

eva<-function(genelist, cell.info, gene.name, workdata, n_components, 
              kmeans.centers, nstart, k){
  
  library("umap")
  library("uwot")
  library("cluster")
  library("factoextra")
  library("funtimes")
  library("mclust")
  library("aricode")
  library("irlba")
  library("cccd")
  
  new_data<-workdata[genelist,]
  
  set.seed(1234)
  seurat.obj<-CreateAssayObject(counts = new_data)
  seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  
  fastpca_data<- as.matrix(seurat.obj@data)
  
  
  pc_df <- as.data.frame(
    prcomp_irlba(fastpca_data, n = n_components,
                 retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
  
  if(sum(1*is.na(pc_df)) >= 1){
    pc_df = t(fastpca_data)
  }else{NA
  }
 
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




### Simulation fine 3456 runs (24 datasets * 9 methods * 4 pcs * 4 k) 



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
    
    for (i in 1:24) {

        

        k.set = c(15,20,25,30)
        pc.set = c(10,15,20,25)
        
        
        
        load(file =paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_fine_", 
                          i, ".RData"))
        workdata=get(paste0("SparSim_fine_", i))
        gene.name<-rownames(workdata)  
        
        genelist.Mcadet<-get(paste0("Mcadet_SparSim_fine_", i))$gene
        genelist.NBDrop<-get(paste0("NBDrop_SparSim_fine_", i))
        genelist.M3Drop<-get(paste0("M3Drop_SparSim_fine_", i))
        genelist.SeuratDisp<-get(paste0("SeuratDisp_SparSim_fine_", i))
        genelist.SeuratMvp<-get(paste0("SeuratMvp_SparSim_fine_", i))
        genelist.SeuratVst<-get(paste0("SeuratVst_SparSim_fine_", i))
        genelist.scryFS<-get(paste0("scryFS_SparSim_fine_", i))
        genelist.BreFS<-get(paste0("BreFS_SparSim_fine_", i))
        genelist.random<-get(paste0("random_SparSim_fine_", i))
        
        
        
        Mcadetlist <- append(Mcadetlist,
                             eva(genelist.Mcadet, cell.info, gene.name, workdata,
                                 n_components=pc.set[p], 
                                 kmeans.centers= 10, 
                                 nstart = 10, 
                                 k = k.set[k]))
        
        NBDroplist <- append(NBDroplist,
                             eva(genelist.NBDrop, cell.info, gene.name, workdata,
                                 n_components=pc.set[p], 
                                 kmeans.centers= 10, 
                                 nstart = 10, 
                                 k = k.set[k]))
        
        M3Droplist <- append(M3Droplist,
                             eva(genelist.M3Drop, cell.info, gene.name, workdata,
                                 n_components=pc.set[p], 
                                 kmeans.centers= 10, 
                                 nstart = 10, 
                                 k = k.set[k]))
        
        
        SeuratDisplist <- append(SeuratDisplist,
                                 eva(genelist.SeuratDisp, cell.info, gene.name, workdata,
                                     n_components=pc.set[p], 
                                     kmeans.centers= 10, 
                                     nstart = 10, 
                                     k = k.set[k]))
        
        
        SeuratVstlist <- append(SeuratVstlist,
                                eva(genelist.SeuratVst, cell.info, gene.name, workdata,
                                    n_components=pc.set[p], 
                                    kmeans.centers= 10, 
                                    nstart = 10, 
                                    k = k.set[k]))
        
        SeuratMvplist <- append(SeuratMvplist,
                                eva(genelist.SeuratMvp, cell.info, gene.name, workdata,
                                    n_components=pc.set[p], 
                                    kmeans.centers= 10, 
                                    nstart = 10, 
                                    k = k.set[k]))
        
        scryFSlist <- append(scryFSlist,
                             eva(genelist.scryFS, cell.info, gene.name, workdata,
                                 n_components=pc.set[p], 
                                 kmeans.centers= 10, 
                                 nstart = 10, 
                                 k = k.set[k]))      
        
        BreFSlist <- append(BreFSlist,
                            eva(genelist.BreFS, cell.info, gene.name, workdata,
                                n_components=pc.set[p], 
                                kmeans.centers= 10, 
                                nstart = 10, 
                                k = k.set[k]))    
        
        randomlist <- append(randomlist,
                             eva(genelist.random, cell.info, gene.name, workdata,
                                 n_components=pc.set[p], 
                                 kmeans.centers= 10, 
                                 nstart = 10, 
                                 k = k.set[k]))      
        
        
        rm(list = paste0("SparSim_fine_", i))
        gc()
        
        
        
        print(c(k,p,i))
        
        
      
    }
    
    assign(paste0("Mcadet_SparSimfineeva_", "PC_", pc.set[p], "_", "k_", k.set[k] ),
           t(matrix(Mcadetlist, nrow = 5, ncol = 24))) 
    
    assign(paste0("NBDrop_SparSimfineeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
           t(matrix(NBDroplist, nrow = 5, ncol = 24))) 
    
    assign(paste0("M3Drop_SparSimfineeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
           t(matrix(M3Droplist, nrow = 5, ncol = 24))) 
    
    
    assign(paste0("SeuratDisp_SparSimfineeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
           t(matrix(SeuratDisplist, nrow = 5, ncol = 24))) 
    
    
    assign(paste0("SeuratMvp_SparSimfineeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
           t(matrix(SeuratMvplist, nrow = 5, ncol = 24))) 
    
    assign(paste0("SeuratVst_SparSimfineeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
           t(matrix(SeuratVstlist, nrow = 5, ncol = 24))) 
    
    assign(paste0("scryFS_SparSimfineeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
           t(matrix(scryFSlist, nrow = 5, ncol = 24))) 
    
    assign(paste0("BreFS_SparSimfineeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
           t(matrix(BreFSlist, nrow = 5, ncol = 24))) 
    
    assign(paste0("random_SparSimfineeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
           t(matrix(randomlist, nrow = 5, ncol = 24))) 
    
    
    
    save(list=paste0("Mcadet_SparSimfineeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/Mcadet_SparSimfineeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("NBDrop_SparSimfineeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/NBDrop_SparSimfineeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    save(list=paste0("M3Drop_SparSimfineeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/M3Drop_SparSimfineeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("SeuratDisp_SparSimfineeva_","PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/SeuratDisp_SparSimfineeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("SeuratMvp_SparSimfineeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/SeuratMvp_SparSimfineeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("SeuratVst_SparSimfineeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/SeuratVst_SparSimfineeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("scryFS_SparSimfineeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/scryFS_SparSimfineeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("BreFS_SparSimfineeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/BreFS_SparSimfineeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("random_SparSimfineeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/random_SparSimfineeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
  }
}


apply(Mcadet_SparSimfineeva_PC_15_k_30, 2, mean)
apply(NBDrop_SparSimfineeva_PC_15_k_30, 2, mean)
apply(M3Drop_SparSimfineeva_PC_15_k_30, 2, mean)
apply(SeuratDisp_SparSimfineeva_PC_15_k_30, 2, mean)
apply(SeuratMvp_SparSimfineeva_PC_15_k_30, 2, mean)
apply(SeuratVst_SparSimfineeva_PC_15_k_30, 2, mean)
apply(BreFS_SparSimfineeva_PC_15_k_30,2,mean)
apply(scryFS_SparSimfineeva_PC_15_k_30, 2, mean)
apply(random_SparSimfineeva_PC_15_k_30, 2, mean)




### SparSim coarse 3456 runs (24 datasets * 9 methods * 4 pcs * 4 k) 



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
    
    for (i in 1:24) {
      
      
      
      k.set = c(15,20,25,30)
      pc.set = c(10,15,20,25)
      
      
      
      load(file =paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_", 
                        i, ".RData"))
      workdata=get(paste0("SparSim_", i))
      gene.name<-rownames(workdata)  
      
      genelist.Mcadet<-get(paste0("Mcadet_SparSim_", i))$gene
      genelist.NBDrop<-get(paste0("NBDrop_SparSim_", i))
      genelist.M3Drop<-get(paste0("M3Drop_SparSim_", i))
      genelist.SeuratDisp<-get(paste0("SeuratDisp_SparSim_", i))
      genelist.SeuratMvp<-get(paste0("SeuratMvp_SparSim_", i))
      genelist.SeuratVst<-get(paste0("SeuratVst_SparSim_", i))
      genelist.scryFS<-get(paste0("scryFS_SparSim_", i))
      genelist.BreFS<-get(paste0("BreFS_SparSim_", i))
      genelist.random<-get(paste0("random_SparSim_", i))
      
      
      
      Mcadetlist <- append(Mcadetlist,
                           eva(genelist.Mcadet, cell.info, gene.name, workdata,
                               n_components=pc.set[p], 
                               kmeans.centers= 10, 
                               nstart = 10, 
                               k = k.set[k]))
      
      NBDroplist <- append(NBDroplist,
                           eva(genelist.NBDrop, cell.info, gene.name, workdata,
                               n_components=pc.set[p], 
                               kmeans.centers= 10, 
                               nstart = 10, 
                               k = k.set[k]))
      
      M3Droplist <- append(M3Droplist,
                           eva(genelist.M3Drop, cell.info, gene.name, workdata,
                               n_components=pc.set[p], 
                               kmeans.centers= 10, 
                               nstart = 10, 
                               k = k.set[k]))
      
      
      SeuratDisplist <- append(SeuratDisplist,
                               eva(genelist.SeuratDisp, cell.info, gene.name, workdata,
                                   n_components=pc.set[p], 
                                   kmeans.centers= 10, 
                                   nstart = 10, 
                                   k = k.set[k]))
      
      
      SeuratVstlist <- append(SeuratVstlist,
                              eva(genelist.SeuratVst, cell.info, gene.name, workdata,
                                  n_components=pc.set[p], 
                                  kmeans.centers= 10, 
                                  nstart = 10, 
                                  k = k.set[k]))
      
      SeuratMvplist <- append(SeuratMvplist,
                              eva(genelist.SeuratMvp, cell.info, gene.name, workdata,
                                  n_components=pc.set[p], 
                                  kmeans.centers= 10, 
                                  nstart = 10, 
                                  k = k.set[k]))
      
      scryFSlist <- append(scryFSlist,
                           eva(genelist.scryFS, cell.info, gene.name, workdata,
                               n_components=pc.set[p], 
                               kmeans.centers= 10, 
                               nstart = 10, 
                               k = k.set[k]))      
      
      BreFSlist <- append(BreFSlist,
                          eva(genelist.BreFS, cell.info, gene.name, workdata,
                              n_components=pc.set[p], 
                              kmeans.centers= 10, 
                              nstart = 10, 
                              k = k.set[k]))    
      
      randomlist <- append(randomlist,
                           eva(genelist.random, cell.info, gene.name, workdata,
                               n_components=pc.set[p], 
                               kmeans.centers= 10, 
                               nstart = 10, 
                               k = k.set[k]))      
      
      
      rm(list = paste0("SparSim_", i))
      gc()
      
      
      
      print(c(k,p,i))
      
      
      
    }
    
    assign(paste0("Mcadet_SparSimcoarseeva_", "PC_", pc.set[p], "_", "k_", k.set[k] ),
           t(matrix(Mcadetlist, nrow = 5, ncol = 24))) 
    
    assign(paste0("NBDrop_SparSimcoarseeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
           t(matrix(NBDroplist, nrow = 5, ncol = 24))) 
    
    assign(paste0("M3Drop_SparSimcoarseeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
           t(matrix(M3Droplist, nrow = 5, ncol = 24))) 
    
    
    assign(paste0("SeuratDisp_SparSimcoarseeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
           t(matrix(SeuratDisplist, nrow = 5, ncol = 24))) 
    
    
    assign(paste0("SeuratMvp_SparSimcoarseeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
           t(matrix(SeuratMvplist, nrow = 5, ncol = 24))) 
    
    assign(paste0("SeuratVst_SparSimcoarseeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
           t(matrix(SeuratVstlist, nrow = 5, ncol = 24))) 
    
    assign(paste0("scryFS_SparSimcoarseeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
           t(matrix(scryFSlist, nrow = 5, ncol = 24))) 
    
    assign(paste0("BreFS_SparSimcoarseeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
           t(matrix(BreFSlist, nrow = 5, ncol = 24))) 
    
    assign(paste0("random_SparSimcoarseeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
           t(matrix(randomlist, nrow = 5, ncol = 24))) 
    
    
    
    save(list=paste0("Mcadet_SparSimcoarseeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/Mcadet_SparSimcoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("NBDrop_SparSimcoarseeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/NBDrop_SparSimcoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    save(list=paste0("M3Drop_SparSimcoarseeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/M3Drop_SparSimcoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("SeuratDisp_SparSimcoarseeva_","PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/SeuratDisp_SparSimcoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("SeuratMvp_SparSimcoarseeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/SeuratMvp_SparSimcoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("SeuratVst_SparSimcoarseeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/SeuratVst_SparSimcoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("scryFS_SparSimcoarseeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/scryFS_SparSimcoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("BreFS_SparSimcoarseeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/BreFS_SparSimcoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
    
    save(list=paste0("random_SparSimcoarseeva_", "PC_", pc.set[p], "_", "k_", k.set[k]),
         file = paste0("C:/Users/13541/Desktop/Thesis/Results/Eva_results/random_SparSimcoarseeva_",
                       "PC_", pc.set[p], "_", "k_", k.set[k], ".RData"))
    
  }
}




## loading data
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/Mcadet_SparSimfineeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/NBDrop_SparSimfineeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/M3Drop_SparSimfineeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/SeuratDisp_SparSimfineeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/SeuratMvp_SparSimfineeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/SeuratVst_SparSimfineeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/scryFS_SparSimfineeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/BreFS_SparSimfineeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/random_SparSimfineeva_PC_15_k_30.RData")

load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/Mcadet_SparSimcoarseeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/NBDrop_SparSimcoarseeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/M3Drop_SparSimcoarseeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/SeuratDisp_SparSimcoarseeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/SeuratMvp_SparSimcoarseeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/SeuratVst_SparSimcoarseeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/scryFS_SparSimcoarseeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/BreFS_SparSimcoarseeva_PC_15_k_30.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Eva_results/random_SparSimcoarseeva_PC_15_k_30.RData")





## boxplot

c(sil.value, purity.value, ari.value, NMI.value, knnretention.value)



generate_figure_simulated <- function(number, input_y, input_title, label, fine,  y_position){
  if(fine){
    Df_simulated = data.frame("Metric" = c(Mcadet_SparSimfineeva_PC_15_k_30[,number],
                                      NBDrop_SparSimfineeva_PC_15_k_30[,number],
                                      M3Drop_SparSimfineeva_PC_15_k_30[,number],
                                      SeuratDisp_SparSimfineeva_PC_15_k_30[,number],
                                      SeuratMvp_SparSimfineeva_PC_15_k_30[,number],
                                      SeuratVst_SparSimfineeva_PC_15_k_30[,number],
                                      scryFS_SparSimfineeva_PC_15_k_30[,number],
                                      BreFS_SparSimfineeva_PC_15_k_30[,number],
                                      random_pSparSimfineeva_PC_15_k_30[,number]),
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
    Df_simulated = data.frame("Metric" = c(Mcadet_SparSimcoarseeva_PC_15_k_30[,number],
                                      NBDrop_SparSimcoarseeva_PC_15_k_30[,number],
                                      M3Drop_SparSimcoarseeva_PC_15_k_30[,number],
                                      SeuratDisp_SparSimcoarseeva_PC_15_k_30[,number],
                                      SeuratMvp_SparSimcoarseeva_PC_15_k_30[,number],
                                      SeuratVst_SparSimcoarseeva_PC_15_k_30[,number],
                                      scryFS_SparSimcoarseeva_PC_15_k_30[,number],
                                      BreFS_SparSimcoarseeva_PC_15_k_30[,number],
                                      random_SparSimcoarseeva_PC_15_k_30[,number]),
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
  
  
  
  
  Df_simulated$Method = factor(Df_simulated$Method,
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
  
  
  df = Df_simulated
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





Box_nhit_Simulated_coarse = generate_figure(fine = F, number = 5,
                                      input_y = "Neighborhood hit",  y_position = 1,
                                      input_title = "Simulated coarse-resolution datasets", label = "C")

Box_nhit_Simulated_fine = generate_figure(fine = T, number = 5,
                                    input_y = "",   y_position = 1,
                                    input_title = "Simulated fine-resolution datasets", label = "D")



(Box_nhit_PBMC_coarse + Box_nhit_PBMC_fine) /
(Box_nhit_Simulated_coarse + Box_nhit_Simulated_fine)  








#### Aggregate all metrics for clustering 



Df_simulated = data.frame("Metric" = c(apply(Mcadet_SparSimfineeva_PC_15_k_30, 1, mean),
                                       apply(NBDrop_SparSimfineeva_PC_15_k_30, 1, mean),
                                       apply(M3Drop_SparSimfineeva_PC_15_k_30, 1, mean),
                                       apply(SeuratDisp_SparSimfineeva_PC_15_k_30, 1, mean),
                                       apply(SeuratMvp_SparSimfineeva_PC_15_k_30, 1, mean),
                                       apply(SeuratVst_SparSimfineeva_PC_15_k_30, 1, mean),
                                       apply(scryFS_SparSimfineeva_PC_15_k_30, 1, mean),
                                       apply(BreFS_SparSimfineeva_PC_15_k_30, 1, mean),
                                       apply(random_SparSimfineeva_PC_15_k_30, 1, mean)),
                          "Method" = c(rep("Mcadet", 24),
                                       rep("NBDrop", 24),
                                       rep("M3Drop", 24),
                                       rep("Brennecke", 24),
                                       rep("Seurat Disp", 24),
                                       rep("Seurat Mvp", 24),
                                       rep("Seurat Vst", 24),
                                       rep("Scry", 24),
                                       rep("Random", 24)))







Df_simulated$Method = factor(Df_simulated$Method,
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


df = Df_simulated
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



Box_ave_simulated_fine <- ggplot(df,aes(Method, Metric, color = Method, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.6, linewidth = 1.1) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 1.3, alpha = 0.9) +
  guides(fill = "none") +
  labs(x = "", y = "" , title = "Simulated fine-resolution datasets" ) +
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
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -1, vjust = 1.5, size = 8, fontface = "bold")+
  geom_text(data = significance_data, aes(x = Method, y = 1, label = Significance), 
            size = 5, color = "black", fontface = "bold") 




##################### Coarse 



Df_simulated = data.frame("Metric" = c(apply(Mcadet_SparSimcoarseeva_PC_15_k_30, 1, mean),
                                       apply(NBDrop_SparSimcoarseeva_PC_15_k_30, 1, mean),
                                       apply(M3Drop_SparSimcoarseeva_PC_15_k_30, 1, mean),
                                       apply(SeuratDisp_SparSimcoarseeva_PC_15_k_30, 1, mean),
                                       apply(SeuratMvp_SparSimcoarseeva_PC_15_k_30, 1, mean),
                                       apply(SeuratVst_SparSimcoarseeva_PC_15_k_30, 1, mean),
                                       apply(scryFS_SparSimcoarseeva_PC_15_k_30, 1, mean),
                                       apply(BreFS_SparSimcoarseeva_PC_15_k_30, 1, mean),
                                       apply(random_SparSimcoarseeva_PC_15_k_30, 1, mean)),
                          "Method" = c(rep("Mcadet", 24),
                                       rep("NBDrop", 24),
                                       rep("M3Drop", 24),
                                       rep("Brennecke", 24),
                                       rep("Seurat Disp", 24),
                                       rep("Seurat Mvp", 24),
                                       rep("Seurat Vst", 24),
                                       rep("Scry", 24),
                                       rep("Random", 24)))



Df_simulated$Method = factor(Df_simulated$Method,
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


df = Df_simulated
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



Box_ave_simulated_coarse <- ggplot(df,aes(Method, Metric, color = Method, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.6, linewidth = 1.1) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 1.3, alpha = 0.9) +
  guides(fill = "none") +
  labs(x = "", y = "Averaged clustering metrics" , title = "Simulated coarse-resolution datasets" ) +
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
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -1, vjust = 1.5, size = 8, fontface = "bold")+
  geom_text(data = significance_data, aes(x = Method, y = 1, label = Significance), 
            size = 5, color = "black", fontface = "bold") 



Box_ave_simulated_coarse + Box_ave_simulated_fine








sum = 0
principal_total = 0
interests_total = 0 
principal_each_year_total_vec = c()
for(y in c(0:29)){
  
principal_each_year = 3400 + 1000 * y 
principal_each_year_FH = 80000*0.07

principal_each_year_total = principal_each_year + principal_each_year_FH


if(principal_each_year_total >= 23000){
  principal_each_year_total = 23000
}
principal_total = principal_total + principal_each_year_total

interests = (principal_each_year_total)*(1.1)^(20-y)
interests_total = interests_total + interests
principal_each_year_total_vec = c(principal_each_year_total_vec, principal_each_year_total )
}



principal_total
interests_total
principal_each_year_total_vec
