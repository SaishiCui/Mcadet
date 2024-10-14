

## PBMC fine data ###

library(Seurat)

for (i in 1:8) {
  for (j in 1:3) {

donor.set<-c(1,2,3,4,5,6,7,8)
time.set<-c(0,3,7)

load(file =paste0("C:/Users/13541/Desktop/Thesis/pbmc/data_fine/pbmc_fine_",  donor.set[i], "_", time.set[j],".RData"))
load(file =paste0("C:/Users/13541/Desktop/Thesis/pbmc/data_fine/pbmc_fine_clabel_", donor.set[i], "_", time.set[j],".RData"))


rawdata<-get(paste0("pbmc_fine_", donor.set[i], "_", time.set[j]))

working.data = rawdata[rowSums(rawdata)!=0,] 

pbmcdata_seurat <- CreateSeuratObject(counts =working.data)
pbmcdata_seurat <- NormalizeData(pbmcdata_seurat, normalization.method = "LogNormalize", scale.factor = 10000)


all.genes <- rownames(rawdata)
pbmcdata_seurat <- ScaleData(pbmcdata_seurat, features = all.genes)
Idents(pbmcdata_seurat)<-factor(get(paste0("pbmc_fine_clabel_", donor.set[i], "_", time.set[j])))


all.markers <- FindAllMarkers(object =pbmcdata_seurat, return.thresh = 1, logfc.threshold = 0)

all.markders.newrank<-all.markers[order(abs(all.markers$avg_log2FC),decreasing = T),]


assign(paste0("rank_df_fine_", donor.set[i], "_", time.set[j]), all.markders.newrank)


save(list = paste0("rank_df_fine_", donor.set[i], "_", time.set[j]),  
     file = paste0("C:/Users/13541/Desktop/Thesis/pbmc/DE_genes/rank_df_fine_", donor.set[i], "_", time.set[j], ".RData"))


rm(list = paste0("pbmc_fine_", donor.set[i], "_", time.set[j]))
gc()
print(c(i,j))
 } 

}



  
## PBMC coarse data ###

library(Seurat)

for (i in 1:8) {
  for (j in 1:3) {
    
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    
    load(file =paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmcdata_",  donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmccell1_", donor.set[i], "_", time.set[j],".RData"))
    
    
    rawdata<-get(paste0("pbmcdata_", donor.set[i], "_", time.set[j]))
    
    working.data = rawdata[rowSums(rawdata)!=0,] 
    
    pbmcdata_seurat <- CreateSeuratObject(counts =working.data)
    pbmcdata_seurat <- NormalizeData(pbmcdata_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
    
    
    all.genes <- rownames(rawdata)
    pbmcdata_seurat <- ScaleData(pbmcdata_seurat, features = all.genes)
    Idents(pbmcdata_seurat)<-factor(get(paste0("pbmccell1_", donor.set[i], "_", time.set[j])))
    
    
    all.markers <- FindAllMarkers(object =pbmcdata_seurat, return.thresh = 1, logfc.threshold = 0)
    
    all.markders.newrank<-all.markers[order(abs(all.markers$avg_log2FC),decreasing = T),]
    
    
    assign(paste0("rank_df_", donor.set[i], "_", time.set[j]), all.markders.newrank)
    
    
    save(list = paste0("rank_df_", donor.set[i], "_", time.set[j]),  
         file = paste0("C:/Users/13541/Desktop/Thesis/pbmc/DE_genes/rank_df_", donor.set[i], "_", time.set[j], ".RData"))
    
    
    rm(list = paste0("pbmcdata_", donor.set[i], "_", time.set[j]))
    gc()
    print(c(i,j))
  } 
  
}









### get DE genes ###

for (i in 1:8) {
  for (j in 1:3) {
    
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    
    rank_df<-get(paste0("rank_df_fine_", donor.set[i], "_", time.set[j]))
    
    num_clusters<-length(unique(rank_df$cluster))
    
    extract_DEgene_list<-NULL

    for (t in 1:num_clusters) {
    cluster_label<-which(rank_df$cluster==unique(rank_df$cluster)[t])
    extract_rank_df<-rank_df[cluster_label,]
    order_by_padjust<-order(extract_rank_df[,"p_val_adj"],decreasing = F)

    len<-sum(1*(extract_rank_df[order_by_padjust,"p_val_adj"]<0.05))
    sig_label<-which(extract_rank_df[order_by_padjust,"p_val_adj"]<0.05)
    
    extract_DEgene<-extract_rank_df[order_by_padjust[sig_label],"gene"]

    
 
    extract_DEgene_list<-append(extract_DEgene_list, extract_DEgene)
    
    }
    
    DEgene_clean=unique(extract_DEgene_list)

    
    assign(paste0("DEgene_fine_", donor.set[i], "_", time.set[j]) , DEgene_clean)
    
    save(list = paste0("DEgene_fine_", donor.set[i], "_", time.set[j]),  
         file = paste0("C:/Users/13541/Desktop/Thesis/pbmc/DE_genes/DEgene_fine_", donor.set[i], "_", time.set[j], ".RData"))
    
    
    
  }
  
}


## coarse 


for (i in 1:8) {
  for (j in 1:3) {
    
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    
    rank_df<-get(paste0("rank_df_", donor.set[i], "_", time.set[j]))
    
    num_clusters<-length(unique(rank_df$cluster))
    
    extract_DEgene_list<-NULL
    
    for (t in 1:num_clusters) {
      cluster_label<-which(rank_df$cluster==unique(rank_df$cluster)[t])
      extract_rank_df<-rank_df[cluster_label,]
      order_by_padjust<-order(extract_rank_df[,"p_val_adj"],decreasing = F)
      
      len<-sum(1*(extract_rank_df[order_by_padjust,"p_val_adj"]<0.05))
      sig_label<-which(extract_rank_df[order_by_padjust,"p_val_adj"]<0.05)
      
      if(len<400){
      extract_DEgene<-extract_rank_df[order_by_padjust[sig_label],"gene"]
      }else{
        extract_DEgene<-extract_rank_df[order_by_padjust[sig_label],"gene"][1:400]
      }
      
      
      extract_DEgene_list<-append(extract_DEgene_list, extract_DEgene)
      
    }
    
    DEgene_clean=unique(extract_DEgene_list)

    
    
    assign(paste0("DEgene_", donor.set[i], "_", time.set[j]) , DEgene_clean)
    
    save(list = paste0("DEgene_", donor.set[i], "_", time.set[j]),  
         file = paste0("C:/Users/13541/Desktop/Thesis/pbmc/DE_genes/DEgene_", donor.set[i], "_", time.set[j], ".RData"))
    
    
    
  }
  
}






