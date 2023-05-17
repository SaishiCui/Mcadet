


### annotated PBMC 
remotes::install_github("mojaveazure/seurat-disk")
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(Matrix)



reference <- LoadH5Seurat("C:/Users/13541/Desktop/Thesis/pbmc_multimodal.h5seurat")



unique(reference[["time"]][,1])
length(unique(reference[["time"]][,1]))

unique(reference[["donor"]][,1])
length(unique(reference[["donor"]][,1]))

unique(reference[["lane"]][,1])
length(unique(reference[["lane"]][,1]))

unique(reference[["orig.ident"]][,1])
length(unique(reference[["orig.ident"]][,1]))

unique(reference[["celltype.l1"]][,1])
length(unique(reference[["celltype.l1"]][,1]))

unique(reference[["celltype.l2"]][,1])
length(unique(reference[["celltype.l2"]][,1]))

unique(reference[["celltype.l3"]][,1])
length(unique(reference[["celltype.l3"]][,1]))




for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    
    sample.name=paste0("P",donor.set[i],"_", time.set[j])
    id<-which(reference[["orig.ident"]][,1]==sample.name)
    sample.data=as.data.frame(reference[['SCT']]@counts[,id])
    cell.info.l1<-reference[["celltype.l1"]][id,1]
    cell.info.l2<-reference[["celltype.l2"]][id,1]
    cell.info.l3<-reference[["celltype.l3"]][id,1]
    
    assign( paste0("pbmcdata_", donor.set[i],"_", time.set[j]),  sample.data)
    assign( paste0("pbmccell1_", donor.set[i],"_", time.set[j]), cell.info.l1) 
    assign( paste0("pbmccell2_", donor.set[i],"_", time.set[j]), cell.info.l2)
    assign( paste0("pbmccell3_", donor.set[i],"_", time.set[j]), cell.info.l3) 
    
    
    save(list= paste0("pbmcdata_", donor.set[i],"_", time.set[j]), 
         file=paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmcdata_", donor.set[i], "_", time.set[j], ".RData" ))
    save(list= paste0("pbmccell1_", donor.set[i],"_", time.set[j]), 
         file=paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmccell1_", donor.set[i], "_", time.set[j], ".RData" ))
    save(list= paste0("pbmccell2_", donor.set[i],"_", time.set[j]), 
         file=paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmccell2_", donor.set[i], "_", time.set[j], ".RData" ))
    save(list= paste0("pbmccell3_", donor.set[i],"_", time.set[j]), 
         file=paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmccell3_", donor.set[i], "_", time.set[j], ".RData" ))
    
    
    rm(list=c(paste0("pbmcdata_", donor.set[i],"_", time.set[j]),
              paste0("pbmccell1_", donor.set[i],"_", time.set[j]),
              paste0("pbmccell2_", donor.set[i],"_", time.set[j]),
              paste0("pbmccell3_", donor.set[i],"_", time.set[j])))
    
    gc()
    print(c(i,j))
  }
  
}


for (i in 1:4) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4)
    time.set<-c(0,3,7)
    load(file = paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmcdata_", donor.set[i], "_", time.set[j], ".RData" ) )
    load(file = paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmccell1_", donor.set[i], "_", time.set[j], ".RData" ) )
    load(file = paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmccell2_", donor.set[i], "_", time.set[j], ".RData" ) )
    load(file = paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmccell3_", donor.set[i], "_", time.set[j], ".RData" ) )

  }
}









