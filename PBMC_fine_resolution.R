

rm(list = ls())
gc()





fine_types<-c("CD4 Naive", "CD8 Naive", "CD4 TEM", "CD8 TEM", "CD4 Proliferating", "CD8 Proliferating",
              "CD4 TCM", "CD8 TCM", "CD4 CTL", "CD14 Mono", "CD16 Mono")

for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =paste0("E:/Thesis/pbmc/pbmccell1_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/pbmccell2_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/pbmcdata_",  donor.set[i], "_", time.set[j],".RData"))
    
    label_fine<-get(paste0("pbmccell2_", donor.set[i], "_", time.set[j]))
    label_coarse<-get(paste0("pbmccell1_", donor.set[i], "_", time.set[j]))
    data<-get(paste0("pbmcdata_",  donor.set[i], "_", time.set[j]))
    
    for (t in 1:11) {
      label_coarse[which(label_fine==fine_types[t])]=fine_types[t]  
    }
    
    del<-c(which(label_coarse=="CD4 Proliferating"),which(label_coarse=="CD8 Proliferating"),which(label_coarse=="CD4 T"))
    
    label_new=label_coarse[-del]
    data_new=data[,-del]
    
    
    
    assign(paste0("pbmc_fine_", donor.set[i], "_", time.set[j]), data_new)
    assign(paste0("pbmc_fine_clabel_", donor.set[i], "_", time.set[j]), label_new )
    
    save(list = paste0("pbmc_fine_", donor.set[i], "_", time.set[j]), 
         file = paste0("E:/Thesis/pbmc/pbmc_fine_", donor.set[i], "_", time.set[j], ".RData") )
    save(list = paste0("pbmc_fine_clabel_", donor.set[i], "_", time.set[j]), 
         file = paste0("E:/Thesis/pbmc/pbmc_fine_clabel_", donor.set[i], "_", time.set[j], ".RData") )
    
    rm(list = paste0("pbmcdata_", donor.set[i], "_", time.set[j]))
    
    print(c(i,j))
  }
  
}

for (i in 1:8) {
  for(j in 1:3){
    
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =paste0("E:/Thesis/pbmc/pbmc_fine_clabel_", donor.set[i], "_", time.set[j],".RData"))
    print(length(unique( get(paste0("pbmc_fine_clabel_",  donor.set[i], "_", time.set[j]  ))  )))
  }
  
}


