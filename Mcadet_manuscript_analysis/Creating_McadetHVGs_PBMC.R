rm(list = ls())
gc()

### Fine resolution 

for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file = 
    paste0("C:/Users/13541/Desktop/Thesis/pbmc/data_fine/pbmc_fine_", 
           donor.set[i], "_", time.set[j], ".RData") )
    result = mcadet(data = get(paste0("pbmc_fine_", donor.set[i], "_", time.set[j])),
           n.comp= 60, run=10, n.feature=NA, nk_percent = 0.010, 
           start_resolution = 0.5,  cell_percent =  0.005,
           MC_iter = 50000, fdr = 0.15, seed = 1234)
    
    assign(paste0("Mcadet_pbmc_fine_",donor.set[i], "_", time.set[j]),
           result)
    
    save(list = paste0("Mcadet_pbmc_fine_",donor.set[i], "_", time.set[j]),
         file = 
paste0("C:/Users/13541/Desktop/Thesis/Generating_results/pbmc_results/Mcadet_pbmc_fine_",
       donor.set[i], "_", time.set[j], ".RData"))
    rm(list =paste0("pbmc_fine_", donor.set[i], "_", time.set[j]) )
    print(c(i,j))
    
  }
}



### 3000 genes

for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file = 
           paste0("C:/Users/13541/Desktop/Thesis/pbmc/data_fine/pbmc_fine_", 
                  donor.set[i], "_", time.set[j], ".RData") )
    result = mcadet(data = get(paste0("pbmc_fine_", donor.set[i], "_", time.set[j])),
                    n.comp= 60, run=10, n.feature=3000, nk_percent = 0.010, 
                    start_resolution = 0.5,  cell_percent =  0.005,
                    MC_iter = 50000, fdr = 0.15, seed = 1234)
    
    assign(paste0("Mcadet_pbmc_fine_3000_",donor.set[i], "_", time.set[j]),
           result)
    
    save(list = paste0("Mcadet_pbmc_fine_3000_",donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Generating_results/pbmc_results/Mcadet_pbmc_fine_3000_",
                  donor.set[i], "_", time.set[j], ".RData"))
    rm(list =paste0("pbmc_fine_", donor.set[i], "_", time.set[j]) )
    print(c(i,j))
    
  }
}








### Coarse resolution 

for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file = 
           paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmcdata_", 
                  donor.set[i], "_", time.set[j], ".RData") )
    result = mcadet(data = get(paste0("pbmcdata_", donor.set[i], "_", time.set[j])),
                    n.comp= 60, run=10, n.feature=NA, nk_percent = 0.010, 
                    start_resolution = 0.5,  cell_percent =  0.005,
                    MC_iter = 50000, fdr = 0.15, seed = 1234)
    
    assign(paste0("Mcadet_pbmc_",donor.set[i], "_", time.set[j]),
           result)
    
    save(list = paste0("Mcadet_pbmc_",donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Generating_results/pbmc_results/Mcadet_pbmc_",
                  donor.set[i], "_", time.set[j], ".RData"))
    rm(list =paste0("pbmcdata_", donor.set[i], "_", time.set[j]) )
    print(c(i,j))
    
  }
}





### 3000 genes

for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file = 
           paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmcdata_", 
                  donor.set[i], "_", time.set[j], ".RData") )
    result = mcadet(data = get(paste0("pbmcdata_", donor.set[i], "_", time.set[j])),
                    n.comp= 60, run=10, n.feature=3000, nk_percent = 0.010, 
                    start_resolution = 0.5,  cell_percent =  0.005,
                    MC_iter = 50000, fdr = 0.15, seed = 1234)
    
    assign(paste0("Mcadet_pbmc_3000_",donor.set[i], "_", time.set[j]),
           result)
    
    save(list = paste0("Mcadet_pbmc_3000_",donor.set[i], "_", time.set[j]),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Generating_results/pbmc_results/Mcadet_pbmc_3000_",
                  donor.set[i], "_", time.set[j], ".RData"))
    rm(list =paste0("pbmcdata_", donor.set[i], "_", time.set[j]) )
    print(c(i,j))
    
  }
}

