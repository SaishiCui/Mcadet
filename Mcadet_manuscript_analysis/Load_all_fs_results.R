rm(list = ls())
gc()


## load SparSim

for (i in 1:24) {
  
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/NBDrop_SparSim_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/NBDrop_SparSim_3000_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/NBDrop_SparSim_fine_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/NBDrop_SparSim_fine_3000_",
                     i , ".RData"))
  
  
  
  
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/M3Drop_SparSim_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/M3Drop_SparSim_3000_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/M3Drop_SparSim_fine_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/M3Drop_SparSim_fine_3000_",
                     i , ".RData"))
  
  
  
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/BreFS_SparSim_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/BreFS_SparSim_fine_",
                     i , ".RData"))
  
  
  
  
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratVst_SparSim_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratVst_SparSim_3000_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratVst_SparSim_fine_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratVst_SparSim_fine_3000_",
                     i , ".RData"))
  
  
  
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratDisp_SparSim_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratDisp_SparSim_3000_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratDisp_SparSim_fine_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratDisp_SparSim_fine_3000_",
                     i , ".RData"))
  
  
  
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratMvp_SparSim_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/SeuratMvp_SparSim_fine_",
                     i , ".RData"))

  
  
  
  
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/scryFS_SparSim_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/scryFS_SparSim_3000_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/scryFS_SparSim_fine_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/scryFS_SparSim_fine_3000_",
                     i , ".RData"))
  
  
  
  
  
  
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/random_SparSim_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/random_SparSim_3000_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/random_SparSim_fine_",
                     i , ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_SparSim/random_SparSim_fine_3000_",
                     i , ".RData"))
  
  
  
  load(file=paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/SparSim_results/Mcadet_SparSim_", i , ".RData"))
  load(file=paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/SparSim_results/Mcadet_SparSim_3000_", i, ".RData" ))
  load(file=paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/SparSim_results/Mcadet_SparSim_fine_", i, ".RData" ))
  load(file=paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/SparSim_results/Mcadet_SparSim_fine_3000_", i, ".RData" ))
  
  

  
  print(i)
  
  
  
  
}











## load PBMC

for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
  
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/NBDrop_PBMC_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/NBDrop_PBMC_3000_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/NBDrop_PBMC_fine_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/NBDrop_PBMC_fine_3000_",
                     donor.set[i], "_", time.set[j], ".RData"))
  
  
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/M3Drop_PBMC_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/M3Drop_PBMC_3000_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/M3Drop_PBMC_fine_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/M3Drop_PBMC_fine_3000_",
                     donor.set[i], "_", time.set[j], ".RData"))
  
  
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/BreFS_PBMC_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/BreFS_PBMC_fine_",
                     donor.set[i], "_", time.set[j], ".RData"))

  
  
  
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/SeuratVst_PBMC_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/SeuratVst_PBMC_3000_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/SeuratVst_PBMC_fine_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/SeuratVst_PBMC_fine_3000_",
                     donor.set[i], "_", time.set[j], ".RData"))
  
  
  
  
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/SeuratDisp_PBMC_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/SeuratDisp_PBMC_3000_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/SeuratDisp_PBMC_fine_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/SeuratDisp_PBMC_fine_3000_",
                     donor.set[i], "_", time.set[j], ".RData"))
  
  
  
  
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/SeuratMvp_PBMC_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/SeuratMvp_PBMC_fine_",
                     donor.set[i], "_", time.set[j], ".RData"))

  
  
  
  
  
  
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/scryFS_PBMC_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/scryFS_PBMC_3000_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/scryFS_PBMC_fine_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/scryFS_PBMC_fine_3000_",
                     donor.set[i], "_", time.set[j], ".RData"))
  
  
  
  
  
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/random_PBMC_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/random_PBMC_3000_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/random_PBMC_fine_",
                     donor.set[i], "_", time.set[j], ".RData"))
  load(file = paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_pbmc/random_PBMC_fine_3000_",
                     donor.set[i], "_", time.set[j], ".RData"))
  
  
  
  load(file=paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/pbmc_results/Mcadet_pbmc_", donor.set[i], "_", time.set[j] , ".RData"))
  load(file=paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/pbmc_results/Mcadet_pbmc_3000_", donor.set[i], "_", time.set[j], ".RData" ))
  load(file=paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/pbmc_results/Mcadet_pbmc_fine_", donor.set[i], "_", time.set[j], ".RData" ))
  load(file=paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/pbmc_results/Mcadet_pbmc_fine_3000_", donor.set[i], "_", time.set[j], ".RData" ))
  
  
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/DE_genes/DEgene_fine_", 
                    donor.set[i], "_", time.set[j], ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/DE_genes/DEgene_", 
                    donor.set[i], "_", time.set[j], ".RData"))
  
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/data_fine/pbmc_fine_clabel_", 
                    donor.set[i], "_", time.set[j], ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/data/pbmccell1_", 
                    donor.set[i], "_", time.set[j], ".RData"))
  print(c(i,j))
  
  
  
  
  }
}








