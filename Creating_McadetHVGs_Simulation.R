rm(list = ls())
gc()

### Fine resolution 

for (i in 1:24) {

    load(file = 
           paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_fine_", 
                  i, ".RData") )
    result = mcadet(data = get(paste0("SparSim_fine_", i)),
                    n.comp= 60, run=10, n.feature=NA, nk_percent = 0.010, 
                    start_resolution = 0.5,  cell_percent =  0.005,
                    MC_iter = 50000, fdr = 0.15, seed = 1234)
    
    assign(paste0("Mcadet_SparSim_fine_",i),
           result)
    
    save(list = paste0("Mcadet_SparSim_fine_",i),
         file = 
           paste0("C:/Users/13541/Desktop/Thesis/Generating_results/SparSim_results/Mcadet_SparSim_fine_",
                  i, ".RData"))
    rm(list =paste0("SparSim_fine_", i) )
    print(i)
    
  
}


 ## 3000 genes


for (i in 1:24) {
  
  load(file = 
         paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_fine_", 
                i, ".RData") )
  result = mcadet(data = get(paste0("SparSim_fine_", i)),
                  n.comp= 60, run=10, n.feature=3000, nk_percent = 0.010, 
                  start_resolution = 0.5,  cell_percent =  0.005,
                  MC_iter = 50000, fdr = 0.15, seed = 1234)
  
  assign(paste0("Mcadet_SparSim_fine_3000_",i),
         result)
  
  save(list = paste0("Mcadet_SparSim_fine_3000_",i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Generating_results/SparSim_results/Mcadet_SparSim_fine_3000_",
                i, ".RData"))
  rm(list =paste0("SparSim_fine_", i) )
  print(i)
  
  
}












### Coarse resolution 


for (i in 1:24) {
  
  load(file = 
         paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_", 
                i, ".RData") )
  result = mcadet(data = get(paste0("SparSim_", i)),
                  n.comp= 60, run=10, n.feature=NA, nk_percent = 0.010, 
                  start_resolution = 0.5,  cell_percent =  0.005,
                  MC_iter = 50000, fdr = 0.15, seed = 1234)
  
  assign(paste0("Mcadet_SparSim_",i),
         result)
  
  save(list = paste0("Mcadet_SparSim_",i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Generating_results/SparSim_results/Mcadet_SparSim_",
                i, ".RData"))
  rm(list =paste0("SparSim_", i) )
  print(i)
  
  
}



## 3000 genes


for (i in 1:24) {
  
  load(file = 
         paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_", 
                i, ".RData") )
  result = mcadet(data = get(paste0("SparSim_", i)),
                  n.comp= 60, run=10, n.feature=3000, nk_percent = 0.010, 
                  start_resolution = 0.5,  cell_percent =  0.005,
                  MC_iter = 50000, fdr = 0.15, seed = 1234)
  
  assign(paste0("Mcadet_SparSim_3000_",i),
         result)
  
  save(list = paste0("Mcadet_SparSim_3000_",i),
       file = 
         paste0("C:/Users/13541/Desktop/Thesis/Generating_results/SparSim_results/Mcadet_SparSim_3000_",
                i, ".RData"))
  rm(list =paste0("SparSim_", i) )
  print(i)
  
  
}



