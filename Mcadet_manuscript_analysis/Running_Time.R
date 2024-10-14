rm(list = ls())
gc()

library(M3Drop)
library(M3DExampleData)
library(dplyr)
library(Seurat)
library(scry)

library("SPARSim")
data(Zheng_param_preset)

sparsim_gene_name<-NULL
for (i in 1:15000) {
  sparsim_gene_name<-append(sparsim_gene_name, paste0("Gene", i ))
}

sparsim_cell_name<-NULL
for (i in 1:20000) {
  sparsim_cell_name<-append(sparsim_cell_name, paste0("Cell", i ))
}





downreg.min.set<- seq(0.1, 0.6, by=0.5/23)
downreg.max.set<- seq(0.3, 0.8, by=0.5/23)
upreg.min.set<- rev(seq(1, 3, by=2/23))
upreg.max.set<- rev(seq(1.5, 3.5, by=2/23))

i = 12
j = 1

set.seed((j-1)*48+1 )
DE_multiplier_1.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2)
set.seed((j-1)*48+2 )
DE_multiplier_1.2 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
set.seed( (j-1)*48+3 )
DE_multiplier_1.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
set.seed( (j-1)*48+4 )
DE_multiplier_1.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
fold_change_multiplier_1 <- c(DE_multiplier_1.1, rep(1, 180), DE_multiplier_1.2, rep(1,400),
                              DE_multiplier_1.3, rep(1,400), DE_multiplier_1.4, rep(1,17936) )


set.seed((j-1)*48+5 )
DE_multiplier_2.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
set.seed((j-1)*48+6 )
DE_multiplier_2.2 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
set.seed( (j-1)*48+7 )
DE_multiplier_2.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
set.seed( (j-1)*48+8 )
DE_multiplier_2.4 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
set.seed( (j-1)*48+9 )
DE_multiplier_2.5 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
fold_change_multiplier_2 <- c(rep(1,20), DE_multiplier_2.1, rep(1,360), DE_multiplier_2.2,
                              rep(1,200), DE_multiplier_2.3, rep(1,200), DE_multiplier_2.4, 
                              rep(1,200), DE_multiplier_2.5, rep(1,17736) )


set.seed((j-1)*48+10 )
DE_multiplier_3.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
set.seed((j-1)*48+11 )
DE_multiplier_3.2 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
set.seed( (j-1)*48+12 )
DE_multiplier_3.3 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
set.seed( (j-1)*48+13 )
DE_multiplier_3.4 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
set.seed( (j-1)*48+14 )
DE_multiplier_3.5 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
fold_change_multiplier_3 <- c(rep(1,40), DE_multiplier_3.1, rep(1,140), DE_multiplier_3.2,
                              DE_multiplier_3.3, rep(1,800), DE_multiplier_3.4, 
                              rep(1,200), DE_multiplier_3.5, rep(1,17536) )



set.seed((j-1)*48+15 )
DE_multiplier_4.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
set.seed((j-1)*48+16 )
DE_multiplier_4.2 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
set.seed( (j-1)*48+17 )
DE_multiplier_4.3 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
set.seed( (j-1)*48+18 )
DE_multiplier_4.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
set.seed( (j-1)*48+19 )
DE_multiplier_4.5 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
set.seed( (j-1)*48+20 )
DE_multiplier_4.6 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
fold_change_multiplier_4 <- c(rep(1,60), DE_multiplier_4.1, rep(1,320), DE_multiplier_4.2, rep(1,400),
                              DE_multiplier_4.3, DE_multiplier_4.4, DE_multiplier_4.5, 
                              rep(1,200), DE_multiplier_4.6, rep(1,17536) )


set.seed((j-1)*48+21 )
DE_multiplier_5.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
set.seed((j-1)*48+22 )
DE_multiplier_5.2 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
set.seed( (j-1)*48+23 )
DE_multiplier_5.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
set.seed( (j-1)*48+24 )
DE_multiplier_5.4 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
set.seed( (j-1)*48+25 )
DE_multiplier_5.5 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
fold_change_multiplier_5 <- c(rep(1,80), DE_multiplier_5.1, rep(1,100), DE_multiplier_5.2,
                              rep(1,200), DE_multiplier_5.3, rep(1,800), DE_multiplier_5.4, 
                              DE_multiplier_5.5, rep(1,17536) )




set.seed((j-1)*48+26 )
DE_multiplier_6.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
set.seed((j-1)*48+27 )
DE_multiplier_6.2 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
set.seed( (j-1)*48+28 )
DE_multiplier_6.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
set.seed( (j-1)*48+29 )
DE_multiplier_6.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
set.seed( (j-1)*48+30 )
DE_multiplier_6.5 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
fold_change_multiplier_6 <- c(rep(1,100), DE_multiplier_6.1, rep(1,280), DE_multiplier_6.2,
                              rep(1,400), DE_multiplier_6.3, rep(1,200), DE_multiplier_5.4, 
                              rep(1,200), DE_multiplier_6.5, rep(1,17536) )



set.seed((j-1)*48+31 )
DE_multiplier_7.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
set.seed((j-1)*48+32 )
DE_multiplier_7.2 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
set.seed( (j-1)*48+33 )
DE_multiplier_7.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
set.seed( (j-1)*48+34 )
DE_multiplier_7.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
fold_change_multiplier_7 <- c(rep(1,120), DE_multiplier_7.1, rep(1,60), DE_multiplier_7.2,
                              rep(1,200), DE_multiplier_7.3, rep(1,400), DE_multiplier_7.4, rep(1,18136) )



set.seed((j-1)*48+35 )
DE_multiplier_8.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
set.seed((j-1)*48+36 )
DE_multiplier_8.2 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
set.seed( (j-1)*48+37 )
DE_multiplier_8.3 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
set.seed( (j-1)*48+38 )
DE_multiplier_8.4 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
fold_change_multiplier_8 <- c(rep(1,140), DE_multiplier_8.1, rep(1,440), DE_multiplier_8.2,
                              rep(1,200), DE_multiplier_8.3, rep(1,400), DE_multiplier_8.4, rep(1,17736) )





set.seed((j-1)*48+39 )
DE_multiplier_9.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
set.seed((j-1)*48+40 )
DE_multiplier_9.2 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
set.seed( (j-1)*48+41 )
DE_multiplier_9.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
set.seed( (j-1)*48+42 )
DE_multiplier_9.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
set.seed( (j-1)*48+43 )
DE_multiplier_9.5 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
fold_change_multiplier_9 <- c(rep(1,160), DE_multiplier_9.1, rep(1,20), DE_multiplier_9.2,
                              rep(1,800), DE_multiplier_9.3, rep(1,200), DE_multiplier_9.4, 
                              DE_multiplier_9.5, rep(1,17536) )




set.seed((j-1)*48+44 )
DE_multiplier_10.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
set.seed((j-1)*48+45 )
DE_multiplier_10.2 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
set.seed( (j-1)*48+46 )
DE_multiplier_10.3 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
set.seed( (j-1)*48+47 )
DE_multiplier_10.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
set.seed( (j-1)*48+48 )
DE_multiplier_10.5 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
fold_change_multiplier_10 <- c(rep(1,180), DE_multiplier_9.1, rep(1,600), DE_multiplier_10.2,
                               DE_multiplier_10.3, rep(1,200), DE_multiplier_10.4, 
                               DE_multiplier_10.5, rep(1,17736) )





cell_type_1 <- SPARSim_create_DE_genes_parameter(
  sim_param =Zheng_param_preset$Zheng_C1, 
  fc_multiplier = fold_change_multiplier_1, 
  N_cells = 200,
  condition_name = "cell_1")


cell_type_2 <- SPARSim_create_DE_genes_parameter(
  sim_param =Zheng_param_preset$Zheng_C1, 
  fc_multiplier = fold_change_multiplier_2, 
  N_cells = 600,
  condition_name = "cell_2")



cell_type_3 <- SPARSim_create_DE_genes_parameter(
  sim_param =Zheng_param_preset$Zheng_C1 , 
  fc_multiplier = fold_change_multiplier_3, 
  N_cells = 800,
  condition_name = "cell_3")


cell_type_4 <- SPARSim_create_DE_genes_parameter(
  sim_param =Zheng_param_preset$Zheng_C1 , 
  fc_multiplier = fold_change_multiplier_4, 
  N_cells = 1200,
  condition_name = "cell_4")


cell_type_5 <- SPARSim_create_DE_genes_parameter(
  sim_param =Zheng_param_preset$Zheng_C1 , 
  fc_multiplier = fold_change_multiplier_5, 
  N_cells = 1400,
  condition_name = "cell_5")


cell_type_6 <- SPARSim_create_DE_genes_parameter(
  sim_param =Zheng_param_preset$Zheng_C1 , 
  fc_multiplier = fold_change_multiplier_6, 
  N_cells = 1600,
  condition_name = "cell_6")


cell_type_7 <- SPARSim_create_DE_genes_parameter(
  sim_param =Zheng_param_preset$Zheng_C1 , 
  fc_multiplier = fold_change_multiplier_7, 
  N_cells = 2000,
  condition_name = "cell_7")


cell_type_8 <- SPARSim_create_DE_genes_parameter(
  sim_param =Zheng_param_preset$Zheng_C1 , 
  fc_multiplier = fold_change_multiplier_8, 
  N_cells = 3200,
  condition_name = "cell_8")



cell_type_9 <- SPARSim_create_DE_genes_parameter(
  sim_param =Zheng_param_preset$Zheng_C1 , 
  fc_multiplier = fold_change_multiplier_9, 
  N_cells = 4000,
  condition_name = "cell_9")




cell_type_10 <- SPARSim_create_DE_genes_parameter(
  sim_param =Zheng_param_preset$Zheng_C1 , 
  fc_multiplier = fold_change_multiplier_10, 
  N_cells = 5000,
  condition_name = "cell_10")


SPARSim_param_with_DE <- list(cell_type_1=cell_type_1,
                              cell_type_2=cell_type_2,
                              cell_type_3=cell_type_3,   
                              cell_type_4=cell_type_4,
                              cell_type_5=cell_type_5,
                              cell_type_6=cell_type_6,
                              cell_type_7=cell_type_7,   
                              cell_type_8=cell_type_8,
                              cell_type_9=cell_type_9,
                              cell_type_10=cell_type_10)

## Can not be seeded (SPARSIM)

SPARSim_result <- SPARSim_simulation(SPARSim_param_with_DE )

set.seed(1234*i*j)
ind<-sample(c(2001:19536), size = 4536 ,replace = F)

sparsim.data=SPARSim_result$count_matrix[-ind,]

rownames(sparsim.data)=sparsim_gene_name
colnames(sparsim.data)=sparsim_cell_name
    






time_cell_mcadet = c()

for(i in 1:20){
  set.seed(i*1234)
  idx = sample(1:20000, size = 1000*i)  
  b = sparsim.data[, idx]
  
  start_time = Sys.time()

  a = mcadet(data = b,
                n.comp= 60, run=10, n.feature=NA, nk_percent = 0.010, 
                start_resolution = 0.5,  cell_percent =  0.005,
                MC_iter = 50000, fdr = 0.15, seed = 1234)
   end_time = Sys.time()
   time_cell_mcadet  = append(time_cell_mcadet , end_time - start_time)
   print(i)
}




time_gene_mcadet = c()
for(i in 1:15){
  set.seed(i*12345)
  idx = sample(1:15000, size = 1000*i)  
  b = sparsim.data[idx, ]
  
  start_time = Sys.time()
  
  a = mcadet(data = b,
             n.comp= 60, run=10, n.feature=NA, nk_percent = 0.010, 
             start_resolution = 0.5,  cell_percent =  0.005,
             MC_iter = 50000, fdr = 0.15, seed = 1234)
  end_time = Sys.time()
  time_gene_mcadet = append(time_gene_mcadet, end_time - start_time)
  print(i)
}





NBDrop_time_cell = c()
for (i in 1:20) {
  set.seed(i*1234)
  idx = sample(1:20000, size = 1000*i)  
  b = sparsim.data[, idx]
  
  start_time = Sys.time()
  
  nb_raw <- NBumiConvertData(b, is.counts=TRUE)
  DANB_fit <- NBumiFitModel(nb_raw)
  NBDropFS<- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thres=0.05, suppress.plot=F)
  end_time = Sys.time()
  NBDrop_time_cell = append(NBDrop_time_cell, end_time - start_time)
  print(i)
  gc()
}



NBDrop_time_gene = c()
for (i in 1:15) {
  set.seed(i*12345)
  idx = sample(1:15000, size = 1000*i)  
  b = sparsim.data[idx,]
  
  start_time = Sys.time()
  
  nb_raw <- NBumiConvertData(b, is.counts=TRUE)
  DANB_fit <- NBumiFitModel(nb_raw)
  NBDropFS<- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thres=0.05, suppress.plot=F)
  end_time = Sys.time()
  NBDrop_time_gene = append(NBDrop_time_gene, end_time - start_time)
  print(i)
  gc()
}



M3Drop_time_cell = c()

for (i in 1:20) {
  set.seed(i*1234)
  idx = sample(1:20000, size = 1000*i)  
  b = sparsim.data[, idx]
  start_time = Sys.time()
  M3_data <- M3DropConvertData(b, is.counts=TRUE)
  M3DropFS<- M3DropFeatureSelection(M3_data ,mt_method="fdr", mt_threshold=0.05)
  end_time = Sys.time()
  M3Drop_time_cell = append(M3Drop_time_cell, end_time - start_time)
  print(i)
  gc()
}




M3Drop_time_gene = c()

for (i in 1:15) {
  set.seed(i*12345)
  idx = sample(1:15000, size = 1000*i)  
  b = sparsim.data[idx,]
  start_time = Sys.time()
  M3_data <- M3DropConvertData(b, is.counts=TRUE)
  M3DropFS<- M3DropFeatureSelection(M3_data ,mt_method="fdr", mt_threshold=0.05)
  end_time = Sys.time()
  M3Drop_time_gene = append(M3Drop_time_gene, end_time - start_time)
  print(i)
  gc()
}






Bre_time_cell = c()

for (i in 1:20) {
  set.seed(i*1234)
  idx = sample(1:20000, size = 1000*i)  
  b = sparsim.data[, idx]
  start_time = Sys.time()
  BreFS = BrenneckeGetVariableGenes(b, suppress.plot = T)
  end_time = Sys.time()
  Bre_time_cell = append(Bre_time_cell, end_time - start_time)
  print(i)
  gc()
}


Bre_time_gene = c()

for (i in 1:15) {
  set.seed(i*12345)
  idx = sample(1:15000, size = 1000*i)  
  b = sparsim.data[idx, ]
  start_time = Sys.time()
  BreFS = BrenneckeGetVariableGenes(b, suppress.plot = T)
  end_time = Sys.time()
  Bre_time_gene = append(Bre_time_gene, end_time - start_time)
  print(i)
  gc()
}








SeuratVst_time_cell = c()

for (i in 1:20) {
  set.seed(i*1234)
  idx = sample(1:20000, size = 1000*i)  
  b = sparsim.data[, idx ]
  start_time = Sys.time()
  seurat.obj<-CreateAssayObject(counts = b)
  seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj.vst <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
  SeuratVstFS<-VariableFeatures(seurat.obj.vst)
  end_time = Sys.time()
  SeuratVst_time_cell = append(SeuratVst_time_cell, end_time - start_time)
  print(i)
  gc()
}



SeuratVst_time_gene = c()

for (i in 1:15) {
  set.seed(i*12345)
  idx = sample(1:15000, size = 1000*i)  
  b = sparsim.data[idx, ]
  start_time = Sys.time()
  seurat.obj<-CreateAssayObject(counts = b)
  seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj.vst <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
  SeuratVstFS<-VariableFeatures(seurat.obj.vst)
  end_time = Sys.time()
  SeuratVst_time_gene = append(SeuratVst_time_gene, end_time - start_time)
  print(i)
  gc()
}





SeuratDisp_time_cell = c()

for (i in 1:20) {
  set.seed(i*1234)
  idx = sample(1:20000, size = 1000*i)  
  b = sparsim.data[, idx ]
  start_time = Sys.time()
  seurat.obj<-CreateAssayObject(counts = b)
  seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj.Disp <- FindVariableFeatures(seurat.obj, selection.method = "disp", nfeatures = 2000)
  SeuratDispFS<-VariableFeatures(seurat.obj.Disp)
  end_time = Sys.time()
  SeuratDisp_time_cell = append(SeuratDisp_time_cell, end_time - start_time)
  print(i)
  gc()
}



SeuratDisp_time_gene = c()

for (i in 1:15) {
  set.seed(i*12345)
  idx = sample(1:15000, size = 1000*i)  
  b = sparsim.data[ idx ,  ]
  start_time = Sys.time()
  seurat.obj<-CreateAssayObject(counts = b)
  seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj.Disp <- FindVariableFeatures(seurat.obj, selection.method = "disp", nfeatures = 2000)
  SeuratDispFS<-VariableFeatures(seurat.obj.Disp)
  end_time = Sys.time()
  SeuratDisp_time_gene = append(SeuratDisp_time_gene, end_time - start_time)
  print(i)
  gc()
}










SeuratMvp_time_cell = c()

for (i in 1:20) {
  set.seed(i*1234)
  idx = sample(1:20000, size = 1000*i)  
  b = sparsim.data[, idx ]
  start_time = Sys.time()
  seurat.obj<-CreateAssayObject(counts = b)
  seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj.Mvp <- FindVariableFeatures(seurat.obj, selection.method = "mvp")
  SeuratMvpFS<-VariableFeatures(seurat.obj.Mvp)
  end_time = Sys.time()
  SeuratMvp_time_cell  = append(SeuratMvp_time_cell , end_time - start_time)
  print(i)
  gc()
}



SeuratMvp_time_gene = c()

for (i in 1:15) {
  set.seed(i*12345)
  idx = sample(1:15000, size = 1000*i)  
  b = sparsim.data[idx, ]
  start_time = Sys.time()
  seurat.obj<-CreateAssayObject(counts = b)
  seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj.Mvp <- FindVariableFeatures(seurat.obj, selection.method = "mvp")
  SeuratMvpFS<-VariableFeatures(seurat.obj.Mvp)
  end_time = Sys.time()
  SeuratMvp_time_gene  = append(SeuratMvp_time_gene , end_time - start_time)
  print(i)
  gc()
}




Scry_time_cell = c()

for (i in 1:20) {
  set.seed(i*1234)
  idx = sample(1:20000, size = 1000*i)  
  b = sparsim.data[, idx ]
  start_time = Sys.time()
  scry=devianceFeatureSelection(b, fam = "poisson")
  scry_vector = names(sort(scry,decreasing = T)[1:2000])
  end_time = Sys.time()
  Scry_time_cell = append(Scry_time_cell, end_time - start_time)
  print(i)
  gc()
}



Scry_time_gene = c()

for (i in 1:15) {
  set.seed(i*12345)
  idx = sample(1:15000, size = 1000*i)  
  b = sparsim.data[idx,]
  start_time = Sys.time()
  scry=devianceFeatureSelection(b, fam = "poisson")
  scry_vector = names(sort(scry,decreasing = T)[1:2000])
  end_time = Sys.time()
  Scry_time_gene = append(Scry_time_gene, end_time - start_time)
  print(i)
  gc()
}




## making the plot 


time_df_cell = data.frame("time" = c(as.numeric(time_cell_mcadet),
                                     as.numeric(NBDrop_time_cell/60),
                                     as.numeric(M3Drop_time_cell/60),
                                     as.numeric(Bre_time_cell/60),
                                     as.numeric(Scry_time_cell/60),
                                     as.numeric(SeuratVst_time_cell/60),
                                     as.numeric(SeuratDisp_time_cell/60),
                                     as.numeric(SeuratMvp_time_cell/60)), 
  "Ncells" = rep(seq(1000,20000,1000), 8),
  "Method" = rep( c("Mcadet", "NBDrop",
                    "M3Drop", "Brennecke",
                    "Scry", "Seurat Vst",
                    "Seurat Disp", "Seurat Mvp")   ,each = 20))



time_df_gene = data.frame("time" = c(as.numeric(time_gene_mcadet),
                                     as.numeric(NBDrop_time_gene/60),
                                     as.numeric(M3Drop_time_gene/60),
                                     as.numeric(Bre_time_gene/60),
                                     as.numeric(Scry_time_gene/60),
                                     as.numeric(SeuratVst_time_gene/60),
                                     as.numeric(SeuratDisp_time_gene/60),
                                     as.numeric(SeuratMvp_time_gene/60)), 
                          "Ngenes" = rep(seq(1000,15000,1000), 8),
                          "Method" = rep( c("Mcadet", "NBDrop",
                                            "M3Drop", "Brennecke",
                                            "Scry", "Seurat Vst",
                                            "Seurat Disp", "Seurat Mvp"), each = 15))

time_df_cell$Method  =  factor(time_df_cell$Method, levels = c("Mcadet", "NBDrop",
                                                               "M3Drop", "Brennecke",
                                                               "Scry", "Seurat Vst",
                                                               "Seurat Disp", "Seurat Mvp"))

time_df_gene$Method  =  factor(time_df_gene$Method, levels = c("Mcadet", "NBDrop",
                                                               "M3Drop", "Brennecke",
                                                               "Scry", "Seurat Vst",
                                                               "Seurat Disp", "Seurat Mvp"))


save(time_df_cell, file = "C:/Users/css22/Desktop/Thesis/Results/Time/time_df_cell.RData")
save(time_df_gene, file = "C:/Users/css22/Desktop/Thesis/Results/Time/time_df_gene.RData")


load(file = "C:/Users/css22/Desktop/Thesis1/Results/Time/time_df_cell.RData")
load(file = "C:/Users/css22/Desktop/Thesis1/Results/Time/time_df_gene.RData")




time_df_cell$Method = factor(time_df_cell$Method, levels = c("Mcadet", "Scry", "Brennecke",
                                                 "Seurat Disp", "Seurat Vst", "Seurat Mvp",
                                                 "M3Drop", "NBDrop"))


library(ggplot2)





linechart_runningtime_cells = ggplot(data =  time_df_cell, 
                    mapping = aes(x = Ncells, y = time, 
                                  colour = Method)) + geom_line(linewidth = 1.2)+
  geom_point(size = 3)+
  labs(x = "", y = "Time (mins)", 
       title = "Running time of different FS methods by number of cells" )+
  theme_bw() + scale_colour_discrete(name = "Method", 
                                     labels = c("Mcadet", "Scry", "Brennecke",
                                                "Seurat Disp", "Seurat Vst", "Seurat Mvp",
                                                "M3Drop", "NBDrop"))+
  scale_x_continuous(limits=c(1000, 20000), breaks=seq(1000, 20000, 2000))+
  theme( panel.grid = element_blank(),  
         panel.grid.major = element_line(color = "grey", linetype = "dashed"),  
         panel.grid.minor = element_line(color = "grey", linetype = "dashed"),
         axis.title = element_text(face ="bold", size = 16),
         plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
         legend.position = "",
         legend.text = element_text(size = 10, face = "bold"),
         legend.title = element_blank(), 
         legend.background = element_rect(color = "black", 
                                          linewidth = 1.5,
                                          fill = "lightblue"),
         axis.text.x = element_text(face = "bold", size = 12, color = "black"), 
         axis.text.y = element_text(face = "bold", size = 12),
         panel.border = element_rect(color = "black", fill = NA, size = 1.5))+
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -1, vjust = 1.5, size = 8, fontface = "bold")











time_df_gene$Method = factor(time_df_gene$Method, levels = c("Mcadet", "Scry", "Brennecke",
                                                             "Seurat Disp", "Seurat Vst", "Seurat Mvp",
                                                             "M3Drop", "NBDrop"))


linechart_runningtime_genes = ggplot(data =  time_df_gene, 
                                     mapping = aes(x = Ngenes, y = time, 
                                                   colour = Method)) + geom_line(linewidth = 1.2)+
  geom_point(size = 3)+
  labs(x = "", y = "Time (mins)", 
       title = "Running time of different FS methods by number of genes" )+
  theme_bw() + scale_colour_discrete(name = "Method", 
                                     labels = c("Mcadet", "Scry", "Brennecke",
                                                "Seurat Disp", "Seurat Vst", "Seurat Mvp",
                                                "M3Drop", "NBDrop"))+
  scale_x_continuous(limits=c(1000, 15000), breaks=seq(1000, 15000, 1000))+
  theme( panel.grid = element_blank(),  
         panel.grid.major = element_line(color = "grey", linetype = "dashed"),  
         panel.grid.minor = element_line(color = "grey", linetype = "dashed"),
         axis.title = element_text(face ="bold", size = 16),
         plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
         legend.position = "bottom",
         legend.text = element_text(size = 10, face = "bold"),
         legend.title = element_blank(), 
         legend.background = element_rect(color = "black", 
                                          linewidth = 1.5,
                                          fill = "lightblue"),
         axis.text.x = element_text(face = "bold", size = 12, color = "black"), 
         axis.text.y = element_text(face = "bold", size = 12),
         panel.border = element_rect(color = "black", fill = NA, size = 1.5))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -1, vjust = 1.5, size = 8, fontface = "bold")


linechart_runningtime_cells / linechart_runningtime_genes
