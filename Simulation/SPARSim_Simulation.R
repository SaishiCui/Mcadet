
library(devtools)
install_gitlab("sysbiobig/sparsim")



library("SPARSim")
data(Zheng_param_preset)

sparsim_gene_name<-NULL
for (i in 1:15000) {
  sparsim_gene_name<-append(sparsim_gene_name, paste0("Gene", i ))
}

sparsim_cell_name<-NULL
for (i in 1:6000) {
  sparsim_cell_name<-append(sparsim_cell_name, paste0("Cell", i ))
}


for (i in 1:24) {
  for (j in 1:2) {
    downreg.min.set<- seq(0.1, 0.6, by=0.5/23)
    downreg.max.set<- seq(0.3, 0.8, by=0.5/23)
    upreg.min.set<- rev(seq(1, 3, by=2/23))
    upreg.max.set<- rev(seq(1.5, 3.5, by=2/23))
    
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
      N_cells = 60,
      condition_name = "cell_1")
    
    
    cell_type_2 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C1, 
      fc_multiplier = fold_change_multiplier_2, 
      N_cells = 180,
      condition_name = "cell_2")
    
    
    
    cell_type_3 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C1 , 
      fc_multiplier = fold_change_multiplier_3, 
      N_cells = 240,
      condition_name = "cell_3")
    
    
    cell_type_4 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C1 , 
      fc_multiplier = fold_change_multiplier_4, 
      N_cells = 360,
      condition_name = "cell_4")
    
    
    cell_type_5 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C1 , 
      fc_multiplier = fold_change_multiplier_5, 
      N_cells = 420,
      condition_name = "cell_5")
    
    
    cell_type_6 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C1 , 
      fc_multiplier = fold_change_multiplier_6, 
      N_cells = 480,
      condition_name = "cell_6")
    
    
    cell_type_7 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C1 , 
      fc_multiplier = fold_change_multiplier_7, 
      N_cells = 600,
      condition_name = "cell_7")
    
    
    cell_type_8 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C1 , 
      fc_multiplier = fold_change_multiplier_8, 
      N_cells = 960,
      condition_name = "cell_8")
    
    
    
    cell_type_9 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C1 , 
      fc_multiplier = fold_change_multiplier_9, 
      N_cells = 1200,
      condition_name = "cell_9")
    
    
    
    
    cell_type_10 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C1 , 
      fc_multiplier = fold_change_multiplier_10, 
      N_cells = 1500,
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
    
    
    
    assign(paste0("SparSim_testing_", (i-1)*2+j), sparsim.data)
    
    save(list= paste0("SparSim_testing_",(i-1)*2+j), 
         file=paste0("E:/Thesis/SPARSim/SparSim_testing_", (i-1)*2+j,  ".RData" ))
    
    a=((i-1)*2+j)/48*100
    print(paste0("Progress: ", a, "%"))
    
    
  }
}








## Generating results simulation datasets

for (i in 1:24) {

    downreg.min.set<- seq(0.4, 0.6, by=0.2/23)
    downreg.max.set<- seq(0.6, 0.8, by=0.2/23)
    upreg.min.set<- rev(seq(1.5, 2, by= 0.5/23))
    upreg.max.set<- rev(seq(2, 2.5, by= 0.5/23))
    
    set.seed( i*10000 + 1 )
    DE_multiplier_1.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2)
    set.seed( i*10000 + 2 )
    DE_multiplier_1.2 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    set.seed( i*10000 + 3 )
    DE_multiplier_1.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    set.seed( i*10000 + 4 )
    DE_multiplier_1.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    fold_change_multiplier_1 <- c(DE_multiplier_1.1, rep(1, 180), DE_multiplier_1.2, rep(1,400),
                                  DE_multiplier_1.3, rep(1,400), DE_multiplier_1.4, rep(1,17936) )
    
    
    set.seed( i*10000 + 5 )
    DE_multiplier_2.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
    set.seed( i*10000 + 6 )
    DE_multiplier_2.2 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
    set.seed( i*10000 + 7 )
    DE_multiplier_2.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    set.seed( i*10000 + 8 )
    DE_multiplier_2.4 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
    set.seed( i*10000 + 9 )
    DE_multiplier_2.5 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
    fold_change_multiplier_2 <- c(rep(1,20), DE_multiplier_2.1, rep(1,360), DE_multiplier_2.2,
                                  rep(1,200), DE_multiplier_2.3, rep(1,200), DE_multiplier_2.4, 
                                  rep(1,200), DE_multiplier_2.5, rep(1,17736) )
    
    
    set.seed( i*10000 + 10 )
    DE_multiplier_3.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
    set.seed(  i*10000 + 11 )
    DE_multiplier_3.2 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    set.seed( i*10000 + 12 )
    DE_multiplier_3.3 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
    set.seed( i*10000 + 13 )
    DE_multiplier_3.4 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
    set.seed( i*10000 + 14 )
    DE_multiplier_3.5 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
    fold_change_multiplier_3 <- c(rep(1,40), DE_multiplier_3.1, rep(1,140), DE_multiplier_3.2,
                                  DE_multiplier_3.3, rep(1,800), DE_multiplier_3.4, 
                                  rep(1,200), DE_multiplier_3.5, rep(1,17536) )
    
    
    
    set.seed( i*10000 + 15 )
    DE_multiplier_4.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
    set.seed( i*10000 + 16 )
    DE_multiplier_4.2 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    set.seed( i*10000 + 17 )
    DE_multiplier_4.3 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
    set.seed( i*10000 + 18 )
    DE_multiplier_4.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    set.seed( i*10000 + 19 )
    DE_multiplier_4.5 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    set.seed( i*10000 + 20 )
    DE_multiplier_4.6 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    fold_change_multiplier_4 <- c(rep(1,60), DE_multiplier_4.1, rep(1,320), DE_multiplier_4.2, rep(1,400),
                                  DE_multiplier_4.3, DE_multiplier_4.4, DE_multiplier_4.5, 
                                  rep(1,200), DE_multiplier_4.6, rep(1,17536) )
    
    
    set.seed( i*10000 + 21 )
    DE_multiplier_5.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
    set.seed( i*10000 + 22 )
    DE_multiplier_5.2 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    set.seed( i*10000 + 23 )
    DE_multiplier_5.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    set.seed( i*10000 + 24 )
    DE_multiplier_5.4 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
    set.seed( i*10000 + 25 )
    DE_multiplier_5.5 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    fold_change_multiplier_5 <- c(rep(1,80), DE_multiplier_5.1, rep(1,100), DE_multiplier_5.2,
                                  rep(1,200), DE_multiplier_5.3, rep(1,800), DE_multiplier_5.4, 
                                  DE_multiplier_5.5, rep(1,17536) )
    
    
    
    
    set.seed( i*10000 + 26 )
    DE_multiplier_6.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
    set.seed( i*10000 + 27 )
    DE_multiplier_6.2 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
    set.seed( i*10000 + 28 )
    DE_multiplier_6.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    set.seed( i*10000 + 29 )
    DE_multiplier_6.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    set.seed( i*10000 + 30 )
    DE_multiplier_6.5 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
    fold_change_multiplier_6 <- c(rep(1,100), DE_multiplier_6.1, rep(1,280), DE_multiplier_6.2,
                                  rep(1,400), DE_multiplier_6.3, rep(1,200), DE_multiplier_5.4, 
                                  rep(1,200), DE_multiplier_6.5, rep(1,17536) )
    
    
    
    set.seed( i*10000 + 31 )
    DE_multiplier_7.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
    set.seed( i*10000 + 32 )
    DE_multiplier_7.2 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
    set.seed( i*10000 + 33 )
    DE_multiplier_7.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    set.seed( i*10000 + 34 )
    DE_multiplier_7.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    fold_change_multiplier_7 <- c(rep(1,120), DE_multiplier_7.1, rep(1,60), DE_multiplier_7.2,
                                  rep(1,200), DE_multiplier_7.3, rep(1,400), DE_multiplier_7.4, rep(1,18136) )
    
    
    
    set.seed( i*10000 + 35 )
    DE_multiplier_8.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
    set.seed( i*10000 + 36 )
    DE_multiplier_8.2 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
    set.seed( i*10000 + 37 )
    DE_multiplier_8.3 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
    set.seed( i*10000 + 38 )
    DE_multiplier_8.4 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
    fold_change_multiplier_8 <- c(rep(1,140), DE_multiplier_8.1, rep(1,440), DE_multiplier_8.2,
                                  rep(1,200), DE_multiplier_8.3, rep(1,400), DE_multiplier_8.4, rep(1,17736) )
    
    
    
    
    
    set.seed( i*10000 + 39 )
    DE_multiplier_9.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
    set.seed( i*10000 + 40 )
    DE_multiplier_9.2 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
    set.seed( i*10000 + 41 )
    DE_multiplier_9.3 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    set.seed( i*10000 + 42 )
    DE_multiplier_9.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    set.seed( i*10000 + 43 )
    DE_multiplier_9.5 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
    fold_change_multiplier_9 <- c(rep(1,160), DE_multiplier_9.1, rep(1,20), DE_multiplier_9.2,
                                  rep(1,800), DE_multiplier_9.3, rep(1,200), DE_multiplier_9.4, 
                                  DE_multiplier_9.5, rep(1,17536) )
    
    
    
    
    set.seed( i*10000 + 44 )
    DE_multiplier_10.1 <- runif(n = 20, min = upreg.min.set[i]*1.2, max = upreg.max.set[i]*1.2 )
    set.seed( i*10000 + 45 )
    DE_multiplier_10.2 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    set.seed( i*10000 + 46 )
    DE_multiplier_10.3 <- runif(n = 200, min = downreg.min.set[i], max = downreg.max.set[i])
    set.seed( i*10000 + 47 )
    DE_multiplier_10.4 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    set.seed( i*10000 + 48 )
    DE_multiplier_10.5 <- runif(n = 200, min = upreg.min.set[i], max = upreg.max.set[i])
    fold_change_multiplier_10 <- c(rep(1,180), DE_multiplier_9.1, rep(1,600), DE_multiplier_10.2,
                                  DE_multiplier_10.3, rep(1,200), DE_multiplier_10.4, 
                                  DE_multiplier_10.5, rep(1,17736) )
    
    
    
  
    
    cell_type_1 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C2, 
      fc_multiplier = fold_change_multiplier_1, 
      N_cells = 60,
      condition_name = "cell_1")
    
    
    cell_type_2 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C2, 
      fc_multiplier = fold_change_multiplier_2, 
      N_cells = 180,
      condition_name = "cell_2")
    
    
    
    cell_type_3 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C2 , 
      fc_multiplier = fold_change_multiplier_3, 
      N_cells = 240,
      condition_name = "cell_3")
    
    
    cell_type_4 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C2 , 
      fc_multiplier = fold_change_multiplier_4, 
      N_cells = 360,
      condition_name = "cell_4")
    
    
    cell_type_5 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C2 , 
      fc_multiplier = fold_change_multiplier_5, 
      N_cells = 420,
      condition_name = "cell_5")
    
    
    cell_type_6 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C2 , 
      fc_multiplier = fold_change_multiplier_6, 
      N_cells = 480,
      condition_name = "cell_6")
    
    
    cell_type_7 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C2 , 
      fc_multiplier = fold_change_multiplier_7, 
      N_cells = 600,
      condition_name = "cell_7")
    
    
    cell_type_8 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C2 , 
      fc_multiplier = fold_change_multiplier_8, 
      N_cells = 960,
      condition_name = "cell_8")
    
    
    
    cell_type_9 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C2 , 
      fc_multiplier = fold_change_multiplier_9, 
      N_cells = 1200,
      condition_name = "cell_9")
    
    
    
    
    cell_type_10 <- SPARSim_create_DE_genes_parameter(
      sim_param =Zheng_param_preset$Zheng_C2 , 
      fc_multiplier = fold_change_multiplier_10, 
      N_cells = 1500,
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
    
    set.seed(1234*i)
    ind<-sample(c(2001:19536), size = 4536 ,replace = F)
    
    sparsim.data=SPARSim_result$count_matrix[-ind,]
    
    rownames(sparsim.data)=sparsim_gene_name
    colnames(sparsim.data)=sparsim_cell_name
    
    
    
    assign(paste0("SparSim_fine_", i), sparsim.data)
    
    save(list= paste0("SparSim_fine_",i), 
         file=paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_fine_",
                     i,  ".RData" ))
    
    rm(list = paste0("SparSim_fine_", i))
    rm(sparsim.data)
    gc()
    
    a=i/24*100
    print(paste0("Progress: ", a, "%"))
    
    
  
}


























