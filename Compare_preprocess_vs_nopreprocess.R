
rm(list=ls())
gc()

mcadet_nopreprocess <- function(data, n.comp=60,  run=10, n.feature=NA, nk_percent = 0.005, start_resolution = 0.5,
                                cell_percent=0.005, MC_iter = 50000, fdr = 0.15, seed = 1234){
  
  
  clean_fucntion = function(x){
    percent = mean(1*(x!=0 ))
  }
  
  percent_vector=apply(data,1,clean_fucntion)            ##  cleaning genes with low expressed cells##
  working.data = data[percent_vector>= cell_percent,]
  gene.name<-rownames(working.data)                       ##  extract gene's name in a list ##
  
  ### Construct corresponding matrix 
  
  P = 1/sum( t(working.data) ) * t(working.data)
  
  row.mass = rowSums(P)  ## calculate row mass 
  col.mass = colSums(P)  ## calculate column mass
  
  
  
  ### Construct standardized Pearson residual matrix 
  
  
  
  expected.matrix = as.matrix(row.mass)%*%t(col.mass)
  sp.res.matrix =  (P-expected.matrix)/sqrt(expected.matrix)
  
  
  rm(working.data)
  rm(expected.matrix)
  gc()
  
  
  
  ## irlba (fast and memory efficient SVD) ##
  
  library(irlba)
  
  irlba.decomp<-irlba(sp.res.matrix, n.comp)
  
  
  U_k = irlba.decomp$u[,c(1:n.comp)]                          ## extract svd U
  V_k = irlba.decomp$v[,c(1:n.comp)]    ## extract svd V (first n.comp PCs)
  
  
  D_alpha_K= diag(irlba.decomp$d[1:n.comp],n.comp)        ## diagonal matrix of singular values 
  
  standard.row.coord  =  sqrt(dim(P)[1])*U_k              ## calculate standard row coordinates (cell coordinates) 
  
  principal.col.coord =  (1/sqrt(col.mass))*V_k%*%D_alpha_K ## calculate principal col coordinates (gene + coordinates) 
  
  
  
  
  
  ### KNN graph and leiden detection
  jaccard = function(a,b){
    numerator = length(intersect(a,b))
    denominator = length(a) + length(b) - numerator
    return( numerator/denominator)
  }
  
  shrunk = function(x){
    which(x==1)
  }
  
  
  library("cccd")
  library("igraph")
  library("RANN")
  
  gene.logR.matrix<-matrix(data=0, nrow = run, ncol = nrow(principal.col.coord))  ## set up an empty matrix
  cluster.list<- c()
  n_cells = nrow(standard.row.coord)
  nn_k = round(n_cells*nk_percent)
  if(nn_k <4){
    nn_k  = 4
  }
  weights_matrix = matrix(0, nrow = n_cells, ncol = nn_k)
  knn_graph_matrix <- nn2(data = standard.row.coord, k = nn_k + 1)$nn.idx[,-1] # +1 to include the point itself
  
  
  for (q in 1:n_cells) {
    for (z in 1:nn_k) {
      neighbor = knn_graph_matrix[q,z]
      weights_matrix[q,z] =jaccard(knn_graph_matrix[q,], knn_graph_matrix[neighbor,])
    }
  }
  
  
  
  # Initialize an empty list to store the edges
  edges <- c()
  
  # Fill the edges list based on the KNN index matrix
  for (i in 1:nrow(knn_graph_matrix)) {
    neighbors <- knn_graph_matrix[i, ] # Exclude the point itself (first column)
    edges <-c(edges, c(rbind(i, neighbors)))
  }
  
  # Create the graph object
  knn_graph <- graph(edges = edges, directed = FALSE)
  
  
  
  for(i in 1:run){
    
    set.seed(seed*i)
    partition = cluster_leiden(graph = knn_graph, 
                               objective_function = c("modularity"),
                               weights  = weights_matrix, 
                               resolution_parameter = start_resolution + 0.1*i)$membership ## Leiden detection
    
    
    gene.rank.matrix<-matrix(NA, ncol = nrow(principal.col.coord), nrow =length(unique(partition)))
    colnames(gene.rank.matrix)<-gene.name
    cluster<-sort(unique(partition))
    cluster.list<-append(cluster.list,length(cluster))
    for (j in 1:length(unique(partition))) {
      
      if(sum(1*(partition==j)) != 1){
        
        gene.rank.matrix[j,]<-rank(sqrt(rowSums((t(matrix( rep(colMeans(standard.row.coord[partition==j,]),   ## Calculate cluster centroid and rank within each cluster
                                                               nrow(principal.col.coord)), ncol =nrow(principal.col.coord)))-principal.col.coord)^2)))
      }else{
        gene.rank.matrix[j,]<-rank(sqrt(rowSums((t(matrix( rep(standard.row.coord[partition==j,],   ## Calculate cluster centroid and rank within each cluster
                                                               nrow(principal.col.coord)), ncol =nrow(principal.col.coord)))-principal.col.coord)^2)))
      }
      
    }
    gene.logR.matrix[i,]<-apply(gene.rank.matrix,2,function(x){log(max(x)/min(x)) })   ##  a matrix with each gene's rank range for each run
    
    
  }
  
  colnames(gene.logR.matrix)<-gene.name
  gene.logR<-colMeans(gene.logR.matrix)     ##  mean log Maximum/Minimum for each gene over all runs ##
  
  
  
  
  if(is.na(n.feature)){
    
    
    logR_matrix<-matrix(NA, nrow = length(cluster.list), ncol = MC_iter)
    for (j in 1:length(cluster.list)) {
      logR_list<-rep(0,MC_iter)
      for (i in 1:MC_iter) {
        sam<-sample( c(1:length(gene.name)), size = cluster.list[j], replace = T)
        logR<- log(max(sam)/min(sam) ) 
        logR_list[i]=logR
      }
      logR_matrix[j,]=logR_list       
    }
    
    logR_hist<-colMeans(logR_matrix)   ## Monte Carlo simulation for log Maximum/Minimum distribution (multiple runs) ##
    
    
    
    ### FDR BH method
    
    pval_list<-rep(NA,length(gene.name))
    for (i in 1:length(gene.name)) {
      p.val=mean(logR_hist>gene.logR[i])    ## upper-tailed p-values calculation ##
      pval_list[i]=p.val
    }
    
    pval_df<-data.frame(gene.name,pval_list)
    sorted_pval<- sort(pval_list)
    sorted_pval<-sorted_pval[sorted_pval<0.9]    # remove last one or two p-value histogram blocks,
    # to make p-values be uniformly distributed ##
    
    
    
    
    position_list<-NULL
    for (i in 1:length(sorted_pval)) {
      if(sorted_pval[i]< fdr*i/length(sorted_pval)){
        position_list<-append(position_list,i)}else{
          break
        }
      
    }
    fdr_thresh<-sorted_pval[length(position_list)]   # find L and corresponding p-value threshold in BH precedure#
    
    
    gene.list<-pval_df[pval_df$pval_list<=fdr_thresh,"gene.name"]
    
    gene.logR.list<-gene.logR[pval_df$pval_list<=fdr_thresh]
    
    p.val.list<-pval_df[pval_df$pval_list<=fdr_thresh,"pval_list"]
    
    df<-data.frame(gene.list, gene.logR.list, p.val.list)  
    colnames(df)<-c("gene", "logR", "p.value")  
    df<-df[order(df$p.value, decreasing = F),]
    
    obj<-df  
    return(obj)
    
  }else{
    DE_label<-order(gene.logR,decreasing = T)[1:n.feature]   # if users know how many genes they want #
    gene.list<-gene.name[DE_label]
    return(gene.list)
    
  }
  
  
}
jaccard <-function(a,b){
  jac = length(intersect(a,b))/(length(union(a,b)))
  return(jac)
}



## fine resolution


### load the data 

for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/DE_genes/DEgene_", 
                      donor.set[p], "_", time.set[k], ".RData"))
  }
}






for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/pbmc_results/Mcadet_pbmc_", 
                      donor.set[p], "_", time.set[k], ".RData"))
  }
}





for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/data/pbmcdata_", 
                      donor.set[p], "_", time.set[k], ".RData"))
  }
}



for (i in 1:24){
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/SparSim_results/Mcadet_SparSim_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/SparSim_results/Mcadet_SparSim_fine_", i, ".RData"))
}


for (i in 1:24){
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/SPARSim/SparSim_Mcadet_data/SparSim_fine_", i, ".RData"))
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/SPARSim/SparSim_Mcadet_data/SparSim_", i, ".RData"))
  
  dta = get(paste0("SparSim_fine_",  i))
  result = mcadet(data = dta)
  assign(paste0("Mcadet_process_SimFine_", i),  result)
  
  
  dta = get(paste0("SparSim_",  i))
  result = mcadet(data = dta)
  assign(paste0("Mcadet_process_SimCoarse", i),  result)
  
  rm(list = paste0("SparSim_fine_",  i))
  rm(list = paste0("SparSim_",  i))
  gc()
  print(i)
}








for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    
    dta = get(paste0("pbmcdata_",  donor.set[p], "_", time.set[k]))
    result = mcadet_nopreprocess(data = dta)
    assign(paste0("Mcadet_noprocess_", donor.set[p], "_", time.set[k]),  result)
    gc()
    
    
  }
}




####### Jaccard 




pbmc_jaccard_macdet = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_", donor.set[p], "_", time.set[k])),
                  get(paste0("Mcadet_pbmc_", donor.set[p], "_", time.set[k]))$gene )
    pbmc_jaccard_macdet = append(pbmc_jaccard_macdet, jac)
  }
}
pbmc_jaccard_macdet





pbmc_jaccard_macdet_noprocess = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_", donor.set[p], "_", time.set[k])),
                  get(paste0("Mcadet_noprocess_", donor.set[p], "_", time.set[k]))$gene )
    pbmc_jaccard_macdet_noprocess = append(pbmc_jaccard_macdet_noprocess, jac)
  }
}
pbmc_jaccard_macdet_noprocess






Sim_jaccard_macdet_noprocess_fine = c()
for (i in 1:24) {
    jac = jaccard(paste0("Gene", c(1:2000)),
                  get(paste0("Mcadet_noprocess_SimFine_", i))$gene)
    Sim_jaccard_macdet_noprocess_fine  = append(Sim_jaccard_macdet_noprocess_fine, jac)
  
}
Sim_jaccard_macdet_noprocess_fine





Sim_jaccard_macdet_process_fine = c()
for (i in 1:24) {
  jac = jaccard(paste0("Gene", c(1:2000)),
                get(paste0("Mcadet_process_SimFine_", i))$gene)
  Sim_jaccard_macdet_process_fine  = append(Sim_jaccard_macdet_process_fine, jac)
  
}
Sim_jaccard_macdet_process_fine




mean(pbmc_jaccard_macdet)
mean(pbmc_jaccard_macdet_noprocess)



############### Plot


pbmc_coarse_process <- c(0.4241026, 0.4465593, 0.4440594, 0.4520214,
                       0.3937210, 0.4754522, 0.3675799, 0.4197609,
                       0.4635993, 0.4582624, 0.4439980, 0.4432071,
                       0.4266350, 0.3898477, 0.4984843, 0.4631937,
                       0.4750831, 0.4119980, 0.3931596, 0.4685694,
                       0.4632325, 0.4768519, 0.4992248, 0.5284810)

pbmc_coarse_noprocess <- c(0.2590341, 0.2562674, 0.2860270, 0.3769612,
                         0.3206669, 0.3540904, 0.2773360, 0.3468697,
                         0.3290960, 0.2752000, 0.2959641, 0.2947695,
                         0.3065442, 0.3046899, 0.2968536, 0.3005698,
                         0.3342522, 0.2481858, 0.3113402, 0.2931422,
                         0.3155130, 0.2817460, 0.2946292, 0.3067606)




pbmc_fine_process <- c(0.3182518, 0.3021219, 0.3302665, 0.3492437,
                         0.3072566, 0.2929751, 0.3413094, 0.3997361, 
                         0.3178476, 0.3376323, 0.3408403, 0.3304862,
                         0.3074992, 0.3258290, 0.3224285, 0.2886717,
                         0.3050751, 0.3243002, 0.3110410, 0.2847073,
                         0.3234855, 0.3289435, 0.3322543, 0.3008487)



pbmc_fine_noprocess <- c(0.3099851, 0.2612586, 0.2868514, 0.2920381,
                           0.2990334, 0.3023889, 0.2759330, 0.2770472,
                           0.2734561, 0.2326241, 0.2816005, 0.2639314,
                           0.2840322, 0.2900336, 0.2980075, 0.2516538,
                           0.2724449, 0.2485453, 0.2588968, 0.2493861,
                           0.3063218, 0.2323120, 0.2370572, 0.2571107)




sim_fine_process <- c(0.3766797, 0.3891085, 0.3311772, 0.3826731, 
                      0.3640275, 0.3374778, 0.3445415, 0.3609520,
                      0.3153915, 0.3082312, 0.3423110, 0.3082969,
                      0.2576153, 0.2368889, 0.2496656, 0.2180923,
                      0.1726714, 0.2390049, 0.2321509, 0.2275556, 
                      0.2087766, 0.1980969, 0.1266726, 0.1313673)



sim_fine_noprocess <- c(0.3777329, 0.3739367, 0.3874126, 0.3708093,
                        0.3728412, 0.4021379, 0.3756432, 0.3950573,
                        0.3673538, 0.3899614, 0.3615568, 0.3831909,
                        0.3417219, 0.3691298, 0.3827116, 0.3518392,
                        0.3514443, 0.3582501, 0.3527597, 0.3583450,
                        0.3463398, 0.3546025, 0.3437715, 0.3407840)




sim_coarse_process <- c(0.4167006, 0.4190169, 0.4328482, 0.4300583,
                        0.4226508, 0.4303639, 0.4093878, 0.4477419,
                        0.4358650, 0.4306661, 0.4249790, 0.4255141,
                        0.4220144, 0.4301920, 0.4247335, 0.4170940,
                        0.4192997, 0.3861053, 0.3954192, 0.3904721,
                        0.4154839, 0.4048027, 0.3992265, 0.3992215)




sim_coarse_noprocess <- c(0.4189048, 0.4055077, 0.4094101, 0.4114187,
                          0.3977770, 0.4248911, 0.4054336, 0.4033207,
                          0.4056407, 0.3866171, 0.4019574, 0.3910235,
                          0.3766105, 0.3928694, 0.4050942, 0.4092044,
                          0.3888697, 0.3834586, 0.3688607, 0.3762242, 
                          0.3841943, 0.3732394, 0.3889650, 0.3706459)






df <- data.frame(y = c(0.4469619, 0.3027712, 0.3217938, 0.2725812, 0.4179107,  0.3950058, 0.2774761, 0.3670556),
                 sd = c(0.03854688, 0.03082228, 0.02388769, 0.0233749, 0.01542097, 0.01554636, 0.07934786, 0.01753022),
                 type = rep(c("W/ Pre-processing", "W/O pre-processing"),4),
                 x = rep(c("PBMC Coarse-resolution", "PBMC Fine-resolution", 
                           "Simulated Coarse-resolution", "Simulated Fine-resolution"), each = 2))



t.test(pbmc_coarse_process, pbmc_coarse_noprocess)
t.test(pbmc_fine_process, pbmc_fine_noprocess)
t.test(sim_coarse_process, sim_coarse_noprocess)
t.test(sim_fine_process, sim_fine_noprocess)




library(ggplot2)
ggplot(df, aes(x = x, y = y, fill = type)) +
  geom_bar(width = 0.7, stat="identity", position = position_dodge(), 
           colour = "black", linewidth = 1) +
  geom_errorbar(aes(ymin = y - 0.5*sd, ymax = y + 0.5*sd), position = position_dodge(0.7), width = 0.5, 
                color = "black", linewidth = 1) +
  labs(title = "",
       x = "",
       y = "Mean Jaccard similarity") +
  scale_fill_manual(values = c("W/ Pre-processing" = "#4575B4" , "W/O pre-processing" = "#D73027")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12), 
    axis.text = element_text(face = "bold", size = 11, color = "black"), 
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    legend.position = "bottom",
    legend.text = element_text(face = "bold", size = 10) 
  ) + 
  geom_hline(yintercept = seq(0, 0.5, by = 0.1), linetype = "dashed", color = "black", size = 0.3) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  guides(fill = guide_legend(title = NULL)) + 

annotate("segment", x = 0.75, xend = 1.25, y = 0.5, yend = 0.5, color = "black", size = 0.8) + 
  annotate("segment", x = 0.75, xend = 0.75, y = 0.5, yend = 0.49, color = "black", size = 0.8) + 
  annotate("segment", x = 1.25, xend = 1.25, y = 0.5, yend = 0.49, color = "black", size = 0.8) +
  annotate("text", x = 1, y = 0.51, label = "***", size = 5, fontface = "bold")  + 
  
  
  annotate("segment", x = 1.75, xend = 2.25, y = 0.5, yend = 0.5, color = "black", size = 0.8) +
  annotate("segment", x = 1.75, xend = 1.75, y = 0.5, yend = 0.49, color = "black", size = 0.8) + 
  annotate("segment", x = 2.25, xend = 2.25, y = 0.5, yend = 0.49, color = "black", size = 0.8) +
  annotate("text", x = 2, y = 0.51, label = "***", size = 5, fontface = "bold") + 
  
  
  annotate("segment", x = 2.75, xend = 3.25, y = 0.5, yend = 0.5, color = "black", size = 0.8) +
  annotate("segment", x = 2.75, xend = 2.75, y = 0.5, yend = 0.49, color = "black", size = 0.8) + 
  annotate("segment", x = 3.25, xend = 3.25, y = 0.5, yend = 0.49, color = "black", size = 0.8) + 
  annotate("text", x = 3, y = 0.51, label = "***", size = 5, fontface = "bold") + 


  
  annotate("segment", x = 3.75, xend = 4.25, y = 0.5, yend = 0.5, color = "black", size = 0.8) +
  annotate("segment", x = 3.75, xend = 3.75, y = 0.5, yend = 0.49, color = "black", size = 0.8) + 
  annotate("segment", x = 4.25, xend = 4.25, y = 0.5, yend = 0.49, color = "black", size = 0.8) + 
  annotate("text", x = 4, y = 0.51, label = "***", size = 5, fontface = "bold") 
