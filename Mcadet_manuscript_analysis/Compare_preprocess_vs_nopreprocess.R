
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
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/DE_genes/DEgene_fine_", 
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
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/pbmc_results/Mcadet_pbmc_fine_", 
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


for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/data_fine/pbmc_fine_", 
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
  result = mcadet_nopreprocess(data = dta)
  assign(paste0("Mcadet_noprocess_SimFine_", i),  result)
  
  
  dta = get(paste0("SparSim_",  i))
  result = mcadet_nopreprocess(data = dta)
  assign(paste0("Mcadet_noprocess_SimCoarse", i),  result)
  
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
    assign(paste0("Mcadet_noprocess_pbmc_coarse_", donor.set[p], "_", time.set[k]),  result)
    gc()
    
    dta = get(paste0("pbmc_fine_",  donor.set[p], "_", time.set[k]))
    result = mcadet_nopreprocess(data = dta)
    assign(paste0("Mcadet_noprocess_pbmc_fine_", donor.set[p], "_", time.set[k]),  result)
    gc()
    
    
  }
}




####### Jaccard 




pbmc_jaccard_macdet_process_coarse = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_", donor.set[p], "_", time.set[k])),
                  get(paste0("Mcadet_pbmc_", donor.set[p], "_", time.set[k]))$gene )
    pbmc_jaccard_macdet_process_coarse = append(pbmc_jaccard_macdet_process_coarse, jac)
  }
}
pbmc_jaccard_macdet_process_coarse





pbmc_jaccard_macdet_noprocess_coarse = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_", donor.set[p], "_", time.set[k])),
                  get(paste0("Mcadet_noprocess_pbmc_coarse_", donor.set[p], "_", time.set[k]))$gene )
    pbmc_jaccard_macdet_noprocess_coarse = append(pbmc_jaccard_macdet_noprocess_coarse, jac)
  }
}
pbmc_jaccard_macdet_noprocess_coarse







pbmc_jaccard_macdet_process_fine = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_fine_", donor.set[p], "_", time.set[k])),
                  get(paste0("Mcadet_pbmc_fine_", donor.set[p], "_", time.set[k]))$gene )
    pbmc_jaccard_macdet_process_fine = append(pbmc_jaccard_macdet_process_fine, jac)
  }
}
pbmc_jaccard_macdet_process_fine





pbmc_jaccard_macdet_noprocess_fine = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_fine_", donor.set[p], "_", time.set[k])),
                  get(paste0("Mcadet_noprocess_pbmc_fine_", donor.set[p], "_", time.set[k]))$gene )
    pbmc_jaccard_macdet_noprocess_fine = append(pbmc_jaccard_macdet_noprocess_fine, jac)
  }
}
pbmc_jaccard_macdet_noprocess_fine






Sim_jaccard_macdet_process_coarse = c()
for (i in 1:24) {
  jac = jaccard(paste0("Gene", c(1:2000)),
                get(paste0("Mcadet_SparSim_", i))$gene)
  Sim_jaccard_macdet_process_coarse  = append(Sim_jaccard_macdet_process_coarse, jac)
  
}
Sim_jaccard_macdet_process_coarse


Sim_jaccard_macdet_process_fine = c()
for (i in 1:24) {
  jac = jaccard(paste0("Gene", c(1:2000)),
                get(paste0("Mcadet_SparSim_fine_", i))$gene)
  Sim_jaccard_macdet_process_fine  = append(Sim_jaccard_macdet_process_fine, jac)
  
}
Sim_jaccard_macdet_process_fine




Sim_jaccard_macdet_noprocess_coarse = c()
for (i in 1:24) {
    jac = jaccard(paste0("Gene", c(1:2000)),
                  get(paste0("Mcadet_noprocess_SimCoarse", i))$gene)
    Sim_jaccard_macdet_noprocess_coarse  = append(Sim_jaccard_macdet_noprocess_coarse, jac)
  
}
Sim_jaccard_macdet_noprocess_coarse





Sim_jaccard_macdet_noprocess_fine = c()
for (i in 1:24) {
  jac = jaccard(paste0("Gene", c(1:2000)),
                get(paste0("Mcadet_noprocess_SimFine_", i))$gene)
  Sim_jaccard_macdet_noprocess_fine  = append(Sim_jaccard_macdet_noprocess_fine, jac)
  
}
Sim_jaccard_macdet_noprocess_fine




mean(pbmc_jaccard_macdet)
mean(pbmc_jaccard_macdet_noprocess)



############### Plot


pbmc_coarse_process <- c(0.3182518, 0.3021219, 0.3302665, 0.3492437,
                         0.3072566, 0.2929751, 0.3413094, 0.3997361, 
                         0.3178476, 0.3376323, 0.3408403, 0.3304862,
                         0.3074992, 0.3258290, 0.3224285, 0.2886717,
                         0.3050751, 0.3243002, 0.3110410, 0.2847073,
                         0.3234855, 0.3289435, 0.3322543, 0.3008487)



pbmc_coarse_noprocess <- c(0.3109573, 0.2620896, 0.2834367, 0.2972386, 0.3008959,
                           0.3024879, 0.2781555, 0.2770865, 0.2808343, 0.2335329,
                           0.2838534, 0.2624925, 0.2850521, 0.2945516, 0.3004940,
                           0.2527106, 0.2779977, 0.2505587, 0.2723960, 0.2551667,
                           0.3068737, 0.2306416, 0.2365294, 0.2620807)




pbmc_fine_process <- c(0.4241026, 0.4465593, 0.4440594, 0.4520214,
                       0.3937210, 0.4754522, 0.3675799, 0.4197609,
                       0.4635993, 0.4582624, 0.4439980, 0.4432071,
                       0.4266350, 0.3898477, 0.4984843, 0.4631937,
                       0.4750831, 0.4119980, 0.3931596, 0.4685694,
                       0.4632325, 0.4768519, 0.4992248, 0.5284810)



pbmc_fine_noprocess <- c(0.2591588, 0.2719424, 0.2865540, 0.3913223, 0.3129594,
                         0.3542435, 0.2623182, 0.3413274, 0.3349562, 0.2780251, 0.2902511, 0.3062815,
                         0.3031268, 0.2988193, 0.2970677, 0.2893284, 0.3231132, 0.2363442,
                         0.2991033, 0.2926554, 0.3167082, 0.2900133, 0.2963155, 0.2997516)




sim_fine_process <- c(0.3766797, 0.3891085, 0.3311772, 0.3826731, 
                      0.3640275, 0.3374778, 0.3445415, 0.3609520,
                      0.3153915, 0.3082312, 0.3423110, 0.3082969,
                      0.2576153, 0.2368889, 0.2496656, 0.2180923,
                      0.1726714, 0.2390049, 0.2321509, 0.2275556, 
                      0.2087766, 0.1980969, 0.1266726, 0.1313673)



sim_fine_noprocess <- c(0.3844581, 0.3788862, 0.3785739, 0.3695946,
                        0.3779070, 0.4009555, 0.3753453, 0.3903915,
                        0.3679595, 0.3924408, 0.3657439, 0.3791519,
                        0.3411648, 0.3641398, 0.3749559, 0.3507463,
                        0.3546592, 0.3569966, 0.3533471, 0.3624017,
                        0.3473137, 0.3506897, 0.3442341, 0.3503801)




sim_coarse_process <- c(0.4167006, 0.4190169, 0.4328482, 0.4300583,
                        0.4226508, 0.4303639, 0.4093878, 0.4477419,
                        0.4358650, 0.4306661, 0.4249790, 0.4255141,
                        0.4220144, 0.4301920, 0.4247335, 0.4170940,
                        0.4192997, 0.3861053, 0.3954192, 0.3904721,
                        0.4154839, 0.4048027, 0.3992265, 0.3992215)




sim_coarse_noprocess <- c(0.4211263, 0.4065858, 0.4125973, 0.4110727,
                          0.3987895, 0.4268026, 0.4087719, 0.4072651,
                          0.4088998, 0.3830214, 0.4074982, 0.3904762,
                          0.3775273, 0.3933356, 0.4034531, 0.4095612,
                          0.3856456, 0.3832765, 0.3690756, 0.3753360,
                          0.3896907, 0.3716667, 0.3977669, 0.3783970)



sd(pbmc_coarse_process)
sd(pbmc_coarse_noprocess)
sd(pbmc_fine_process)
sd(pbmc_fine_noprocess)
sd(sim_coarse_process)
sd(sim_coarse_noprocess)
sd(sim_fine_process)
sd(sim_fine_noprocess)





df <- data.frame(y = c(0.3217938 , 0.2749214, 0.4469619, 0.3013203, 0.4179107,  0.3950058, 0.2774761, 0.3670556),
                 sd = c(0.02388769, 0.02351793, 0.03854688, 0.03242294, 0.01542097, 0.01554636, 0.07934786, 0.01753022),
                 type = rep(c("W/ Pre-processing", "W/O pre-processing"),4),
                 x = rep(c("PBMC Coarse-resolution", "PBMC Fine-resolution", 
                           "Simulated Coarse-resolution", "Simulated Fine-resolution"), each = 2))



t.test(pbmc_coarse_process, pbmc_coarse_noprocess, alternative = "greater")
t.test(pbmc_fine_process, pbmc_fine_noprocess, alternative = "greater")
t.test(sim_coarse_process, sim_coarse_noprocess, alternative = "greater")
t.test(sim_fine_process, sim_fine_noprocess, alternative = "greater")




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
  annotate("text", x = 4, y = 0.51, label = "NS", size = 5, fontface = "bold") 
