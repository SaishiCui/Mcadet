
library("igraph")
library("cccd")
library("irlba")






knn_eva<-function(workdata, genelist, k){
  
  gene.name<-rownames(workdata)  
  
  label<-NULL
  for (i in 1:length(genelist)) {
    label[i]=which(gene.name==genelist[i]) 
  }
  
  new_data<-workdata[label,]
  
  
  fastpca<-prcomp_irlba(log(new_data+1), n = 50, retx = TRUE, center = TRUE, scale. = F)
  position<-fastpca$rotation
  
  knn.graph.cells=nng(position, mutual = F, k= k)   ## NNgraph algorithm based on Euclidean distance
  adjacency.matrix <- igraph::as_adjacency_matrix(knn.graph.cells)
  
  knn_ret_list<-NULL
  for (n in 1:length(cell.info)) {
    knn_ret<-mean(cell.info[n]==cell.info[which(adjacency.matrix[n,]==1)])
    knn_ret_list<-append(knn_ret_list, knn_ret)
  }
  return(mean(knn_ret_list))
}







knn_retention<-NULL
for(i in 1:8){
  for (j in 1:3) {

    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    
    
  load(file =paste0("E:/Thesis/pbmc/pbmc_fine_clabel_", donor.set[i], "_", time.set[j],".RData"))
  load(file =paste0("E:/Thesis/pbmc/mcadet_", donor.set[i], "_", time.set[j],".RData"))
  load(file =paste0("E:/Thesis/pbmc/mcadet_fine_", donor.set[i], "_", time.set[j],".RData"))
  load(file =paste0("E:/Thesis/pbmc/NB_fine_All_",     donor.set[i], "_", time.set[j],".RData"))
  load(file =paste0("E:/Thesis/pbmc/Bre_fine_All_",    donor.set[i], "_", time.set[j],".RData"))
  load(file =paste0("E:/Thesis/pbmc/Seurat_fine_All_", donor.set[i], "_", time.set[j],".RData"))
  load(file =paste0("E:/Thesis/pbmc/disp_fine_All_", donor.set[i], "_", time.set[j],".RData"))
  load(file =paste0("E:/Thesis/pbmc/random_fine_All_", donor.set[i], "_", time.set[j],".RData"))
  
  load(file= paste0("E:/Thesis/pbmc/pbmc_fine_",  donor.set[i], "_", time.set[j],".RData"))
  workdata=get(paste0("pbmc_fine_", donor.set[i], "_", time.set[j]))
  
  
  cell.info<-get(paste0("pbmc_fine_clabel_", donor.set[i], "_", time.set[j]))
  
  knn.matrix<-matrix(NA, nrow = 4, ncol = 7)
  for (m in 1:4) {
    start.time<-Sys.time()
    k.set<-c(20,30,40,50)
    
    genelist_mcadet_fine<-get(paste0("mcadet_fine_", donor.set[i], "_", time.set[j]))$genelist[1:100]
    genelist_mcadet<-get(paste0("mcadet_", donor.set[i], "_", time.set[j]))$genelist[1:100]
    genelist_NB<-get(paste0("NB_fine_All_",donor.set[i], "_", time.set[j]))[1:100]
    genelist_Bre<-get(paste0("Bre_fine_All_",donor.set[i], "_", time.set[j]))[1:100]
    genelist_Seurat<-get(paste0("Seurat_fine_All_",donor.set[i], "_", time.set[j]))[1:100]
    genelist_disp<-get(paste0("disp_fine_All_",donor.set[i], "_", time.set[j]))[1:100]
    genelist_random<-get(paste0("random_fine_All_",donor.set[i], "_", time.set[j]))[1:100]
    
    one.time.run<-c(knn_eva(workdata, genelist_mcadet_fine, k.set[m]),
                    knn_eva(workdata, genelist_mcadet, k.set[m]),
                    knn_eva(workdata, genelist_NB, k.set[m]),
                    knn_eva(workdata, genelist_Bre, k.set[m]),
                    knn_eva(workdata, genelist_Seurat, k.set[m]),
                    knn_eva(workdata, genelist_disp, k.set[m]),
                    knn_eva(workdata, genelist_random, k.set[m]))
    
    knn.matrix[m,]<-one.time.run
    end.time<-Sys.time()
    print(c(i,j,m))
    print(end.time-start.time)
  }
  
  knn_retention<- append(knn_retention , colMeans(knn.matrix))
  rm(list = paste0("pbmc_fine_",  donor.set[i], "_", time.set[j]))
  gc()

  }
}



knn_retention_df_pbmc<-data.frame(
  "knnret"=knn_retention,
  "method"=
    as.factor(rep(c("Mcadet_fine","Mcadet","NBdrop","Brennecke", "Seurat", "Disp","Random"),24))
)


knn_retention_df_pbmc$method<-
  factor(knn_retention_df_pbmc$method, 
         levels = c("Mcadet_fine","Mcadet", "Seurat", "NBdrop", "Disp","Brennecke","Random"))


knn_retention_df_pbmc_100_fine<-knn_retention_df_pbmc

save(knn_retention_df_pbmc_100_fine, file = "E:/Thesis/pbmc/knn_retention_df_pbmc_100_fine.RData")
save(knn_retention_df_pbmc_300_fine, file = "E:/Thesis/pbmc/knn_retention_df_pbmc_300_fine.RData")
save(knn_retention_df_pbmc_600_fine, file = "E:/Thesis/pbmc/knn_retention_df_pbmc_600_fine.RData")
save(knn_retention_df_pbmc_1000_fine, file = "E:/Thesis/pbmc/knn_retention_df_pbmc_1000_fine.RData")
save(knn_retention_df_pbmc_1500_fine, file = "E:/Thesis/pbmc/knn_retention_df_pbmc_1500_fine.RData")
save(knn_retention_df_pbmc_2000_fine, file = "E:/Thesis/pbmc/knn_retention_df_pbmc_2000_fine.RData")


load(file = "E:/Thesis/pbmc/knn_retention_df_pbmc_100_fine.RData")
load(file = "E:/Thesis/pbmc/knn_retention_df_pbmc_2000_fine.RData")
load(file = "E:/Thesis/pbmc/knn_retention_df_pbmc_1500_fine.RData")
load(file = "E:/Thesis/pbmc/knn_retention_df_pbmc_1000_fine.RData")
load(file = "E:/Thesis/pbmc/knn_retention_df_pbmc_600_fine.RData")
load(file = "E:/Thesis/pbmc/knn_retention_df_pbmc_300_fine.RData")









library(ggpubr)



### boxplot
B<-ggboxplot(knn_retention_df_pbmc_2000, x="method", y="knnret",
             color="method", add="jitter") +labs(x="", y="NNgraph retention rate")+
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(knn_retention_df_pbmc_2000$knnret),
                                         linetype=2) +
  stat_compare_means(method = "kruskal.test",label.y = 1.2, label.x = 1.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 1.1,
                     ref.group = ".all.")+theme(legend.position="none")             # Pairwise comparison against all


knn_retention_df_pbmc_2000$method<-factor(knn_retention_df_pbmc_2000$method, 
                  levels = c("Mcadet","Seurat", "NBdrop", "Disp","Brennecke","Random"))
B<-ggboxplot(knn_retention_df_pbmc_2000, x="method", y="knnret",
             color="method", add="jitter") +labs(x="", y="NNgraph retention rate")+
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(knn_retention_df_pbmc_2000$knnret),
                                         linetype=2) + 
  stat_compare_means(method = "kruskal.test",label.y = 1.15)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 1.0,
                     ref.group = ".all.")+theme(legend.position="none")             # Pairwise comparison against all





### Linechart KNN retention ###

knnmean<-NULL
knnsd<-NULL
for (i in 1:6) {
  number.set<-c(100,300,600,1000,1500,2000)
  df<-get(paste0("knn_retention_df_pbmc_", number.set[i], "_fine"))
  knnmean<-append(knnmean,tapply(df$knnret,  factor(df$method), mean))
  knnsd<-append(knnsd,tapply(df$knnret,  factor(df$method), sd))
}

line.chart.knn.df<-data.frame(knnmean, knnsd, "method"=rep(c( "Mcadet_fine","Mcadet", "NBdrop",
                                                             "Brennecke", "Seurat", "Disp", "Random"),6),
                              "n_genes"=rep(c(100,300,600,1000,1500,2000), each=7))


line.chart.knn.df$method<-factor(line.chart.knn.df$method, 
                                 levels = c("Mcadet_fine","Mcadet","Seurat", "NBdrop","Disp", "Brennecke", "Random"))





library(ggplot2)



knn.only=ggplot(data = line.chart.knn.df, mapping = aes(x = n_genes, y =knnmean , 
                                                                             colour = method, shape=method,linetype = method))+
  geom_line(size=0.8)+geom_point()+theme_bw()+theme(legend.position="none")+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Mean NNgraph retention rate")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(panel.background = element_blank(),axis.line = element_blank(),
        panel.grid.major = element_blank())+
  scale_x_continuous(limits=c(100,2000), breaks = c(100,300,600,1000,1500,2000))









