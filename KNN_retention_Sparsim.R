
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
  
  
  fastpca<-prcomp_irlba(log(new_data+1), n = 50, retx = TRUE, center = TRUE, scale. = FALSE)
  position<-fastpca$rotation
  
  knn.graph.cells=nng(position, mutual = F, k= k)   ## NNgraph algorithm based on Euclidean distance
  adjacency.matrix <- igraph::as_adjacency_matrix(knn.graph.cells)
  
  knn_ret_list<-NULL
  for (i in 1:length(cell.info)) {
    knn_ret<-mean(cell.info[i]==cell.info[which(adjacency.matrix[i,]==1)])
    knn_ret_list<-append(knn_ret_list, knn_ret)
  }
  return(mean(knn_ret_list))
}







knn_retention<-NULL
for(i in 1:16){
  load(file =paste0("E:/Thesis/SPARSim/mcadet_Sparsim_", i, ".RData"))
  load(file =paste0("E:/Thesis/SPARSim/mcadet_Sparsim_up", i, ".RData"))
  load(file =paste0("E:/Thesis/SPARSim/mcadet_Sparsim_mid", i, ".RData"))
  load(file =paste0("E:/Thesis/SPARSim/mcadet_Sparsim_down", i, ".RData"))
  load(file =paste0("E:/Thesis/SPARSim/NB_Sparsim_",        i,".RData"))
  load(file =paste0("E:/Thesis/SPARSim/Bre_Sparsim_",       i,".RData")) 
  load(file =paste0("E:/Thesis/SPARSim/Seurat_Sparsim_",    i,".RData"))
  load(file =paste0("E:/Thesis/SPARSim/disp_Sparsim_",      i,".RData"))
  load(file =paste0("E:/Thesis/SPARSim/random_Sparsim_",    i,".RData"))
  
  load(file =paste0("E:/Thesis/SPARSim/SparSim_", i, ".RData"))
  workdata=get(paste0("SparSim_", i))
  
  
  cell.info<-c(rep("Type 1",  60),
               rep("Type 2",  180),
               rep("Type 3",  240),
               rep("Type 4",  360),
               rep("Type 5",  420),
               rep("Type 6",  480),
               rep("Type 7",  600),
               rep("Type 8",  960),
               rep("Type 9",  1200),
               rep("Type 10", 1500))
  
  knn.matrix<-matrix(NA, nrow = 4, ncol = 9)
  for (j in 1:4) {
    start.time<-Sys.time()
    k.set<-c(20,30,40,50)
    

  genelist_mcadet<-get(paste0("mcadet_Sparsim_",i))$genelist[1:1500]
  genelist_mcadet_up<-get(paste0("mcadet_Sparsim_up",i))$genelist[1:1500]
  genelist_mcadet_mid<-get(paste0("mcadet_Sparsim_mid",i))$genelist[1:1500]
  genelist_mcadet_down<-get(paste0("mcadet_Sparsim_down",i))$genelist[1:1500]
  genelist_NB<-get(paste0("NB_Sparsim_",i))[1:1500]
  genelist_Bre<-get(paste0("Bre_Sparsim_",i))[1:1500]
  genelist_Seurat<-get(paste0("Seurat_Sparsim_",i))[1:1500]
  genelist_disp<-get(paste0("disp_Sparsim_",i))[1:1500]
  genelist_random<-get(paste0("random_Sparsim_",i))[1:1500]

  one.time.run<-c(knn_eva(workdata, genelist_mcadet, k.set[j]),
  knn_eva(workdata, genelist_mcadet_up, k.set[j]),
  knn_eva(workdata, genelist_mcadet_mid, k.set[j]),
  knn_eva(workdata, genelist_mcadet_down, k.set[j]),
  knn_eva(workdata, genelist_NB, k.set[j]),
  knn_eva(workdata, genelist_Bre, k.set[j]),
  knn_eva(workdata, genelist_Seurat, k.set[j]),
  knn_eva(workdata, genelist_disp, k.set[j]),
  knn_eva(workdata, genelist_random, k.set[j]))
  
  knn.matrix[j,]<-one.time.run
  
  

  end.time<-Sys.time()
  print(c(i,j))
  print(end.time-start.time)
  }
  
  knn_retention<- append(knn_retention , colMeans(knn.matrix))
  rm(list = paste0("SparSim_", i))
  gc()
}




knn_retention_df_Sparsim<-data.frame(
  "knnret"=knn_retention,
  "method"=
    as.factor(rep(c("Mcadet","Mcadet.up","Mcadet.mid", "Mcadet.down",
                    "NBdrop","Brennecke", "Seurat", "Disp","Random"),16))
)


knn_retention_df_Sparsim$method<-
  factor(knn_retention_df_Sparsim$method, 
         levels = c("Mcadet.down","Mcadet.mid","Mcadet.up", "Mcadet", "Seurat", "Brennecke", "NBdrop", "Disp","Random"))

knn_retention_df_Sparsim_1500<-knn_retention_df_Sparsim

save(knn_retention_df_Sparsim_1500, file = "E:/Thesis/SPARSim/knn_retention_df_Sparsim_1500.RData")



load(file = "E:/Thesis/SPARSim/knn_retention_df_Sparsim_2000.RData")
load(file = "E:/Thesis/SPARSim/knn_retention_df_Sparsim_1500.RData")
load(file = "E:/Thesis/SPARSim/knn_retention_df_Sparsim_1000.RData")
load(file = "E:/Thesis/SPARSim/knn_retention_df_Sparsim_600.RData")
load(file = "E:/Thesis/SPARSim/knn_retention_df_Sparsim_300.RData")
load(file = "E:/Thesis/SPARSim/knn_retention_df_Sparsim_100.RData")






knn_retention_df_Sparsim<-knn_retention_df_Sparsim_2000






library(ggpubr)


remove.list.knn=c(which(knn_retention_df_Sparsim$method=="Mcadet.up"),
                  which(knn_retention_df_Sparsim$method=="Mcadet.mid"),
                  which(knn_retention_df_Sparsim$method=="Mcadet.down"))

### boxplot
B<-ggboxplot(knn_retention_df_Sparsim, x="method", y="knnret",
             color="method", add="jitter") +labs(x="", y="NNgraph retention rate")+
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(knn_retention_df_Sparsim$knnret),
                                         linetype=2) + 
  stat_compare_means(method = "anova",label.y = 0.9, label.x = 1.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",label.y = 0.8,
                     ref.group = ".all.")+theme(legend.position="none")             # Pairwise comparison against all

B<-ggboxplot(knn_retention_df_Sparsim[-remove.list.knn,], x="method", y="knnret",
             color="method", add="jitter") +labs(x="", y="NNgraph retention rate")+
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(knn_retention_df_Sparsim[-remove.list.knn,]$knnret),
                                         linetype=2) + 
  stat_compare_means(method = "anova",label.y = 0.9)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",label.y = 0.8,
                     ref.group = ".all.")+theme(legend.position="none")             # Pairwise comparison against all


B<-ggboxplot(knn_retention_df_Sparsim, x="method", y="knnret",
             color="method", add="jitter") +labs(x="", y="NNgraph retention rate")+
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(knn_retention_df_Sparsim$knnret),
                                         linetype=2) + 
  stat_compare_means(method = "kruskal.test",label.y = 0.9)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 0.8,
                     ref.group = ".all.")+theme(legend.position="none")             # Pairwise comparison against all

B<-ggboxplot(knn_retention_df_Sparsim[-remove.list.knn,], x="method", y="knnret",
             color="method", add="jitter") +labs(x="", y="NNgraph retention rate")+
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(knn_retention_df_Sparsim[-remove.list.knn,]$knnret),
                                         linetype=2) + 
  stat_compare_means(method = "kruskal.test",label.y = 0.9)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 0.8,
                     ref.group = ".all.")+theme(legend.position="none")             # Pairwise comparison against all





### Linechart KNN retention ###

knnmean<-NULL
knnsd<-NULL
for (i in 1:6) {
  number.set<-c(100,300,600,1000,1500,2000)
  df<-get(paste0("knn_retention_df_Sparsim_", number.set[i]))
  knnmean<-append(knnmean,tapply(df$knnret,  factor(df$method), mean))
  knnsd<-append(knnsd,tapply(df$knnret,  factor(df$method), sd))
}

line.chart.knn.df<-data.frame(knnmean, knnsd, "method"=rep(c("Mcadet.down", "Mcadet.mid",
                                                             "Mcadet.up", "Mcadet", "NBdrop",
                                                             "Brennecke", "Seurat", "Disp", "Random"),6),
                              "n_genes"=rep(c(100,300,600,1000,1500,2000), each=9))


line.chart.knn.df$method<-factor(line.chart.knn.df$method, 
                                 levels = c("Mcadet","Mcadet.up", "Mcadet.mid","Mcadet.down",
                                            "Seurat", "NBdrop", "Brennecke", "Disp","Random"))



method.remove.knn<-c( which(line.chart.knn.df$method=="Mcadet.down" ),
                      which(line.chart.knn.df$method=="Mcadet.mid" ), 
                      which(line.chart.knn.df$method=="Mcadet.up" ))


library(ggplot2)

knn.all=ggplot(data = line.chart.knn.df, mapping = aes(x = n_genes, y =knnmean , 
                                                 colour = method, shape=method,linetype = method))+
  geom_line(size=0.8)+geom_point()+theme_bw()+theme(legend.position="none")+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Mean KNN retention rate")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(panel.background = element_blank(), axis.line = element_blank(),
        panel.grid.major = element_blank())+
  scale_x_continuous(limits=c(100,2000), breaks = c(100,300,600,1000,1500,2000))



knn.only=ggplot(data = line.chart.knn.df[-method.remove.knn,], mapping = aes(x = n_genes, y =knnmean , 
                                                      colour = method, shape=method,linetype = method))+
  geom_line(size=0.8)+geom_point()+theme_bw()+theme(legend.position="none")+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Mean NNgraph retention rate")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(panel.background = element_blank(),axis.line = element_blank(),
        panel.grid.major = element_blank())+
  scale_x_continuous(limits=c(100,2000), breaks = c(100,300,600,1000,1500,2000))



library(patchwork)
(i+j)/(knn.only+a)/(c+b)





