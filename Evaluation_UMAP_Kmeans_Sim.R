
rm(list = ls())
gc()

library(patchwork)
library(dplyr)
library(forcats)










## evaluation function function 1 ###

eva<-function(genelist, cell.info, gene.name, workdata, n_components, kmeans.centers, nstart){
  
  library("umap")
  library("uwot")
  library("cluster")
  library("factoextra")
  library("funtimes")
  library("mclust")
  library("aricode")
  
  label<-NULL
  for (i in 1:length(genelist)) {
    label[i]=which(gene.name==genelist[i]) 
  }
  
  new_data<-t(workdata[label,])
  
  set.seed(1234)
  
  umap<-uwot::umap(log(new_data+1), n_components=n_components)
  
  umap.df=as.data.frame(umap)
  
  partition.umap<-kmeans(umap.df,centers = kmeans.centers, nstart = nstart)
  

  sil=silhouette(partition.umap$cluster, dist(umap.df ))
  
  sil.mean.cluster<-NULL
  for (m in 1:kmeans.centers) {
    sil.mean.cluster[m]<-mean(sil[sil[,1]==m,3])
  }
  
  sil.value<-mean(sil.mean.cluster)
  purity.value<-purity(partition.umap$cluster, as.factor(cell.info))$pur
  ari.value<-adjustedRandIndex(partition.umap$cluster,as.factor(cell.info))
  NMI.value<- NMI(partition.umap$cluster,as.factor(cell.info))
  
  
  return(c(sil.value, purity.value, ari.value, NMI.value))
}




 




### simulation Sparsim ######### run 1728 times (8 datasets * 9 methods * 4 metrics * 6 pcs) ####

rm(list = ls())
gc()

for (i in 1:16) {
  load(file =paste0("E:/Thesis/SPARSim/SparSim_", i, ".RData"))
  load(file =paste0("E:/Thesis/SPARSim/mcadet_Sparsim_", i, ".RData"))
  load(file =paste0("E:/Thesis/SPARSim/mcadet_Sparsim_up", i, ".RData"))
  load(file =paste0("E:/Thesis/SPARSim/mcadet_Sparsim_mid", i, ".RData"))
  load(file =paste0("E:/Thesis/SPARSim/mcadet_Sparsim_down", i, ".RData"))
  load(file =paste0("E:/Thesis/SPARSim/NB_Sparsim_",        i,".RData"))
  load(file =paste0("E:/Thesis/SPARSim/Bre_Sparsim_",       i,".RData")) 
  load(file =paste0("E:/Thesis/SPARSim/Seurat_Sparsim_",    i,".RData"))
  load(file =paste0("E:/Thesis/SPARSim/disp_Sparsim_",      i,".RData"))
  load(file =paste0("E:/Thesis/SPARSim/random_Sparsim_",    i,".RData"))
  
}





final.eva.matrix<-matrix(data=NA, nrow = 6, ncol = 576)

for (k in 5:10) {
  
  final.eva<-NULL
  for (i in 1:16) {

      start.time<-Sys.time()
      
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
      
      
      workdata=get(paste0("SparSim_", i))
      gene.name<-rownames(workdata)  


      working.data = workdata[rowSums(workdata)!=0,]          ##  keep genes with positive expression (delete empty genes)##
      clean=T
      if(clean){
        mean.expression<-apply(working.data,1,mean)
        working.data=working.data[mean.expression>=0.005,]
      }
      gene.name.2<-rownames(working.data)               ##  extract gene's name in a list ## 
      n.feature.1=4000
      n.feature.2=2000
      
      gene.range<-colSums(get(paste0("mcadet_Sparsim_", i))$gene.range.matrix )       ## take the sum of rank range for each gene for all runs (equivalent to take the average) 
      gene.min<-colSums(get(paste0("mcadet_Sparsim_", i))$gene.min.matrix )           ## take the sum of minimum rank for each gene for all runs (equivalent to take the average) 
      
      
      order1<-order(gene.range,decreasing = T)[1:n.feature.1]   ## first round of filtering genes
      gene.name.order1<-gene.name.2[order1]
      
      filter.list.1<-gene.min[order1]
      
      up.prob=1
      if(up.prob!=1){
        order2<-order(filter.list.1,decreasing = F)[1:round((n.feature.2*up.prob))]  ## second round of filtering genes
        order3<-order(filter.list.1,decreasing = T)[1:round((n.feature.2*(1-up.prob)))]
        gene.list<-gene.name.order1[c(order2,order3) ]                   ## final features list
      }else{
        order2<-order(filter.list.1,decreasing = F)[1:round((n.feature.2*up.prob))]  ## second round of filtering genes
        genelist.mca<-gene.name.order1[c(order2) ]                   ## final features list
      }
      
      up.prob=0.8
      if(up.prob!=1){
        order2<-order(filter.list.1,decreasing = F)[1:round((n.feature.2*up.prob))]  ## second round of filtering genes
        order3<-order(filter.list.1,decreasing = T)[1:round((n.feature.2*(1-up.prob)))]
        genelist.mca.up<-gene.name.order1[c(order2,order3) ]                   ## final features list
      }else{
        order2<-order(filter.list.1,decreasing = F)[1:round((n.feature.2*up.prob))]  ## second round of filtering genes
        gene.list<-gene.name.order1[c(order2) ]                   ## final features list
      }
      
      
      up.prob=0.5
      if(up.prob!=1){
        order2<-order(filter.list.1,decreasing = F)[1:round((n.feature.2*up.prob))]  ## second round of filtering genes
        order3<-order(filter.list.1,decreasing = T)[1:round((n.feature.2*(1-up.prob)))]
        genelist.mca.mid<-gene.name.order1[c(order2,order3) ]                   ## final features list
      }else{
        order2<-order(filter.list.1,decreasing = F)[1:round((n.feature.2*up.prob))]  ## second round of filtering genes
        gene.list<-gene.name.order1[c(order2) ]                   ## final features list
      }
      
      
      up.prob=0.2
      if(up.prob!=1){
        order2<-order(filter.list.1,decreasing = F)[1:round((n.feature.2*up.prob))]  ## second round of filtering genes
        order3<-order(filter.list.1,decreasing = T)[1:round((n.feature.2*(1-up.prob)))]
        genelist.mca.down<-gene.name.order1[c(order2,order3) ]                   ## final features list
      }else{
        order2<-order(filter.list.1,decreasing = F)[1:round((n.feature.2*up.prob))]  ## second round of filtering genes
        gene.list<-gene.name.order1[c(order2) ]                   ## final features list
      }
      
      
      genelist.NB<-get(paste0("NB_Sparsim_", i))[1:1500]
      genelist.Bre<-get(paste0("Bre_Sparsim_", i))[1:1500]
      genelist.Seurat<-get(paste0("Seurat_Sparsim_", i))[1:1500] 
      
      genelist.disp<-get(paste0("disp_Sparsim_", i))[1:1500]
      genelist.random<-get(paste0("random_Sparsim_", i))[1:1500]
      
      
      one.time.run<-c(eva(genelist=genelist.mca.down, cell.info = cell.info, gene.name = gene.name,
                          workdata=workdata, n_components=k, kmeans.centers=10, nstart=50),
                      eva(genelist=genelist.mca.mid, cell.info = cell.info, gene.name = gene.name,
                          workdata=workdata, n_components=k, kmeans.centers=10, nstart=50),
                      eva(genelist=genelist.mca.up, cell.info = cell.info, gene.name = gene.name,
                          workdata=workdata, n_components=k, kmeans.centers=10, nstart=50),
                      eva(genelist=genelist.mca, cell.info = cell.info, gene.name = gene.name,
                          workdata=workdata, n_components=k, kmeans.centers=10, nstart=50),
                      eva(genelist=genelist.NB, cell.info = cell.info, gene.name = gene.name,
                          workdata=workdata,  n_components=k, kmeans.centers=10, nstart=50),
                      eva(genelist=genelist.Bre, cell.info = cell.info, gene.name = gene.name,
                          workdata=workdata,  n_components=k, kmeans.centers=10, nstart=50),
                      eva(genelist=genelist.Seurat, cell.info = cell.info, gene.name = gene.name,
                          workdata=workdata,  n_components=k, kmeans.centers=10, nstart=50),
                      eva(genelist=genelist.disp, cell.info = cell.info, gene.name = gene.name,
                          workdata=workdata, n_components=k, kmeans.centers=10, nstart=50),
                      eva(genelist=genelist.random, cell.info = cell.info, gene.name = gene.name,
                          workdata=workdata,  n_components=k, kmeans.centers=10, nstart=50))
      
      
      final.eva=append(final.eva,one.time.run)
      end.time<-Sys.time()
      print(c(k,i))
      print(end.time-start.time)
      
      
    
  }
  
  final.eva.matrix[k-4,]<-final.eva
  
  
}

final.umap.Sparsim_1500<-apply(final.eva.matrix,2,mean)



save(final.umap.Sparsim_1500, file = "E:/Thesis/SPARSim/final.umap.Sparsim_1500.RData")



load(file = "E:/Thesis/SPARSim/final.umap.Sparsim_2000.RData")
load(file = "E:/Thesis/SPARSim/final.umap.Sparsim_100.RData")
load(file = "E:/Thesis/SPARSim/final.umap.Sparsim_300.RData")
load(file = "E:/Thesis/SPARSim/final.umap.Sparsim_600.RData")
load(file = "E:/Thesis/SPARSim/final.umap.Sparsim_1000.RData")
load(file = "E:/Thesis/SPARSim/final.umap.Sparsim_1500.RData")



for (j in 1:6) {
  for (i in 1:4) {
  metric.set= c("Silhoutte", "Purity", "ARI", "NMI")
  number.set<-c(100,300,600,1000,1500,2000)
  final<-get(paste0("final.umap.Sparsim_", number.set[j]))
  
  assign(paste0(metric.set[i], "_Sparsim_", number.set[j], ".df"),data.frame(final[seq(i,4*9*16,4)],
                                                 "method"=
                                                   as.factor(rep(c("Mcadet.down","Mcadet.mid","Mcadet.up",
                                                                "Mcadet","NBdrop","Brennecke","Seurat","Disp","Random"),16)) ))
  df<-get(paste0(metric.set[i], "_Sparsim_", number.set[j], ".df"))
  df$method<-factor(df$method, 
                    levels = c("Mcadet.down","Mcadet.mid","Mcadet.up","Mcadet","Seurat", "Brennecke","NBdrop", "Disp","Random"))
  colnames(df)[1]="Y"   
  
  assign(paste0(metric.set[i], "_Sparsim_", number.set[j]  , ".df"),df)
  }
 }


library(ggpubr)

 ### boxplot ###

remove.list=c(which(NMI_Sparsim_2000.df$method=="Mcadet.up"),
              which(NMI_Sparsim_2000.df$method=="Mcadet.mid"),
              which(NMI_Sparsim_2000.df$method=="Mcadet.down"))


d<-ggboxplot(NMI_Sparsim_2000.df[-remove.list,], x="method", y="Y",
             color="method", add="jitter",legend="none") + 
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(NMI_Sparsim_2000.df[-remove.list,]$Y),
                                         linetype=2) + 
  stat_compare_means(method = "anova",label.y = 1.2, label.x = 1.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",label.y = 1.0,
                     ref.group = ".all.")+labs(y="NMI")+labs(x="")      # Pairwise comparison against all



c<-ggboxplot(ARI_Sparsim_2000.df[-remove.list,], x="method", y="Y",
             color="method", add="jitter",legend="none") + 
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(ARI_Sparsim_2000.df[-remove.list,]$Y),
                                         linetype=2) + 
  stat_compare_means(method = "anova",label.y = 1.2, label.x = 1.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",label.y = 1.0,
                     ref.group = ".all.")+labs(y="ARI")+labs(x="")     # Pairwise comparison against all



b<-ggboxplot(Purity_Sparsim_2000.df[-remove.list,], x="method", y="Y",
             color="method", add="jitter",legend="none") + 
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(Purity_Sparsim_2000.df[-remove.list,]$Y),
                                         linetype=2) + 
  stat_compare_means(method = "anova",label.y = 1.2 , label.x = 1.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",label.y = 1.0,
                     ref.group = ".all.")+labs(y="Purity")+labs(x="")      # Pairwise comparison against all




a<-ggboxplot(Silhoutte_Sparsim_2000.df[-remove.list,], x="method", y="Y",
             color="method", add="jitter",legend="none") + 
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(Silhoutte_Sparsim_2000.df[-remove.list,]$Y),
                                         linetype=2) + 
  stat_compare_means(method = "anova",label.y = 1.2 , label.x = 1.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",label.y = 1.0,
                     ref.group = ".all.")+labs(y="Silhoutte score")+labs(x="")      # Pairwise comparison against all



library(patchwork)


(C+B+a)/(b+c+d)


### Linechart of ARI NMI Purity Silhoutte ###


## Silhoutte

silmean<-NULL
silsd<-NULL
for (i in 1:6) {
  number.set<-c(100,300,600,1000,1500,2000)
  df<-get(paste0("Silhoutte", "_Sparsim_", number.set[i]  , ".df"))
  silmean<-append(silmean,tapply(df$Y,  factor(df$method), mean))
  silsd<-append(silsd,tapply(df$Y,  factor(df$method), sd))
}

line.chart.sil.df<-data.frame(silmean, silsd, "method"=rep(c("Mcadet.down", "Mcadet.mid",
                                                            "Mcadet.up", "Mcadet", "NBdrop",
                                                            "Brennecke", "Seurat", "Disp", "Random"),6),
                             "n_genes"=rep(c(100,300,600,1000,1500,2000), each=9))


line.chart.sil.df$method<-factor(line.chart.sil.df$method, 
                                levels = c("Mcadet","Mcadet.up", "Mcadet.mid","Mcadet.down",
                                           "Seurat", "NBdrop", "Brennecke", "Disp","Random"))



method.remove.sil<-c( which(line.chart.sil.df$method=="Mcadet.down" ),
                      which(line.chart.sil.df$method=="Mcadet.mid" ), 
                      which(line.chart.sil.df$method=="Mcadet.up" ))

library(ggplot2)

a=ggplot(data = line.chart.sil.df[-method.remove.sil,], mapping = aes(x = n_genes, y =silmean , 
                                                 colour = method, shape=method,linetype = method))+
  geom_line(size=0.8)+geom_point()+theme_bw()+theme(legend.position="none")+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Mean silhoutte score")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(panel.background = element_blank(),axis.line = element_blank(),
        panel.grid.major = element_blank())+
  scale_x_continuous(limits=c(100,2000), breaks = c(100,300,600,1000,1500,2000))








## Purity

puritymean<-NULL
puritysd<-NULL
for (i in 1:6) {
  number.set<-c(100,300,600,1000,1500,2000)
  df<-get(paste0("Purity", "_Sparsim_", number.set[i]  , ".df"))
  puritymean<-append(puritymean,tapply(df$Y,  factor(df$method), mean))
  puritysd<-append(puritysd,tapply(df$Y,  factor(df$method), sd))
}

line.chart.purity.df<-data.frame(puritymean, puritysd, "method"=rep(c("Mcadet.down", "Mcadet.mid",
                                                                      "Mcadet.up", "Mcadet", "NBdrop",
                                                                      "Brennecke", "Seurat", "Disp", "Random"),6),
                                 "n_genes"=rep(c(100,300,600,1000,1500,2000), each=9))


line.chart.purity.df$method<-factor(line.chart.purity.df$method, 
                                    levels = c("Mcadet","Mcadet.up", "Mcadet.mid","Mcadet.down",
                                               "Seurat", "NBdrop", "Brennecke", "Disp","Random"))



method.remove.purity<-c( which(line.chart.purity.df$method=="Mcadet.down" ),
                         which(line.chart.purity.df$method=="Mcadet.mid" ), 
                         which(line.chart.purity.df$method=="Mcadet.up" ))

library(ggplot2)

b=ggplot(data = line.chart.purity.df[-method.remove.purity,], mapping = aes(x = n_genes, y =puritymean , 
                                                    colour = method, shape=method,linetype = method))+
  geom_line(size=0.8)+geom_point()+theme_bw()+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Mean purity")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(panel.background = element_blank(),axis.line = element_blank(),
        panel.grid.major = element_blank())+
  scale_x_continuous(limits=c(100,2000), breaks = c(100,300,600,1000,1500,2000))









## ARI

arimean<-NULL
arisd<-NULL
for (i in 1:6) {
  number.set<-c(100,300,600,1000,1500,2000)
  df<-get(paste0("ARI", "_Sparsim_", number.set[i]  , ".df"))
  arimean<-append(arimean,tapply(df$Y,  factor(df$method), mean))
  arisd<-append(arisd,tapply(df$Y,  factor(df$method), sd))
}

line.chart.ari.df<-data.frame(arimean, arisd, "method"=rep(c("Mcadet.down", "Mcadet.mid",
                                                             "Mcadet.up", "Mcadet", "NBdrop",
                                                             "Brennecke", "Seurat", "Disp", "Random"),6),
                              "n_genes"=rep(c(100,300,600,1000,1500,2000), each=9))


line.chart.ari.df$method<-factor(line.chart.ari.df$method, 
                                 levels = c("Mcadet","Mcadet.up", "Mcadet.mid","Mcadet.down",
                                            "Seurat", "NBdrop", "Brennecke", "Disp","Random"))



method.remove.ari<-c( which(line.chart.ari.df$method=="Mcadet.down" ),
                      which(line.chart.ari.df$method=="Mcadet.mid" ), 
                      which(line.chart.ari.df$method=="Mcadet.up" ))

library(ggplot2)

c=ggplot(data = line.chart.ari.df[-method.remove.ari,], mapping = aes(x = n_genes, y =arimean , 
                                                 colour = method, shape=method,linetype = method))+
  geom_line(size=0.8)+geom_point()+theme_bw()+theme(legend.position="none")+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Mean ARI")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(panel.background = element_blank(),axis.line = element_blank(),
        panel.grid.major = element_blank())+
  scale_x_continuous(limits=c(100,2000), breaks = c(100,300,600,1000,1500,2000))









## NMI

nmimean<-NULL
nmisd<-NULL
for (i in 1:6) {
  number.set<-c(100,300,600,1000,1500,2000)
  df<-get(paste0("NMI", "_Sparsim_", number.set[i]  , ".df"))
  nmimean<-append(nmimean,tapply(df$Y,  factor(df$method), mean))
  nmisd<-append(nmisd,tapply(df$Y,  factor(df$method), sd))
}

line.chart.nmi.df<-data.frame(nmimean, nmisd, "method"=rep(c("Mcadet.down", "Mcadet.mid",
                                                             "Mcadet.up", "Mcadet", "NBdrop",
                                                             "Brennecke", "Seurat", "Disp", "Random"),6),
                              "n_genes"=rep(c(100,300,600,1000,1500,2000), each=9))


line.chart.nmi.df$method<-factor(line.chart.nmi.df$method, 
                                 levels = c("Mcadet","Mcadet.up", "Mcadet.mid","Mcadet.down",
                                            "Seurat", "NBdrop", "Brennecke", "Disp","Random"))



method.remove.nmi<-c( which(line.chart.nmi.df$method=="Mcadet.down" ),
                      which(line.chart.nmi.df$method=="Mcadet.mid" ), 
                      which(line.chart.nmi.df$method=="Mcadet.up" ))

library(ggplot2)

d=ggplot(data = line.chart.nmi.df[-method.remove.nmi,], mapping = aes(x = n_genes, y =nmimean , 
                                                 colour = method, shape=method,linetype = method))+
  geom_line(size=0.8)+geom_point()+theme_bw()+theme(legend.position="none")+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Mean NMI")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(panel.background = element_blank(),axis.line = element_blank(),
        panel.grid.major = element_blank())+
  scale_x_continuous(limits=c(100,2000), breaks = c(100,300,600,1000,1500,2000))


h=ggplot(data = line.chart.nmi.df[-method.remove.nmi,], mapping = aes(x = n_genes, y =nmimean, 
                                                                      colour = method, shape=method,linetype = method))+
  geom_line(size=0.8)+geom_point()+theme_bw()+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Mean NMI")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(panel.background = element_blank(),axis.line = element_blank(),
        panel.grid.major = element_blank())+
  scale_x_continuous(limits=c(100,2000), breaks = c(100,300,600,1000,1500,2000))






library(patchwork)

(linechart_Sparsim+knn.only)/(a+b)/(c+d)
















