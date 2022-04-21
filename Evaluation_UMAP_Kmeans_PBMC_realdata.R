rm(list = ls())
gc()

library(patchwork)
library(dplyr)
library(forcats)

### load needed datasets ###

for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =paste0("E:/Thesis/pbmc/pbmc_fine_clabel_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/pbmccell2_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/pbmccell3_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/mcadet_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/mcadet_fine_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/NB_",     donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/Bre_",    donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/Seurat_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/disp_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/random_", donor.set[i], "_", time.set[j],".RData"))
  }
}






## evaluation function function  ###

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




### real data pbmc 3456 runs (24 datasets * 7 methods * 4 metrics * 6 pcs) 

final.eva.matrix<-matrix(data=NA, nrow = 6, ncol = 672)
for (k in 5:10) {
  final.eva<-NULL
  for (i in 1:8) {
    for (j in 1:3) {
      start.time<-Sys.time()
      donor.set<-c(1,2,3,4,5,6,7,8)
      time.set<-c(0,3,7)
      load(file =paste0("E:/Thesis/pbmc/pbmc_fine_",  donor.set[i], "_", time.set[j],".RData"))
      
      cell.info<-get(paste0("pbmc_fine_clabel_", donor.set[i], "_", time.set[j]))
      workdata=get(paste0("pbmc_fine_", donor.set[i], "_", time.set[j]))
      gene.name<-rownames(workdata)  
      
      genelist.mca_fine<-get(paste0("mcadet_fine_", donor.set[i], "_", time.set[j]))$genelist[1:100]
      genelist.mca<-get(paste0("mcadet_", donor.set[i], "_", time.set[j]))$genelist[1:100]
      genelist.NB<-get(paste0("NB_", donor.set[i], "_", time.set[j]))[1:100]
      genelist.Bre<-get(paste0("Bre_", donor.set[i], "_", time.set[j]))[1:100] 
      genelist.Seurat<-get(paste0("Seurat_", donor.set[i], "_",time.set[j]))[1:100]
      genelist.disp<-get(paste0("disp_", donor.set[i], "_",time.set[j]))[1:100]
      genelist.random<-get(paste0("random_", donor.set[i], "_",time.set[j]))[1:100]
      
      one.time.run<-c(eva(genelist=genelist.mca_fine, cell.info = cell.info, gene.name = gene.name,
                          workdata=workdata, n_components=k, kmeans.centers=14, nstart=50),
                      eva(genelist=genelist.mca, cell.info = cell.info, gene.name = gene.name,
                          workdata=workdata, n_components=k, kmeans.centers=14, nstart=50),
                      eva(genelist=genelist.NB, cell.info = cell.info, gene.name = gene.name,
                          workdata=workdata,  n_components=k, kmeans.centers=14, nstart=50),
                      eva(genelist=genelist.Bre, cell.info = cell.info, gene.name = gene.name,
                          workdata=workdata,  n_components=k, kmeans.centers=14, nstart=50),
                      eva(genelist=genelist.Seurat, cell.info = cell.info, gene.name = gene.name,
                          workdata=workdata,  n_components=k, kmeans.centers=14, nstart=50),
                      eva(genelist=genelist.disp, cell.info = cell.info, gene.name = gene.name,
                          workdata=workdata, n_components=k, kmeans.centers=14, nstart=50),
                      eva(genelist=genelist.random, cell.info = cell.info, gene.name = gene.name,
                          workdata=workdata,  n_components=k, kmeans.centers=14, nstart=50))
      
      rm(list = paste0("pbmc_fine_", donor.set[i], "_", time.set[j]))
      gc()
      
      
      final.eva<-append(final.eva, one.time.run)
      end.time<-Sys.time()
      print(c(k,i,j))
      print(end.time-start.time)
      
    }
  }
  final.eva.matrix[k-4,]<-final.eva
  
}


final.umap.pbmc_100_fine<-apply(final.eva.matrix,2,mean)


save(final.umap.pbmc_100_fine, file="E:/Thesis/pbmc/final.umap.pbmc_100_fine.RData")


load(file="E:/Thesis/pbmc/final.umap.pbmc_2000_fine.RData")


load(file="E:/Thesis/pbmc/final.umap.pbmc_100.RData" )
load(file="E:/Thesis/pbmc/final.umap.pbmc_300.RData" )
load(file="E:/Thesis/pbmc/final.umap.pbmc_600.RData" )
load(file="E:/Thesis/pbmc/final.umap.pbmc_1000.RData" )
load(file="E:/Thesis/pbmc/final.umap.pbmc_1500.RData" )
load(file="E:/Thesis/pbmc/final.umap.pbmc_2000.RData" )

load(file="E:/Thesis/pbmc/final.umap.pbmc_100_fine.RData" )
load(file="E:/Thesis/pbmc/final.umap.pbmc_300_fine.RData" )
load(file="E:/Thesis/pbmc/final.umap.pbmc_600_fine.RData" )
load(file="E:/Thesis/pbmc/final.umap.pbmc_1000_fine.RData" )
load(file="E:/Thesis/pbmc/final.umap.pbmc_1500_fine.RData" )
load(file="E:/Thesis/pbmc/final.umap.pbmc_2000_fine.RData" )



    ### split into 4 metrics datasets for each number of genes ###


for (j in 1:6) {
  
  for (i in 1:4) {
    metric.set= c("Silhoutte", "Purity", "ARI", "NMI")
    number.set<-c(100,300,600,1000,1500,2000)
    final<-get(paste0("final.umap.pbmc_", number.set[j], "_fine" ) )
    assign(paste0(metric.set[i], "_pbmc_", number.set[j], ".df"),
           data.frame(final[seq(i, 4*7*8*3,4)],
                      "method"=as.factor(rep(c("Mcadet_fine","Mcadet","NBdrop","Brennecke","Seurat","Disp","Random"),24)),
                      "donor"=c(rep(c(rep("donor 1",7), rep("donor 2",7), rep("donor 3",7),
                                      rep("donor 4",7), rep("donor 5",7), rep("donor 7",7), 
                                      rep("donor 6",7), rep("donor 8",7)),3)),
                      "time"=c(rep(0,56),rep(3,56),rep(7,56))))
    
    df<-get(paste0(metric.set[i], "_pbmc_", number.set[j], ".df"))
    df$method<-factor(df$method, 
                      levels = c("Mcadet_fine","Mcadet","Seurat", "NBdrop", "Disp","Brennecke","Random"))
    colnames(df)[1]="Y"   
    assign(paste0(metric.set[i], "_pbmc_", number.set[j], ".df"),df)
  }
  
}





   #### boxplot tailored for 2,000 gene list ###


library(ggpubr)


d<-ggboxplot(NMI_pbmc_2000.df[-which(NMI_pbmc_2000.df$method=="Mcadet"),] , x="method", y="Y",
             color="method", add="jitter",legend="none") + 
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(NMI_pbmc_2000.df[-which(NMI_pbmc_2000.df$method=="Mcadet"),]$Y),
                                         linetype=2) + 
  stat_compare_means(method = "kruskal.test",label.y = 1.15)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 0.97,
                     ref.group = ".all.")+labs(x="", y="NMI")     # Pairwise comparison against all





c<-ggboxplot(ARI_pbmc_2000.df[-which(NMI_pbmc_2000.df$method=="Mcadet"),], x="method", y="Y",
             color="method", add="jitter",legend="none") + 
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(ARI_pbmc_2000.df[-which(NMI_pbmc_2000.df$method=="Mcadet"),]$Y),
                                         linetype=2) + 
  stat_compare_means(method = "kruskal.test",label.y = 1.05)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 0.92,
                     ref.group = ".all.")+labs(x="",y="ARI")     # Pairwise comparison against all



b<-ggboxplot(Purity_pbmc_2000.df[-which(NMI_pbmc_2000.df$method=="Mcadet"),], x="method", y="Y",
             color="method", add="jitter",legend="none") + 
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(Purity_pbmc_2000.df[-which(NMI_pbmc_2000.df$method=="Mcadet"),]$Y),
                                         linetype=2) + 
  stat_compare_means(method = "kruskal.test",label.y = 1.05)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 0.92,
                     ref.group = ".all.")+labs(x="",y="Purity")     # Pairwise comparison against all




a<-ggboxplot(Silhoutte_pbmc_2000.df, x="method", y="Y",
             color="method", add="jitter",legend="none") + 
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(Silhoutte_pbmc_2000.df$Y),
                                         linetype=2) + 
  stat_compare_means(method = "kruskal.test",label.y = 1)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 0.9,
                     ref.group = ".all.")+labs(x="", y="Silhoutte score")     # Pairwise comparison against all




d<-ggboxplot(NMI_pbmc_2000.df , x="method", y="Y",
             color="method", add="jitter",legend="none") + 
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(NMI_pbmc_2000.df$Y),
                                         linetype=2) +
  stat_compare_means(method = "kruskal.test",label.y = 1.15)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 0.97,
                     ref.group = ".all.")+labs(x="", y="NMI")     # Pairwise comparison against all





c<-ggboxplot(ARI_pbmc_2000.df, x="method", y="Y",
             color="method", add="jitter",legend="none") + 
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(ARI_pbmc_2000.df$Y),
                                         linetype=2) + 
  stat_compare_means(method = "kruskal.test",label.y = 1.05)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 0.92,
                     ref.group = ".all.")+labs(x="",y="ARI")     # Pairwise comparison against all



b<-ggboxplot(Purity_pbmc_2000.df, x="method", y="Y",
             color="method", add="jitter",legend="none") + 
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(Purity_pbmc_2000.df$Y),
                                         linetype=2) +
  stat_compare_means(method = "kruskal.test",label.y = 1.05)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 0.92,
                     ref.group = ".all.")+labs(x="",y="Purity")     # Pairwise comparison against all




a<-ggboxplot(Silhoutte_pbmc_2000.df, x="method", y="Y",
             color="method", add="jitter",legend="none") + 
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(Silhoutte_pbmc_2000.df$Y),
                                         linetype=2) + 
  stat_compare_means(method = "kruskal.test",label.y = 1.05)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 0.95,
                     ref.group = ".all.")+labs(x="", y="Silhoutte score")     # Pairwise comparison against all







library(patchwork)
(C+B+a)/(b+c+d)

(C+B+a)/(b+c+d)







### PBMC Linechart of ARI NMI Purity Silhoutte ###


  ## Silhoutte

silmean<-NULL
silsd<-NULL
for (i in 1:6) {
  number.set<-c(100,300,600,1000,1500,2000)
  df<-get(paste0("Silhoutte", "_pbmc_", number.set[i]  , ".df"))
  silmean<-append(silmean,tapply(df$Y,  factor(df$method), mean))
  silsd<-append(silsd,tapply(df$Y,  factor(df$method), sd))
}

line.chart.sil.df<-data.frame(silmean, silsd, "method"=rep(c("Mcadet_fine", "Mcadet", "NBdrop",
                                                             "Brennecke", "Seurat", "Disp", "Random"),6),
                              "n_genes"=rep(c(100,300,600,1000,1500,2000), each=7))


line.chart.sil.df$method<-factor(line.chart.sil.df$method, 
                                 levels = c("Mcadet_fine","Mcadet", "Seurat", "NBdrop", "Disp", "Brennecke", "Random"))


library(ggplot2)

  
a=ggplot(data = line.chart.sil.df, mapping = aes(x = n_genes, y =silmean , 
                                                 colour = method, shape=method,linetype = method))+
  geom_line(size=0.8)+geom_point()+theme_bw()+theme(legend.position="none")+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Mean silhoutte score")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(panel.grid.major = element_blank(),
        panel.background = element_blank(),axis.line = element_blank())+
  scale_x_continuous(limits=c(100,2000), breaks = c(100,300,600,1000,1500,2000))











  ## Purity

puritymean<-NULL
puritysd<-NULL
for (i in 1:6) {
  number.set<-c(100,300,600,1000,1500,2000)
  df<-get(paste0("Purity", "_pbmc_", number.set[i]  , ".df"))
  puritymean<-append(puritymean,tapply(df$Y,  factor(df$method), mean))
  puritysd<-append(puritysd,tapply(df$Y,  factor(df$method), sd))
}

line.chart.purity.df<-data.frame(puritymean, puritysd, "method"=rep(c("Mcadet_fine","Mcadet", "NBdrop",
                                                                      "Brennecke", "Seurat", "Disp", "Random"),6),
                                 "n_genes"=rep(c(100,300,600,1000,1500,2000), each=7))


line.chart.purity.df$method<-factor(line.chart.purity.df$method, 
                                    levels = c("Mcadet_fine","Mcadet","Seurat", "NBdrop", "Disp","Brennecke","Random"))



library(ggplot2)

b=ggplot(data = line.chart.purity.df, mapping = aes(x = n_genes, y =puritymean , 
                                                    colour = method, shape=method,linetype = method))+
  geom_line(size=0.8)+geom_point()+theme_bw()+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Mean purity")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(axis.line = element_blank(),panel.background = element_blank(),
        panel.grid.major = element_blank())+
  scale_x_continuous(limits=c(100,2000), breaks = c(100,300,600,1000,1500,2000))






  ## ARI

arimean<-NULL
arisd<-NULL
for (i in 1:6) {
  number.set<-c(100,300,600,1000,1500,2000)
  df<-get(paste0("ARI", "_pbmc_", number.set[i]  , ".df"))
  arimean<-append(arimean,tapply(df$Y,  factor(df$method), mean))
  arisd<-append(arisd,tapply(df$Y,  factor(df$method), sd))
}

line.chart.ari.df<-data.frame(arimean, arisd, "method"=rep(c("Mcadet_fine","Mcadet", "NBdrop",
                                                             "Brennecke", "Seurat", "Disp", "Random"),6),
                              "n_genes"=rep(c(100,300,600,1000,1500,2000), each=7))


line.chart.ari.df$method<-factor(line.chart.ari.df$method, 
                                 levels = c("Mcadet_fine","Mcadet","Seurat", "NBdrop",  "Disp","Brennecke","Random"))




library(ggplot2)

c=ggplot(data = line.chart.ari.df, mapping = aes(x = n_genes, y =arimean , 
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
  df<-get(paste0("NMI", "_pbmc_", number.set[i]  , ".df"))
  nmimean<-append(nmimean,tapply(df$Y,  factor(df$method), mean))
  nmisd<-append(nmisd,tapply(df$Y,  factor(df$method), sd))
}

line.chart.nmi.df<-data.frame(nmimean, nmisd, "method"=rep(c( "Mcadet_fine","Mcadet", "NBdrop",
                                                              "Brennecke", "Seurat", "Disp", "Random"),6),
                              "n_genes"=rep(c(100,300,600,1000,1500,2000), each=7))


line.chart.nmi.df$method<-factor(line.chart.nmi.df$method, 
                                 levels = c("Mcadet_fine","Mcadet","Seurat", "NBdrop", "Disp", "Brennecke","Random"))



library(ggplot2)

d=ggplot(data = line.chart.nmi.df, mapping = aes(x = n_genes, y =nmimean , 
                                                 colour = method, shape=method,linetype = method))+
  geom_line(size=0.8)+geom_point()+theme_bw()+theme(legend.position="none")+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Mean NMI")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(panel.background = element_blank(), axis.line = element_blank(),
        panel.grid.major = element_blank())+
  scale_x_continuous(limits=c(100,2000), breaks = c(100,300,600,1000,1500,2000))




library('patchwork')

(jac_pbmc_linechart_p+knn.only)/(a+b)/(c+d)








