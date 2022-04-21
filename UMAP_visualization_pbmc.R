



library("umap")
library("uwot")
library("cluster")
library("factoextra")
library("funtimes")
library("mclust")
library("aricode")

load(file =paste0("E:/Thesis/pbmc/pbmcdata_",  2, "_", 7,".RData"))

load(file =paste0("E:/Thesis/pbmc/pbmccell1_",  2, "_", 7,".RData"))
load(file =paste0("E:/Thesis/pbmc/pbmccell2_",  2, "_", 7,".RData"))
load(file =paste0("E:/Thesis/pbmc/mcadet_", 2, "_", 7,".RData"))
load(file =paste0("E:/Thesis/pbmc/NB_",     2, "_", 7,".RData"))
load(file =paste0("E:/Thesis/pbmc/Bre_",    2, "_", 7,".RData"))
load(file =paste0("E:/Thesis/pbmc/Seurat_", 2, "_", 7,".RData"))
load(file =paste0("E:/Thesis/pbmc/disp_", 2, "_", 7,".RData"))
load(file =paste0("E:/Thesis/pbmc/random_", 2, "_", 7,".RData"))

library(RColorBrewer)
display.brewer.all() 
allcolour<-c(brewer.pal(9, "Set1")[1], "blue", "green",
             brewer.pal(9, "Set1")[4:8])


cell.info<-get(paste0("pbmccell1_", 2, "_", 7))
workdata=get(paste0("pbmcdata_", 2, "_", 7))


  ### randomly select 2000 ##
genelist<-random_2_7
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
random_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position = "top",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(name="", values = allcolour)+
  annotate("text",x=-3, y=20, label="Random",fontface='bold', size=5)






  ### mcadet select 2000 ###

genelist<-mcadet_2_7$genelist
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
mcadet_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(name="",values = allcolour)+
  annotate("text",x=-3, y=18, label="Mcadet",fontface='bold', size=5)





   ### Seurat select 2000 ###

genelist<-Seurat_2_7
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
seurat_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+
  theme(legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(values = allcolour)+
  annotate("text",x=2.5, y=15, label="Seurat",fontface='bold', size=5)



### NBdrop select 2000 ###

genelist<-NB_2_7
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
NB_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(values = allcolour)+
  annotate("text",x=0, y=20, label="NBdrop",fontface='bold', size=5)




### Brennecke select 2000 ###

genelist<-Bre_2_7
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
Bre_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(values = allcolour)+
  annotate("text",x=10, y=25, label="Brennecke",fontface='bold', size=5)



### Disp select 2000 ###

genelist<-disp_2_7
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
Disp_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+
  theme(legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(values = allcolour)+
  annotate("text",x=3, y=20, label="Disp", fontface='bold', size=5)



 ### all genes ###



new_data<-t(workdata)

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
nonfs_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+
  theme(legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(values = allcolour)+
  annotate("text",x=2, y=15, label="All genes", fontface='bold', size=5)


  ### umap plot###

library("patchwork")
(nonfs_umap+random_umap+Bre_umap)/(Disp_umap+NB_umap+mcadet_umap)







## only keep T cells ###


T_cell_label<-c(which(cell.info=="CD4 T"),which(cell.info=="CD8 T"), which(cell.info=="other T"))
cell.info<-cell.info[T_cell_label]

  ## All genes

new_data<-t(workdata[,T_cell_label])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
nonfs_umap_T<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+
  theme(legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(values = allcolour)+
  annotate("text",x=1.25, y=6, label="All genes", fontface='bold', size=5)




### randomly select 2000 ##
genelist<-random_2_7
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,T_cell_label])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
random_umap_T<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position = "top",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(name="",values = allcolour)+
  annotate("text",x=-0.8, y=4, label="Random",fontface='bold', size=5)






### Brennecke select 2000 ##

genelist<-Bre_2_7
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,T_cell_label])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
Bre_umap_T<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(values = allcolour)+
  annotate("text",x=-2, y=35, label="Brennecke",fontface='bold', size=5)




### Disp select 2000 ##

genelist<-disp_2_7
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,T_cell_label])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
disp_umap_T<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(values = allcolour)+
  annotate("text",x=1, y=10, label="Disp",fontface='bold', size=5)




### NBdrop select 2000 ##

genelist<-NB_2_7
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,T_cell_label])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
NB_umap_T<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(values = allcolour)+
  annotate("text",x=1, y=6, label="NBdrop",fontface='bold', size=5)




### Mcadet select 2000 ##

genelist<-mcadet_2_7$genelist
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,T_cell_label])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
mcadet_umap_T<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(values = allcolour)+
  annotate("text",x=-2, y=3.5, label="Mcadet",fontface='bold', size=5)



library("patchwork")
(nonfs_umap_T+random_umap_T+Bre_umap_T)/(disp_umap_T+NB_umap_T+mcadet_umap_T)








####### T cells finer resolution #########

load(file =paste0("E:/Thesis/pbmc/pbmc_fine_",  2, "_", 7,".RData"))
load(file =paste0("E:/Thesis/pbmc/pbmc_fine_clabel_",  2, "_", 7,".RData"))

cell.info<-get(paste0("pbmc_fine_clabel_", 2, "_", 7))
workdata=get(paste0("pbmc_fine_", 2, "_", 7))

types<-unique(pbmccell2_2_7)[c(3,4,5,7,10,11,21)]
display.brewer.all() 
allcolour<-c(brewer.pal(9, "Set1")[1], "blue",
             brewer.pal(9, "Set1")[4], brewer.pal(9, "Set1")[5],  "black", brewer.pal(9, "Set1")[7], "green")


T_cell_label<-c(which(cell.info==types[1]),
                which(cell.info==types[2]), 
                which(cell.info==types[3]),
                which(cell.info==types[4]),
                which(cell.info==types[5]), 
                which(cell.info==types[6]),
                which(cell.info==types[7]))

cell.info<-cell.info[T_cell_label]


## All genes

new_data<-t(workdata[,T_cell_label])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"
umap.df$group=factor(umap.df$group, levels = c("CD4 Naive", "CD8 Naive", "CD4 TCM", "CD8 TCM", "CD4 TEM",
                                               "CD8 TEM", "CD4 CTL"))


library(ggplot2)
nonfs_umap_T<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+
  theme(legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(values = allcolour)+
  annotate("text",x=1.25, y=4, label="All genes", fontface='bold', size=5)




### randomly select 2000 ##

genelist<-random_fine_All_2_7[1:2000]
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,T_cell_label])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"
umap.df$group=factor(umap.df$group, levels = c("CD4 Naive", "CD8 Naive", "CD4 TCM", "CD8 TCM", "CD4 TEM",
                                               "CD8 TEM", "CD4 CTL"))



library(ggplot2)
random_umap_T<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position = "top",legend.text=element_text(size=12), legend.title =element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(size=2)))+scale_color_manual(name="", values = allcolour)+
  annotate("text",x=0.6, y=4, label="Random",fontface='bold', size=5)





### Brennecke select 2000 ##

genelist<-Bre_fine_All_2_7[1:2000]
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,T_cell_label])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"
umap.df$group=factor(umap.df$group, levels = c("CD4 Naive", "CD8 Naive", "CD4 TCM", "CD8 TCM", "CD4 TEM",
                                               "CD8 TEM", "CD4 CTL"))




library(ggplot2)
Bre_umap_T<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(values = allcolour)+
  annotate("text",x=-2, y=5, label="Brennecke",fontface='bold', size=5)




### Disp select 2000 ##

genelist<-disp_fine_All_2_7[1:2000]
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,T_cell_label])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"
umap.df$group=factor(umap.df$group, levels = c("CD4 Naive", "CD8 Naive", "CD4 TCM", "CD8 TCM", "CD4 TEM",
                                               "CD8 TEM", "CD4 CTL"))




library(ggplot2)
disp_umap_T<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(values = allcolour)+
  annotate("text",x=1, y=5, label="Disp",fontface='bold', size=5)




### NBdrop select 2000 ##

genelist<-NB_fine_All_2_7[1:2000]
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,T_cell_label])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"
umap.df$group=factor(umap.df$group, levels = c("CD4 Naive", "CD8 Naive", "CD4 TCM", "CD8 TCM", "CD4 TEM",
                                               "CD8 TEM", "CD4 CTL"))



library(ggplot2)
NB_umap_T<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position = "none",legend.text=element_text(size=12), legend.title =element_text(size=12))+
  guides(colour = guide_legend(override.aes = list(size=2)))+scale_color_manual(values = allcolour)+
  annotate("text",x=1.5, y=6, label="NBdrop",fontface='bold', size=5)




### Mcadet select 2000 ##

genelist<-mcadet_fine_2_7$genelist
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,T_cell_label])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"
umap.df$group=factor(umap.df$group, levels = c("CD4 Naive", "CD8 Naive", "CD4 TCM", "CD8 TCM", "CD4 TEM",
                                               "CD8 TEM", "CD4 CTL"))




library(ggplot2)
mcadet_umap_T<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(values = allcolour)+
  annotate("text",x=1.5, y=4, label="Mcadet fine",fontface='bold', size=5)






### Seurat select 2000 ##

genelist<-Seurat_fine_All_2_7[1:2000]
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,T_cell_label])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2, pca=50,  pca_method = "irlba")

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"
umap.df$group=factor(umap.df$group, levels = c("CD4 Naive", "CD8 Naive", "CD4 TCM", "CD8 TCM", "CD4 TEM",
                                               "CD8 TEM", "CD4 CTL"))




library(ggplot2)
Seurat_umap_T<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+scale_color_manual(values = allcolour)+
  annotate("text",x=1.5, y=6, label="Seurat",fontface='bold', size=5)


















library("patchwork")
library(gridExtra)
(nonfs_umap_T+random_umap_T+Bre_umap_T)/(disp_umap_T+NB_umap_T+mcadet_umap_T)
(nonfs_umap_T+random_umap_T+Bre_umap_T)/(disp_umap_T+Seurat_umap_T+mcadet_umap_T)
grid.arrange(nonfs_umap_T,random_umap_T,Bre_umap_T,disp_umap_T,Seurat_umap_T, NB_umap_T, mcadet_umap_T)

(Seurat_umap_T)







############################################### Violin plot pbmc ##########################################################


exp.rank_high<- rownames(pbmcdata_1_0)[order(apply(pbmcdata_1_0,1,mean), decreasing = T)][1:5000]
exp.rank_low<- rownames(pbmcdata_1_0)[order(apply(pbmcdata_1_0,1,mean), decreasing = F)][1:8000]

all.select<-intersect(intersect(intersect(NB_8_0, Bre_8_0),Seurat_8_0),mcadet_8_0$genelist)

all.select<-intersect(intersect(NB_8_0,Seurat_8_0),mcadet_8_0$genelist)



intersect(all.select,exp.rank_low)
intersect(all.select,exp.rank_high)


union1<-union(NB_1_0, Seurat_1_0)

mcadet.only<-setdiff(mcadet_1_0$genelist,union1)

intersect(mcadet.only,exp.rank_high)
intersect(mcadet.only,exp.rank_low)


notselmca<-setdiff(rownames(pbmcdata_8_0),mcadet_8_0$genelist)
notselnb<-setdiff(rownames(pbmcdata_8_0), NB_8_0)
notselseu<-setdiff(rownames(pbmcdata_8_0), Seurat_8_0)


notsel<-intersect(intersect(notselmca,notselnb), notselseu)

intersect(notsel,exp.rank_high)[1:50]
intersect(notsel,exp.rank_low)[1:50]


selother<-intersect(NB_8_0, Seurat_8_0 )

intersect(intersect(selother, notselmca),exp.rank_low )
intersect(intersect(selother, notselmca),exp.rank_high )






### all select high expression
mean(as.numeric(pbmcdata_8_0["GNLY",]))


df.violin.GNLY<-data.frame(exp=t(pbmcdata_8_0["GNLY",]), group= factor(pbmccell1_8_0))
colnames(df.violin.GNLY)[1]="exp"
p.GNLY<-ggplot(df.violin.GNLY, aes(x = group, y = log(exp+1) ))+geom_violin(aes(fill = group), trim = F)+
  labs(title="GNLY", y="Log gene expression", x="")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 12, face = "bold"))+
  annotate("text",x=4.5, y=2, label="mean=1.09", size=4)






### all select low expression
mean(as.numeric(pbmcdata_8_0["CLNK",]))


df.violin.CLNK<-data.frame(exp=t(pbmcdata_8_0["CLNK",]), group= factor(pbmccell1_8_0))
colnames(df.violin.CLNK)[1]="exp"
p.CLNK<-ggplot(df.violin.CLNK, aes(x = group, y = log(exp+1) ))+geom_violin(aes(fill = group), trim = F)+
  labs(title="CLNK", y="Log gene expression", x="")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 12, face = "bold"))+  annotate("text",x=4.5, y=1, label="mean=0.005", size=4)





### only selected by Mcadet high

intersect(mcadet.only,exp.rank_high)
intersect(mcadet.only,exp.rank_low)


mean(as.numeric(pbmcdata_1_0["ILK",]))

df.violin.C15orf39<-data.frame(exp=t(pbmcdata_1_0["C15orf39",]), group=pbmccell1_1_0)
colnames(df.violin.C15orf39)[1]="exp"
p.C15orf39<-ggplot(df.violin.C15orf39, aes(x = group, y = log(exp +1) ))+geom_violin(aes(fill = group), trim = FALSE)+
  labs(title="C15orf39", y="Log Gene Expression", x="")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 12, face = "bold"))

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df,t(pbmcdata_1_0["C15orf39",]) )
colnames(umap.df)[3]<-"Expression"

umap_C15orf39<-ggplot(umap.df,aes(V1,V2, color=Expression))+
  geom_point(size=1.3, alpha=1.2)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+theme(plot.title = element_text(hjust = 0.5))+
  scale_color_gradient(low = "cyan",high = "red")+labs(title = "C15orf39")





df.violin.THEM4<-data.frame(exp=t(pbmcdata_1_0["THEM4",]), group=pbmccell1_1_0)
colnames(df.violin.THEM4)[1]="exp"
p.THEM4<-ggplot(df.violin.THEM4, aes(x = group, y = log(exp +1) ))+geom_violin(aes(fill = group), trim = FALSE)+
  labs(title="THEM4", y="Log Gene Expression", x="")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 12, face = "bold"))


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df,t(pbmcdata_1_0["THEM4",]) )
colnames(umap.df)[3]<-"Expression"

umap_THEM4<-ggplot(umap.df,aes(V1,V2, color=Expression))+
  geom_point(size=1.3, alpha=1.2)+theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x="UMAP 1", y="UMAP 2")+ scale_color_gradient(low = "cyan",high = "red")+labs(title = "THEM4")


library(patchwork)
(p.C15orf39+umap_C15orf39)/(p.THEM4+umap_THEM4)












### only selected by Mcadet low

mean(as.numeric(pbmcdata_8_0["TNS1",]))

df.violin.TNS1<-data.frame(exp=t(pbmcdata_8_0["TNS1",]), group=pbmccell1_8_0)
colnames(df.violin.TNS1)[1]="exp"
p.TNS1<-ggplot(df.violin.TNS1, aes(x = group, y = log(exp +1) ))+geom_violin(aes(fill = group), trim = FALSE)+
  labs(title="TNS1", y="Log Gene Expression", x="")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 12, face = "bold"))+annotate("text",x=4.5, y=1, label="mean=0.005", size=4)




### not selected by Mcadet but selected by others high
mean(as.numeric(pbmcdata_8_0["IL1B",]))

df.violin.IL1B<-data.frame(exp=t(pbmcdata_8_0["IL1B",]), group=pbmccell1_8_0)
colnames(df.violin.IL1B)[1]="exp"
p.IL1B<-ggplot(df.violin.IL1B, aes(x = group, y = log(exp +1) ))+geom_violin(aes(fill = group), trim = FALSE)+
  labs(title="IL1B", y="Log Gene Expression", x="")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 12, face = "bold"))+annotate("text",x=4.5, y=2.5, label="mean=1.39", size=4)







### not selected by Mcadet but selected by others low
mean(as.numeric(pbmcdata_8_0["CCL8",]))

df.violin.CCL8<-data.frame(exp=t(pbmcdata_8_0["CCL8",]), group=pbmccell1_8_0)
colnames(df.violin.CCL8)[1]="exp"
p.CCL8<-ggplot(df.violin.CCL8, aes(x = group, y = log(exp +1) ))+geom_violin(aes(fill = group), trim = FALSE)+
  labs(title="CCL8", y="Log Gene Expression", x="")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 12, face = "bold"))+annotate("text",x=4.5, y=1.5, label="mean=0.001", size=4)









### not selected by all high
mean(as.numeric(pbmcdata_8_0["ISG15",]))

df.violin.ISG15<-data.frame(exp=t(pbmcdata_8_0["ISG15",]), group=pbmccell1_8_0)
colnames(df.violin.ISG15)[1]="exp"
p.ISG15<-ggplot(df.violin.ISG15, aes(x = group, y = log(exp +1) ))+geom_violin(aes(fill = group), trim = FALSE)+
  labs(title="ISG15", y="Log Gene Expression", x="")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 12, face = "bold"))+annotate("text",x=4.5, y=2.0, label="mean=0.33", size=4)





### not selected by all low
mean(as.numeric(pbmcdata_8_0["AL645608.8",]))

df.violin.AL645608.8<-data.frame(exp=t(pbmcdata_8_0["AL645608.8",]), group=pbmccell1_8_0)
colnames(df.violin.AL645608.8)[1]="exp"
p.AL645608.8<-ggplot(df.violin.AL645608.8, aes(x = group, y = log(exp +1) ))+geom_violin(aes(fill = group), trim = FALSE)+
  labs(title="AL645608.8", y="Log Gene Expression", x="")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 12, face = "bold"))+annotate("text",x=4.5, y=1.2, label="mean=0.0004", size=4)










library(patchwork)
(p.GNLY+p.CLNK)/(p.ISG15+p.AL645608.8)/(p.SUN2+p.TNS1)/(p.IL1B+p.CCL8)
















