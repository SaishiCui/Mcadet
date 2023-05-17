
rm(list = ls())
gc()

library("umap")
library("uwot")
library("cluster")
library("factoextra")
library("funtimes")
library("mclust")
library("aricode")
library("Seurat")
library("irlba")


donor.set<-c(1,2,3,4,5,6,7,8)
time.set<-c(0,3,7)
p = 7
k= 3

load(file =paste0("C:/Users/css22/Desktop/Thesis/pbmc/data/pbmcdata_", 
                  donor.set[p], "_", time.set[k], ".RData"))
load(file =paste0("C:/Users/css22/Desktop/Thesis/pbmc/data_fine/pbmc_fine_", 
                  donor.set[p], "_", time.set[k], ".RData"))

load(file =paste0("C:/Users/css22/Desktop/Thesis/pbmc/data/pbmccell1_", 
                  donor.set[p], "_", time.set[k], ".RData"))

load(file =paste0("C:/Users/css22/Desktop/Thesis/pbmc/data_fine/pbmc_fine_clabel_", 
                  donor.set[p], "_", time.set[k], ".RData"))

load(file =paste0("C:/Users/css22/Desktop/Thesis/pbmc/DE_genes/DEgene_", 
                  donor.set[p], "_", time.set[k], ".RData"))

load(file =paste0("C:/Users/css22/Desktop/Thesis/pbmc/DE_genes/DEgene_fine_", 
                  donor.set[p], "_", time.set[k], ".RData"))








## coarse  resolution


cell.info<- get(paste0("pbmccell1_", donor.set[p], "_", time.set[k]))
cell.info <- factor(cell.info, levels = c("B", "CD4 T", "CD8 T", "DC", "Mono", "NK", "other", "other T"))
workdata=get(paste0("pbmcdata_", donor.set[p], "_", time.set[k]))
cbPalette <- c( "turquoise3", "blue3", "springgreen4", "#E69F00", "#ef1828", "darkorchid2", "#3f60aa", "#9ec417")


### All genes  ##

new_data<-workdata

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
all_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "All genes")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

all_umap









### True  ##

true_pbmc_vector = get(paste0("DEgene_", donor.set[p], "_", time.set[k]))

genelist<- true_pbmc_vector
gene.name<- rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-workdata[label,]



set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
true_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  scale_color_manual(values = cbPalette)+
  theme_classic()+labs(x="", y="", title = "True HVGs")+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))


true_umap








### random default 2000 ##

genelist<-get(paste0("random_PBMC_", donor.set[p], "_", time.set[k]))

gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
random_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Random selection")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

random_umap




### mcadet default ###

genelist<- get(paste0("Mcadet_pbmc_",donor.set[p], "_", time.set[k] ))$gene
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}


new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)


m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
mcadet_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Mcadet")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

mcadet_umap




### disp select 2000 ###

genelist<-  get(paste0("SeuratDisp_PBMC_",donor.set[p], "_", time.set[k] ))
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
disp_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  scale_color_manual(values = cbPalette)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2", title = "Seurat Disp")+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))


disp_umap




### NBdrop  ###

genelist<- get(paste0("NBDrop_PBMC_",donor.set[p], "_", time.set[k] ))
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
NB_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "NBDrop")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

NB_umap






### Brennecke  ###

genelist<- get(paste0("BreFS_PBMC_",donor.set[p], "_", time.set[k] ))
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
Bre_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Brennecke")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

Bre_umap





### M3Drop ###


genelist<- get(paste0("M3Drop_PBMC_",donor.set[p], "_", time.set[k] ))
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}


new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
M3Drop_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "M3Drop")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

M3Drop_umap




### Seurat Vst ###

genelist<- get(paste0("SeuratVst_PBMC_",donor.set[p], "_", time.set[k] ))
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
Vst_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Seurat Vst")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))


Vst_umap



### Seurat Mvp ###

genelist<- get(paste0("SeuratMvp_PBMC_",donor.set[p], "_", time.set[k] ))

gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)


m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
Mvp_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Seurat Mvp")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

Mvp_umap




### Scry ###


genelist<- get(paste0("scryFS_PBMC_",donor.set[p], "_", time.set[k] ))
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
scry_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha= 1)+
  theme_classic()+labs(x="", y="", title = "Scry")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

scry_umap





legend = ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=2, alpha=1)+
  theme_classic()+labs(x="", y="", title = "")+
  scale_color_manual(values = cbPalette, name = "Cell type")+
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title = element_text(face ="bold", size = 12),
    legend.position = "right",
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.background = element_rect(color = "black", 
                                     linewidth = 1.5,
                                     fill = "lightblue"))+
  guides(colour = guide_legend(override.aes = list(size=3)))



library(patchwork)

legend


random_umap + all_umap + true_umap + 
  mcadet_umap + NB_umap + M3Drop_umap + scry_umap + 
  Bre_umap+ disp_umap + Vst_umap + Mvp_umap 




#### fine resolution 


cell.info<- get(paste0("pbmc_fine_clabel_", donor.set[p], "_", time.set[k]))
cell.info = factor(cell.info, levels = c("CD4 CTL", "CD4 Naive", "CD4 TCM", "CD4 TEM", "CD8 Naive", "CD8 TCM",
                                                       "CD8 TEM", "dnT", "gdT", "MAIT", "Treg"))
workdata=get(paste0("pbmc_fine_", donor.set[p], "_", time.set[k]))


cbPalette <- c("springgreen4", "#3f60aa","#ef1828", "blue3", "turquoise3", "darkorchid2", "#E69F00", "#9ec417","darkgrey","#9b3a74", "black")
                                          
### All genes  ##

new_data<-workdata

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
all_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "All genes")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

all_umap









### True  ##

true_pbmc_vector = get(paste0("DEgene_fine_", donor.set[p], "_", time.set[k]))

genelist<- true_pbmc_vector
gene.name<- rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-workdata[label,]



set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
true_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  scale_color_manual(values = cbPalette)+
  theme_classic()+labs(x="", y="", title = "True HVGs")+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))


true_umap








### random default 2000 ##

genelist<-get(paste0("random_PBMC_fine_", donor.set[p], "_", time.set[k]))

gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
random_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Random selection")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

random_umap




### mcadet default ###

genelist<- get(paste0("Mcadet_pbmc_fine_",donor.set[p], "_", time.set[k] ))$gene
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}


new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)


m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
mcadet_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Mcadet")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

mcadet_umap




### disp select 2000 ###

genelist<-  get(paste0("SeuratDisp_PBMC_fine_",donor.set[p], "_", time.set[k] ))
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
disp_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  scale_color_manual(values = cbPalette)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2", title = "Seurat Disp")+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))


disp_umap




### NBdrop  ###

genelist<- get(paste0("NBDrop_PBMC_fine_",donor.set[p], "_", time.set[k] ))
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
NB_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "NBDrop")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

NB_umap






### Brennecke  ###

genelist<- get(paste0("BreFS_PBMC_fine_",donor.set[p], "_", time.set[k] ))
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
Bre_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Brennecke")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

Bre_umap





### M3Drop ###


genelist<- get(paste0("M3Drop_PBMC_fine_",donor.set[p], "_", time.set[k] ))
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}


new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
M3Drop_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "M3Drop")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

M3Drop_umap




### Seurat Vst ###

genelist<- get(paste0("SeuratVst_PBMC_fine_",donor.set[p], "_", time.set[k] ))
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
Vst_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Seurat Vst")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))


Vst_umap



### Seurat Mvp ###

genelist<- get(paste0("SeuratMvp_PBMC_fine_",donor.set[p], "_", time.set[k] ))

gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)


m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
Mvp_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Seurat Mvp")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

Mvp_umap




### Scry ###


genelist<- get(paste0("scryFS_PBMC_fine_",donor.set[p], "_", time.set[k] ))
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-workdata[label,]

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
scry_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha= 1)+
  theme_classic()+labs(x="", y="", title = "Scry")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

scry_umap





legend = ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=2, alpha=1)+
  theme_classic()+labs(x="", y="", title = "")+
  scale_color_manual(values = cbPalette, name = "Cell type")+
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title = element_text(face ="bold", size = 12),
    legend.position = "right",
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.background = element_rect(color = "black", 
                                     linewidth = 1.5,
                                     fill = "lightblue"))+
  guides(colour = guide_legend(override.aes = list(size=3)))



library(patchwork)

legend


random_umap + all_umap + true_umap + 
  mcadet_umap + NB_umap + M3Drop_umap + scry_umap + 
  Bre_umap+ disp_umap + Vst_umap + Mvp_umap 



