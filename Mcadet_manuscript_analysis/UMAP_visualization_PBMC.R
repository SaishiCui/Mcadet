
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
library("vegan")
library("patchwork")

donor.set<-c(1,2,3,4,5,6,7,8)
time.set<-c(0,3,7)
p = 7
k= 3

load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/data/pbmcdata_", 
                  donor.set[p], "_", time.set[k], ".RData"))
load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/data_fine/pbmc_fine_", 
                  donor.set[p], "_", time.set[k], ".RData"))

load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/data/pbmccell1_", 
                  donor.set[p], "_", time.set[k], ".RData"))

load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/data_fine/pbmc_fine_clabel_", 
                  donor.set[p], "_", time.set[k], ".RData"))

load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/DE_genes/DEgene_", 
                  donor.set[p], "_", time.set[k], ".RData"))

load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/M3Drop_PBMC_fine_", 
                  donor.set[p], "_", time.set[k], ".RData"))


load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/NBDrop_PBMC_fine_", 
                  donor.set[p], "_", time.set[k], ".RData"))


load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/BreFS_PBMC_fine_", 
                  donor.set[p], "_", time.set[k], ".RData"))

load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/scryFS_PBMC_fine_", 
                  donor.set[p], "_", time.set[k], ".RData"))

load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/SeuratDisp_PBMC_fine_", 
                  donor.set[p], "_", time.set[k], ".RData"))

load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/SeuratMvp_PBMC_fine_", 
                  donor.set[p], "_", time.set[k], ".RData"))


load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/SeuratVst_PBMC_fine_", 
                  donor.set[p], "_", time.set[k], ".RData"))


load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/DE_genes/DEgene_fine_", 
                  donor.set[p], "_", time.set[k], ".RData"))


load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/random_PBMC_fine_",
                  donor.set[p], "_", time.set[k], ".RData"))


load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/random_PBMC_fine_",
                  donor.set[p], "_", time.set[k], ".RData"))

load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/pbmc_results/Mcadet_pbmc_fine_",
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
















### random default 2000 ##






# Assuming you have multiple UMAP embeddings (umap1 and umap2) to align
# Here umap1.df and umap2.df are data frames containing UMAP embeddings





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
cluster_results_true <- kmeans(pc_df, 11)

umap<-uwot::umap(pc_df , n_components=2, seed = 123456)

umap.df_random<-as.data.frame(umap)
umap.df_random<-cbind(umap.df_random, cell.info, factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                                   "6", "7", "8", "9", "10", "11")))
colnames(umap.df_random)[3]<-"group"
colnames(umap.df_random)[4]<-"group_cluster"






# Perform Procrustes alignment
proc <- procrustes(umap.df[, 1:2], umap.df_random[, 1:2])

# Extract the aligned coordinates
aligned_umap_random <- proc$Yrot

# Convert the aligned coordinates back to a data frame
aligned_umap.df_random <- as.data.frame(aligned_umap_random)
aligned_umap.df_random <- cbind(aligned_umap.df_random, cell.info, 
                                factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                      "6", "7", "8", "9", "10", "11")))
colnames(aligned_umap.df_random)[3] <- "group"
colnames(aligned_umap.df_random)[4] <- "group_cluster"




library(ggplot2)
random_umap<-ggplot(aligned_umap.df_random, aes(V1,V2, color=group))+
  geom_point(size=1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Random selection")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")



random_umap_cluster <-ggplot(aligned_umap.df_random, aes(V1,V2, color=group_cluster))+
  geom_point(size=1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")







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
cluster_results_true <- kmeans(pc_df, 11)
umap<-uwot::umap(pc_df , n_components=2, seed = 123456)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info, factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                                   "6", "7", "8", "9", "10", "11")))
colnames(umap.df)[3]<-"group"
colnames(umap.df)[4]<-"group_cluster"



true_umap <- ggplot(umap.df, aes(V1, V2, color = group)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_manual(values = cbPalette) +
  theme_classic() +
  labs(x = "", y = "", title = "Semi-true HVGs") +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "left",
    legend.text = element_text(size = 15, face = "bold"),
    legend.title = element_blank(), 
    legend.background = element_rect(fill = "lightblue", colour = "black", linewidth = 0.5, linetype = "solid"),
    axis.title.x = element_text(hjust = 0.15, size = 12, face = "bold"),  
    axis.title.y = element_text(hjust = 0.2, angle = 90, vjust = 0, size = 12, face = "bold"))+
  guides(colour = guide_legend(override.aes = list(size=3))) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")










true_umap_cluster <-ggplot(umap.df,aes(V1,V2, color=group_cluster))+
  geom_point(size = 1, alpha = 1) +
  scale_color_manual(values = cbPalette) +
  theme_classic() +
  labs(x = "UMAP1", y = "UMAP2", title = "Semi-true HVGs") +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "left",
    legend.text = element_text(size = 15, face = "bold"),
    legend.title = element_blank(), 
    legend.background = element_rect(fill = "lightblue", colour = "black", linewidth = 0.5, linetype = "solid"),
    axis.title.x = element_text(hjust = 0.15, size = 12, face = "bold"),  
    axis.title.y = element_text(hjust = 0.2, angle = 90, vjust = 0, size = 12, face = "bold"))+
  guides(colour = guide_legend(override.aes = list(size = 3)))+
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")








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
cluster_results_true <- kmeans(pc_df, 11)


umap<-uwot::umap(pc_df , n_components=2, seed = 123456)


umap.df_mcadet<-as.data.frame(umap)
umap.df_mcadet<-cbind(umap.df_mcadet, cell.info, factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                                                 "6", "7", "8", "9", "10", "11")))
colnames(umap.df_mcadet)[3]<-"group"
colnames(umap.df_mcadet)[4]<-"group_cluster"




# Perform Procrustes alignment
proc <- procrustes(umap.df[, 1:2], umap.df_mcadet[, 1:2])

# Extract the aligned coordinates
aligned_umap_mcadet <- proc$Yrot

# Convert the aligned coordinates back to a data frame
aligned_umap.df_mcadet <- as.data.frame(aligned_umap_mcadet)
aligned_umap.df_mcadet <- cbind(aligned_umap.df_mcadet, cell.info,
                                factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                                "6", "7", "8", "9", "10", "11")))
colnames(aligned_umap.df_mcadet)[3] <- "group"
colnames(aligned_umap.df_mcadet)[4]<-"group_cluster"




mcadet_umap<-ggplot(aligned_umap.df_mcadet, aes(V1,V2, color=group))+
  geom_point(size=1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Mcadet")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")




mcadet_umap_cluster<-ggplot(aligned_umap.df_mcadet, aes(V1,V2, color=group_cluster))+
  geom_point(size=1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Mcadet")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "F", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")





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

cluster_results_true <- kmeans(pc_df, 11)

umap<-uwot::umap(pc_df , n_components=2)



umap.df_SeuratDisp<-as.data.frame(umap)
umap.df_SeuratDisp<-cbind(umap.df_SeuratDisp, cell.info, factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                                                 "6", "7", "8", "9", "10", "11")))
colnames(umap.df_SeuratDisp)[3]<-"group"
colnames(umap.df_SeuratDisp)[4]<-"group_cluster"




# Perform Procrustes alignment
proc <- procrustes(umap.df[, 1:2], umap.df_SeuratDisp[, 1:2])

# Extract the aligned coordinates
aligned_umap_SeuratDisp <- proc$Yrot

# Convert the aligned coordinates back to a data frame
aligned_umap.df_SeuratDisp <- as.data.frame(aligned_umap_SeuratDisp)
aligned_umap.df_SeuratDisp <- cbind(aligned_umap.df_SeuratDisp, cell.info,
                                factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                                "6", "7", "8", "9", "10", "11")))
colnames(aligned_umap.df_SeuratDisp)[3] <- "group"
colnames(aligned_umap.df_SeuratDisp)[4]<-"group_cluster"







SeuratDisp_umap<-ggplot(aligned_umap.df_SeuratDisp, aes(V1,V2, color=group))+
  geom_point(size=1, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2", title = "Seurat Disp")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15),
    axis.title.x = element_text(hjust = 0.15, size = 12, face = "bold"),  
    axis.title.y = element_text(hjust = 0.2, angle = 90, vjust = 0, size = 12, face = "bold"))+ 
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "G", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")







SeuratDisp_umap_cluster<-ggplot(aligned_umap.df_SeuratDisp, aes(V1,V2, color=group_cluster))+
  geom_point(size=1, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2", title = "Seurat Disp")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15),
    ,
    axis.title.x = element_text(hjust = 0.15, size = 12, face = "bold"),  
    axis.title.y = element_text(hjust = 0.2, angle = 90, vjust = 0, size = 12, face = "bold")  
   )+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "G", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")








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


cluster_results_true <- kmeans(pc_df, 11)

umap<-uwot::umap(pc_df , n_components=2)



umap.df_NBDrop<-as.data.frame(umap)
umap.df_NBDrop<-cbind(umap.df_NBDrop, cell.info, factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                                                         "6", "7", "8", "9", "10", "11")))
colnames(umap.df_NBDrop)[3]<-"group"
colnames(umap.df_NBDrop)[4]<-"group_cluster"




# Perform Procrustes alignment
proc <- procrustes(umap.df[, 1:2], umap.df_NBDrop[, 1:2])

# Extract the aligned coordinates
aligned_umap_NBDrop <- proc$Yrot

# Convert the aligned coordinates back to a data frame
aligned_umap.df_NBDrop <- as.data.frame(aligned_umap_NBDrop)
aligned_umap.df_NBDrop <- cbind(aligned_umap.df_NBDrop, cell.info,
                                    factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                                    "6", "7", "8", "9", "10", "11")))
colnames(aligned_umap.df_NBDrop)[3] <- "group"
colnames(aligned_umap.df_NBDrop)[4]<-"group_cluster"







NBDrop_umap<-ggplot(aligned_umap.df_NBDrop, aes(V1,V2, color=group))+
  geom_point(size=1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "NBDrop")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "left",
    legend.text = element_text(size = 15, face = "bold"),
    legend.title = element_blank(), 
    legend.background = element_rect(fill = "lightblue", colour = "black", linewidth = 0.5, linetype = "solid"))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")



NBDrop_umap_cluster<-ggplot(aligned_umap.df_NBDrop, aes(V1,V2, color=group_cluster))+
  geom_point(size=1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "NBDrop")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "left",
    legend.text = element_text(size = 15, face = "bold"),
    legend.title = element_blank(), 
    legend.background = element_rect(fill = "lightblue", colour = "black", linewidth = 0.5, linetype = "solid"))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")














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

cluster_results_true <- kmeans(pc_df, 11)

umap<-uwot::umap(pc_df , n_components=2)



umap.df_BreFS<-as.data.frame(umap)
umap.df_BreFS<-cbind(umap.df_BreFS, cell.info, factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                                                 "6", "7", "8", "9", "10", "11")))
colnames(umap.df_BreFS)[3]<-"group"
colnames(umap.df_BreFS)[4]<-"group_cluster"




# Perform Procrustes alignment
proc <- procrustes(umap.df[, 1:2], umap.df_BreFS[, 1:2])

# Extract the aligned coordinates
aligned_umap_BreFS <- proc$Yrot

# Convert the aligned coordinates back to a data frame
aligned_umap.df_BreFS <- as.data.frame(aligned_umap_BreFS)
aligned_umap.df_BreFS <- cbind(aligned_umap.df_BreFS, cell.info,
                                factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                                "6", "7", "8", "9", "10", "11")))
colnames(aligned_umap.df_BreFS)[3] <- "group"
colnames(aligned_umap.df_BreFS)[4]<-"group_cluster"







BreFS_umap<-ggplot(aligned_umap.df_BreFS, aes(V1,V2, color=group))+
  geom_point(size=1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Brennecke")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")



BreFS_umap_cluster<-ggplot(aligned_umap.df_BreFS, aes(V1,V2, color=group_cluster))+
  geom_point(size=1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Brennecke")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")






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


cluster_results_true <- kmeans(pc_df, 11)

umap<-uwot::umap(pc_df , n_components=2)



umap.df_M3Drop<-as.data.frame(umap)
umap.df_M3Drop<-cbind(umap.df_M3Drop, cell.info, factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                                               "6", "7", "8", "9", "10", "11")))
colnames(umap.df_M3Drop)[3]<-"group"
colnames(umap.df_M3Drop)[4]<-"group_cluster"




# Perform Procrustes alignment
proc <- procrustes(umap.df[, 1:2], umap.df_M3Drop[, 1:2])

# Extract the aligned coordinates
aligned_umap_M3Drop <- proc$Yrot

# Convert the aligned coordinates back to a data frame
aligned_umap.df_M3Drop <- as.data.frame(aligned_umap_M3Drop)
aligned_umap.df_M3Drop <- cbind(aligned_umap.df_M3Drop, cell.info,
                               factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                               "6", "7", "8", "9", "10", "11")))
colnames(aligned_umap.df_M3Drop)[3] <- "group"
colnames(aligned_umap.df_M3Drop)[4]<-"group_cluster"







M3Drop_umap<-ggplot(aligned_umap.df_M3Drop, aes(V1,V2, color=group))+
  geom_point(size=1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "M3Drop")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "F", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")



M3Drop_umap_cluster<-ggplot(aligned_umap.df_M3Drop, aes(V1,V2, color=group_cluster))+
  geom_point(size=1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "M3Drop")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "F", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")






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



cluster_results_true <- kmeans(pc_df, 11)

umap<-uwot::umap(pc_df , n_components=2)



umap.df_SeuratVst <-as.data.frame(umap)
umap.df_SeuratVst <-cbind(umap.df_SeuratVst, cell.info, factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                                                 "6", "7", "8", "9", "10", "11")))
colnames(umap.df_SeuratVst)[3]<-"group"
colnames(umap.df_SeuratVst)[4]<-"group_cluster"




# Perform Procrustes alignment
proc <- procrustes(umap.df[, 1:2], umap.df_SeuratVst[, 1:2])

# Extract the aligned coordinates
aligned_umap_SeuratVst <- proc$Yrot

# Convert the aligned coordinates back to a data frame
aligned_umap.df_SeuratVst <- as.data.frame(aligned_umap_SeuratVst)
aligned_umap.df_SeuratVst <- cbind(aligned_umap.df_SeuratVst, cell.info,
                                factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                                "6", "7", "8", "9", "10", "11")))
colnames(aligned_umap.df_SeuratVst)[3] <- "group"
colnames(aligned_umap.df_SeuratVst)[4]<-"group_cluster"







SeuratVst_umap<-ggplot(aligned_umap.df_SeuratVst, aes(V1,V2, color=group))+
  geom_point(size=1, alpha=1)+ 
  theme_classic()+labs(x="", y="", title = "Seurat Vst")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "H", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")



SeuratVst_umap_cluster<-ggplot(aligned_umap.df_SeuratVst, aes(V1,V2, color=group_cluster))+
  geom_point(size=1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Seurat Vst")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "H", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")






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


cluster_results_true <- kmeans(pc_df, 11)

umap<-uwot::umap(pc_df , n_components=2)



umap.df_SeuratMvp <-as.data.frame(umap)
umap.df_SeuratMvp <-cbind(umap.df_SeuratMvp, cell.info, factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                                                        "6", "7", "8", "9", "10", "11")))
colnames(umap.df_SeuratMvp)[3]<-"group"
colnames(umap.df_SeuratMvp)[4]<-"group_cluster"




# Perform Procrustes alignment
proc <- procrustes(umap.df[, 1:2], umap.df_SeuratMvp[, 1:2])

# Extract the aligned coordinates
aligned_umap_SeuratMvp <- proc$Yrot

# Convert the aligned coordinates back to a data frame
aligned_umap.df_SeuratMvp <- as.data.frame(aligned_umap_SeuratMvp)
aligned_umap.df_SeuratMvp <- cbind(aligned_umap.df_SeuratMvp, cell.info,
                                   factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                                   "6", "7", "8", "9", "10", "11")))
colnames(aligned_umap.df_SeuratMvp)[3] <- "group"
colnames(aligned_umap.df_SeuratMvp)[4]<-"group_cluster"







SeuratMvp_umap <-ggplot(aligned_umap.df_SeuratMvp, aes(V1,V2, color=group))+
  geom_point(size=1, alpha=1)+ 
  theme_classic()+labs(x="", y="", title = "Seurat Mvp")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "I", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")



SeuratMvp_umap_cluster<-ggplot(aligned_umap.df_SeuratMvp, aes(V1,V2, color=group_cluster))+
  geom_point(size=1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Seurat Mvp")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "I", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")




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


cluster_results_true <- kmeans(pc_df, 11)

umap<-uwot::umap(pc_df , n_components=2)



umap.df_Scry <-as.data.frame(umap)
umap.df_Scry <-cbind(umap.df_Scry, cell.info, factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                                                        "6", "7", "8", "9", "10", "11")))
colnames(umap.df_Scry)[3]<-"group"
colnames(umap.df_Scry)[4]<-"group_cluster"




# Perform Procrustes alignment
proc <- procrustes(umap.df[, 1:2], umap.df_Scry[, 1:2])

# Extract the aligned coordinates
aligned_umap_Scry <- proc$Yrot

# Convert the aligned coordinates back to a data frame
aligned_umap.df_Scry <- as.data.frame(aligned_umap_Scry)
aligned_umap.df_Scry <- cbind(aligned_umap.df_Scry, cell.info,
                                   factor(cluster_results_true$cluster, levels = c("1", "2", "3", "4", "5",
                                                                                   "6", "7", "8", "9", "10", "11")))
colnames(aligned_umap.df_Scry)[3] <- "group"
colnames(aligned_umap.df_Scry)[4]<-"group_cluster"







Scry_umap <-ggplot(aligned_umap.df_Scry, aes(V1,V2, color=group))+
  geom_point(size=1, alpha=1)+ 
  theme_classic()+labs(x="", y="", title = "Scry")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")



Scry_umap_cluster<-ggplot(aligned_umap.df_Scry, aes(V1,V2, color=group_cluster))+
  geom_point(size=1, alpha=1)+
  theme_classic()+labs(x="", y="", title = "Scry")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.1, vjust = 1.1, size = 9, fontface = "bold")





(true_umap + random_umap + mcadet_umap) / (true_umap_cluster + random_umap_cluster + mcadet_umap_cluster)
