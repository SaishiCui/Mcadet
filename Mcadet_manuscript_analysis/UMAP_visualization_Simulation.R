
rm(list = ls())
gc()

true_SparSim_vector = c()
for(i in 1:2000){
  true_SparSim_vector = append(true_SparSim_vector, paste0("Gene", i))
  
}

library("umap")
library("uwot")
library("cluster")
library("factoextra")
library("funtimes")
library("mclust")
library("aricode")
library("Seurat")

j=12

load(file =paste0("C:/Users/css22/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_fine_", j, ".RData"))

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
cell.info = factor(cell.info, levels = c("Type 1",
                                         "Type 2",
                                         "Type 3",
                                         "Type 4",
                                         "Type 5",
                                         "Type 6",
                                         "Type 7",
                                         "Type 8",
                                         "Type 9",
                                         "Type 10"))

workdata=get(paste0("SparSim_fine_", j))

cbPalette <- c("springgreen4", "#3f60aa","#9b3a74", "blue3", "turquoise3", "darkorchid2", "#E69F00", "#9ec417","darkgrey", "#ef1828")


### All genes  ##

new_data<-workdata

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)


m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]

pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
               retx = TRUE, center = TRUE, scale. = TRUE)$rotation)



umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
all_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "All genes")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

all_umap








### True 2000 ##

genelist<- true_SparSim_vector
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

pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)



umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
true_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "True HVGs")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))


true_umap







### random default 2000 ##

genelist<-get(paste0("random_SparSim_fine_", j))
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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
random_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "Random selection")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

random_umap




### mcadet default ###

genelist<- get(paste0("Mcadet_SparSim_fine_", j))$gene
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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
mcadet_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "Mcadet")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))
mcadet_umap



### disp select 2000 ###

genelist<- get(paste0("SeuratDisp_SparSim_fine_", j))
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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
disp_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2", title = "Seurat Disp")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title = element_text(face ="bold", size = 12),
    legend.position = "none",
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.background = element_rect(color = "black", 
                                     linewidth = 1.5,
                                     fill = "lightblue"))+
  guides(colour = guide_legend(override.aes = list(size=3)))

disp_umap




### NBdrop  ###

genelist<- get(paste0("NBDrop_SparSim_fine_", j))
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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)



umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
NB_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "NBDrop")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

NB_umap






### Brennecke  ###

genelist<- get(paste0("BreFS_SparSim_fine_", j))
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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
Bre_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "Brennecke")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

Bre_umap





### M3Drop ###


genelist<- get(paste0("M3Drop_SparSim_fine_", j))
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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
M3Drop_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "M3Drop")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

M3Drop_umap




### Seurat Vst ###

genelist<- get(paste0("SeuratVst_SparSim_fine_", j))
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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
Vst_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "Seurat Vst")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

Vst_umap



### Seurat Mvp ###

genelist<- get(paste0("SeuratMvp_SparSim_fine_", j))

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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
Mvp_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "Seurat Mvp")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

Mvp_umap




### Scry ###


genelist<- get(paste0("scryFS_SparSim_fine_", j))
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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
scry_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "Scry")+
  scale_color_manual(values = cbPalette)+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))
scry_umap



legend = ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1, alpha=1)+
  scale_color_manual(values = cbPalette, name = "Cell type")+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2", title = "")+
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




###########  SparSim coarse 

j=12

load(file =paste0("C:/Users/css22/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_", j, ".RData"))

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

cell.info = factor(cell.info, levels = c("Type 1",
                                         "Type 2",
                                         "Type 3",
                                         "Type 4",
                                         "Type 5",
                                         "Type 6",
                                         "Type 7",
                                         "Type 8",
                                         "Type 9",
                                         "Type 10"))
workdata=get(paste0("SparSim_", j))

cbPalette <- c("springgreen4", "#3f60aa","#9b3a74", "blue3", "turquoise3", "darkorchid2", "#E69F00", "#9ec417","darkgrey", "#ef1828")



### All genes  ##

new_data<-workdata

set.seed(1234)
seurat.obj<-CreateAssayObject(counts = new_data)
seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)


m= as.matrix(seurat.obj@data)[which(apply(as.matrix(seurat.obj@data), 1, sum)  != 0),]

pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)



umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
all_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "All genes")+
  scale_color_manual(values = cbPalette, name = "Cell type")+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

all_umap








### True 2000 ##

genelist<- true_SparSim_vector
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

pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)



umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
true_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "True HVGs")+
  scale_color_manual(values = cbPalette, name = "Cell type")+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))


true_umap







### random default 2000 ##

genelist<-get(paste0("random_SparSim_", j))
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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
random_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "Random selection")+
  scale_color_manual(values = cbPalette, name = "Cell type")+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))
random_umap




### mcadet default ###

genelist<- get(paste0("Mcadet_SparSim_", j))$gene
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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
mcadet_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "Mcadet")+
  scale_color_manual(values = cbPalette, name = "Cell type")+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))
mcadet_umap



### disp select 2000 ###

genelist<- get(paste0("SeuratDisp_SparSim_", j))
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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
disp_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2", title = "Seurat Disp")+
  scale_color_manual(values = cbPalette, name = "Cell type")+
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title = element_text(face ="bold", size = 12),
    legend.position = "none",
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.background = element_rect(color = "black", 
                                     linewidth = 1.5,
                                     fill = "lightblue"))+
  guides(colour = guide_legend(override.aes = list(size=3)))
disp_umap




### NBdrop  ###

genelist<- get(paste0("NBDrop_SparSim_", j))
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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)



umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
NB_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "NBDrop")+
  scale_color_manual(values = cbPalette, name = "Cell type")+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

NB_umap






### Brennecke  ###

genelist<- get(paste0("BreFS_SparSim_", j))
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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
Bre_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "Brennecke")+
  scale_color_manual(values = cbPalette, name = "Cell type")+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))
Bre_umap





### M3Drop ###


genelist<- get(paste0("M3Drop_SparSim_", j))
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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
M3Drop_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "M3Drop")+
  scale_color_manual(values = cbPalette, name = "Cell type")+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

M3Drop_umap




### Seurat Vst ###

genelist<- get(paste0("SeuratVst_SparSim_", j))
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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
Vst_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "Seurat Vst")+
  scale_color_manual(values = cbPalette, name = "Cell type")+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

Vst_umap



### Seurat Mvp ###

genelist<- get(paste0("SeuratMvp_SparSim_", j))

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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
Mvp_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "Seurat Mvp")+
  scale_color_manual(values = cbPalette, name = "Cell type")+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

Mvp_umap




### Scry ###


genelist<- get(paste0("scryFS_SparSim_", j))
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
pc_df <- as.data.frame(prcomp_irlba(m, n = 15,
                                    retx = TRUE, center = TRUE, scale. = TRUE)$rotation)
umap<-uwot::umap(pc_df , n_components=2)


umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
scry_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="", y="", title = "Scry")+
  scale_color_manual(values = cbPalette, name = "Cell type")+
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "none",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))

scry_umap



legend = ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=0.1, alpha=0.8)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2", title = "")+
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

