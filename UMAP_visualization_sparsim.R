

library("umap")
library("uwot")
library("cluster")
library("factoextra")
library("funtimes")
library("mclust")
library("aricode")

load(file =paste0("E:/Thesis/SPARSim/SparSim_", 16, ".RData"))
load(file =paste0("E:/Thesis/SPARSim/mcadet_Sparsim_", 16, ".RData"))
load(file =paste0("E:/Thesis/SPARSim/NB_Sparsim_",     16,".RData"))
load(file =paste0("E:/Thesis/SPARSim/Bre_Sparsim_",    16,".RData")) 
load(file =paste0("E:/Thesis/SPARSim/Seurat_Sparsim_", 16,".RData"))
load(file =paste0("E:/Thesis/SPARSim/disp_Sparsim_",   16,".RData"))
load(file =paste0("E:/Thesis/SPARSim/random_Sparsim_", 16,".RData"))

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

workdata=get(paste0("SparSim_",16))



### randomly select 2000 ##
genelist<-random_Sparsim_16
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
random_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+
  theme(legend.position = "top",legend.text=element_text(size=15), legend.title =element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  scale_color_discrete(breaks = c('Type 1','Type 2','Type 3','Type 4',
                                  "Type 5", "Type 6", "Type 7", "Type 8",
                                  "Type 9", "Type 10"))+annotate("text",x=0, y=2.5, label="Random",fontface='bold', size=5)








### mcadet select 2000 ###
genelist<-mcadet_Sparsim_16$genelist
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
mcadet_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position="none")+annotate("text",x=-2, y=5, label="Mcadet",fontface='bold', size=5)



### disp select 2000 ###
genelist<-disp_Sparsim_16
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
disp_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+theme_classic()+
  theme(legend.position="none")+annotate("text",x=-1, y=8, label="Disp",fontface='bold', size=5)




### NBdrop select 2000 ###

genelist<-NB_Sparsim_16
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
NB_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position="none")+annotate("text",x=0, y=8, label="NBdrop",fontface='bold', size=5)







### Brennecke select 2000 ###

genelist<-Bre_Sparsim_16
gene.name<-rownames(workdata)  

label<-NULL
for (i in 1:length(genelist)) {
  label[i]=which(gene.name==genelist[i]) 
}

new_data<-t(workdata[label,])

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
Bre_umap<-ggplot(umap.df,aes(V1,V2, color=group))+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+ 
  theme(legend.position="none")+
  annotate("text",x=-1, y=4, label="Brennecke",fontface='bold', size=5)







### all genes ###



new_data<-t(workdata)

set.seed(1234)

umap<-uwot::umap(log(new_data+1), n_components=2)

umap.df<-as.data.frame(umap)
umap.df<-cbind(umap.df, cell.info)
colnames(umap.df)[3]<-"group"

library(ggplot2)
nonfs_umap<-ggplot(umap.df,aes(V1,V2, color=group),legend="none")+
  geom_point(size=1.3, alpha=1)+
  theme_classic()+labs(x="UMAP 1", y="UMAP 2")+
  theme(legend.position="none")+
 annotate("text",x=0.5, y=6, label="All genes",fontface='bold', size=5)



### umap ###

library("patchwork")
(nonfs_umap+random_umap+Bre_umap)/(disp_umap+NB_umap+mcadet_umap)






