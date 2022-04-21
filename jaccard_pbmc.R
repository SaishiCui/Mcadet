

library(Seurat)

for (i in 1:8) {
  for (j in 1:3) {

donor.set<-c(1,2,3,4,5,6,7,8)
time.set<-c(0,3,7)


load(file =paste0("E:/Thesis/pbmc/pbmc_fine_",  donor.set[i], "_", time.set[j],".RData"))
load(file =paste0("E:/Thesis/pbmc/pbmc_fine_clabel_", donor.set[i], "_", time.set[j],".RData"))

rawdata<-get(paste0("pbmc_fine_", donor.set[i], "_", time.set[j]))

working.data = rawdata[rowSums(rawdata)!=0,] 

pbmcdata_seurat <- CreateSeuratObject(counts =working.data)
pbmcdata_seurat <- NormalizeData(pbmcdata_seurat, normalization.method = "LogNormalize", scale.factor = 10000)


all.genes <- rownames(rawdata)
pbmcdata_seurat <- ScaleData(pbmcdata_seurat, features = all.genes)
Idents(pbmcdata_seurat)<-factor(get(paste0("pbmc_fine_clabel_", donor.set[i], "_", time.set[j])))


all.markers <- FindAllMarkers(object =pbmcdata_seurat, return.thresh = 1, logfc.threshold = 0)

all.markders.newrank<-all.markers[order(abs(all.markers$avg_log2FC),decreasing = T),]


assign(paste0("rank_df_fine_", donor.set[i], "_", time.set[j]), all.markders.newrank)


save(list = paste0("rank_df_fine_", donor.set[i], "_", time.set[j]),  
     file = paste0("E:/Thesis/pbmc/rank_df_fine_", donor.set[i], "_", time.set[j], ".RData"))


rm(list = paste0("pbmcdata_", donor.set[i], "_", time.set[j]))
gc()
print(c(i,j))
 } 

}



for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file = paste0("E:/Thesis/pbmc/rank_df_fine_", donor.set[i], "_", time.set[j], ".RData"))
    load(file =paste0("E:/Thesis/pbmc/mcadet_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/mcadet_fine_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/NB_fine_All_",     donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/Bre_fine_All_",    donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/Seurat_fine_All_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/disp_fine_All_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/random_fine_All_", donor.set[i], "_", time.set[j],".RData"))
  }
}





for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file = paste0("E:/Thesis/pbmc/rank_df_", donor.set[i], "_", time.set[j], ".RData"))
    load(file =paste0("E:/Thesis/pbmc/mcadet_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/NB_all_",     donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/Bre_all_",    donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/Seurat_all_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/disp_all_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/random_all_", donor.set[i], "_", time.set[j],".RData"))
  }
}


### get DE genes ###

for (i in 1:8) {
  for (j in 1:3) {
    
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    
    rank_df<-get(paste0("rank_df_", donor.set[i], "_", time.set[j]))
    
    num_clusters<-length(unique(rank_df$cluster))
    extract_DEgene_list<-NULL

    for (t in 1:num_clusters) {
    cluster_label<-which(rank_df$cluster==unique(rank_df$cluster)[t])
    extract_rank_df<-rank_df[cluster_label,]
    order_by_padjust<-order(extract_rank_df[,"p_val_adj"],decreasing = F)

    len<-sum(1*(extract_rank_df[order_by_padjust,"p_val_adj"]<0.05))
    sig_label<-which(extract_rank_df[order_by_padjust,"p_val_adj"]<0.05)
    if(len<400){
      extract_DEgene<-extract_rank_df[order_by_padjust[sig_label],"gene"]
      
    }else{
      extract_DEgene<-extract_rank_df[order_by_padjust,"gene"][1:400]
      
    }
    
 
    extract_DEgene_list<-append(extract_DEgene_list, extract_DEgene)
    
    }
    
    DEgene_clean_coarse=unique(extract_DEgene_list)

    
    assign(paste0("DEgene_coarse", donor.set[i], "_", time.set[j]) , DEgene_clean_coarse)
    
    save(list = paste0("DEgene_coarse", donor.set[i], "_", time.set[j]),  
         file = paste0("E:/Thesis/pbmc/DEgene_coarse", donor.set[i], "_", time.set[j], ".RData"))
    
    
    
  }
  
}




for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file = paste0("E:/Thesis/pbmc/DEgene_fine_", donor.set[i], "_", time.set[j], ".RData"))
    load(file =paste0("E:/Thesis/pbmc/mcadet_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/mcadet_fine_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/NB_fine_All_",     donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/Bre_fine_All_",    donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/Seurat_fine_All_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/disp_fine_All_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/random_fine_All_", donor.set[i], "_", time.set[j],".RData"))
  }
}





jaccard=function(a,b){
  return(c(length(intersect(a,b)), length(intersect(a,b))/(length(a)+length(b)-length(intersect(a,b)))))}


jac_pbmc<-NULL
jacmean=NULL
jacsd=NULL
for (k in 1:20) {
  num.set=seq(100,2000,100)

for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    
    one.time.run<-c(
    jaccard(get(paste0("mcadet_fine_", donor.set[i], "_", time.set[j]))$genelist[1:num.set[k]], 
              get(paste0("DEgene_fine_", donor.set[i], "_", time.set[j])))[2],
    jaccard(get(paste0("mcadet_", donor.set[i], "_", time.set[j]))$genelist[1:num.set[k]], 
            get(paste0("DEgene_fine_", donor.set[i], "_", time.set[j])))[2],
    jaccard(get(paste0("Seurat_fine_All_", donor.set[i], "_", time.set[j]))[1:num.set[k]], 
            get(paste0("DEgene_fine_", donor.set[i], "_", time.set[j])))[2],
    jaccard(get(paste0("NB_fine_All_", donor.set[i], "_", time.set[j]))[1:num.set[k]], 
            get(paste0("DEgene_fine_", donor.set[i], "_", time.set[j])))[2],
    jaccard(get(paste0("Bre_fine_All_", donor.set[i], "_", time.set[j]))[1:num.set[k]], 
            get(paste0("DEgene_fine_", donor.set[i], "_", time.set[j])))[2],
    jaccard(get(paste0("disp_fine_All_", donor.set[i], "_", time.set[j]))[1:num.set[k]], 
            get(paste0("DEgene_fine_", donor.set[i], "_", time.set[j])))[2],
    jaccard(get(paste0("random_fine_All_", donor.set[i], "_", time.set[j]))[1:num.set[k]], 
            get(paste0("DEgene_fine_", donor.set[i], "_", time.set[j])))[2])
    jac_pbmc<-append(jac_pbmc, one.time.run)
    
  }
  
 }
  df<-data.frame(jac_pbmc, method=c("Mcadet_fine", "Mcadet", "Seurat",
                                     "NBdrop", "Brennecke", "Disp", "Random"))
  colnames(df)[1]="Y"
  df$method=factor(df$method, 
                   levels = c("Mcadet_fine","Mcadet","Seurat", "NBdrop", "Disp", "Brennecke","Random"))
  
  
  jacmean<-append(jacmean,tapply(df$Y,  df$method, mean))
  jacsd<-append(jacsd,tapply(df$Y,  factor(df$method), sd))
  
}


jac_pbmc_linechart_df<-as.data.frame( cbind(jacmean, jacsd, "n_genes"=rep(seq(100,2000,100),each=7 ),
                                            "method"=rep(c("Mcadet_fine", "Mcadet", "Seurat",
                                                           "NBdrop", "Brennecke", "Disp", "Random"),20 ) )  )
jac_pbmc_linechart_df$n_genes=as.numeric(jac_pbmc_linechart_df$n_genes)
jac_pbmc_linechart_df$jacmean=as.numeric(jac_pbmc_linechart_df$jacmean)
jac_pbmc_linechart_df$jacsd=as.numeric(jac_pbmc_linechart_df$jacsd)
jac_pbmc_linechart_df$method= factor(jac_pbmc_linechart_df$method, 
                                     levels = c("Mcadet_fine","Mcadet","Seurat", "NBdrop", "Disp", "Brennecke","Random"))


### Linechart ###


jac_pbmc_linechart_p<- ggplot(data =jac_pbmc_linechart_df, mapping = aes(x = n_genes, y = jacmean, 
                                               colour = method, shape=method,linetype = method))+
  geom_line(size=0.8)+geom_point()+theme_bw()+theme(legend.position="none")+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Mean Jaccard similarity index")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(panel.background = element_blank(),axis.line = element_blank(),
        panel.grid.major = element_blank())+
  scale_x_continuous(limits=c(200,2000), breaks = seq(200,2000,200))






 ### Jaccard boxplot ###

library(ggpubr)

C<-ggboxplot(jac_pbmc_df_coarse, x="method", y="Jaccard",
          color="method", add="jitter") +labs(x="", y="Jaccard Similarity Index")+
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(jac_pbmc_df_coarse$Jaccard),
                                         linetype=2) + # ??????base mean????????????
  stat_compare_means(method = "kruskal.test",label.y =0.7)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",label.y =0.6,
                     ref.group = ".all.")+theme(legend.position="none")             # Pairwise comparison against all






C<-ggboxplot(jac_pbmc_df[-which(jac_pbmc_df$method=="Mcadet_fine"),], x="method", y="Jaccard",
             color="method", add="jitter") +labs(x="", y="Jaccard Similarity Index")+
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(jac_pbmc_df[- which(jac_pbmc_df$method=="Mcadet_fine"),]$Jaccard),
                                         linetype=2) + # ??????base mean????????????
  stat_compare_means(method = "kruskal.test",label.y =0.5)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",label.y =0.4,
                     ref.group = ".all.")+theme(legend.position="none")             # Pairwise comparison against all









   ### stable analysis



jaccard=function(a,b){
  return(c(length(intersect(a,b)), length(intersect(a,b))/(length(a)+length(b)-length(intersect(a,b)))))}

mcadet_jac_fine<-NULL
mcadet_jac<-NULL
NB_jac<-NULL
Bre_jac<-NULL
Seurat_jac<-NULL
disp_jac<-NULL
random_jac<-NULL

for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    list_mca_fine<-get(paste0("mcadet_fine_", donor.set[i], "_", time.set[j]))
    list_mca<-get(paste0("mcadet_", donor.set[i], "_", time.set[j]))
    list_NB<-get(paste0("NB_fine_All_", donor.set[i], "_", time.set[j]))
    list_Bre<-get(paste0("Bre_fine_All_", donor.set[i], "_", time.set[j]))
    list_Seurat<-get(paste0("Seurat_fine_All_", donor.set[i], "_", time.set[j]))
    list_disp<-get(paste0("disp_fine_All_", donor.set[i], "_", time.set[j]))
    list_random<-get(paste0("random_fine_All_", donor.set[i], "_", time.set[j]))
    
    
    jac_mca_fine<-jaccard(mcadet_fine_1_0$genelist, list_mca_fine$genelist)
    jac_mca<-jaccard(mcadet_1_0$genelist, list_mca$genelist)
    jac_NB<-jaccard(NB_fine_All_1_0[1:2000], list_NB[1:2000])
    jac_Bre<-jaccard(Bre_fine_All_1_0[1:2000], list_Bre[1:2000])
    jac_Seurat<-jaccard(Seurat_fine_All_1_0[1:2000], list_Seurat[1:2000])    
    jac_disp<-jaccard(disp_fine_All_1_0[1:2000], list_disp[1:2000])
    jac_random<-jaccard(random_fine_All_1_0[1:2000], list_random[1:2000])   
    
    
    mcadet_jac_fine<-append(mcadet_jac_fine,jac_mca_fine)
    mcadet_jac<-append(mcadet_jac,jac_mca)
    NB_jac<-append(NB_jac,jac_NB)
    Bre_jac<-append(Bre_jac,jac_Bre)
    Seurat_jac<-append(Seurat_jac,jac_Seurat)
    disp_jac<-append(disp_jac,jac_disp)
    random_jac<-append(random_jac,jac_random)
    
  }
}

mcadet_jac_fine<-mcadet_jac_fine[-c(1,2)]
mcadet_jac_fine<-mcadet_jac_fine[seq(2,46,2)]
mcadet_jac<-mcadet_jac[-c(1,2)]
mcadet_jac<-mcadet_jac[seq(2,46,2)]
NB_jac<-NB_jac[-c(1,2)]
NB_jac<-NB_jac[seq(2,46,2)]
Bre_jac<-Bre_jac[-c(1,2)]
Bre_jac<-Bre_jac[seq(2,46,2)]
Seurat_jac<-Seurat_jac[-c(1,2)]
Seurat_jac<-Seurat_jac[seq(2,46,2)]
disp_jac<-disp_jac[-c(1,2)]
disp_jac<-disp_jac[seq(2,46,2)]
random_jac<-random_jac[-c(1,2)]
random_jac<-random_jac[seq(2,46,2)]



jac.df<-data.frame("Jaccard"<-c(mcadet_jac_fine, mcadet_jac,NB_jac,Bre_jac,Seurat_jac, disp_jac, random_jac), 
                   "Method"<-c(rep("Mcadet_fine",23),rep("Mcadet",23), rep("NBdrop",23), rep("Brennecke", 23),
                               rep("Seurat", 23), rep("disp", 23), rep("random", 23)),
                   "list"<- rep(c(1:23),7))

colnames(jac.df)<-c("Jaccard", "Method", "list")
jac.df$Method=factor(jac.df$Method, levels = c("disp", "NBdrop", "Mcadet_fine","Mcadet", "Seurat", "Brennecke", "random"))


library(ggplot2)
f<-ggplot(data = jac.df, mapping = aes(x = list, y = Jaccard, linetype = Method, colour = Method, shape=Method,
                                       fill = Method))+theme_bw()+
  geom_line(size=0.8) + geom_point() +labs(x="", y="Jaccard similarity index")
+theme(legend.position="none")





### Percentage of CD4 and CD8 driver genes ###

CD4<-rank_df_1_0[rank_df_1_0$cluster=="CD4 T",]
CD4_makers<-CD4$gene[1:1000]
CD8<-rank_df_1_0[rank_df_1_0$cluster=="CD8 T",]
CD8_makers<-CD8$gene[1:1000]


mean(CD4_makers%in%mcadet_1_0_fine_30$genelist)
mean(CD8_makers%in%mcadet_1_0$genelist)

mean(CD4_makers%in%NB_1_0)
mean(CD8_makers%in%NB_1_0)

mean(CD4_makers%in%Seurat_1_0)
mean(CD8_makers%in%Seurat_1_0)

