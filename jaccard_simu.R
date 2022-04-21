
rm(list = ls())
gc()





for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)

    load(file =paste0("E:/Thesis/pbmc/mcadet_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/NB_",     donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/Bre_",    donor.set[i], "_", time.set[j],".RData")) ### pbmc real data
    load(file =paste0("E:/Thesis/pbmc/Seurat_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/disp_", donor.set[i], "_", time.set[j],".RData"))
    load(file =paste0("E:/Thesis/pbmc/random_", donor.set[i], "_", time.set[j],".RData"))
  }
}





jaccard=function(a,b){
  return(c(length(intersect(a,b)), length(intersect(a,b))/(length(a)+length(b)-length(intersect(a,b)))))}




#### Sparsim simulation #########




load(file =paste0("E:/Thesis/SPARSim/SparSim_", i, ".RData"))

for (i in 1:16) {
  

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


load(file =paste0("E:/Thesis/SPARSim/SparSim_", 1, ".RData"))
for (i in 1:16) {
  list<-get(paste0("mcadet_Sparsim_try2_", i))
  print(jaccard(rownames(SparSim_1)[1:2000], list$genelist))
  }

for (i in 1:16){ 
  list<-get(paste0("mcadet_Sparsim_", i))
  print(jaccard(rownames(SparSim_1)[1:2000], list$genelist))
}

for (i in 1:16){ 
  list<-get(paste0("NB_Sparsim_", i))
  print(jaccard(rownames(SparSim_1)[1:2000], list))
}

for (i in 1:16){ 
  list<-get(paste0("Seurat_Sparsim_", i))
  print(jaccard(rownames(SparSim_1)[1:2000], list))
}








final.jaccard<-NULL
for (i in 1:16) {

    jaccard=function(a,b){
    return(c(length(intersect(a,b)), length(intersect(a,b))/(length(a)+length(b)-length(intersect(a,b)))))}
  
    
    dt.list<-rownames(SparSim_1)[1:2000]
    
    

    
    genelist.mca<-get(paste0("mcadet_Sparsim_try2_", i))$genelist
    genelist.mca.up<-get(paste0("mcadet_Sparsim_", i))$genelist
    genelist.mca.mid<-get(paste0("mcadet_Sparsim_mid", i))$genelist
    genelist.mca.down<-get(paste0("mcadet_Sparsim_down", i))$genelist
    
    genelist.NB<-get(paste0("NB_Sparsim_", i))
    genelist.Bre<-get(paste0("Bre_Sparsim_", i))  
    genelist.Seurat<-get(paste0("Seurat_Sparsim_", i)) 
    
    genelist.disp<-get(paste0("disp_Sparsim_", i))
    genelist.random<-get(paste0("random_Sparsim_", i))
    
    one.time.run<-c(  jaccard(dt.list, genelist.mca.down)[2],
                      jaccard(dt.list, genelist.mca.mid)[2],
                      jaccard(dt.list, genelist.mca.up)[2],
                      jaccard(dt.list, genelist.mca)[2],
                      jaccard(dt.list, genelist.NB)[2],
                      jaccard(dt.list, genelist.Bre)[2],
                      jaccard(dt.list, genelist.Seurat)[2],
                      jaccard(dt.list, genelist.disp)[2],
                      jaccard(dt.list, genelist.random)[2] )
    
    
    final.jaccard<-append(final.jaccard, one.time.run)
  }




jaccard_df_Sparsim<-data.frame(
  "Jaccard"=final.jaccard,
  "method"=
    as.factor(rep(c("Mcadet.down","Mcadet.mid","Mcadet.ori", "Mcadet","NBdrop","Brennecke", "Seurat", "Disp","Random"),16))
)


jaccard_df_Sparsim$method<-
  factor(jaccard_df_Sparsim$method, 
         levels = c("Mcadet.down","Mcadet.mid","Mcadet.ori", "Mcadet", "Seurat", "Brennecke", "NBdrop", "Disp","Random"))

save(jaccard_df_Sparsim, file = "E:/Thesis/SPARSim/jaccard_df_Sparsim.RData")



load(file = "E:/Thesis/SPARSim/jaccard_df_Sparsim.RData")


library(ggpubr)


remove.list.jac=c(
              which(jaccard_df_Sparsim$method=="Mcadet.mid"),
              which(jaccard_df_Sparsim$method=="Mcadet.down"))

 ### boxplot
C<-ggboxplot(jaccard_df_Sparsim[-remove.list.jac,], x="method", y="Jaccard",
          color="method", add="jitter") +labs(x="", y="Jaccard similarity index")+
  rotate_x_text(angle = 45)+  geom_hline(yintercept = mean(jaccard_df_Sparsim[-remove.list.jac,]$Jaccard),
                                         linetype=2) + # ??????base mean????????????
  stat_compare_means(method = "anova",label.y = 0.7)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",label.y = 0.6,
                     ref.group = ".all.")+theme(legend.position="none")             # Pairwise comparison against all


### barplot

Jaccard<-0.5*c(jaccard_df_Sparsim$Jaccard[10:18]+jaccard_df_Sparsim$Jaccard[1:9],
jaccard_df_Sparsim$Jaccard[19:27]+jaccard_df_Sparsim$Jaccard[28:36],
jaccard_df_Sparsim$Jaccard[37:45]+jaccard_df_Sparsim$Jaccard[46:54],
jaccard_df_Sparsim$Jaccard[55:63]+jaccard_df_Sparsim$Jaccard[64:72],
jaccard_df_Sparsim$Jaccard[73:81]+jaccard_df_Sparsim$Jaccard[82:90],
jaccard_df_Sparsim$Jaccard[91:99]+jaccard_df_Sparsim$Jaccard[100:108],
jaccard_df_Sparsim$Jaccard[109:117]+jaccard_df_Sparsim$Jaccard[118:126],
jaccard_df_Sparsim$Jaccard[127:135]+jaccard_df_Sparsim$Jaccard[136:144])


method<-rep(c("Mcadet.down","Mcadet.mid","Mcadet.up", "Mcadet", "Seurat", "Brennecke", "NBdrop", "Disp","Random"),8)
level<-c(rep(8,9), rep(7,9), rep(6, 9), rep(5, 9), rep(4, 9), rep(3, 9), rep(2, 9), rep(1, 9) )


jaccard_df_Sparsim_bar<-data.frame("Jaccard"=as.numeric(Jaccard),
                                   "method"=factor(method, 
                          levels = c("Mcadet","Mcadet.up","Mcadet.mid", "Mcadet.down", "NBdrop","Seurat", "Brennecke",  "Disp","Random")) ,
                                   "level"=level)



library(ggplot2)
library(RColorBrewer)
nb.cols <- 9
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

remove.list.jac=c(which(jaccard_df_Sparsim$method=="Mcadet.up"),
                  which(jaccard_df_Sparsim$method=="Mcadet.mid"),
                  which(jaccard_df_Sparsim$method=="Mcadet.down"))

jaccard_df_Sparsim_bar_half1<- jaccard_df_Sparsim_bar[1:36, ]
jaccard_df_Sparsim_bar_half1_clean<-jaccard_df_Sparsim_bar_half1[-remove.list.jac,]
jaccard_df_Sparsim_bar_half2<- jaccard_df_Sparsim_bar[37:72, ]
jaccard_df_Sparsim_bar_half2_clean<-jaccard_df_Sparsim_bar_half2[-remove.list.jac,]

level1_4<-ggplot(data = jaccard_df_Sparsim_bar_half1_clean, 
       mapping = aes(x = factor(level), 
                     y = Jaccard, fill =method )) + 
  geom_bar(stat = 'identity', position ="dodge", width = 0.5)+
  labs(x="", y="Jaccard similarity index")+ylim(0,0.5)+scale_fill_manual(values = mycolors) 


level5_8<-ggplot(data = jaccard_df_Sparsim_bar_half2_clean, 
                 mapping = aes(x = factor(level), 
                               y = Jaccard, fill =method )) + 
  geom_bar(stat = 'identity', position ="dodge", width = 0.5)+
  labs(x="DE level", y="Jaccard Similarity Index")+ylim(0,0.5)+scale_fill_manual(values = mycolors)+
  theme(legend.position="none")


library(patchwork)
level1_4/level5_8










############################################### Violin plot ##########################################################

dt.list<-rownames(gene.info_0.03_0.4)[order(apply(gene.info_0.03_0.4[,5:14],1,sd),decreasing = T)][1:2000]
exp.rank_high<- rownames(simdata_0.03_0.4)[order(apply(simdata_0.03_0.4,1,mean), decreasing = T)][1:2000]
exp.rank_low<- rownames(simdata_0.03_0.4)[order(apply(simdata_0.03_0.4,1,mean), decreasing = F)][1:5000]

all.select<-intersect(intersect(intersect(intersect(NB_0.03_0.4, Bre_0.03_0.4),Seurat_0.03_0.4),mcadet_0.03_0.4$genelist),
                      dt.list)


intersect(all.select,exp.rank_low)
intersect(all.select,exp.rank_high)


union1<-union(NB_0.03_0.4,Bre_0.03_0.4)
union2<-union(union1,Seurat_0.03_0.4)
set3<-setdiff(mcadet_0.03_0.4$genelist,union2)
truth<-intersect(set3,dt.list)


intersect(truth,exp.rank_high)
intersect(truth,exp.rank_low)


notselmca<-setdiff(rownames(simdata_0.03_0.4),mcadet_0.03_0.4$genelist)
notselnb<-setdiff(rownames(simdata_0.03_0.4),NB_0.03_0.4)
notselbre<-setdiff(rownames(simdata_0.03_0.4),Bre_0.03_0.4)
notselseu<-setdiff(rownames(simdata_0.03_0.4),Seurat_0.03_0.4)
notseldt<-setdiff(rownames(simdata_0.03_0.4),dt.list)

notsel<-intersect(notselmca,intersect(notselnb, intersect(notselbre, intersect(notselseu, notseldt))))

intersect(notsel,exp.rank_high)
intersect(notsel,exp.rank_low)[1:20]


selother<-intersect(NB_0.03_0.4, intersect(Bre_0.03_0.4, intersect(Seurat_0.03_0.4, dt.list)) )

intersect(intersect(selother, notselmca),exp.rank_low )
intersect(intersect(selother, notselmca),exp.rank_high )






### all select high expression
mean(as.numeric(simdata_0.03_0.4[10439,]))
gene.info_0.03_0.4[10439,]

df.violin.10439<-data.frame(exp=t(simdata_0.03_0.4[10439,]), group=cell.info_0.03_0.4$Group)
colnames(df.violin.10439)[1]="exp"
p.10439<-ggplot(df.violin.10439, aes(x = group, y = exp))+geom_violin(aes(fill = group), trim = FALSE)+
  annotate("text", x=c(0.6,c(1:10),5), y=c(rep(-10,11),50), label=c("df:",2.64,1,1,1,1,1,1,1,1,1,"mean gene expression: 8.66"), size=4)+
  labs(title="Gene10439", y="Gene Expression", x="")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 12, face = "bold"))






### all select low expression
mean(as.numeric(simdata_0.03_0.4[5298,]))
gene.info_0.03_0.4[5298,]

df.violin.5298<-data.frame(exp=t(simdata_0.03_0.4[5298,]), group=cell.info_0.03_0.4$Group)
colnames(df.violin.5298)[1]="exp"
p.5298<-ggplot(df.violin.5298, aes(x = group, y = exp))+geom_violin(aes(fill = group), trim = FALSE)+
  annotate("text", x=c(0.6,c(1:10),5), y=c(rep(-1.5,11),8.75), label=c("df:",1,1,1,1,1,1,1,1.78,1.64,1,"mean gene expression: 0.13"), size=4)+
  labs(title="Gene5298", y="Gene Expression", x="")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 12, face = "bold"))




### only selected by Mcadet high
mean(as.numeric(simdata_0.03_0.4[3047,]))
gene.info_0.03_0.4[3047,]

df.violin.3047<-data.frame(exp=t(simdata_0.03_0.4[3047,]), group=cell.info_0.03_0.4$Group)
colnames(df.violin.3047)[1]="exp"
p.3047<-ggplot(df.violin.3047, aes(x = group, y = exp))+geom_violin(aes(fill = group), trim = FALSE)+
  annotate("text", x=c(0.6,c(1:10),5), y=c(rep(-10,11),45), label=c("df:",0.78,1,1,1.58,1,1,1,1,1,1,"mean gene expression: 8.31"), size=4)+
  labs(title="Gene3047", y="Gene Expression", x="")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 12, face = "bold"))


### only selected by Mcadet low
mean(as.numeric(simdata_0.03_0.4[10173,]))
gene.info_0.03_0.4[10173,]

df.violin.10173<-data.frame(exp=t(simdata_0.03_0.4[10173,]), group=cell.info_0.03_0.4$Group)
colnames(df.violin.10173)[1]="exp"
p.10173<-ggplot(df.violin.10173, aes(x = group, y = exp))+geom_violin(aes(fill = group), trim = FALSE)+
  annotate("text", x=c(0.6,c(1:10),5), y=c(rep(-1,11),6), label=c("df:",1,1,1,3.06,1,1,1,1,1,1,"mean gene expression: 0.10"), size=4)+
  labs(title="Gene10173", y="Gene Expression", x="")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 12, face = "bold"))







### not selected by Mcadet but selected by others high
mean(as.numeric(simdata_0.03_0.4[97,]))
gene.info_0.03_0.4[97,]

df.violin.97<-data.frame(exp=t(simdata_0.03_0.4[97,]), group=cell.info_0.03_0.4$Group)
colnames(df.violin.97)[1]="exp"
p.97<-ggplot(df.violin.97, aes(x = group, y = exp))+geom_violin(aes(fill = group), trim = FALSE)+
  annotate("text", x=c(0.6,c(1:10),5), y=c(rep(-10,11),50), label=c("df:",1,1,0.5,1,1,1,1,1,1,1,"mean gene expression: 7.53"), size=4)+
  labs(title="Gene97", y="Gene Expression", x="")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 12, face = "bold"))


### not selected by Mcadet but selected by others low
mean(as.numeric(simdata_0.03_0.4[1552,]))
gene.info_0.03_0.4[1552,]


df.violin.1552<-data.frame(exp=t(simdata_0.03_0.4[1552,]), group=cell.info_0.03_0.4$Group)
colnames(df.violin.1552)[1]="exp"
p.1552<-ggplot(df.violin.1552, aes(x = group, y = exp))+geom_violin(aes(fill = group), trim = FALSE)+
  annotate("text", x=c(0.6,c(1:10),5), y=c(rep(-2,11),15), label=c("df:",1,1,1,1,1,1,1,1,0.51,1,"mean gene expression: 0.30"), size=4)+
  labs(title="Gene1552", y="Gene Expression", x="")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 12, face = "bold"))








### not selected by all high
mean(as.numeric(simdata_0.03_0.4[7,]))
gene.info_0.03_0.4[7,]

df.violin.7<-data.frame(exp=t(simdata_0.03_0.4[7,]), group=cell.info_0.03_0.4$Group)
colnames(df.violin.7)[1]="exp"
p.7<-ggplot(df.violin.7, aes(x = group, y = exp))+geom_violin(aes(fill = group), trim = FALSE)+
  annotate("text", x=c(0.6,c(1:10),5), y=c(rep(-10,11),100), label=c("df:",1,1,1,1,1,1,1,1,1,1,"mean gene expression: 34.6"), size=4)+
  labs(title="Gene7", y="Gene Expression", x="")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 12, face = "bold"))


### not selected by all low
mean(as.numeric(simdata_0.03_0.4[80,]))
gene.info_0.03_0.4[80,]

df.violin.80<-data.frame(exp=t(simdata_0.03_0.4[80,]), group=cell.info_0.03_0.4$Group)
colnames(df.violin.80)[1]="exp"
p.80<-ggplot(df.violin.80, aes(x = group, y = exp))+geom_violin(aes(fill = group), trim = FALSE)+
  annotate("text", x=c(0.6,c(1:10),5), y=c(rep(-1,11),8), label=c("df:",1,1,1,1,1,1,1,1,1,1,"mean gene expression: 0.25"), size=4)+
  labs(title="Gene48", y="Gene Expression", x="")+theme(plot.title = element_text(hjust = 0.5))+theme(legend.position="none")+
  theme(axis.text.x = element_text(size = 12, face = "bold"))




library(patchwork)
(p.10439+p.5298)/(p.7+p.80)/(p.3047+p.10173)/(p.97+p.1552)







#### Line chart of jaccard ####

df_1<-NULL
df_2<-NULL
for (k in 1:20) {
  pt.matrix<-matrix(NA,nrow = 16, ncol = 9)
  for (j in 1:16) {
    method.set<-c("mcadet_down","mcadet_mid","mcadet_up",
                  "mcadet","NB","Bre","Seurat","disp","random")
    number.list<-seq(100,2000,100)
    n.feature.1=2*number.list[k]
    n.feature.2=number.list[k]
    
    data=get(paste0("SparSim_",j))
    working.data = data[rowSums(data)!=0,]          ##  keep genes with positive expression (delete empty genes)##
    clean=T
    if(clean){
      mean.expression<-apply(working.data,1,mean)
      working.data=working.data[mean.expression>=0.005,]
    }
    gene.name<-rownames(working.data)               ##  extract gene's name in a list ## 
    
    gene.range<-colSums(get(paste0("mcadet_Sparsim_", j))$gene.range.matrix )       ## take the sum of rank range for each gene for all runs (equivalent to take the average) 
    gene.min<-colSums(get(paste0("mcadet_Sparsim_", j))$gene.min.matrix )           ## take the sum of minimum rank for each gene for all runs (equivalent to take the average) 
    
    
    order1<-order(gene.range,decreasing = T)[1:n.feature.1]   ## first round of filtering genes
    gene.name.order1<-gene.name[order1]
    
    filter.list.1<-gene.min[order1]
    
    up.prob=1
    if(up.prob!=1){
      order2<-order(filter.list.1,decreasing = F)[1:round((n.feature.2*up.prob))]  ## second round of filtering genes
      order3<-order(filter.list.1,decreasing = T)[1:round((n.feature.2*(1-up.prob)))]
      gene.list<-gene.name.order1[c(order2,order3) ]                   ## final features list
    }else{
      order2<-order(filter.list.1,decreasing = F)[1:round((n.feature.2*up.prob))]  ## second round of filtering genes
      gene.list<-gene.name.order1[c(order2) ]                   ## final features list
    }
    
    up.prob=0.8
    if(up.prob!=1){
      order2<-order(filter.list.1,decreasing = F)[1:round((n.feature.2*up.prob))]  ## second round of filtering genes
      order3<-order(filter.list.1,decreasing = T)[1:round((n.feature.2*(1-up.prob)))]
      gene.list.up<-gene.name.order1[c(order2,order3) ]                   ## final features list
    }else{
      order2<-order(filter.list.1,decreasing = F)[1:round((n.feature.2*up.prob))]  ## second round of filtering genes
      gene.list<-gene.name.order1[c(order2) ]                   ## final features list
    }
    
    
    up.prob=0.5
    if(up.prob!=1){
      order2<-order(filter.list.1,decreasing = F)[1:round((n.feature.2*up.prob))]  ## second round of filtering genes
      order3<-order(filter.list.1,decreasing = T)[1:round((n.feature.2*(1-up.prob)))]
      gene.list.mid<-gene.name.order1[c(order2,order3) ]                   ## final features list
    }else{
      order2<-order(filter.list.1,decreasing = F)[1:round((n.feature.2*up.prob))]  ## second round of filtering genes
      gene.list<-gene.name.order1[c(order2) ]                   ## final features list
    }
    
    
    up.prob=0.2
    if(up.prob!=1){
      order2<-order(filter.list.1,decreasing = F)[1:round((n.feature.2*up.prob))]  ## second round of filtering genes
      order3<-order(filter.list.1,decreasing = T)[1:round((n.feature.2*(1-up.prob)))]
      gene.list.down<-gene.name.order1[c(order2,order3) ]                   ## final features list
    }else{
      order2<-order(filter.list.1,decreasing = F)[1:round((n.feature.2*up.prob))]  ## second round of filtering genes
      gene.list<-gene.name.order1[c(order2) ]                   ## final features list
    }
    
    
    assign(paste0("mcadet_down","_Sparsim_", j, "_", number.list[k]), gene.list.down)
    assign(paste0("mcadet_mid","_Sparsim_", j, "_", number.list[k]), gene.list.mid)
    assign(paste0("mcadet_up","_Sparsim_", j, "_", number.list[k]), gene.list.up)
    assign(paste0("mcadet","_Sparsim_", j, "_", number.list[k]), gene.list)
    assign(paste0("NB","_Sparsim_", j, "_", number.list[k]), get(paste0("NB_Sparsim_", j))[1:number.list[k]] )
    assign(paste0("Bre","_Sparsim_", j, "_", number.list[k]), get(paste0("Bre_Sparsim_", j))[1:number.list[k]] )
    assign(paste0("Seurat","_Sparsim_", j, "_", number.list[k]), get(paste0("Seurat_Sparsim_", j))[1:number.list[k]] )
    assign(paste0("disp","_Sparsim_", j, "_", number.list[k]), get(paste0("disp_Sparsim_", j))[1:number.list[k]] )
    assign(paste0("random","_Sparsim_", j, "_", number.list[k]), get(paste0("random_Sparsim_", j))[1:number.list[k]] )
    
    one.time.run<-c(mean(gene.list.down%in%rownames(SparSim_1)[1:2000]),
                    mean(gene.list.mid%in%rownames(SparSim_1)[1:2000]),
                    mean(gene.list.up%in%rownames(SparSim_1)[1:2000]),
                    mean(gene.list%in%rownames(SparSim_1)[1:2000]),
                    mean(get(paste0("NB_Sparsim_", j))[1:number.list[k]]%in%rownames(SparSim_1)[1:2000]),
                    mean(get(paste0("Bre_Sparsim_", j))[1:number.list[k]]%in%rownames(SparSim_1)[1:2000]),
                    mean(get(paste0("Seurat_Sparsim_", j))[1:number.list[k]]%in%rownames(SparSim_1)[1:2000]),
                    mean(get(paste0("disp_Sparsim_", j))[1:number.list[k]]%in%rownames(SparSim_1)[1:2000]),
                    mean(get(paste0("random_Sparsim_", j))[1:number.list[k]]%in%rownames(SparSim_1)[1:2000]))
    
    pt.matrix[j,]<-one.time.run
    a=((k-1)*16+j)/320*100
    print(paste0("Progress: ", a, "%"))
    
  }
  
  df_1<-append(df_1, apply(pt.matrix,2,mean) )
  df_2<-append(df_2, apply(pt.matrix,2,sd))
  
  
}




line.chart.pt_df<-data.frame("pt"=df_1, "sd"=df_2,
                             "n_genes"=rep(seq(100,2000,100),each=9), "method"=rep(c("Mcadet.down", "Mcadet.mid",
                                                                                     "Mcadet.up", "Mcadet", "NBdrop",
                                                                                     "Brennecke", "Seurat", "Disp", "Random"),20) )



line.chart.pt_df$method<-factor(line.chart.pt_df$method, 
                                levels = c("Mcadet","Mcadet.up", "Mcadet.mid","Mcadet.down",
                                           "Seurat", "NBdrop", "Brennecke", "Disp","Random"))

save(line.chart.pt_df, file = "E:/Thesis/SPARSim/line.chart.pt_df.RData")

load(file = "E:/Thesis/SPARSim/line.chart.jac_df.RData")
load(file = "E:/Thesis/SPARSim/line.chart.pt_df.RData")




method.remove.pt<-c( which(line.chart.pt_df$method=="Mcadet.down" ),
                     which(line.chart.pt_df$method=="Mcadet.mid" ), 
                     which(line.chart.pt_df$method=="Mcadet.up" ))

method.remove.jac<-c( which(line.chart.jac_df$method=="Mcadet.down" ),
                      which(line.chart.jac_df$method=="Mcadet.mid" ), 
                      which(line.chart.jac_df$method=="Mcadet.up" ))

library(ggplot2)

linechart_Sparsim=ggplot(data = line.chart.jac_df[-method.remove.jac,], mapping = aes(x = n_genes, y = Mean_Jaccard_similarity_index, 
                                                 colour = method, shape=method,linetype = method))+
  geom_line(size=0.8)+geom_point()+theme_bw()+theme(legend.position="none")+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Mean Jaccard similarity index")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(panel.background = element_blank(),axis.line = element_blank(),
        panel.grid.major = element_blank())+
  scale_x_continuous(limits=c(200,2000), breaks = seq(200,2000,200))


geom_errorbar(aes(ymin=Mean_Jaccard_similarity_index-sd, ymax=Mean_Jaccard_similarity_index+sd), 
              width=0.1,colour="black") 


j=ggplot(data = line.chart.pt_df, mapping = aes(x = n_genes, y = pt, 
                                                colour = method, shape=method,linetype = method))+
  geom_errorbar(aes(ymin=pt-sd, ymax=pt+sd), 
                width=0.1,colour="black")+
  geom_line(size=0.8)+geom_point()+theme_bw()+
  scale_shape_manual(values=seq(0,15))+labs(x="Number of genes selected", y="Percentage of true IGs")+
  theme(axis.title.x = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 10,  face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(panel.background = element_blank(),axis.line = element_blank(),
        panel.grid.major = element_blank())+
  scale_x_continuous(limits=c(200,2000), breaks = seq(200,2000,200))




library(patchwork)  


(b+a)/(d+c)
