


 donor.set<-c(1,2,3,4,5,6,7,8)
 time.set<-c(0,3,7)
 i = 1
 j= 1
 
 load(file =paste0("C:/Users/css22/Desktop/Thesis/pbmc/data/pbmcdata_", 
                   donor.set[i], "_", time.set[j], ".RData"))
 
 
 
### PBMC coarse 

Mcadet_hist = c()
NBDrop_hist = c()
M3Drop_hist = c()
Bre_hist = c()
Scry_hist = c()
SeuratVst_hist = c()
SeuratDisp_hist = c()
SeuratMvp_hist = c()
random_hist = c()
TrueHVG_hist = c()



for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/data/pbmcdata_", 
                      donor.set[i], "_", time.set[j], ".RData"))
    work_data = get(paste0("pbmcdata_", donor.set[i], "_", time.set[j]))
    
    Mcadet_hist =  append( Mcadet_hist, log(apply(work_data[get(paste0("Mcadet_pbmc_", donor.set[i], "_", time.set[j]))$gene,], 1, mean)))
    NBDrop_hist =  append( NBDrop_hist, log(apply(work_data[get(paste0("NBDrop_PBMC_", donor.set[i], "_", time.set[j])),], 1, mean)))
    M3Drop_hist =  append( M3Drop_hist, log(apply(work_data[get(paste0("M3Drop_PBMC_", donor.set[i], "_", time.set[j])),], 1, mean)))
    Bre_hist =  append( Bre_hist, log(apply(work_data[get(paste0("BreFS_PBMC_", donor.set[i], "_", time.set[j])),], 1, mean)))
    Scry_hist =  append( Scry_hist, log(apply(work_data[get(paste0("scryFS_PBMC_", donor.set[i], "_", time.set[j])),], 1, mean)))
    SeuratVst_hist =  append( SeuratVst_hist, log(apply(work_data[get(paste0("SeuratVst_PBMC_", donor.set[i], "_", time.set[j])),], 1, mean)))
    SeuratDisp_hist =  append( SeuratDisp_hist, log(apply(work_data[get(paste0("SeuratDisp_PBMC_", donor.set[i], "_", time.set[j])),], 1, mean)))
    SeuratMvp_hist =  append( SeuratMvp_hist, log(apply(work_data[get(paste0("SeuratMvp_PBMC_", donor.set[i], "_", time.set[j])),], 1, mean)))
    random_hist =  append( random_hist, log(apply(work_data[get(paste0("random_PBMC_", donor.set[i], "_", time.set[j])),], 1, mean)))
    TrueHVG_hist =  append( TrueHVG_hist, log(apply(work_data[get(paste0("DEgene_", donor.set[i], "_", time.set[j])),], 1, mean)))
    
    
    rm(list = paste0("pbmcdata_", donor.set[i], "_", time.set[j]))
    gc()
    print(c(i,j))
  }
  
}


random_hist = random_hist[random_hist != -Inf]


histdf <- data.frame("log_mean_expression" = c(Mcadet_hist, NBDrop_hist, M3Drop_hist, Bre_hist,
                                               Scry_hist, SeuratVst_hist, SeuratDisp_hist, SeuratMvp_hist, random_hist),
                     "Method" = factor(c(rep("Mcadet", length(Mcadet_hist)), rep("NBDrop", length(NBDrop_hist)), rep("M3Drop", length(M3Drop_hist)),
                                rep("Brennecke", length(Bre_hist)), rep("Scry", length(Scry_hist)), rep("Seurat Vst", length(SeuratVst_hist)),
                                rep("Seurat Disp", length(SeuratDisp_hist)), rep("Seurat Mvp", length(SeuratMvp_hist)),
                                rep("random", length(random_hist))), levels = c("Mcadet", "NBDrop", "M3Drop", "Brennecke", "Scry", 
                                                                                "Seurat Vst", "Seurat Disp", "Seurat Mvp", "random")))

histmeandf <- aggregate(histdf$log_mean_expression, by=list(Method=histdf$Method), mean)
colnames(histmeandf) = c("Method", "xintercept")



library(ggplot2)
library(viridis)
library(hrbrthemes)
library(RColorBrewer)
display.brewer.all()






true_hvg_df <- data.frame("log_mean_expression" = rep(TrueHVG_hist, 9),
                          "Method" = factor(rep(levels(histdf$Method), each = length(TrueHVG_hist)), levels = levels(histdf$Method)))




ggplot() +
  geom_density(data = histdf, aes(x = log_mean_expression, fill = Method), alpha = 0.8) +  # 原有方法的密度
  geom_density(data = true_hvg_df, aes(x = log_mean_expression), fill = "#CCCCFF", alpha = 0.8) +  # TrueHVG 的密度
  facet_wrap(~ Method, scales = "free") +
  theme_bw() +
  labs(x = "log mean expression", y = "Density", 
       title = "" ) +  
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(colour = "red", fill = "#CCCCFF"),
    axis.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 12, color = "black"), 
    axis.text.y = element_text(face = "bold", size = 12)
  ) + 
  scale_x_continuous(breaks = seq(-8, 8, 2)) +
  geom_vline(data = histmeandf, aes(xintercept = xintercept), linetype = "dashed", color = "blue", linewidth = 0.8)  # 添加均值垂线













ggplot(histdf, aes(x=log_mean_expression, fill = Method))+geom_density(position = "fill")+
  labs(x = "log mean expression", y = "Density", 
       title = "PBMC coarse datasets" )+ theme(
    legend.position="none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.text= element_text(size=10, face="bold"),
    strip.background = element_rect(colour="red", fill="#CCCCFF"),
    axis.title = element_text(face ="bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 12, color = "black"), 
    axis.text.y = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.background = element_rect(color = "black", 
                                     linewidth = 1.5,
                                     fill = "lightblue"))+ 
  scale_x_continuous(breaks=seq(-8, 8, 2))










### PBMC fine 

Mcadet_hist = c()
NBDrop_hist = c()
M3Drop_hist = c()
Bre_hist = c()
Scry_hist = c()
SeuratVst_hist = c()
SeuratDisp_hist = c()
SeuratMvp_hist = c()
random_hist = c()
TrueHVG_hist = c()


for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/data_fine/pbmc_fine_", 
                      donor.set[i], "_", time.set[j], ".RData"))
    work_data = get(paste0("pbmc_fine_", donor.set[i], "_", time.set[j]))
    
    Mcadet_hist =  append( Mcadet_hist, log(apply(work_data[get(paste0("Mcadet_pbmc_fine_", donor.set[i], "_", time.set[j]))$gene,], 1, mean)))
    NBDrop_hist =  append( NBDrop_hist, log(apply(work_data[get(paste0("NBDrop_PBMC_fine_", donor.set[i], "_", time.set[j])),], 1, mean)))
    M3Drop_hist =  append( M3Drop_hist, log(apply(work_data[get(paste0("M3Drop_PBMC_fine_", donor.set[i], "_", time.set[j])),], 1, mean)))
    Bre_hist =  append( Bre_hist, log(apply(work_data[get(paste0("BreFS_PBMC_fine_", donor.set[i], "_", time.set[j])),], 1, mean)))
    Scry_hist =  append( Scry_hist, log(apply(work_data[get(paste0("scryFS_PBMC_fine_", donor.set[i], "_", time.set[j])),], 1, mean)))
    SeuratVst_hist =  append( SeuratVst_hist, log(apply(work_data[get(paste0("SeuratVst_PBMC_fine_", donor.set[i], "_", time.set[j])),], 1, mean)))
    SeuratDisp_hist =  append( SeuratDisp_hist, log(apply(work_data[get(paste0("SeuratDisp_PBMC_fine_", donor.set[i], "_", time.set[j])),], 1, mean)))
    SeuratMvp_hist =  append( SeuratMvp_hist, log(apply(work_data[get(paste0("SeuratMvp_PBMC_fine_", donor.set[i], "_", time.set[j])),], 1, mean)))
    random_hist =  append( random_hist, log(apply(work_data[get(paste0("random_PBMC_fine_", donor.set[i], "_", time.set[j])),], 1, mean)))
    TrueHVG_hist =  append( TrueHVG_hist, log(apply(work_data[get(paste0("DEgene_fine_", donor.set[i], "_", time.set[j])),], 1, mean)))
    
    
    rm(list = paste0("pbmc_fine_", donor.set[i], "_", time.set[j]))
    gc()
    print(c(i,j))
  }
  
}


random_hist = random_hist[random_hist != -Inf]


histdf <- data.frame("log_mean_expression" = c(Mcadet_hist, NBDrop_hist, M3Drop_hist, Bre_hist,
                                               Scry_hist, SeuratVst_hist, SeuratDisp_hist, SeuratMvp_hist, random_hist),
                     "Method" = factor(c(rep("Mcadet", length(Mcadet_hist)), rep("NBDrop", length(NBDrop_hist)), rep("M3Drop", length(M3Drop_hist)),
                                         rep("Brennecke", length(Bre_hist)), rep("Scry", length(Scry_hist)), rep("Seurat Vst", length(SeuratVst_hist)),
                                         rep("Seurat Disp", length(SeuratDisp_hist)), rep("Seurat Mvp", length(SeuratMvp_hist)),
                                         rep("random", length(random_hist))), levels = c("Mcadet", "NBDrop", "M3Drop", "Brennecke", "Scry", 
                                                                                         "Seurat Vst", "Seurat Disp", "Seurat Mvp", "random")))


histmeandf <- aggregate(histdf$log_mean_expression, by=list(Method=histdf$Method), mean)
colnames(histmeandf) = c("Method", "xintercept")



library(ggplot2)
library(viridis)
library(hrbrthemes)
library(RColorBrewer)
display.brewer.all()


true_hvg_df <- data.frame("log_mean_expression" = rep(TrueHVG_hist, 9),
                          "Method" = factor(rep(levels(histdf$Method), each = length(TrueHVG_hist)), levels = levels(histdf$Method)))




ggplot() +
  geom_density(data = histdf, aes(x = log_mean_expression, fill = Method), alpha = 0.8) +  # 原有方法的密度
  geom_density(data = true_hvg_df, aes(x = log_mean_expression), fill = "#CCCCFF", alpha = 0.8) +  # TrueHVG 的密度
  facet_wrap(~ Method, scales = "free") +
  theme_bw() +
  labs(x = "log mean expression", y = "Density", 
       title = "" ) +  
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(colour = "red", fill = "#CCCCFF"),
    axis.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 12, color = "black"), 
    axis.text.y = element_text(face = "bold", size = 12)
  ) + 
  scale_x_continuous(breaks = seq(-8, 8, 2)) +
  geom_vline(data = histmeandf, aes(xintercept = xintercept), linetype = "dashed", color = "blue", linewidth = 0.8)  # 添加均值垂线






ggplot(histdf, aes(x=log_mean_expression, fill = Method))+geom_density(position = "fill")+
  labs(x = "log mean expression", y = "Density", 
       title = "PBMC fine datasets" )+ theme(
         legend.position="none",
         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
         strip.text= element_text(size=10, face="bold"),
         strip.background = element_rect(colour="red", fill="#CCCCFF"),
         axis.title = element_text(face ="bold", size = 13),
         axis.text.x = element_text(face = "bold", size = 12, color = "black"), 
         axis.text.y = element_text(face = "bold", size = 12),
         legend.text = element_text(size = 10, face = "bold"),
         legend.title = element_text(size = 12, face = "bold"),
         legend.background = element_rect(color = "black", 
                                          linewidth = 1.5,
                                          fill = "lightblue"))+ 
  scale_x_continuous(breaks=seq(-8, 8, 2))








### SparSim coarse 

Mcadet_hist = c()
NBDrop_hist = c()
M3Drop_hist = c()
Bre_hist = c()
Scry_hist = c()
SeuratVst_hist = c()
SeuratDisp_hist = c()
SeuratMvp_hist = c()
random_hist = c()
TrueHVG_hist = c()



for (i in 1:24) {


    load(file =paste0("C:/Users/css22/Desktop/Thesis1/SPARSim/SparSim_Mcadet_data/SparSim_", i, ".RData"))
    work_data = get(paste0("SparSim_", i))
    
    Mcadet_hist =  append( Mcadet_hist, log(apply(work_data[get(paste0("Mcadet_SparSim_", i))$gene,], 1, mean)))
    NBDrop_hist =  append( NBDrop_hist, log(apply(work_data[get(paste0("NBDrop_SparSim_", i)),], 1, mean)))
    M3Drop_hist =  append( M3Drop_hist, log(apply(work_data[get(paste0("M3Drop_SparSim_", i)),], 1, mean)))
    Bre_hist =  append( Bre_hist, log(apply(work_data[get(paste0("BreFS_SparSim_", i)),], 1, mean)))
    Scry_hist =  append( Scry_hist, log(apply(work_data[get(paste0("scryFS_SparSim_", i)),], 1, mean)))
    SeuratVst_hist =  append( SeuratVst_hist, log(apply(work_data[get(paste0("SeuratVst_SparSim_", i)),], 1, mean)))
    SeuratDisp_hist =  append( SeuratDisp_hist, log(apply(work_data[get(paste0("SeuratDisp_SparSim_", i)),], 1, mean)))
    SeuratMvp_hist =  append( SeuratMvp_hist, log(apply(work_data[get(paste0("SeuratMvp_SparSim_", i)),], 1, mean)))
    random_hist =  append( random_hist, log(apply(work_data[get(paste0("random_SparSim_", i)),], 1, mean)))
    TrueHVG_hist =  append( TrueHVG_hist, log(apply(work_data[paste0("Gene", c(1:2000)),], 1, mean)))
    
    
    rm(list = paste0("SparSim_", i))
    gc()
    print(i)
  
  
}


random_hist = random_hist[random_hist != -Inf]


histdf <- data.frame("log_mean_expression" = c(Mcadet_hist, NBDrop_hist, M3Drop_hist, Bre_hist,
                                               Scry_hist, SeuratVst_hist, SeuratDisp_hist, SeuratMvp_hist, random_hist),
                     "Method" = factor(c(rep("Mcadet", length(Mcadet_hist)), rep("NBDrop", length(NBDrop_hist)), rep("M3Drop", length(M3Drop_hist)),
                                         rep("Brennecke", length(Bre_hist)), rep("Scry", length(Scry_hist)), rep("Seurat Vst", length(SeuratVst_hist)),
                                         rep("Seurat Disp", length(SeuratDisp_hist)), rep("Seurat Mvp", length(SeuratMvp_hist)),
                                         rep("random", length(random_hist))), levels = c("Mcadet", "NBDrop", "M3Drop", "Brennecke", "Scry", 
                                                                                         "Seurat Vst", "Seurat Disp", "Seurat Mvp", "random")))

histmeandf <- aggregate(histdf$log_mean_expression, by=list(Method=histdf$Method), mean)
colnames(histmeandf) = c("Method", "xintercept")



library(ggplot2)
library(viridis)
library(hrbrthemes)
library(RColorBrewer)





true_hvg_df <- data.frame("log_mean_expression" = rep(TrueHVG_hist, 9),
                          "Method" = factor(rep(levels(histdf$Method), each = length(TrueHVG_hist)), levels = levels(histdf$Method)))




ggplot() +
  geom_density(data = histdf, aes(x = log_mean_expression, fill = Method), alpha = 0.8) +  # 原有方法的密度
  geom_density(data = true_hvg_df, aes(x = log_mean_expression), fill = "#CCCCFF", alpha = 0.8) +  # TrueHVG 的密度
  facet_wrap(~ Method, scales = "free") +
  theme_bw() +
  labs(x = "log mean expression", y = "Density", 
       title = "" ) +  
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(colour = "red", fill = "#CCCCFF"),
    axis.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 12, color = "black"), 
    axis.text.y = element_text(face = "bold", size = 12)
  ) + 
  scale_x_continuous(breaks = seq(-8, 8, 2)) +
  geom_vline(data = histmeandf, aes(xintercept = xintercept), linetype = "dashed", color = "blue", linewidth = 0.8)  # 添加均值垂线





ggplot(histdf, aes(x=log_mean_expression, fill = Method))+geom_density(position = "fill")+
  labs(x = "log mean expression", y = "Density", 
       title = "Simulation coarse datasets" )+ theme(
         legend.position="none",
         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
         strip.text= element_text(size=10, face="bold"),
         strip.background = element_rect(colour="red", fill="#CCCCFF"),
         axis.title = element_text(face ="bold", size = 13),
         axis.text.x = element_text(face = "bold", size = 12, color = "black"), 
         axis.text.y = element_text(face = "bold", size = 12),
         legend.text = element_text(size = 10, face = "bold"),
         legend.title = element_text(size = 12, face = "bold"),
         legend.background = element_rect(color = "black", 
                                          linewidth = 1.5,
                                          fill = "lightblue"))+ 
  scale_x_continuous(breaks=seq(-8, 8, 2))






### SparSim  fine 

Mcadet_hist = c()
NBDrop_hist = c()
M3Drop_hist = c()
Bre_hist = c()
Scry_hist = c()
SeuratVst_hist = c()
SeuratDisp_hist = c()
SeuratMvp_hist = c()
random_hist = c()
TrueHVG_hist = c()


for (i in 1:24) {
  
  
  load(file =paste0("C:/Users/css22/Desktop/Thesis1/SPARSim/SparSim_Mcadet_data/SparSim_fine_", i, ".RData"))
  work_data = get(paste0("SparSim_fine_", i))
  
  Mcadet_hist =  append( Mcadet_hist, log(apply(work_data[get(paste0("Mcadet_SparSim_fine_", i))$gene,], 1, mean)))
  NBDrop_hist =  append( NBDrop_hist, log(apply(work_data[get(paste0("NBDrop_SparSim_fine_", i)),], 1, mean)))
  M3Drop_hist =  append( M3Drop_hist, log(apply(work_data[get(paste0("M3Drop_SparSim_fine_", i)),], 1, mean)))
  Bre_hist =  append( Bre_hist, log(apply(work_data[get(paste0("BreFS_SparSim_fine_", i)),], 1, mean)))
  Scry_hist =  append( Scry_hist, log(apply(work_data[get(paste0("scryFS_SparSim_fine_", i)),], 1, mean)))
  SeuratVst_hist =  append( SeuratVst_hist, log(apply(work_data[get(paste0("SeuratVst_SparSim_fine_", i)),], 1, mean)))
  SeuratDisp_hist =  append( SeuratDisp_hist, log(apply(work_data[get(paste0("SeuratDisp_SparSim_fine_", i)),], 1, mean)))
  SeuratMvp_hist =  append( SeuratMvp_hist, log(apply(work_data[get(paste0("SeuratMvp_SparSim_fine_", i)),], 1, mean)))
  random_hist =  append( random_hist, log(apply(work_data[get(paste0("random_SparSim_fine_", i)),], 1, mean)))
  TrueHVG_hist =  append( TrueHVG_hist, log(apply(work_data[paste0("Gene", c(1:2000)),], 1, mean)))
  
  rm(list = paste0("SparSim_fine_", i))
  gc()
  print(i)
  
  
}


random_hist = random_hist[random_hist != -Inf]


histdf <- data.frame("log_mean_expression" = c(Mcadet_hist, NBDrop_hist, M3Drop_hist, Bre_hist,
                                               Scry_hist, SeuratVst_hist, SeuratDisp_hist, SeuratMvp_hist, random_hist),
                     "Method" = factor(c(rep("Mcadet", length(Mcadet_hist)), rep("NBDrop", length(NBDrop_hist)), rep("M3Drop", length(M3Drop_hist)),
                                         rep("Brennecke", length(Bre_hist)), rep("Scry", length(Scry_hist)), rep("Seurat Vst", length(SeuratVst_hist)),
                                         rep("Seurat Disp", length(SeuratDisp_hist)), rep("Seurat Mvp", length(SeuratMvp_hist)),
                                         rep("random", length(random_hist))), levels = c("Mcadet", "NBDrop", "M3Drop", "Brennecke", "Scry", 
                                                                                         "Seurat Vst", "Seurat Disp", "Seurat Mvp", "random")))

histmeandf <- aggregate(histdf$log_mean_expression, by=list(Method=histdf$Method), mean)
colnames(histmeandf) = c("Method", "xintercept")



library(ggplot2)
library(viridis)
library(hrbrthemes)
library(RColorBrewer)





true_hvg_df <- data.frame("log_mean_expression" = rep(TrueHVG_hist, 9),
                          "Method" = factor(rep(levels(histdf$Method), each = length(TrueHVG_hist)), levels = levels(histdf$Method)))




ggplot() +
  geom_density(data = histdf, aes(x = log_mean_expression, fill = Method), alpha = 0.8) +  # 原有方法的密度
  geom_density(data = true_hvg_df, aes(x = log_mean_expression), fill = "#CCCCFF", alpha = 0.8) +  # TrueHVG 的密度
  facet_wrap(~ Method, scales = "free") +
  theme_bw() +
  labs(x = "log mean expression", y = "Density", 
       title = "" ) +  
  theme(
    legend.position = "None",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(colour = "red", fill = "#CCCCFF"),
    axis.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", size = 12, color = "black"), 
    axis.text.y = element_text(face = "bold", size = 12)
  ) + 
  scale_x_continuous(breaks = seq(-8, 8, 2)) +
  geom_vline(data = histmeandf, aes(xintercept = xintercept), linetype = "dashed", color = "blue", linewidth = 0.8)  # 添加均值垂线






ggplot(histdf, aes(x=log_mean_expression, fill = Method))+geom_density(position = "fill")+
  labs(x = "log mean expression", y = "Density", 
       title = "Simulation fine datasets" )+ theme(
         legend.position="none",
         plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
         strip.text= element_text(size=10, face="bold"),
         strip.background = element_rect(colour="red", fill="#CCCCFF"),
         axis.title = element_text(face ="bold", size = 13),
         axis.text.x = element_text(face = "bold", size = 12, color = "black"), 
         axis.text.y = element_text(face = "bold", size = 12),
         legend.text = element_text(size = 10, face = "bold"),
         legend.title = element_text(size = 12, face = "bold"),
         legend.background = element_rect(color = "black", 
                                          linewidth = 1.5,
                                          fill = "lightblue"))+ 
  scale_x_continuous(breaks=seq(-8, 8, 2))


