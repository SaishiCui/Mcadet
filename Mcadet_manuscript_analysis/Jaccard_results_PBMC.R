
jaccard <-function(a,b){
  jac = length(intersect(a,b))/(length(union(a,b)))
  return(jac)
}


## fine resolution


### load the data 

for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/DE_genes/DEgene_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
  }
}


for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/pbmc_results/Mcadet_pbmc_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
  }
}


#######################################







pbmc_fine_jaccard_macdet = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
jac = jaccard(get(paste0("DEgene_fine_", donor.set[p], "_", time.set[k])),
        get(paste0("Mcadet_pbmc_fine_", donor.set[p], "_", time.set[k]))$gene )
pbmc_fine_jaccard_macdet = append(pbmc_fine_jaccard_macdet, jac)
  }
}
pbmc_fine_jaccard_macdet
mean(pbmc_fine_jaccard_macdet)






pbmc_fine_jaccard_NBDrop = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_fine_", donor.set[p], "_", time.set[k])),
                  get(paste0("NBDrop_PBMC_fine_", donor.set[p], "_", time.set[k])) )
    pbmc_fine_jaccard_NBDrop = append(pbmc_fine_jaccard_NBDrop, jac)
  }
}
pbmc_fine_jaccard_NBDrop
mean(pbmc_fine_jaccard_NBDrop)




pbmc_fine_jaccard_M3Drop = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_fine_", donor.set[p], "_", time.set[k])),
                  get(paste0("M3Drop_PBMC_fine_", donor.set[p], "_", time.set[k])) )
    pbmc_fine_jaccard_M3Drop = append(pbmc_fine_jaccard_M3Drop, jac)
  }
}
pbmc_fine_jaccard_M3Drop
mean(pbmc_fine_jaccard_M3Drop)





pbmc_fine_jaccard_SeuratDisp = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_fine_", donor.set[p], "_", time.set[k])),
                  get(paste0("SeuratDisp_PBMC_fine_", donor.set[p], "_", time.set[k])) )
    pbmc_fine_jaccard_SeuratDisp = append(pbmc_fine_jaccard_SeuratDisp, jac)
  }
}
pbmc_fine_jaccard_SeuratDisp
mean(pbmc_fine_jaccard_SeuratDisp)





pbmc_fine_jaccard_SeuratMvp = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_fine_", donor.set[p], "_", time.set[k])),
                  get(paste0("SeuratMvp_PBMC_fine_", donor.set[p], "_", time.set[k])) )
    pbmc_fine_jaccard_SeuratMvp = append(pbmc_fine_jaccard_SeuratMvp, jac)
  }
}
pbmc_fine_jaccard_SeuratMvp
mean(pbmc_fine_jaccard_SeuratMvp)




pbmc_fine_jaccard_SeuratVst = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_fine_", donor.set[p], "_", time.set[k])),
                  get(paste0("SeuratVst_PBMC_fine_", donor.set[p], "_", time.set[k])) )
    pbmc_fine_jaccard_SeuratVst = append(pbmc_fine_jaccard_SeuratVst, jac)
  }
}
pbmc_fine_jaccard_SeuratVst
mean(pbmc_fine_jaccard_SeuratVst)







pbmc_fine_jaccard_scryFS = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_fine_", donor.set[p], "_", time.set[k])),
                  get(paste0("scryFS_PBMC_fine_", donor.set[p], "_", time.set[k])) )
    pbmc_fine_jaccard_scryFS = append(pbmc_fine_jaccard_scryFS, jac)
  }
}
pbmc_fine_jaccard_scryFS
mean(pbmc_fine_jaccard_scryFS)




pbmc_fine_jaccard_BreFS = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_fine_", donor.set[p], "_", time.set[k])),
                  get(paste0("BreFS_PBMC_fine_", donor.set[p], "_", time.set[k])) )
    pbmc_fine_jaccard_BreFS = append(pbmc_fine_jaccard_BreFS, jac)
  }
}
pbmc_fine_jaccard_BreFS
mean(pbmc_fine_jaccard_BreFS)




pbmc_fine_jaccard_random = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_fine_", donor.set[p], "_", time.set[k])),
                  get(paste0("random_PBMC_fine_", donor.set[p], "_", time.set[k])) )
    pbmc_fine_jaccard_random = append(pbmc_fine_jaccard_random, jac)
  }
}
pbmc_fine_jaccard_random
mean(pbmc_fine_jaccard_random)








Df_jaccard_pbmc_fine = data.frame("Jaccard" = c(pbmc_fine_jaccard_macdet,
                                                pbmc_fine_jaccard_NBDrop,
                                                pbmc_fine_jaccard_M3Drop,
                                                pbmc_fine_jaccard_BreFS,
                                                pbmc_fine_jaccard_SeuratDisp,
                                                pbmc_fine_jaccard_SeuratMvp,
                                                pbmc_fine_jaccard_SeuratVst,
                                                pbmc_fine_jaccard_scryFS,
                                                pbmc_fine_jaccard_random),
                                     
                                     "Method" = c(rep("Mcadet", 24),
                                                  rep("NBDrop", 24),
                                                  rep("M3Drop", 24),
                                                  rep("Brennecke", 24),
                                                  rep("Seurat Disp", 24),
                                                  rep("Seurat Mvp", 24),
                                                  rep("Seurat Vst", 24),
                                                  rep("Scry", 24),
                                                  rep("Random", 24)))



Df_jaccard_pbmc_fine$Method = factor(Df_jaccard_pbmc_fine$Method,
                                        levels = c("Mcadet",
                                                   "Seurat Vst",
                                                   "Seurat Disp",
                                                   "Seurat Mvp",
                                                   "NBDrop",
                                                   "M3Drop",
                                                   "Brennecke",
                                                   "Scry",
                                                   "Random"))

library(ggplot2)
library(ggpubr)
library(patchwork)


cbPalette <- c("#cc340c", "#CC79A7", "#56B4E9",
               "#009E73", "#E69F00", "#F0E442", 
               "#999999","#0072B2","#D55E00")








df = Df_jaccard_pbmc_fine
for (method in c("Seurat Vst",
                 "Seurat Disp",
                 "Seurat Mvp",
                 "NBDrop",
                 "M3Drop",
                 "Brennecke",
                 "Scry",
                 "Random")) {
  group1 <- df[df$Method == "Mcadet",]$Jaccard
  group2 <- df[df$Method == method,]$Jaccard
  
  test_result <- t.test(group1, group2, alternative = "greater")
  
  if (test_result$p.value >= 0.05) {
    df$Significance[df$Method == method] <- "NS"
  } else if( test_result$p.value < 0.05 & test_result$p.value >= 0.01 ){
    df$Significance[df$Method == method] <- "*"
  }else if( test_result$p.value < 0.01 & test_result$p.value >= 0.001 ){
    df$Significance[df$Method == method] <- "**"
  }else if( test_result$p.value < 0.001){
    df$Significance[df$Method == method] <- "***"
  }
  
  
}



significance_data <- unique(df[, c("Method", "Significance")])
significance_data$y <- 0.6
significance_data$Significance[1] = ""





box_jaccard_PBMC_fine <- ggplot(df, 
                                  aes(Method, Jaccard, color = Method, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.6, linewidth = 1.1) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 1.3, alpha = 0.9) +
  guides(fill = "none") +
  labs(x = "", y = "", title = "PBMC fine-resolution datasets") +
  theme_bw() + 
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette)+
  geom_hline(yintercept = mean(df$Jaccard), linetype = 2, linewidth = 1) +
  theme(
    panel.grid = element_blank(),  
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.position = "none",
    axis.text.x = element_text(face = "bold", size = 13, color = "black"), 
    axis.text.y = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", color = "black", size = 15), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.9) 
  )+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -1, vjust = 1.5, size = 8, fontface = "bold")+
  geom_text(data = significance_data, aes(x = Method, y = y, label = Significance), 
            size = 5, color = "black", fontface = "bold") 







## Coarse resolution




for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/DE_genes/DEgene_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    
  }
}


pbmc_jaccard_macdet = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_", donor.set[p], "_", time.set[k])),
                  get(paste0("Mcadet_pbmc_", donor.set[p], "_", time.set[k]))$gene )
    pbmc_jaccard_macdet = append(pbmc_jaccard_macdet, jac)
  }
}
pbmc_jaccard_macdet
mean(pbmc_jaccard_macdet)






pbmc_jaccard_NBDrop = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_", donor.set[p], "_", time.set[k])),
                  get(paste0("NBDrop_PBMC_", donor.set[p], "_", time.set[k])) )
    pbmc_jaccard_NBDrop = append(pbmc_jaccard_NBDrop, jac)
  }
}
pbmc_jaccard_NBDrop
mean(pbmc_jaccard_NBDrop)




pbmc_jaccard_M3Drop = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_", donor.set[p], "_", time.set[k])),
                  get(paste0("M3Drop_PBMC_", donor.set[p], "_", time.set[k])) )
    pbmc_jaccard_M3Drop = append(pbmc_jaccard_M3Drop, jac)
  }
}
pbmc_jaccard_M3Drop
mean(pbmc_jaccard_M3Drop)





pbmc_jaccard_SeuratDisp = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_", donor.set[p], "_", time.set[k])),
                  get(paste0("SeuratDisp_PBMC_", donor.set[p], "_", time.set[k])) )
    pbmc_jaccard_SeuratDisp = append(pbmc_jaccard_SeuratDisp, jac)
  }
}
pbmc_jaccard_SeuratDisp
mean(pbmc_jaccard_SeuratDisp)





pbmc_jaccard_SeuratMvp = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_", donor.set[p], "_", time.set[k])),
                  get(paste0("SeuratMvp_PBMC_", donor.set[p], "_", time.set[k])) )
    pbmc_jaccard_SeuratMvp = append(pbmc_jaccard_SeuratMvp, jac)
  }
}
pbmc_jaccard_SeuratMvp
mean(pbmc_jaccard_SeuratMvp)




pbmc_jaccard_SeuratVst = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_", donor.set[p], "_", time.set[k])),
                  get(paste0("SeuratVst_PBMC_", donor.set[p], "_", time.set[k])) )
    pbmc_jaccard_SeuratVst = append(pbmc_jaccard_SeuratVst, jac)
  }
}
pbmc_jaccard_SeuratVst
mean(pbmc_jaccard_SeuratVst)







pbmc_jaccard_scryFS = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_", donor.set[p], "_", time.set[k])),
                  get(paste0("scryFS_PBMC_", donor.set[p], "_", time.set[k])) )
    pbmc_jaccard_scryFS = append(pbmc_jaccard_scryFS, jac)
  }
}
pbmc_jaccard_scryFS
mean(pbmc_jaccard_scryFS)




pbmc_jaccard_BreFS = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_", donor.set[p], "_", time.set[k])),
                  get(paste0("BreFS_PBMC_", donor.set[p], "_", time.set[k])) )
    pbmc_jaccard_BreFS = append(pbmc_jaccard_BreFS, jac)
  }
}
pbmc_jaccard_BreFS
mean(pbmc_jaccard_BreFS)




pbmc_jaccard_random = c()
for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    jac = jaccard(get(paste0("DEgene_", donor.set[p], "_", time.set[k])),
                  get(paste0("random_PBMC_", donor.set[p], "_", time.set[k])) )
    pbmc_jaccard_random = append(pbmc_jaccard_random, jac)
  }
}
pbmc_jaccard_random
mean(pbmc_jaccard_random)







Df_jaccard_pbmc_coarse = data.frame("Jaccard" = c(pbmc_jaccard_macdet,
                                           pbmc_jaccard_NBDrop,
                                           pbmc_jaccard_M3Drop,
                                           pbmc_jaccard_BreFS,
                                           pbmc_jaccard_SeuratDisp,
                                           pbmc_jaccard_SeuratMvp,
                                           pbmc_jaccard_SeuratVst,
                                           pbmc_jaccard_scryFS,
                                           pbmc_jaccard_random),
                                  
                                  "Method" = c(rep("Mcadet", 24),
                                               rep("NBDrop", 24),
                                               rep("M3Drop", 24),
                                               rep("Brennecke", 24),
                                               rep("Seurat Disp", 24),
                                               rep("Seurat Mvp", 24),
                                               rep("Seurat Vst", 24),
                                               rep("Scry", 24),
                                               rep("Random", 24)))



Df_jaccard_pbmc_coarse$Method = factor(Df_jaccard_pbmc_coarse$Method,
                                     levels = c("Mcadet",
                                                "Seurat Vst",
                                                "Seurat Disp",
                                                "Seurat Mvp",
                                                "NBDrop",
                                                "M3Drop",
                                                "Brennecke",
                                                "Scry",
                                                "Random"))


cbPalette <- c("#cc340c", "#CC79A7", "#56B4E9",
               "#009E73", "#E69F00", "#F0E442", 
               "#999999","#0072B2","#D55E00")




df = Df_jaccard_pbmc_coarse
for (method in c("Seurat Vst",
                 "Seurat Disp",
                 "Seurat Mvp",
                 "NBDrop",
                 "M3Drop",
                 "Brennecke",
                 "Scry",
                 "Random")) {
  group1 <- df[df$Method == "Mcadet",]$Jaccard
  group2 <- df[df$Method == method,]$Jaccard
  
  test_result <- t.test(group1, group2, alternative = "greater")
  
  if (test_result$p.value >= 0.05) {
    df$Significance[df$Method == method] <- "NS"
  } else if( test_result$p.value < 0.05 & test_result$p.value >= 0.01 ){
    df$Significance[df$Method == method] <- "*"
  }else if( test_result$p.value < 0.01 & test_result$p.value >= 0.001 ){
    df$Significance[df$Method == method] <- "**"
  }else if( test_result$p.value < 0.001){
    df$Significance[df$Method == method] <- "***"
  }
  
  
}



significance_data <- unique(df[, c("Method", "Significance")])
significance_data$y <- 0.6
significance_data$Significance[1] = ""



box_jaccard_PBMC_coarse <- ggplot(df, 
                            aes(Method, Jaccard, color = Method, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.6, linewidth = 1.1) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 1.3, alpha = 0.9) +
  guides(fill = "none") +
  labs(x = "", y = "Jaccard Similarity", title = "PBMC coarse-resolution datasets") +
  theme_bw() + 
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette)+
  geom_hline(yintercept = mean(df$Jaccard), linetype = 2, linewidth = 1) +
  theme(
    panel.grid = element_blank(),  
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.position = "none",
    axis.text.x = element_text(face = "bold", size = 13, color = "black"), 
    axis.text.y = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", color = "black", size = 15),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.9) 
  )+
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -1, vjust = 1.5, size = 8, fontface = "bold")+
  geom_text(data = significance_data, aes(x = Method, y = y, label = Significance), 
            size = 5, color = "black", fontface = "bold") 





box_jaccard_PBMC_coarse + box_jaccard_PBMC_fine 




