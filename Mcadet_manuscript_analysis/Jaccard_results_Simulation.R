



true_SparSim_vector = c()
for(i in 1:2000){
  true_SparSim_vector = append(true_SparSim_vector, paste0("Gene", i))
  
}


## fine resolution


SparSim_fine_jaccard_macdet = c()
for (i in 1:24){
    jac = jaccard(true_SparSim_vector,
                  get(paste0("Mcadet_SparSim_fine_", i))$gene)
    SparSim_fine_jaccard_macdet = append(SparSim_fine_jaccard_macdet, jac)
  
}
SparSim_fine_jaccard_macdet
mean(SparSim_fine_jaccard_macdet)



SparSim_fine_jaccard_NBDrop = c()
for (i in 1:24){
  jac = jaccard(true_SparSim_vector,
                get(paste0("NBDrop_SparSim_fine_", i)))
  SparSim_fine_jaccard_NBDrop = append(SparSim_fine_jaccard_NBDrop, jac)
  
}
SparSim_fine_jaccard_NBDrop
mean(SparSim_fine_jaccard_NBDrop)




SparSim_fine_jaccard_M3Drop = c()
for (i in 1:24){
  jac = jaccard(true_SparSim_vector,
                get(paste0("M3Drop_SparSim_fine_", i)))
  SparSim_fine_jaccard_M3Drop = append(SparSim_fine_jaccard_M3Drop, jac)
  
}
SparSim_fine_jaccard_M3Drop
mean(SparSim_fine_jaccard_M3Drop)




SparSim_fine_jaccard_SeuratDisp = c()
for (i in 1:24){
  jac = jaccard(true_SparSim_vector,
                get(paste0("SeuratDisp_SparSim_fine_", i)))
  SparSim_fine_jaccard_SeuratDisp = append(SparSim_fine_jaccard_SeuratDisp, jac)
  
}
SparSim_fine_jaccard_SeuratDisp
mean(SparSim_fine_jaccard_SeuratDisp)




SparSim_fine_jaccard_SeuratMvp = c()
for (i in 1:24){
  jac = jaccard(true_SparSim_vector,
                get(paste0("SeuratMvp_SparSim_fine_", i)))
  SparSim_fine_jaccard_SeuratMvp = append(SparSim_fine_jaccard_SeuratMvp, jac)
  
}
SparSim_fine_jaccard_SeuratMvp
mean(SparSim_fine_jaccard_SeuratMvp)




SparSim_fine_jaccard_SeuratVst = c()
for (i in 1:24){
  jac = jaccard(true_SparSim_vector,
                get(paste0("SeuratVst_SparSim_fine_", i)))
  SparSim_fine_jaccard_SeuratVst = append(SparSim_fine_jaccard_SeuratVst, jac)
  
}
SparSim_fine_jaccard_SeuratVst
mean(SparSim_fine_jaccard_SeuratVst)




SparSim_fine_jaccard_scryFS = c()
for (i in 1:24){
  jac = jaccard(true_SparSim_vector,
                get(paste0("scryFS_SparSim_fine_", i)))
  SparSim_fine_jaccard_scryFS = append(SparSim_fine_jaccard_scryFS, jac)
  
}
SparSim_fine_jaccard_scryFS
mean(SparSim_fine_jaccard_scryFS)



SparSim_fine_jaccard_BreFS = c()
for (i in 1:24){
  jac = jaccard(true_SparSim_vector,
                get(paste0("BreFS_SparSim_fine_", i)))
  SparSim_fine_jaccard_BreFS = append(SparSim_fine_jaccard_BreFS, jac)
  
}
SparSim_fine_jaccard_BreFS
mean(SparSim_fine_jaccard_BreFS)




SparSim_fine_jaccard_random = c()
for (i in 1:24){
  jac = jaccard(true_SparSim_vector,
                get(paste0("random_SparSim_fine_", i)))
  SparSim_fine_jaccard_random = append(SparSim_fine_jaccard_random, jac)
  
}
SparSim_fine_jaccard_random
mean(SparSim_fine_jaccard_random)



Df_jaccard_sparsim_fine = data.frame("Jaccard" = c(SparSim_fine_jaccard_macdet,
SparSim_fine_jaccard_NBDrop,
SparSim_fine_jaccard_M3Drop,
SparSim_fine_jaccard_BreFS,
SparSim_fine_jaccard_SeuratDisp,
SparSim_fine_jaccard_SeuratMvp,
SparSim_fine_jaccard_SeuratVst,
SparSim_fine_jaccard_scryFS,
SparSim_fine_jaccard_random),

"Method" = c(rep("Mcadet", 24),
  rep("NBDrop", 24),
  rep("M3Drop", 24),
  rep("Brennecke", 24),
  rep("Seurat Disp", 24),
  rep("Seurat Mvp", 24),
  rep("Seurat Vst", 24),
  rep("Scry", 24),
  rep("Random", 24)))

Df_jaccard_sparsim_fine$Method = factor(Df_jaccard_sparsim_fine$Method,
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

cbPalette <- c("#cc340c", "#CC79A7", "#56B4E9",
               "#009E73", "#E69F00", "#F0E442", 
               "#999999","#0072B2","#D55E00")



df = Df_jaccard_sparsim_fine
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



box_jaccard_simulated_fine <- ggplot(df, 
                                       aes(Method, Jaccard, color = Method, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.6, linewidth = 1.1) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 1.3, alpha = 0.9) +
  guides(fill = "none") +
  labs(x = "", y = "", title = "Simulated fine-resolution datasets") +
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
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -1, vjust = 1.5, size = 8, fontface = "bold")+
  geom_text(data = significance_data, aes(x = Method, y = y, label = Significance), 
            size = 5, color = "black", fontface = "bold") 




## Coarse resolution


SparSim_jaccard_macdet = c()
for (i in 1:24){
  jac = jaccard(true_SparSim_vector,
                get(paste0("Mcadet_SparSim_", i))$gene)
  SparSim_jaccard_macdet = append(SparSim_jaccard_macdet, jac)
  
}
SparSim_jaccard_macdet
mean(SparSim_jaccard_macdet)



SparSim_jaccard_NBDrop = c()
for (i in 1:24){
  jac = jaccard(true_SparSim_vector,
                get(paste0("NBDrop_SparSim_", i)))
  SparSim_jaccard_NBDrop = append(SparSim_jaccard_NBDrop, jac)
  
}
SparSim_jaccard_NBDrop
mean(SparSim_jaccard_NBDrop)




SparSim_jaccard_M3Drop = c()
for (i in 1:24){
  jac = jaccard(true_SparSim_vector,
                get(paste0("M3Drop_SparSim_", i)))
  SparSim_jaccard_M3Drop = append(SparSim_jaccard_M3Drop, jac)
  
}
SparSim_jaccard_M3Drop
mean(SparSim_jaccard_M3Drop)




SparSim_jaccard_SeuratDisp = c()
for (i in 1:24){
  jac = jaccard(true_SparSim_vector,
                get(paste0("SeuratDisp_SparSim_", i)))
  SparSim_jaccard_SeuratDisp = append(SparSim_jaccard_SeuratDisp, jac)
  
}
SparSim_jaccard_SeuratDisp
mean(SparSim_jaccard_SeuratDisp)




SparSim_jaccard_SeuratMvp = c()
for (i in 1:24){
  jac = jaccard(true_SparSim_vector,
                get(paste0("SeuratMvp_SparSim_", i)))
  SparSim_jaccard_SeuratMvp = append(SparSim_jaccard_SeuratMvp, jac)
  
}
SparSim_jaccard_SeuratMvp
mean(SparSim_jaccard_SeuratMvp)




SparSim_jaccard_SeuratVst = c()
for (i in 1:24){
  jac = jaccard(true_SparSim_vector,
                get(paste0("SeuratVst_SparSim_", i)))
  SparSim_jaccard_SeuratVst = append(SparSim_jaccard_SeuratVst, jac)
  
}
SparSim_jaccard_SeuratVst
mean(SparSim_jaccard_SeuratVst)




SparSim_jaccard_scryFS = c()
for (i in 1:24){
  jac = jaccard(true_SparSim_vector,
                get(paste0("scryFS_SparSim_", i)))
  SparSim_jaccard_scryFS = append(SparSim_jaccard_scryFS, jac)
  
}
SparSim_jaccard_scryFS
mean(SparSim_jaccard_scryFS)




SparSim_jaccard_BreFS = c()
for (i in 1:24){
  jac = jaccard(true_SparSim_vector,
                get(paste0("BreFS_SparSim_", i)))
  SparSim_jaccard_BreFS = append(SparSim_jaccard_BreFS, jac)
  
}
SparSim_jaccard_BreFS
mean(SparSim_jaccard_BreFS)




SparSim_jaccard_random = c()
for (i in 1:24){
  jac = jaccard(true_SparSim_vector,
                get(paste0("random_SparSim_", i)))
  SparSim_jaccard_random = append(SparSim_jaccard_random, jac)
  
}
SparSim_jaccard_random
mean(SparSim_jaccard_random)











Df_jaccard_sparsim = data.frame("Jaccard" = c(SparSim_jaccard_macdet,
                                                   SparSim_jaccard_NBDrop,
                                                   SparSim_jaccard_M3Drop,
                                                   SparSim_jaccard_BreFS,
                                                   SparSim_jaccard_SeuratDisp,
                                                   SparSim_jaccard_SeuratMvp,
                                                   SparSim_jaccard_SeuratVst,
                                                   SparSim_jaccard_scryFS,
                                                   SparSim_jaccard_random),
                                     
                                     "Method" = c(rep("Mcadet", 24),
                                                  rep("NBDrop", 24),
                                                  rep("M3Drop", 24),
                                                  rep("Brennecke", 24),
                                                  rep("Seurat Disp", 24),
                                                  rep("Seurat Mvp", 24),
                                                  rep("Seurat Vst", 24),
                                                  rep("Scry", 24),
                                                  rep("Random", 24)))




Df_jaccard_sparsim$Method = factor(Df_jaccard_sparsim$Method,
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

cbPalette <- c("#cc340c", "#CC79A7", "#56B4E9",
               "#009E73", "#E69F00", "#F0E442", 
               "#999999","#0072B2","#D55E00")


df = Df_jaccard_sparsim
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



box_jaccard_simulated_coarse <- ggplot(df, 
                                  aes(Method, Jaccard, color = Method, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.6, linewidth = 1.1) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 1.3, alpha = 0.9) +
  guides(fill = "none") +
  labs(x = "", y = "Jaccard Similarity", title = "Simulated coarse-resolution datasets") +
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
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -1, vjust = 1.5, size = 8, fontface = "bold")+
  geom_text(data = significance_data, aes(x = Method, y = y, label = Significance), 
            size = 5, color = "black", fontface = "bold") 





(box_jaccard_PBMC_coarse + box_jaccard_PBMC_fine) / (box_jaccard_simulated_coarse + box_jaccard_simulated_fine)







