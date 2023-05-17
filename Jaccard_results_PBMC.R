## fine resolution

for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =paste0("C:/Users/13541/Desktop/Thesis/pbmc/DE_genes/DEgene_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
  }
}


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


Df_jaccard_pbmc_fine$"normalized" =(Df_jaccard_pbmc_fine$Jaccard -
                                      mean(Df_jaccard_pbmc_fine$Jaccard))/
  sd(Df_jaccard_pbmc_fine$Jaccard)


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

cbPalette <- c("#cc340c", "#CC79A7", "#56B4E9",
               "#009E73", "#E69F00", "#F0E442", 
               "#999999","#0072B2","#D55E00")

ggplot(Df_jaccard_pbmc_fine, 
       aes(Method, normalized, color = Method )) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  guides(fill = "none") +
  labs(x = "", y = "Normalized Jaccard Similarity", 
       title = "PBMC fine datasets")+theme_bw()+ 
  scale_color_manual(values = cbPalette)+
  stat_compare_means(label = "p.signif", method = "t.test",label.y = 4.2,
                     ref.group = ".all.")+
  geom_hline(yintercept = mean(Df_jaccard_pbmc_fine$normalized),
             linetype=2)+
  theme( panel.grid = element_blank(),  axis.title = element_text(face ="bold"),
         plot.title = element_text(hjust = 0.5, face= "bold", size = 13),
         legend.position = "none",
         axis.text.x = element_text(face = "bold", size = 13, color="black"), 
         axis.text.y = element_text(face = "bold", size = 13))








## Coarse resolution

for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =paste0("C:/Users/13541/Desktop/Thesis/pbmc/DE_genes/DEgene_", 
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







Df_jaccard_pbmc = data.frame("Jaccard" = c(pbmc_jaccard_macdet,
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


Df_jaccard_pbmc$"normalized" =(Df_jaccard_pbmc$Jaccard -
                               mean(Df_jaccard_pbmc$Jaccard))/
  sd(Df_jaccard_pbmc$Jaccard)


Df_jaccard_pbmc$Method = factor(Df_jaccard_pbmc$Method,
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

ggplot(Df_jaccard_pbmc, 
       aes(Method, normalized, color = Method )) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  guides(fill = "none") +
  labs(x = "", y = "Normalized Jaccard Similarity", 
       title = "PBMC coarse datasets")+theme_bw()+ 
  scale_color_manual(values = cbPalette)+
  stat_compare_means(label = "p.signif", method = "t.test",label.y = 4.2,
                     ref.group = ".all.")+
  geom_hline(yintercept = mean(Df_jaccard_pbmc$normalized),
             linetype=2)+
  theme( panel.grid = element_blank(),  axis.title = element_text(face ="bold"),
         plot.title = element_text(hjust = 0.5, face= "bold", size = 13),
         legend.position = "none",
         axis.text.x = element_text(face = "bold", size = 13, color="black"), 
         axis.text.y = element_text(face = "bold", size = 13))











