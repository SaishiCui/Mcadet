rm(list = ls())
gc()






for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/pbmc/DE_genes/DEgene_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/pbmc_results/Mcadet_pbmc_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/BreFS_PBMC_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/M3Drop_PBMC_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/NBDrop_PBMC_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/scryFS_PBMC_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/SeuratDisp_PBMC_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/SeuratMvp_PBMC_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/css22/Desktop/Thesis1/Generating_results/OhterFS/OtherFS_results_PBMC/SeuratVst_PBMC_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    
    
    
    
    load(file = paste0("C:/Users/css22/Desktop/Thesis1/pbmc/data_fine/pbmc_fine_clabel_", donor.set[p], "_", time.set[k], ".RData"))
    load(file = paste0("C:/Users/css22/Desktop/Thesis1/pbmc/data_fine/pbmc_fine_", donor.set[p], "_", time.set[k], ".RData"))
    
  }
}

#####################################

donor.set<-c(1,2,3,4,5,6,7,8)
time.set<-c(0,3,7)

new_discover_list = c()
for (p in 1:8) {
  for (k in 1:3) {
a = intersect(get(paste0("DEgene_fine_", donor.set[p] , "_", time.set[k])), 
              get(paste0("Mcadet_pbmc_fine_", donor.set[p] , "_", time.set[k]))$gene)


b = Reduce(union, list(get(paste0("BreFS_PBMC_fine_", donor.set[p] , "_", time.set[k])), 
                       get(paste0("M3Drop_PBMC_fine_", donor.set[p] , "_", time.set[k])), 
                       get(paste0("NBDrop_PBMC_fine_", donor.set[p] , "_", time.set[k])), 
                       get(paste0("scryFS_PBMC_fine_", donor.set[p] , "_", time.set[k])), 
                       get(paste0("SeuratDisp_PBMC_fine_", donor.set[p] , "_", time.set[k])), 
                       get(paste0("SeuratMvp_PBMC_fine_", donor.set[p] , "_", time.set[k])), 
                       get(paste0("SeuratVst_PBMC_fine_", donor.set[p] , "_", time.set[k]))))



new_discover_list = c(new_discover_list, setdiff(a,b))
   }
  }

sort( table(new_discover_list), decreasing = T)[1:20]

as.data.frame(sort( table(new_discover_list), decreasing = T)[1:20])


df <- data.frame(
  new_discover_list = c("SPINT2", "ADPRM", "CHMP7", "BCAS4", "CITED4", "CLN5", "EPHX2", "PRKCA", 
                        "ATP6V0E2", "CRLF3", "GYG1", "LYRM4", "SYPL1", "TMIGD2", "APBA2", "IMPDH2", 
                        "MAN1C1", "PRNP", "RAB9A", "SUSD3"),
  Freq = c(19, 18, 18, 17, 17, 16, 16, 16, 14, 14, 14, 14, 14, 14, 13, 13, 13, 13, 13, 13)
)


ggplot(df, aes(x = reorder(new_discover_list, Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", width = 0.7,
           colour = "black", linewidth = 1) +  
  coord_flip() +
  geom_text(aes(label = Freq), hjust = -0.3, size = 5, fontface = "bold") + 
  labs(title = "Frequency of informative genes discovered by Mcadet exclusively",
       x = "Gene", 
       y = "Frequency") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(), 
    panel.grid.major.y = element_line(color = "grey", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1.2),  
    axis.title.y = element_blank(), 
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),  
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),  
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),  
    legend.position = "none",  
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5) 
  )

####################### 



mean_exp_matrix = matrix(0, nrow = 24, ncol = 11)
donor.set<-c(1,2,3,4,5,6,7,8)
time.set<-c(0,3,7)

for (p in 1:8) {
  for (k in 1:3) {
    mean_exp = c()
    assign("dta", get(paste0("pbmc_fine_",  donor.set[p], "_",  time.set[k] )))
    assign("cl", get(paste0("pbmc_fine_clabel_",  donor.set[p], "_",  time.set[k] )))
    for (name in names(table(pbmc_fine_clabel_1_0))) {
      mexp <- mean(as.numeric(dta[which(rownames(dta) == "BCAS4"), which(cl == name)]))
      mean_exp <- c(mean_exp, mexp)
    }
    mean_exp_matrix[(p-1)*3+k,] = mean_exp
  }
}




df <- data.frame(mexp =apply(mean_exp_matrix, 2, mean),
                 ct = names(table(pbmc_fine_clabel_1_0))
                 )




ggplot(df, aes(x = ct, y = mexp, fill = ct)) +
  geom_bar(stat = "identity", color = "black", width = 0.6, linewidth = 1) + 
  labs(y = "Mean gene expression", title = "BCAS4 expression by fine-resolution cell types") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey", linetype = "dashed"),  
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1.2),  
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),  
    axis.text.x = element_text(size = 12, face = "bold", color = "black", angle = 45, hjust = 1),  
    axis.text.y = element_text(size = 12, face = "bold", color = "black"), 
    legend.position = "none",  
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  
  )+
  geom_hline(yintercept = mean(df$mexp), linetype = "dashed", color = "red", size = 1) 
