library("devtools")
devtools::install_github('YosefLab/SymSim')

library("SymSim")



## Simulate one population
## First, we simulate the case where there is one population.


ngenes <- 500
true_counts_res <- SimulateTrueCounts(ncells_total=300, ngenes=ngenes, evf_type="one.population", Sigma=0.4, randseed=0)
tsne_true_counts <- PlotTsne(meta=true_counts_res[[3]], data=log2(true_counts_res[[1]]+1), evf_type="one.population", n_pc=20, label='pop', saving = F, plotname="one.population")
tsne_true_counts[[2]]

tsne_true_counts

  

## Simulate multiple discrete populations
## When there are multiple populations, users need to provide a tree. 
## A tree with five leaves (five populations) can be generated as follows:

phyla1 <- Phyla5()
  
## or read from a file with the tree in Newick format:

phyla2 <- read.tree(system.file("extdata", "Newick_ABCDE.txt", package = "SymSim"))
par(mfrow=c(1,2))
plot(phyla1)
plot(phyla2)



true_counts_res <- SimulateTrueCounts(ncells_total=300, min_popsize=50, i_minpop=2, ngenes=ngenes, nevf=10, evf_type="discrete", n_de_evf=9, vary="s", Sigma=0.4, phyla=Phyla5(), randseed=0)
true_counts_res_dis <- true_counts_res
tsne_true_counts <- PlotTsne(meta=true_counts_res[[3]], data=log2(true_counts_res[[1]]+1), evf_type="discrete", n_pc=20, label='pop', saving = F, plotname="discrete populations (true counts)")
tsne_true_counts[[2]]


data(gene_len_pool)
gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], meta_cell=true_counts_res[[3]], protocol="nonUMI", alpha_mean=0.1, alpha_sd=0.05, gene_len=gene_len, depth_mean=1e5, depth_sd=3e3)
tsne_nonUMI_counts <- PlotTsne(meta=observed_counts[[2]], data=log2(observed_counts[[1]]+1), evf_type="discrete", n_pc=20, label='pop', saving = F, plotname="observed counts nonUMI")
tsne_nonUMI_counts[[2]]

observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], meta_cell=true_counts_res[[3]], protocol="UMI", alpha_mean=0.05, alpha_sd=0.02, gene_len=gene_len, depth_mean=5e4, depth_sd=3e3)
tsne_UMI_counts <- PlotTsne(meta=observed_counts[[2]], data=log2(observed_counts[[1]]+1), evf_type="discrete", n_pc=20, label='pop', saving = F, plotname="observed counts UMI")
tsne_UMI_counts[[2]]



## Add batch effects 
## We can divide the data we simulated using the previous steps into multiple batches and add batch effects to each batch. 

observed_counts_2batches <- DivideBatches(observed_counts_res = observed_counts, nbatch = 2, batch_effect_size = 1)
tsne_batches <- PlotTsne(meta=observed_counts_2batches[[2]], data=log2(observed_counts_2batches[[1]]+1), evf_type="discrete", n_pc=20, label='batch', saving = F, plotname="observed counts in batches")
tsne_batches[[2]]


## Simulate continuous populations
## We use the same tree as used for the simulation of discrete populations above. 
## To visualize the continuous populations, we color the cells by the edges they belong to on the tree. 
## We then label both the internal and tip nodes of the tree.

plot(phyla1, show.tip.label=F)
nodelabels()
tiplabels()



true_counts_res <- SimulateTrueCounts(ncells_total=500, ngenes=ngenes, nevf=20, evf_type="continuous", n_de_evf=12, vary="s", Sigma=0.4, phyla=Phyla5(), randseed=1)
tsne_true_counts <- PlotTsne(meta=true_counts_res[[3]], data=log2(true_counts_res[[1]]+1), evf_type="continuous", n_pc=20, label='pop', saving = F, plotname="continuous populations (true counts)")
tsne_true_counts[[2]]














library(ape)

# 生成一个随机树
set.seed(1234)
tree <- rtree(n = 6)  # 生成一个具有10个叶节点的随机树

# 可视化树
plot(tree, main = "Random Tree")




symsim_generate <- function(ncells_total = 2000, ngenes = 5000, ce = 0.1, Sigma = 0.4,
                            nevf_para = 0.015, n_de_evf_para = 0.01, phyla,
                            random_state){
  #### simulates single-cell RNA sequencing data with SymSim: https://www.nature.com/articles/s41467-019-10500-w

  true_counts_res <- SimulateTrueCounts(ncells_total=ncells_total, ngenes=ngenes, 
                                        nevf = nevf_para*ngenes, n_de_evf = n_de_evf_para*ngenes, 
                                        evf_type="continuous", vary="s", Sigma=Sigma, phyla= phyla, randseed=random_state)
  
  data(gene_len_pool)
  gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
  observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]],
                                         meta_cell=true_counts_res[[3]], protocol="UMI", 
                                         alpha_mean=ce,
                                         alpha_sd=0.02, gene_len=gene_len, depth_mean=5e5, depth_sd=3e4)
  
  rownames(observed_counts$counts) = paste0("gene_", c(1:5000))
  colnames(observed_counts$counts) = paste0("cell_", c(1:2000))
  observed_count_table = as.data.frame(observed_counts$counts)
  
  
  TrajInfo <- getTrajectoryGenes(observed_counts$cell_meta)
  branches = as.character(TrajInfo["branch"][,1])
  pseudotime = TrajInfo["pseudotime"]
  gene_pseudotime_cor <- matrix(NA, nrow = length(unique(branches)), ncol = ngenes)
  rownames(gene_pseudotime_cor) = unique(branches)
  colnames(gene_pseudotime_cor) = rownames(observed_counts$counts)
  
  
  
  
  for (gene in rownames(observed_count_table)) {
    for (branch in unique(branches)) {
      branch_cells <- which(branches == branch)
      gene_expr <- as.numeric(observed_count_table[gene, branch_cells])
      pseudotime_branch <- as.numeric(pseudotime[branch_cells,])
      
      
      if(length(unique(gene_expr)) > 1 && length(unique(pseudotime_branch )) > 1){
        gene_pseudotime_cor[branch, gene]  =  abs(cor(gene_expr, pseudotime_branch, method = "spearman"))
      } else {
        gene_pseudotime_cor[branch, gene] = NA # Avoid correlation calculation when there is no variance
      }
      
    }
    
    
  }
  
  max_gene_pseudotime_cor <- apply(gene_pseudotime_cor, 2, function(x){max(x, na.rm =T)})

  
  output = list()
  output$X = t(observed_count_table)
  output$groups = branches
  output$pseudotime =  as.numeric(pseudotime[,1])
  output$corr = max_gene_pseudotime_cor
  
  
  return(output)
}



for (i in 1:50) {
  
  set.seed(123*i)
  tree <- rtree(n = (i-1) %/%10 + 2 ) 
  
 result <- symsim_generate(random_state = 123*i, nevf_para = 0.015, n_de_evf_para = 0.01, Sigma = 0.4, phyla = tree)
 assign(paste0("Traj_data_", i), result)
 saveRDS( get(paste0("Traj_data_", i)), file = paste0("C:/Users/css22/Desktop/Thesis1/Traj/Traj_data_", i, ".RDS"))
 print(i)
}





for (i in 1:50) {
  raw_dta=readRDS(file = paste0("C:/Users/css22/Desktop/Thesis1/Traj/Traj_data_", i, ".RDS"))
  dta = t(raw_dta$X)
  result = mcadet(data = dta)
  assign(paste0("Mcadet_Traj_", i),  result)
  rm(list = paste0("Traj_data_",  i))
  gc()
  print(i)
}







######### Mcadet selected genes 

total_pseudo_cor_vector_Mcadetgenes = c()

for(i in 1:50){
  
  
  assign(paste0("Traj_data_", i), readRDS(file = paste0("C:/Users/css22/Desktop/Thesis1/Traj/Traj_data_", i, ".RDS") ))
  



gene_names <- get(paste0("Mcadet_Traj_", i))$gene
branches = get(paste0("Traj_data_", i))$groups
pseudotime = get(paste0("Traj_data_", i))$pseudotime
dta = t(get(paste0("Traj_data_", i))$X)
gene_pseudotime_cor <- matrix(NA, nrow = length(unique(branches)), ncol = length(gene_names))
rownames(gene_pseudotime_cor) = unique(branches)
colnames(gene_pseudotime_cor) = gene_names 





for (gene in gene_names) {
  for (branch in unique(branches)) {
    branch_cells <- which(branches == branch)
    gene_expr <- as.numeric(dta[gene, branch_cells])
    pseudotime_branch <- as.numeric(pseudotime[branch_cells])
    
    
    if(length(unique(gene_expr)) > 1 && length(unique(pseudotime_branch )) > 1){
      gene_pseudotime_cor[branch, gene]  =  abs(cor(gene_expr, pseudotime_branch, method = "spearman"))
    } else {
      gene_pseudotime_cor[branch, gene] = NA # Avoid correlation calculation when there is no variance
    }
    
  }
  
  
}


each_vector <- apply(gene_pseudotime_cor, 2, max)

total_pseudo_cor_vector_Mcadetgenes <- c(total_pseudo_cor_vector_Mcadetgenes, each_vector)

print(i)

}













######## All genes

total_pseudo_cor_vector_allgenes = c()

for(i in 1:50){


assign(paste0("Traj_data_", i), readRDS(file = paste0("C:/Users/css22/Desktop/Thesis1/Traj/Traj_data_", i, ".RDS") ))

  
branches = get(paste0("Traj_data_", i))$groups
pseudotime = get(paste0("Traj_data_", i))$pseudotime
dta = t(get(paste0("Traj_data_", i))$X)
gene_names = rownames(dta)

gene_pseudotime_cor <- matrix(NA, nrow = length(unique(branches)), ncol = length(gene_names))
rownames(gene_pseudotime_cor) = unique(branches)
colnames(gene_pseudotime_cor) = gene_names 



for (gene in gene_names) {
  for (branch in unique(branches)) {
    branch_cells <- which(branches == branch)
    gene_expr <- as.numeric(dta[gene, branch_cells])
    pseudotime_branch <- as.numeric(pseudotime[branch_cells])
    
    
    if(length(unique(gene_expr)) > 1 && length(unique(pseudotime_branch )) > 1){
      gene_pseudotime_cor[branch, gene]  =  abs(cor(gene_expr, pseudotime_branch, method = "spearman"))
    } else {
      gene_pseudotime_cor[branch, gene] = NA # Avoid correlation calculation when there is no variance
    }
    
  }
  
  
}


each_vector <- apply(gene_pseudotime_cor, 2, max)

total_pseudo_cor_vector_allgenes <- c(total_pseudo_cor_vector_allgenes, each_vector)

print(i)
}




hist(total_pseudo_cor_vector_Mcadetgenes, breaks = 100)
hist(total_pseudo_cor_vector_allgenes, breaks = 100)









library(ggplot2)
library(dplyr)


# total_pseudo_cor_vector_Mcadetgenes <- ...
# total_pseudo_cor_vector_allgenes <- ...


data <- data.frame(
  value = c(total_pseudo_cor_vector_Mcadetgenes, total_pseudo_cor_vector_allgenes),
  group = factor(rep(c("Mcadet Genes", "All Genes"), c(length(total_pseudo_cor_vector_Mcadetgenes), length(total_pseudo_cor_vector_allgenes))))
)


ggplot(data, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  labs(title = "",
       x = "Analytical Spearman Correlation",
       y = "Density") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12), 
    axis.text = element_text(face = "bold", size = 11), 
    axis.line = element_line(color = "black"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    legend.text = element_text(face = "bold", size = 10) 
  ) +
  guides(fill = guide_legend(title = NULL)) + 
  scale_fill_manual(values = c("Mcadet Genes" = "blue", "All Genes" = "red")) +
  geom_vline(xintercept = seq(from = 0, to = 1, by = 0.25), linetype = "dashed", color = "black", linewidth = 0.3) + 
  geom_hline(yintercept = seq(0, 3.5, by = 0.5), linetype = "dashed", color = "black", linewidth = 0.3) 






########## Other FS

for (i in 1:50) {
  dta = t(get(paste0("Traj_data_",  i))$X)
  BreFS = BrenneckeGetVariableGenes(dta, suppress.plot = T, fdr=0.1, minBiolDisp=0.5)
  assign(paste0("BreFS_Traj_",i), BreFS$Gene)
  print(i)
}


for (i in 1:50) {
  dta = t(get(paste0("Traj_data_",  i))$X)
  M3_data <- M3DropConvertData(dta, is.counts=TRUE)
  M3DropFS<- M3DropFeatureSelection(M3_data ,mt_method="fdr", mt_threshold=0.05)
  assign(paste0("M3Drop_Traj_",i), M3DropFS$Gene)
  print(i)
}



for (i in 1:50) {
  dta = t(get(paste0("Traj_data_",  i))$X)
  nb_raw <- NBumiConvertData(dta, is.counts=TRUE)
  DANB_fit <- NBumiFitModel(nb_raw)
  NBDropFS<- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thres=0.05, suppress.plot=F)
  assign(paste0("NBDrop_Traj_",i), NBDropFS$Gene)
  print(i)
}





for (i in 1:50) {
  dta = t(get(paste0("Traj_data_",  i))$X)
  nb_raw <- NBumiConvertData(dta, is.counts=TRUE)
  DANB_fit <- NBumiFitModel(nb_raw)
  NBDropFS<- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thres=0.05, suppress.plot=F)
  assign(paste0("NBDrop_Traj_",i), NBDropFS$Gene)
  print(i)
}



for (i in 1:50) {
  dta = t(get(paste0("Traj_data_",  i))$X)
  set.seed(12345*1+2222)
  rand_genes = sample(rownames(dta), 1000)
  assign(paste0("Random_Traj_",i), rand_genes)
  print(i)
}





for (i in 1:50) {
  dta = t(get(paste0("Traj_data_",  i))$X)
  scry=devianceFeatureSelection(as.matrix(dta), fam = "poisson")
  scryFS = names(sort(scry,decreasing = T)[1:1000])
  assign(paste0("scryFS_Traj_",i), scryFS)
  print(i)
}




for (i in 1:50) {
  dta = t(get(paste0("Traj_data_",  i))$X)
  seurat.obj<-CreateAssayObject(counts =  dta )
  seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj.disp <- FindVariableFeatures(seurat.obj, selection.method = "disp", nfeatures = 1000)
  SeuratDispFS<-VariableFeatures(seurat.obj.disp)
  assign(paste0("SeuratDisp_Traj_",i), SeuratDispFS)
  print(i)
}




for (i in 1:50) {
  dta = t(get(paste0("Traj_data_",  i))$X)
  seurat.obj<-CreateAssayObject(counts =  dta)
  seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj.mvp <- FindVariableFeatures(seurat.obj, selection.method = "mvp")
  SeuratMvpFS<-VariableFeatures(seurat.obj.mvp)
  assign(paste0("SeuratMvp_Traj_",i), SeuratMvpFS)
  print(i)
}




for (i in 1:50) {
  dta = t(get(paste0("Traj_data_",  i))$X)
  seurat.obj<-CreateAssayObject(counts = dta)
  seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj.vst <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 1000)
  SeuratVstFS<-VariableFeatures(seurat.obj.vst)
  assign(paste0("SeuratVst_Traj_",i), SeuratVstFS)
  print(i)
}





ASC_fun<-function(method_name){
total_pseudo_cor_vector_FSMethodgenes = c()
for(i in 1:50){
  
  
  assign(paste0("Traj_data_", i), readRDS(file = paste0("C:/Users/css22/Desktop/Thesis1/Traj/Traj_data_", i, ".RDS") ))
  
  
  
  
  gene_names <- get(paste0(method_name, "_Traj_", i))
  
  gene_names  <- gsub("-", "_", gene_names)
  
  branches = get(paste0("Traj_data_", i))$groups
  pseudotime = get(paste0("Traj_data_", i))$pseudotime
  dta = t(get(paste0("Traj_data_", i))$X)
  gene_pseudotime_cor <- matrix(NA, nrow = length(unique(branches)), ncol = length(gene_names))
  rownames(gene_pseudotime_cor) = unique(branches)
  colnames(gene_pseudotime_cor) = gene_names 
  
  
  
  
  
  for (gene in gene_names) {
    for (branch in unique(branches)) {
      branch_cells <- which(branches == branch)
      gene_expr <- as.numeric(dta[gene, branch_cells])
      pseudotime_branch <- as.numeric(pseudotime[branch_cells])
      
      
      if(length(unique(gene_expr)) > 1 && length(unique(pseudotime_branch )) > 1){
        gene_pseudotime_cor[branch, gene]  =  abs(cor(gene_expr, pseudotime_branch, method = "spearman"))
      } else {
        gene_pseudotime_cor[branch, gene] = NA # Avoid correlation calculation when there is no variance
      }
      
    }
    
    
  }
  
  
  each_vector <- apply(gene_pseudotime_cor, 2, max)
  
  total_pseudo_cor_vector_FSMethodgenes <- c(total_pseudo_cor_vector_FSMethodgenes, each_vector)
  
  print(i)
  
}
return(total_pseudo_cor_vector_FSMethodgenes)
}

total_pseudo_cor_vector_BreFSgenes = ASC_fun(method_name = "BreFS")
total_pseudo_cor_vector_M3Dropgenes = ASC_fun(method_name = "M3Drop")
total_pseudo_cor_vector_NBDropgenes = ASC_fun(method_name = "NBDrop")
total_pseudo_cor_vector_Randomgenes = ASC_fun(method_name = "Random")

total_pseudo_cor_vector_scryFSgenes = ASC_fun(method_name = "scryFS")
total_pseudo_cor_vector_SeuratDispgenes = ASC_fun(method_name = "SeuratDisp")
total_pseudo_cor_vector_SeuratMvpgenes = ASC_fun(method_name = "SeuratMvp")
total_pseudo_cor_vector_SeuratVstgenes = ASC_fun(method_name = "SeuratVst")




mean(total_pseudo_cor_vector_allgenes, na.rm = T)
mean(total_pseudo_cor_vector_Randomgenes, na.rm = T)
mean(total_pseudo_cor_vector_Mcadetgenes, na.rm = T)

mean(total_pseudo_cor_vector_BreFSgenes, na.rm = T)
mean(total_pseudo_cor_vector_M3Dropgenes, na.rm = T)
mean(total_pseudo_cor_vector_NBDropgenes, na.rm = T)
mean(total_pseudo_cor_vector_scryFSgenes, na.rm = T)
mean(total_pseudo_cor_vector_SeuratDispgenes, na.rm = T)
mean(total_pseudo_cor_vector_SeuratMvpgenes, na.rm = T)
mean(total_pseudo_cor_vector_SeuratVstgenes, na.rm = T)






df <- data.frame(y=c(0.252625,0.2525581, 0.3614572, 0.2039174, 0.1799792,
                     0.2308134, 0.3595574, 0.3451776, 0.2587468, 0.2644396),
                 x = factor(c("All Genes", "Random",
                              "Mcadet", "Brennecke", "M3Drop",
                              "NBDrop", "Scry", "Seurat Disp",
                              "Seurat Mvp", "Seurat Vst"), levels = c("Mcadet", "All Genes", "Random",
                                                                    "Brennecke", "M3Drop",
                                                                    "NBDrop", "Scry", "Seurat Disp",
                                                                    "Seurat Mvp", "Seurat Vst")),
                 sd = c(0.1581, 0.1581758, 0.1774935, 0.1407827,
                        0.1176975,  0.1470688, 0.1915522, 0.1960411,
                        0.1699013,0.1758486)
                 )




method_colors <- c("All Genes" = "#D73027", "Random" = "#4575B4", "Mcadet" = "#91BFDB", 
                   "Brennecke" = "#E0F3F8", "M3Drop" = "#FFFFBF", "NBDrop" = "#FDAE61", 
                   "Scry" = "#F46D43", "Seurat Disp" = "#A6D96A", "Seurat Mvp" = "#66BD63",
                   "Seurat Vst" = "#1A9850")


df$Significance = c("***", "***", "", "***", "***", "***", "NS", "***", "***", "***")
df$y_position = 0.5


ggplot(df, aes(x = x, y = y, fill = x)) +
  geom_bar(stat = "identity", width = 0.5, color = "black", alpha = 0.8, linewidth = 1) +
  geom_errorbar(aes(ymin = y - 0.25 * sd, ymax = y + 0.25 * sd), width = 0.2, color = "black", linewidth = 1) +
  labs(title = "",
       x = "",
       y = "Mean Analytical Spearman Correlation") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
    axis.title = element_text(face = "bold", size = 13), 
    axis.text = element_text(face = "bold", size = 12), 
    axis.line = element_line(color = "black"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    legend.title = element_blank()  
  ) + 
  geom_hline(yintercept = seq(0, 0.5, by = 0.1), linetype = "dashed", color = "black", size = 0.3) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_fill_manual(values = method_colors) +
  geom_text(aes(x = x, y = y_position, label = Significance), 
            size = 5, color = "black", fontface = "bold") 



t.test(total_pseudo_cor_vector_Mcadetgenes, total_pseudo_cor_vector_allgenes, alternative = "greater")
t.test(total_pseudo_cor_vector_Mcadetgenes, total_pseudo_cor_vector_Randomgenes, alternative = "greater")
t.test(total_pseudo_cor_vector_Mcadetgenes, total_pseudo_cor_vector_BreFSgenes, alternative = "greater")
t.test(total_pseudo_cor_vector_Mcadetgenes, total_pseudo_cor_vector_M3Dropgenes, alternative = "greater")
t.test(total_pseudo_cor_vector_Mcadetgenes, total_pseudo_cor_vector_NBDropgenes, alternative = "greater")
t.test(total_pseudo_cor_vector_Mcadetgenes, total_pseudo_cor_vector_scryFSgenes, alternative = "greater")
t.test(total_pseudo_cor_vector_Mcadetgenes, total_pseudo_cor_vector_SeuratDispgenes, alternative = "greater")
t.test(total_pseudo_cor_vector_Mcadetgenes, total_pseudo_cor_vector_SeuratMvpgenes, alternative = "greater")


