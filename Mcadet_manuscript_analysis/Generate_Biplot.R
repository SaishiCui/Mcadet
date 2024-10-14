
rm(list = ls())
gc()


load(file = paste0("C:/Users/css22/Desktop/Thesis1/SPARSim/SparSim_Mcadet_data/SparSim_", 12, ".RData") )
load(file = paste0("C:/Users/css22/Desktop/Thesis1/pbmc/data/pbmccell1_1_0", ".RData") )
load(file = paste0("C:/Users/css22/Desktop/Thesis1/pbmc/data/pbmcdata_1_0", ".RData") )
load(file = paste0("C:/Users/css22/Desktop/Thesis1/pbmc/DE_genes/DEgene_1_0", ".RData") )



data = pbmcdata_1_0
n.comp=60 
run=10
n.feature=NA
nk_percent = 0.005
start_resolution = 0.5
cell_percent=0.005 
MC_iter = 50000
fdr = 0.15
seed = 1234
  
  
clean_fucntion = function(x){
    percent = mean(1*(x!=0 ))
  }
  
percent_vector=apply(data,1,clean_fucntion)            ##  cleaning genes with low expressed cells##
working.data = data[percent_vector>= cell_percent,]
gene.name<-rownames(working.data)                       ##  extract gene's name in a list ##
  
dim(working.data)
  
## fuzzy coding step
  
fuzzy_func = function(x){
    fuzzy_trans  = (x-min(x))/(max(x) - min(x))
    return(fuzzy_trans)
  }
  
  
fuzzy.pos=apply(working.data,1,fuzzy_func)
  
fuzzy.neg=1-fuzzy.pos                                             ## matrix of gene -
  
  
Aug= matrix(data=0, nrow = dim(fuzzy.pos)[1], ncol = 2*dim(fuzzy.pos)[2])      ## define a empty augmented matrix (n x 2p)
  
  
Aug[,seq(1,2*dim(working.data)[1],by=2)]=fuzzy.pos[,1:dim(fuzzy.pos)[2]]            ## fill the odd rows in gene +
Aug[,seq(2,2*dim(working.data)[1],by=2)]=fuzzy.neg[,1:dim(fuzzy.neg)[2]]            ## fill the even rows in gene -  
  
  
  
  ### Construct corresponding matrix 
  
P = 1/sum(Aug) * Aug
  
row.mass = rowSums(P)  ## calculate row mass 
col.mass = colSums(P)  ## calculate column mass
  
  
  
  ### Construct standardized Pearson residual matrix 
  
  
  
expected.matrix = as.matrix(row.mass)%*%t(col.mass)
sp.res.matrix =  (P-expected.matrix)/sqrt(expected.matrix)
  
  
rm(working.data)
rm(Aug)
rm(expected.matrix)
gc()
  
  
  
  ## irlba (fast and memory efficient SVD) ##
  
library(irlba)
  
irlba.decomp<-irlba(sp.res.matrix, n.comp)
  
  
U_k = irlba.decomp$u[,c(1:n.comp)]                          ## extract svd U
V_k = irlba.decomp$v[seq(1,dim(P)[2], by=2),c(1:n.comp)]    ## extract svd V (first n.comp PCs)
  
  
D_alpha_K= diag(irlba.decomp$d[1:n.comp],n.comp)        ## diagonal matrix of singular values 
  
standard.row.coord  =  sqrt(dim(P)[1])*U_k              ## calculate standard row coordinates (cell coordinates) 
  
principal.col.coord =  (1/sqrt(col.mass[seq(1,dim(P)[2], by=2)]))*V_k%*%D_alpha_K ## calculate principal col coordinates (gene + coordinates) 


dim(standard.row.coord)
dim(principal.col.coord)


### Caculate distance 

standard.row.coord[which(pbmccell1_1_0 == "B"), 1:60]
B_centroid <- apply(standard.row.coord[which(pbmccell1_1_0 == "B"), 1:60], 2, mean)
CD4_centroid <- apply(standard.row.coord[which(pbmccell1_1_0 == "CD4 T"), 1:60], 2, mean)
CD8_centroid <- apply(standard.row.coord[which(pbmccell1_1_0 == "CD8 T"), 1:60], 2, mean)
NK_centroid <- apply(standard.row.coord[which(pbmccell1_1_0 == "NK"), 1:60], 2, mean)
Mono_centroid <- apply(standard.row.coord[which(pbmccell1_1_0 == "Mono"), 1:60], 2, mean)
DC_centroid <- apply(standard.row.coord[which(pbmccell1_1_0 == "DC"), 1:60], 2, mean)


library(ggplot2)
B_labels = which(pbmccell1_1_0 == "B")[1:240]
CD4T_labels = which(pbmccell1_1_0 == "CD4 T")[1:360]
CD8T_labels = which(pbmccell1_1_0 == "CD8 T")[1:360]
NK_labels = which(pbmccell1_1_0 == "NK")[1:360]
Mono_labels = which(pbmccell1_1_0 == "Mono")[1:360]
DC_labels = which(pbmccell1_1_0 == "DC")[1:120]


CD4T_label = which(gene.name == "IL7R")
CD8T_label = which(gene.name == "CD8A")
Mono_label = which(gene.name == "CD68")
NK_label = which(gene.name == "NCR1")
DC_label   = which(gene.name == "CLEC9A")
B_label = which(gene.name == "MS4A1")

dist_df <- data.frame(cell_type = rep(c("B", "CD4 T", "CD8 T", "NK", "Mono", "DC"), 6),
           marker_gene = rep(c("MS4A1 (Marker gene for B cells)",
                               "IL7R (Marker gene for CD4 T cells)",
                               "CCL5 (Marker gene for CD8 T cells)",
                               "NCR1 (Marker gene for NK cells)",
                               "CD68 (Marker gene for Mono cells)",
                               "CLEC9A (Marker gene for DC cells)"),each = 6),
           
           dis = c(sqrt(sum((principal.col.coord[B_label, 1:60] - B_centroid)^2)),
                   sqrt(sum((principal.col.coord[B_label, 1:60] - CD4_centroid)^2)),
                   sqrt(sum((principal.col.coord[B_label, 1:60] - CD8_centroid)^2)),
                   sqrt(sum((principal.col.coord[B_label, 1:60] - NK_centroid)^2)),
                   sqrt(sum((principal.col.coord[B_label, 1:60] - Mono_centroid)^2)),
                   sqrt(sum((principal.col.coord[B_label, 1:60] - DC_centroid)^2)),
                   
                   sqrt(sum((principal.col.coord[CD4T_label, 1:60] - B_centroid)^2)),
                   sqrt(sum((principal.col.coord[CD4T_label, 1:60] - CD4_centroid)^2)),
                   sqrt(sum((principal.col.coord[CD4T_label, 1:60] - CD8_centroid)^2)),
                   sqrt(sum((principal.col.coord[CD4T_label, 1:60] - NK_centroid)^2)),
                   sqrt(sum((principal.col.coord[CD4T_label, 1:60] - Mono_centroid)^2)),
                   sqrt(sum((principal.col.coord[CD4T_label, 1:60] - DC_centroid)^2)),
                   
                   sqrt(sum((principal.col.coord[CD8T_label, 1:60] - B_centroid)^2)),
                   sqrt(sum((principal.col.coord[CD8T_label, 1:60] - CD4_centroid)^2)),
                   sqrt(sum((principal.col.coord[CD8T_label, 1:60] - CD8_centroid)^2)),
                   sqrt(sum((principal.col.coord[CD8T_label, 1:60] - NK_centroid)^2)),
                   sqrt(sum((principal.col.coord[CD8T_label, 1:60] - Mono_centroid)^2)),
                   sqrt(sum((principal.col.coord[CD8T_label, 1:60] - DC_centroid)^2)),
                   
                   
                   sqrt(sum((principal.col.coord[NK_label, 1:60] - B_centroid)^2)),
                   sqrt(sum((principal.col.coord[NK_label, 1:60] - CD4_centroid)^2)),
                   sqrt(sum((principal.col.coord[NK_label, 1:60] - CD8_centroid)^2)),
                   sqrt(sum((principal.col.coord[NK_label, 1:60] - NK_centroid)^2)),
                   sqrt(sum((principal.col.coord[NK_label, 1:60] - Mono_centroid)^2)),
                   sqrt(sum((principal.col.coord[NK_label, 1:60] - DC_centroid)^2)),
                   
                   sqrt(sum((principal.col.coord[Mono_label, 1:60] - B_centroid)^2)),
                   sqrt(sum((principal.col.coord[Mono_label, 1:60] - CD4_centroid)^2)),
                   sqrt(sum((principal.col.coord[Mono_label, 1:60] - CD8_centroid)^2)),
                   sqrt(sum((principal.col.coord[Mono_label, 1:60] - NK_centroid)^2)),
                   sqrt(sum((principal.col.coord[Mono_label, 1:60] - Mono_centroid)^2)),
                   sqrt(sum((principal.col.coord[Mono_label, 1:60] - DC_centroid)^2)),
                   
                   sqrt(sum((principal.col.coord[DC_label, 1:60] - B_centroid)^2)),
                   sqrt(sum((principal.col.coord[DC_label, 1:60] - CD4_centroid)^2)),
                   sqrt(sum((principal.col.coord[DC_label, 1:60] - CD8_centroid)^2)),
                   sqrt(sum((principal.col.coord[DC_label, 1:60] - NK_centroid)^2)),
                   sqrt(sum((principal.col.coord[DC_label, 1:60] - Mono_centroid)^2)),
                   sqrt(sum((principal.col.coord[DC_label, 1:60] - DC_centroid)^2)))
           )



ggplot(dist_df, aes(x = cell_type, y = dis)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", fill = "skyblue", linewidth = 1) +
  facet_wrap(~ marker_gene, scales = "free_y") +
  labs(x = "", y = "Euclidean Distance (top 10 PCs)", title = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black"),
    axis.title.y = element_text(face = "bold", size = 14, color = "black"),
    axis.text.x = element_text(face = "bold", size = 12, color = "black"),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    strip.text = element_text(face = "bold", size = 12, color = "black")
  ) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Euclidean Distance (top 60 PCs)")



######### PBMC data #########





df_cells = data.frame( standard.row.coord[c( CD4T_labels, CD8T_labels, NK_labels, Mono_labels, B_labels, DC_labels), 1:2],
                       group = c(rep("CD4 T", 360), rep("CD8 T", 360), rep("NK", 360), rep("Mono", 360), rep("B", 240),
                                 rep("DC", 120)))


df_genes =data.frame(principal.col.coord[c(CD4T_label, Mono_label, NK_label, DC_label, B_label, CD8T_label),c(1,2)])

arrow_data <- data.frame(
  x_start = rep(df_genes$X1[4], each = 5),
  y_start = rep(df_genes$X2[4], each = 5),
  x_end = c(df_genes$X1[4]-0.3, df_genes$X1[4]-0.5, df_genes$X1[4]+0.7, df_genes$X1[4]+1.3, df_genes$X1[4]+0.85),
  y_end = c(df_genes$X2[4]-0, df_genes$X2[4]+2, df_genes$X2[4]+2, df_genes$X2[4]+2.3, df_genes$X2[4]+1.5)
)


ggplot(df_cells, aes(x = X1, y = X2, color = group)) +
  geom_point(size = 2.5) +
  labs(title = "",
       x = "PC 1",
       y = "PC 2",
       color = "group") +
  theme_minimal() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth  = 1.2),
        axis.title = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 12))+
  annotate("point", x = df_genes$X1[1], y = df_genes$X2[1], color = "black", shape = 24, size = 5, fill = "black") +  
  annotate("point", x = df_genes$X1[2], y = df_genes$X2[2], color = "black", shape = 23, size = 5, fill = "black" ) + 
  annotate("point", x = df_genes$X1[3], y = df_genes$X2[3], color = "black", shape = 22, size = 5, fill = "black" ) +
  annotate("point", x = df_genes$X1[4], y = df_genes$X2[4], color = "black", shape = 21, size = 5, fill = "black" ) +
  annotate("point", x = df_genes$X1[5], y = df_genes$X2[5], color = "black", shape = 25, size = 5, fill = "black" ) +
  annotate("point", x = df_genes$X1[6], y = df_genes$X2[6], color = "black", shape = 24, size = 5, fill = "blue" ) +
  
annotate("text", x = df_genes$X1[1], y = df_genes$X2[1], 
         label = "IL7R", vjust = -1, hjust = 0.2, size = 6, fontface = "bold")+
  annotate("text", x = df_genes$X1[2], y = df_genes$X2[2], 
           label = "CD68", vjust = -1, hjust = 0.5, size = 6, fontface = "bold")+
  annotate("text", x = df_genes$X1[3], y = df_genes$X2[3], 
           label = "NCR1", vjust = -1.0, hjust = 0.5, size = 6, fontface = "bold")+
  annotate("text", x = df_genes$X1[4], y = df_genes$X2[4], 
           label = "CLEC9A", vjust = 1.5, hjust = 0.5, size = 6, fontface = "bold")+
  annotate("text", x = df_genes$X1[5], y = df_genes$X2[5], 
           label = "MS4A1", vjust = -0.5, hjust = 0.5, size = 6, fontface = "bold")+
  annotate("text", x = df_genes$X1[6], y = df_genes$X2[6], 
           label = "CCL5", vjust = 1.5,  hjust = 0.6, size = 6, fontface = "bold")+
  
  
  geom_segment(data = arrow_data, aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth =1)














#######################


library(ggplot2)
df_cells = data.frame( standard.row.coord[1:480, 1:2], group = c(rep("Cell type 1", 60), rep("Cell type 2", 180), rep("Cell type 3", 240)) )
df_genes =data.frame(principal.col.coord[c(21,752,996),c(1,2)])
arrow_data <- data.frame(
  x_start = rep(df_genes$X1[1], each = 3),
  y_start = rep(df_genes$X2[1], each = 3),
  x_end = c(df_genes$X1[1]+0.1, df_genes$X1[1]-0, df_genes$X1[1]-0.2),
  y_end = c(df_genes$X2[1]+0.5, df_genes$X2[1]-0.4, df_genes$X2[1]+0)
)


ggplot(df_cells, aes(x = X1, y = X2, color = group)) +
  geom_point(size = 3) +
  labs(title = "",
       x = "PC 1",
       y = "PC 2",
       color = "group") +
  theme_minimal() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth  = 1.2),
        axis.title = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text = element_text(face = "bold", size = 12),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 12))+
  annotate("point", x = df_genes$X1[1], y = df_genes$X2[1], color = "black", shape = 24, size = 5, fill = "black") +  
  annotate("point", x = df_genes$X1[2], y = df_genes$X2[2], color = "black", shape = 23, size = 5, fill = "black" ) + 
  annotate("point", x = df_genes$X1[3], y = df_genes$X2[3], color = "black", shape = 22, size = 5, fill = "black" ) + 
  
 annotate("text", x = df_genes$X1[1], y = df_genes$X2[1], 
           label = "HVG of cell type 2", vjust = -2, hjust = 0.5, size = 4, fontface = "bold")+
  annotate("text", x = df_genes$X1[2], y = df_genes$X2[2], 
           label = "HVG of cell type 3", vjust = -2, hjust = 0.5, size = 4, fontface = "bold")+
  annotate("text", x = df_genes$X1[3], y = df_genes$X2[3], 
           label = "HVG of cell type 1", vjust = -2, hjust = 0.5, size = 4, fontface = "bold")+
  geom_segment(data = arrow_data, aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.2, "cm")), color = "red", linewidth =1)


###
## Gene 1801
## Gene 31
## Gene 1532


plot(x = standard.row.coord[1:60,1],
     y = standard.row.coord[1:60,2])

  
  
