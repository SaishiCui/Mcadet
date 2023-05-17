
## loading required package
library("Seurat")
library("cccd")
library("leidenAlg")
library("irlba")




### This function has 7 parameters in total:

### 1). "data" is the input data, a raw count matrix
###     genes on columns and cells on rows, the matrix should have 
###     row names which represent genes' name.

### 2). "n.comp" is the number of PCs selected from MCA decomposition, default is 60.

### 3). "run" is the number of iterations you want to run to get a more robust result, default is 10.

### 4). "n.feature" is the number of genes to select, default is NA, which means the proposed statistical testing method will be used

### 5). "start_resolution" is the starting resolution parameter for Leiden algorithm, increasing by 0.1*run 

### 6). "cell_percent": If a gene is only expressed in "cell_percent" (i.e., 0.5%) of cells, we 
## preclude these noninformative genes. Default is 0.005

## 7). MC_iter: Monte-Carlo iteration times, default is 50,000

## 8). fdr: False discovery rate for the proposed statistical testing rate

## 9). seed: random seed for Leiden algorithm




### In order to run this function successfully, you need to install these packages: "irlba", "igraph", "cccd", 
### 'leidenalg', and "leiden"




########### function ###########


mcadet<-function(data, n.comp=60,  run=10, n.feature=NA, nk_percent = 0.005, start_resolution = 0.5,
                 cell_percent=0.005, MC_iter = 50000, fdr = 0.15, seed = NULL){
  
  
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
  
  
  
  
  
  ### KNN graph and leiden detection
  jaccard = function(a,b){
    numerator = length(intersect(a,b))
    denominator = length(a) + length(b) - numerator
    return( numerator/denominator)
  }
  
  shrunk = function(x){
    which(x==1)
  }
  
  library("leidenAlg")
  library("cccd")
  library("igraph")
  
  
  gene.logR.matrix<-matrix(data=0, nrow = run, ncol = nrow(principal.col.coord))  ## set up an empty matrix
  cluster.list<- c()
  n_cells = nrow(standard.row.coord)
  nn_k = round(n_cells*nk_percent)
  if(nn_k <4){
    nn_k  = 4
  }
  weights_matrix = matrix(0, nrow = n_cells, ncol = nn_k)
  knn.graph.cells=nng(standard.row.coord, mutual = F, k= nn_k, algorithm = "kd_tree")
  adj_matrix = igraph::as_adjacency_matrix(knn.graph.cells)
  
  knn_graph=t(apply(adj_matrix,1,shrunk))
  
  
  for (q in 1:n_cells) {
    for (z in 1:nn_k) {
      neighbor = knn_graph[q,z]
      weights_matrix[q,z] =jaccard(knn_graph[q,],knn_graph[neighbor,])
    }
  }
  weights_vector = as.vector(t(weights_matrix))
  
  
  
  
  for(i in 1:run){
    
    set.seed(seed*i)
    partition = find_partition(knn.graph.cells, edge_weights = weights_vector, resolution = start_resolution + 0.1*i, niter = 2L)+1 ## Leiden detection
    
    gene.rank.matrix<-matrix(NA, ncol = nrow(principal.col.coord), nrow =length(unique(partition)))
    colnames(gene.rank.matrix)<-gene.name
    cluster<-sort(unique(partition))
    cluster.list<-append(cluster.list,length(cluster))
    for (j in 1:length(unique(partition))) {
      
      if(sum(1*(partition==j)) != 1){
        
        gene.rank.matrix[j,]<-rank(sqrt(rowSums((t(matrix( rep(colMeans(standard.row.coord[partition==j,]),   ## Calculate cluster centroid and rank within each cluster
                                                               nrow(principal.col.coord)), ncol =nrow(principal.col.coord)))-principal.col.coord)^2)))
      }else{
        gene.rank.matrix[j,]<-rank(sqrt(rowSums((t(matrix( rep(standard.row.coord[partition==j,],   ## Calculate cluster centroid and rank within each cluster
                                                               nrow(principal.col.coord)), ncol =nrow(principal.col.coord)))-principal.col.coord)^2)))
      }
      
    }
    gene.logR.matrix[i,]<-apply(gene.rank.matrix,2,function(x){log(max(x)/min(x)) })   ##  a matrix with each gene's rank range for each run
    
    
  }
  
  colnames(gene.logR.matrix)<-gene.name
  gene.logR<-colMeans(gene.logR.matrix)     ##  mean log Maximum/Minimum for each gene over all runs ##
  
  
  
  
  if(is.na(n.feature)){
    
    
    logR_matrix<-matrix(NA, nrow = length(cluster.list), ncol = MC_iter)
    for (j in 1:length(cluster.list)) {
      logR_list<-rep(0,MC_iter)
      for (i in 1:MC_iter) {
        sam<-sample( c(1:length(gene.name)), size = cluster.list[j], replace = T)
        logR<- log(max(sam)/min(sam) ) 
        logR_list[i]=logR
      }
      logR_matrix[j,]=logR_list       
    }
    
    logR_hist<-colMeans(logR_matrix)   ## Monte Carlo simulation for log Maximum/Minimum distribution (multiple runs) ##
    
    
    
    ### FDR BH method
    
    pval_list<-rep(NA,length(gene.name))
    for (i in 1:length(gene.name)) {
      p.val=mean(logR_hist>gene.logR[i])    ## upper-tailed p-values calculation ##
      pval_list[i]=p.val
    }
    
    pval_df<-data.frame(gene.name,pval_list)
    sorted_pval<- sort(pval_list)
    sorted_pval<-sorted_pval[sorted_pval<0.9]    # remove last one or two p-value histogram blocks,
    # to make p-values be uniformly distributed ##
    
    
    
    
    position_list<-NULL
    for (i in 1:length(sorted_pval)) {
      if(sorted_pval[i]< fdr*i/length(sorted_pval)){
        position_list<-append(position_list,i)}else{
          break
        }
      
    }
    fdr_thresh<-sorted_pval[length(position_list)]   # find L and corresponding p-value threshold in BH precedure#
    
    
    gene.list<-pval_df[pval_df$pval_list<=fdr_thresh,"gene.name"]
    
    gene.logR.list<-gene.logR[pval_df$pval_list<=fdr_thresh]
    
    p.val.list<-pval_df[pval_df$pval_list<=fdr_thresh,"pval_list"]
    
    df<-data.frame(gene.list, gene.logR.list, p.val.list)  
    colnames(df)<-c("gene", "logR", "p.value")  
    df<-df[order(df$p.value, decreasing = F),]
    
    obj<-df  
    return(obj)
    
  }else{
    DE_label<-order(gene.logR,decreasing = T)[1:n.feature]   # if users know how many genes they want #
    gene.list<-gene.name[DE_label]
    return(gene.list)

  }
  
  
}


