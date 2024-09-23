#' Mcadet: Feature Selection for scRNA-seq Data
#'
#' This function performs feature selection for scRNA-seq (count-based) data. It integrates Multiple Correspondence Analysis (MCA), graph-based community detection, and a novel statistical testing approach to select the most significant features (genes).
#'
#' @param data A raw count matrix where rows represent cells and columns represent genes.
#' @param n.comp The number of principal components (PCs) selected from MCA decomposition. Default is 60.
#' @param run The number of iterations to run for a more robust result. Default is 10.
#' @param n.feature The number of genes to select. Default is NA, which means the proposed statistical testing method will be used.
#' @param start_resolution The starting resolution parameter for the Leiden algorithm, increasing by 0.1 for each run. Default is 0.5.
#' @param cell_percent The percentage threshold for filtering low-expressed genes. Genes expressed in fewer cells than this threshold will be excluded. Default is 0.005.
#' @param MC_iter Monte-Carlo iteration times. Default is 50,000.
#' @param fdr The false discovery rate (FDR) threshold for the statistical testing. Default is 0.15.
#' @param seed The random seed for the Leiden algorithm. Default is 1234.
#'
#' @return A list of selected genes, or a data frame with p-values and log ratios if `n.feature` is set to NA.
#'
#' @examples
#' \dontrun{
#' data <- matrix(rnbinom(1000,2,0.5), nrow = 100, ncol = 10)
#' selected_genes <- mcadet(data)
#' }
#'
#' @export


mcadet<-function(data, n.comp=60,  run=10, n.feature=NA, nk_percent = 0.005, start_resolution = 0.5,
                 cell_percent=0.005, MC_iter = 50000, fdr = 0.15, seed = 1234){
  
  ## loading required package
  library("Seurat")
  library("cccd")
  library("igraph")
  library("irlba")
  
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
  

  library("cccd")
  library("igraph")
  library("RANN")
  
  gene.logR.matrix<-matrix(data=0, nrow = run, ncol = nrow(principal.col.coord))  ## set up an empty matrix
  cluster.list<- c()
  n_cells = nrow(standard.row.coord)
  nn_k = round(n_cells*nk_percent)
  if(nn_k <4){
    nn_k  = 4
  }
  weights_matrix = matrix(0, nrow = n_cells, ncol = nn_k)
  knn_graph_matrix <- nn2(data = standard.row.coord, k = nn_k + 1)$nn.idx[,-1] # +1 to include the point itself
  
  
  for (q in 1:n_cells) {
    for (z in 1:nn_k) {
      neighbor = knn_graph_matrix[q,z]
      weights_matrix[q,z] =jaccard(knn_graph_matrix[q,], knn_graph_matrix[neighbor,])
    }
  }

  
  
  # Initialize an empty list to store the edges
  edges <- c()
  
  # Fill the edges list based on the KNN index matrix
  for (i in 1:nrow(knn_graph_matrix)) {
    neighbors <- knn_graph_matrix[i, ] # Exclude the point itself (first column)
    edges <-c(edges, c(rbind(i, neighbors)))
  }
  
  # Create the graph object
  knn_graph <- graph(edges = edges, directed = FALSE)
  

  
  for(i in 1:run){
    
    set.seed(seed*i)
    partition = cluster_leiden(graph = knn_graph, 
                               objective_function = c("modularity"),
                               weights  = weights_matrix, 
                               resolution_parameter = start_resolution + 0.1*i)$membership ## Leiden detection
 
    
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





















