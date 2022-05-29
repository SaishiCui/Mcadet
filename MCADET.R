

### This function has 7 parameters in total:

### 1). "data" is the input data, a raw count matrix
###     genes on columns and cells on rows, the matrix should have 
###     row names which represent genes' name.

### 2). "n.comp" is the number of PCs selected from MCA decomposition, default is 100.

### 3). "k.list" is a list containing the number of k used to draw K-NN graph, 
###  default is 0.05*(ncols of your input data set), namely, 0.05* the number of the cells.

### 4). "run" is the number of iterations you want to run to get a more robust result, default is 10.

### 5). "step" should be specified when you want to do iterations, default is 2.

### 6). "seed" is the seed for random generator when perform leiden algorithm, the dault is NULL. 

### 7). "n.feature" is the number of genes to be selected, default is NA (for users who have domain knowledge ).

### 8). "p" is the type I error for feature selection, default is 0.05.

### 9). "fdr" is the false discovery rate for BH method to correct multiple testing issue, default is 0.15.

### 10). "clean" is the parameter when =True means clean the genes with mean gene expression less than 0.005.



### In order to run this function successfully, you need to install these packages: "irlba", "igraph", "cccd", 
### 'leidenalg', and "leiden"





mcadet<-function(data, n.comp=100, k.list=NULL, step=2, run=10, seed=NULL, n.feature=NA, p=0.05, clean.thresh=0.005,
                 fdr=0.15){
  
  
  mean.expression<-apply(data,1,mean)
  working.data=data[mean.expression>=clean.thresh,]      ##  keep genes with positive expression (delete empty and low expression genes)##
  

  
  
  gene.name<-rownames(working.data)               ##  extract gene's name in a list ## 
  
  if(is.null(k.list)){
    step=step
    k= round(ncol(working.data)*0.01)
    k.set = seq(k-step*round(run/2), k+step*round(run/2)-1, step) 
    
  }else{
    k.set=k.list
  }
  
  
  
  
  ## fuzzy coding step
  
  
  min.working=apply(working.data,1,min)         ## calculate the minimum expression of each gene
  
  max.working=apply(working.data,1,max)         ## calculate the maximum expression of each gene
  
  range=max.working-min.working                 ## calculate the difference between the max and min
  
  range.matrix=matrix(rep(range,dim(working.data)[2] ), ncol=dim(working.data)[2])    ## make the range as a matrix 
  
  min.matrix = matrix(rep(min.working, dim(working.data)[2]), nrow=dim(working.data)[1]) ## make the minimum as a matrix
  
  
  fuzzy.pos=as.matrix((working.data-min.working)/range.matrix)      ## matrix of gene +
  
  
  
  fuzzy.neg=1-fuzzy.pos                                             ## matrix of gene -
  
  
  Aug= matrix(data=0, nrow = 2*dim(fuzzy.pos)[1], ncol = dim(fuzzy.pos)[2])      ## define a empty augmented matrix (2p x n)
  
  
  Aug[seq(1,2*dim(working.data)[1],by=2),]=fuzzy.pos[1:dim(fuzzy.pos)[1],]            ## fill the odd rows in gene +
  Aug[seq(2,2*dim(working.data)[1],by=2),]=fuzzy.neg[1:dim(fuzzy.neg)[1],]            ## fill the even rows in gene -  
  
  
  
  Aug=t(Aug)    ## transpose matrix as an n x 2p matrix
  
  
  
  
  
  ### Construct corresponding matrix 
  
  P = 1/sum(Aug) * Aug
  
  row.mass = rowSums(P)  ## calculate row mass 
  col.mass = colSums(P)  ## calculate column mass
  
  
  
  ### Construct standardized Pearson residual matrix 
  rm(working.data)
  rm(Aug)
  rm(min.working)
  rm(max.working)
  rm(range)
  rm(range.matrix)
  rm(min.matrix)
  gc()
  
  expected.matrix = as.matrix(row.mass)%*%t(col.mass)
  sp.res.matrix =  (P-expected.matrix)/sqrt(expected.matrix)
  
  rm(expected.matrix)
  gc()
  
  ## irlba (fast and memory efficient SVD) ##
  
  
  library(irlba)
  
  irlba.decomp<-irlba(sp.res.matrix, n.comp)
  
  
  U_k = irlba.decomp$u[,c(1:n.comp)]                          ## extract svd U
  V_k = irlba.decomp$v[seq(1,dim(P)[2], by=2),c(1:n.comp)]    ## extract svd V (first n.comp PCs)
  
  Dc.matrix=matrix(rep(1/sqrt(col.mass[seq(1,dim(P)[2] ,by=2)]), n.comp), nrow = dim(P)[2]/2) ## diagonal matrix of column mass (gene +)
  
  
  D_alpha_K= diag(irlba.decomp$d[1:n.comp],n.comp)        ## diagonal matrix of singular values 
  
  
  
  standard.row.coord  =  sqrt(dim(P)[1])*U_k              ## calculate standard row coordinates (cell coordinates) 
  
  
  principal.col.coord =  Dc.matrix*V_k%*%D_alpha_K        ## calculate principal row coordinates (gene + coordinates) 
  
  
  
  
  
  
  ### KNN graph and leiden detection
  
  library("igraph")
  library("cccd")
  library("leiden")
  

  run=length(k.set)
  gene.logR.matrix<-matrix(data=0, nrow = run, ncol = nrow(principal.col.coord))  ## set up an empty matrix
  cluster.list<-NULL
  for(i in 1:run){
    
    knn.graph.cells=nng(standard.row.coord, mutual = F, k= k.set[i])   ## NNgraph algorithm based on Euclidean distance
    adjacency.matrix <- igraph::as_adjacency_matrix(knn.graph.cells)
    partition <- leiden(adjacency.matrix, seed=seed)   ## Leiden detection
    
    
    
    
    gene.rank.matrix<-matrix(NA, ncol = nrow(principal.col.coord), nrow =length(unique(partition)))
    colnames(gene.rank.matrix)<-gene.name
    cluster<-sort(unique(partition))
    cluster.list<-append(cluster.list,length(cluster))
    for (j in 1:length(unique(partition))) {
      gene.rank.matrix[j,]<-rank(sqrt(rowSums((t(matrix(rep(colSums(standard.row.coord[partition==j,])/n.comp,   ## Calculate cluster centroid and rank within each cluster
                                                            nrow(principal.col.coord)), ncol =nrow(principal.col.coord)))-principal.col.coord)^2)))
    
      }
    gene.logR.matrix[i,]<-apply(gene.rank.matrix,2,function(x){log(max(x)/min(x)) })   ##  a matrix with each gene's rank range for each run


  }
  
  colnames(gene.logR.matrix)<-gene.name
  gene.logR<-colMeans(gene.logR.matrix)     ##  mean log Maximum/Minimum for each gene over all runs ##
    


  
  if(is.na(n.feature)){
  
    
    logR_matrix<-matrix(NA, nrow = length(cluster.list), ncol = 50000)
    for (j in 1:length(cluster.list)) {
      logR_list<-rep(0,50000)
      for (i in 1:50000) {
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
      position_list<-append(position_list,i)}
      
    }
    fdr_thresh<-sorted_pval[length(position_list)]   # find L and corresponding p-value threshold in BH precedure#
    
    
    gene.list<-pval_df[pval_df$pval_list<=fdr_thresh,"gene.name"]
    
    gene.logR.list<-gene.logR[pval_df$pval_list<=fdr_thresh]
    
    p.val.list<-pval_df[pval_df$pval_list<=fdr_thresh,"pval_list"]
    
    df<-data.frame(gene.list, gene.logR.list, p.val.list)  
    colnames(df)<-c("gene", "logR", "p.value")  
    df<-df[order(df$p.value, decreasing = F),]
    
    obj<-df  
    
    
  }else{
   DE_label<-order(gene.logR,decreasing = T)[1:n.feature]   # if users know how many genes they want #
   gene.list<-gene.name[DE_label]
   
   p.val.list<-NULL
   for (i in 1:length(gene.name)) {
     p.val=mean(logR_hist>gene.logR[i]) 
     p.val.list<-append(p.val.list,p.val)
     }
   
   df<-data.frame(gene.list,gene.logR.list[DE_label], p.val.list)  
   colnames(df)<-c("gene", "logR", "p.value")  
   df<-df[order(df$p.value, decreasing = F),]
   obj<-df  
  }
  
  
  
  
  return(obj)
  
}


















