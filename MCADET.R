rm(list = ls())

library(reticulate)
py_available()
py_module_available("numpy")
py_config()


library("Seurat")
library("uwot")
library("igraph")
library("cccd")
library("leiden")




# Setup the Seurat Object


# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "C:/Users/cuisa/Desktop/BST_PhD/research/Dissertation/pbmc3k/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)




## QC and selecting cells for further analysis ##

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


## We filter cells that have unique feature counts over 2,500 or less than 200 and 
## We filter cells that have >5% mitochondrial counts



pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

raw_data=as.data.frame(pbmc[['RNA']]@data)





mca.list.comp200<-mcadet(raw_data, n.comp = 200, run = 10, step = 2,seed = 1)
mca.list.comp150<-mcadet(data=raw_data, n.comp = 150, run = 10, step = 2, seed = 1)
mca.list.comp100<-mcadet(raw_data, n.comp = 100, run = 10, step = 2, seed = 1)
mca.list.comp70<-mcadet(raw_data, n.comp = 70, run =10, step = 2, seed = 1)
mca.list.comp50<-mcadet(raw_data, n.comp = 50, run = 10, step = 2, seed = 1)
mca.list.comp30<-mcadet(raw_data, n.comp = 30, run = 10, step = 2, seed = 1)



### This function has 7 parameters in total:

### 1). "data" is the input data, a raw count matrix
###     genes on columns and cells on rows, the matrix should have 
###     row names which represent genes' name.

### 2). "n.comp" is the number of PCs selected from MCA decomposition, default is 100.

### 3). "k" is number of K used to draw K-NN graph, default is 5.

### 4). "run" is the number of iterations you want to run to get a more robust result, default is 1.

### 5). "step" should be specified when you want to do iterations, if your k=5, 
###      run = 10, and step = 5, it means that you do 10 iterations on k={5,10,15,...50}, respectively.
###      default is 0.

### 6). "n.feature.1" is the number of genes to select in the first round, default is 3,000.

### 7). "n.feature.2" is the number of genes to select in the second round, default is 2,000.


### In order to run this function successfully, you need to install these packages: "irlba", "igraph", "cccd", 
### 'leidenalg', and "leiden"

load(file = "C:/Users/cuisa/Desktop/BST_PhD/research/Dissertation/MCA/sim.data.0.08.RData")
load(file = "C:/Users/cuisa/Desktop/BST_PhD/research/Dissertation/MCA/sim.data.0.05.RData")
load(file = "C:/Users/cuisa/Desktop/BST_PhD/research/Dissertation/MCA/sim.data.0.03.RData")
load(file = "C:/Users/cuisa/Desktop/BST_PhD/research/Dissertation/MCA/sim.data.0.3.RData")
load(file = "C:/Users/cuisa/Desktop/BST_PhD/research/Dissertation/MCA/sim.data.0.2.RData")
load(file = "C:/Users/cuisa/Desktop/BST_PhD/research/Dissertation/MCA/sim.data.0.1.RData")




mca.sim.0.2<-mcadet(data=sim.data.0.2, n.comp=100, k=NA, run=10, step=2, seed=1, n.feature.1=3000, 
                 n.feature.2=2000, clean=TRUE)
mca.sim.0.1<-mcadet(data=sim.data.0.1, n.comp=100, k=NA, run=10, step=2, seed=1, n.feature.1=3000, 
                    n.feature.2=2000, clean=TRUE)
mca.sim.0.3<-mcadet(data=sim.data.0.3, n.comp=100, k=NA, run=10, step=2, seed=1, n.feature.1=3000, 
                    n.feature.2=2000, clean=TRUE)



mca.sim.0.08<-mcadet(data=sim.data.0.08, n.comp=100, k=NA, run=10, step=2, seed=1, n.feature.1=1500, 
                    n.feature.2=1000, clean=TRUE)
mca.sim.0.05<-mcadet(data=sim.data.0.05, n.comp=100, k=NA, run=10, step=2, seed=1, n.feature.1=1500, 
                    n.feature.2=1000, clean=TRUE)
mca.sim.0.03<-mcadet(data=sim.data.0.03, n.comp=100, k=NA, run=10, step=2, seed=1, n.feature.1=1500, 
                    n.feature.2=1000, clean=TRUE)




save(mca.sim.0.2,file = "C:/Users/cuisa/Desktop/BST_PhD/research/Dissertation/MCA/mca.sim.0.2.RData")
save(mca.sim.0.1,file = "C:/Users/cuisa/Desktop/BST_PhD/research/Dissertation/MCA/mca.sim.0.1.RData")
save(mca.sim.0.3,file = "C:/Users/cuisa/Desktop/BST_PhD/research/Dissertation/MCA/mca.sim.0.3.RData")
save(mca.sim.0.08,file = "C:/Users/cuisa/Desktop/BST_PhD/research/Dissertation/MCA/mca.sim.0.08.RData")
save(mca.sim.0.05,file = "C:/Users/cuisa/Desktop/BST_PhD/research/Dissertation/MCA/mca.sim.0.05.RData")
save(mca.sim.0.03,file = "C:/Users/cuisa/Desktop/BST_PhD/research/Dissertation/MCA/mca.sim.0.03.RData")



mcadet<-function(data, n.comp=100, k=NA, run=10, step=2, seed=NULL, n.feature.1=3000, n.feature.2=2000,
                 clean=TRUE){



working.data = data[rowSums(data)!=0,]          ##  keep genes with positive expression (delete empty genes)##

if(clean){
mean.expression<-apply(working.data,1,mean)
working.data=working.data[mean.expression>=0.005,]
}
gene.name<-rownames(working.data)               ##  extract gene's name in a list ## 

if(is.na(k)){
  k=round(ncol(working.data)*0.01)
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


k.set = seq(k-step*round(run/2), k+step*(round(run/2)-1), step)             ## set up a list of k for knn (if run = 1, then this list only cotains one value)
gene.range.matrix<-matrix(data=0, nrow = run, ncol = nrow(principal.col.coord))  ## set up an empty matrix
gene.min.matrix<-matrix(data=0, nrow = run, ncol = nrow(principal.col.coord))    ## set up an empty matrix
for(i in 1:run){

   knn.graph.cells=nng(standard.row.coord, mutual = F, k= k.set[i])   ## NNgraph algorithm based on Euclidean distance
   adjacency.matrix <- igraph::as_adjacency_matrix(knn.graph.cells)
   partition <- leiden(adjacency.matrix, seed=seed)   ## Leiden detection


   

   gene.rank.matrix<-matrix(NA, ncol = nrow(principal.col.coord), nrow =length(unique(partition)))
   colnames(gene.rank.matrix)<-gene.name
   cluster<-sort(unique(partition))
   for (j in 1:length(unique(partition))) {
        gene.rank.matrix[j,]<-rank(sqrt(rowSums((t(matrix(rep(colSums(standard.row.coord[partition==j,])/n.comp,   ## Calculate cluster centroid and rank within each cluster
              nrow(principal.col.coord)), ncol =nrow(principal.col.coord)))-principal.col.coord)^2)))
   }
   
   
   range.func<-function(x){
       diff(range(x))      ## define a range function
       }

  gene.range.matrix[i,]<-apply(gene.rank.matrix,2,range.func)   ##  a matrix with each gene's rank range for each run
  gene.min.matrix[i,]<-apply(gene.rank.matrix, 2, min)          ##  a matrix with each gene's minimum rank for each run

  
  }


colnames(gene.range.matrix)<-gene.name
colnames(gene.min.matrix)<-gene.name

gene.range<-colSums(gene.range.matrix)       ## take the sum of rank range for each gene for all runs (equivalent to take the average) 
gene.min<-colSums(gene.min.matrix)           ## take the sum of minimum rank for each gene for all runs (equivalent to take the average) 


order1<-order(gene.range,decreasing = T)[1:n.feature.1]   ## first round of filtering genes
gene.name.order1<-gene.name[order1]

filter.list.1<-gene.min[order1]

order2<-order(filter.list.1,decreasing = F)[1:n.feature.2]  ## second round of filtering genes


gene.list<-gene.name.order1[order2]                   ## final features list





obj<-list(genelist=gene.list, gene.range.matrix=gene.range.matrix, gene.min.matrix=gene.min.matrix)  ## the final list to return




return(obj)

}

















