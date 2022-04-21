
rm(list = ls())
gc()
library(M3Drop)
library(M3DExampleData)
library(dplyr)
library(Seurat)
library(patchwork)


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






# M3Drop Michaelis-Menten (should normalize the data)

norm <- M3DropConvertData(raw_data, is.counts=TRUE)
norm <- M3DropConvertData(log2(norm+1), is.log=TRUE, pseudocount=1)


M3Drop_genes_2000<- M3DropFeatureSelection(norm,mt_method="fdr", mt_threshold=1)

M3Drop_genes<- M3DropFeatureSelection(norm,mt_method="fdr", mt_threshold=0.05)



M3drop=M3Drop_genes$Gene
M3drop_2000=M3Drop_genes_2000$Gene[1:2000]


#################  Splatter simulation  #####################

rm(list = ls())
gc()
library(M3Drop)
for (i in 1:4) {
  for (j in 1:4) {
    deprob.set<-c(0.02,0.025,0.03,0.035)
    defacScale.set<-c(0.3,0.4,0.5,0.6)
    load(file =paste0("E:/Thesis/simdata_",deprob.set[i], "_", defacScale.set[j],".RData"))
    
  }
}










library(M3Drop)
for (i in 1:4) {
  for (j in 1:4) {
    deprob.set<-c(0.02,0.025,0.03,0.035)
    defacScale.set<-c(0.3,0.4,0.5,0.6)
    
    raw_data=get(paste0("simdata_",deprob.set[i],"_", defacScale.set[j]))
    
    nb_raw <- NBumiConvertData(raw_data, is.counts=TRUE)
    DANB_fit <- NBumiFitModel(nb_raw)
    
    
    NBDropFS<- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thres=1, suppress.plot=TRUE)
    
    assign(paste0("NB_",deprob.set[i],"_", defacScale.set[j]), NBDropFS$Gene[1:2000])
    assign(paste0("Bre_",deprob.set[i],"_", defacScale.set[j]),
           BrenneckeGetVariableGenes(log(raw_data/10000+1), fdr = 2)[,1][1:2000])
    
    save(list = paste0("NB_",deprob.set[i],"_", defacScale.set[j]),  
         file = paste0("E:/Thesis/NB_",deprob.set[i], "_", defacScale.set[j],".RData"))

    save(list = paste0("Bre_",deprob.set[i],"_", defacScale.set[j]),  
         file = paste0("E:/Thesis/Bre_",deprob.set[i], "_", defacScale.set[j],".RData"))
    
  }
}



library(Seurat)

for(i in 1:4){
  for(j in 1:4){
    deprob.set<-c(0.02,0.025,0.03,0.035)
    defacScale.set<-c(0.3,0.4,0.5,0.6)
    
    seurat.obj<-CreateAssayObject(counts = get(paste0("simdata_",deprob.set[i],"_", defacScale.set[j])))
    seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat.obj.vst <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
    seurat.list<-VariableFeatures(seurat.obj.vst)
    assign(paste0("Seurat_",deprob.set[i],"_", defacScale.set[j]), seurat.list)
    save(list = paste0("Seurat_",deprob.set[i],"_", defacScale.set[j]),  
         file = paste0("E:/Thesis/Seurat_",deprob.set[i], "_", defacScale.set[j],".RData"))
    
  }
}






###########  SparSim simulation  #######

rm(list = ls())
gc()




library(Seurat)
library(M3Drop)


for (i in 1:16) {
     
    load(file = paste0("E:/Thesis/SPARSim/SparSim_", i,".RData") )
    raw_data=get(paste0("SparSim_",i ))
    
    nb_raw <- NBumiConvertData(raw_data, is.counts=TRUE)
    DANB_fit <- NBumiFitModel(nb_raw)
    
    
    NBDropFS<- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thres=1, suppress.plot=TRUE)
    
    assign(paste0("NB_Sparsim_", i), NBDropFS$Gene[1:2000])
    assign(paste0("Bre_Sparsim_", i), BrenneckeGetVariableGenes(log(raw_data/10000+1), fdr = 2)[,1][1:2000])
    
    
    seurat.obj<-CreateAssayObject(counts = raw_data)
    seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat.obj.vst <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
    seurat.list<-VariableFeatures(seurat.obj.vst)
    assign(paste0("Seurat_Sparsim_", i), seurat.list)

    gene.name<-rownames(raw_data)  
    
    disp.value<-apply(log(raw_data/10000+1),1,sd)
    genelist.disp<-names(disp.value[order(disp.value,decreasing = T)[1:2000]])
    assign(paste0("disp_Sparsim_", i), genelist.disp)

    
    set.seed(i*100)
    genelist.random<- gene.name[sample(1:nrow(raw_data),2000)]
    assign(paste0("random_Sparsim_", i), genelist.random)

    
    save(list = paste0("NB_Sparsim_", i),  
         file = paste0("E:/Thesis/SPARSim/NB_Sparsim_", i, ".RData"))
    
    save(list = paste0("Bre_Sparsim_", i),  
         file = paste0("E:/Thesis/SPARSim/Bre_Sparsim_", i, ".RData"))
    
    save(list = paste0("Seurat_Sparsim_", i),  
         file = paste0("E:/Thesis/SPARSim/Seurat_Sparsim_", i, ".RData"))
    
    save(list = paste0("disp_Sparsim_", i), 
         file = paste0("E:/Thesis/SPARSim/disp_Sparsim_", i, ".RData" ))
    
    save(list = paste0("random_Sparsim_", i),
         file = paste0("E:/Thesis/SPARSim/random_Sparsim_", i, ".RData" ))
    
    rm(list = paste0("SparSim_",i))
    gc()
    
}











########### real data pbmc ####




library(M3Drop)
library(Seurat)
for (i in 1:8) {
  for (j in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file = paste0("E:/Thesis/pbmc/pbmc_fine_", donor.set[i], "_", time.set[j], ".RData" ) )
    raw_data=get(paste0("pbmc_fine_", donor.set[i], "_", time.set[j]))
    
    gene.name <- rownames(raw_data)
    
    numb<-sum(1*(rowSums(raw_data)!=0))
    
    nb_raw <- NBumiConvertData(raw_data, is.counts=TRUE)
    DANB_fit <- NBumiFitModel(nb_raw)
    NBDropFS<- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thres=1, suppress.plot=TRUE)
    assign(paste0("NB_fine_All_",donor.set[i], "_", time.set[j]), NBDropFS$Gene)
    assign(paste0("Bre_fine_All_",donor.set[i],"_", time.set[j]), 
           BrenneckeGetVariableGenes(log(raw_data/10000+1), fdr = 2)[,1])
    
    
    seurat.obj<-CreateAssayObject(counts =raw_data )
    seurat.obj<- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat.obj.vst <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = nrow(raw_data))
    seurat.list<-VariableFeatures(seurat.obj.vst)
    assign(paste0("Seurat_fine_All_", donor.set[i], "_", time.set[j]), seurat.list)
    
    disp.value<-apply(log(raw_data/10000+1),1,sd)
    genelist.disp<-names(disp.value[order(disp.value,decreasing = T)[1:numb]])
    assign(paste0("disp_fine_All_",donor.set[i], "_", time.set[j]), genelist.disp)

    set.seed((i-1)*3+j)
    genelist.random<- gene.name[sample(1:nrow(raw_data), nrow(raw_data) , replace = F)]
    assign(paste0("random_fine_All_",donor.set[i], "_", time.set[j]), genelist.random)

 
    save(list = paste0("NB_fine_All_",donor.set[i],"_", time.set[j]),  
         file = paste0("E:/Thesis/pbmc/NB_fine_All_",donor.set[i], "_",  time.set[j],".RData"))
    
    save(list = paste0("Bre_fine_All_", donor.set[i],"_", time.set[j]),  
         file = paste0("E:/Thesis/pbmc/Bre_fine_All_", donor.set[i], "_", time.set[j],".RData"))
    
    save(list = paste0("Seurat_fine_All_", donor.set[i],"_", time.set[j]),  
         file = paste0("E:/Thesis/pbmc/Seurat_fine_All_", donor.set[i], "_", time.set[j],".RData"))
    
    save(list = paste0("disp_fine_All_",donor.set[i], "_", time.set[j]),
         file = paste0("E:/Thesis/pbmc/disp_fine_All_", donor.set[i], "_", time.set[j], ".RData" ))
   
    save(list = paste0("random_fine_All_",donor.set[i], "_", time.set[j]),
         file = paste0("E:/Thesis/pbmc/random_fine_All_", donor.set[i], "_", time.set[j], ".RData" ))
     
    rm(list =paste0("pbmc_fine_", donor.set[i], "_", time.set[j]) )
    gc()
    
    print(c(i,j))
  }
}








