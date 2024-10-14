
jaccard = function(a,b){
  numerator = length(intersect(a,b))
  denominator = length(a) + length(b) - numerator
  return( numerator/denominator)
}


true_SparSim_vector = c()
for(i in 1:2000){
  true_SparSim_vector = append(true_SparSim_vector, paste0("Gene", i))
  
}


for (p in 1:8) {
  for (k in 1:3) {
    donor.set<-c(1,2,3,4,5,6,7,8)
    time.set<-c(0,3,7)
    load(file =paste0("C:/Users/13541/Desktop/Thesis/pbmc/DE_genes/DEgene_fine_", 
                      donor.set[p], "_", time.set[k], ".RData"))
    load(file =paste0("C:/Users/13541/Desktop/Thesis/pbmc/DE_genes/DEgene_", 
                      donor.set[p], "_", time.set[k], ".RData"))
  }
}



########## Linegraph jaccard similarity #############


generate_linegraph_data_Jac = function(Method, dataset, resolution){

  mean_vec = c()
  sd_vec = c()
  number_vec = c()
  method_vec = c()

for (j in 1:15 ) {
  number_set = seq(200,3000,200)
  jaccard_list = c()
  if(dataset == "SparSim"){
  for (i in 1:24) {
    if(resolution == "coarse"){  
      if(Method == "Mcadet"){
        gene = get(paste0(Method, "_", dataset, "_3000_", i))[1:number_set[j]]
      }else{
        gene = get(paste0(Method, "_", dataset, "_3000_", i))[1:number_set[j]]
      }}else{
        if(Method == "Mcadet"){
          gene = get(paste0(Method, "_", dataset, "_fine_3000_", i))[1:number_set[j]]
        }else{
          gene = get(paste0(Method, "_", dataset, "_fine_3000_", i))[1:number_set[j]]
        }
      }
    jac = jaccard(gene, true_SparSim_vector)
    jaccard_list = append(jaccard_list, jac)  
    }
}else{
    for (p in 1:8) {
      for (k in 1:3) {
        donor.set<-c(1,2,3,4,5,6,7,8)
        time.set<-c(0,3,7)
        if(resolution == "coarse"){  
          if(Method == "Mcadet"){
            gene = get(paste0(Method, "_", "pbmc", "_3000_", donor.set[p], "_", time.set[k]))[1:number_set[j]]
          }else{
            gene = get(paste0(Method, "_", dataset, "_3000_", donor.set[p], "_", time.set[k]))[1:number_set[j]]
          }}else{
            if(Method == "Mcadet"){
              gene = get(paste0(Method, "_", "pbmc", "_fine_3000_", donor.set[p], "_", time.set[k]))[1:number_set[j]]
            }else{
              gene = get(paste0(Method, "_", dataset, "_fine_3000_", donor.set[p], "_", time.set[k]))[1:number_set[j]]
            }
          }
        if(resolution == "coarse"){
          true_PBMC_vector = get(paste0("DEgene_", donor.set[p], "_", time.set[k]))
        }else{
          true_PBMC_vector = get(paste0("DEgene_fine_", donor.set[p], "_", time.set[k]))
        }
        
        jac = jaccard(gene, true_PBMC_vector)
        jaccard_list = append(jaccard_list, jac)
      }
      }
    }
  mean_vec = append(mean_vec, mean(jaccard_list))
  sd_vec = append(sd_vec, sd(jaccard_list)) 
  number_vec = append(number_vec, number_set[j])
  
  if(Method == "Mcadet"){
    method_vec = append(method_vec, "Mcadet")
  }else{
    method_vec = append(method_vec, Method)
  }
  
}
  
  obj = list(mean_vec, sd_vec, number_vec, method_vec)
  return(obj)
}









generate_linegraph_df_Jac <- function(datatype, resolution){

mean_vec = c()
sd_vec = c()
number_vec = c()
method_vec = c()

for (i in 1:7) {
  method_list = c("Mcadet", "scryFS", "NBDrop", "M3Drop",
                  "SeuratDisp", "SeuratVst", "random")
  
  mean_vec = append(mean_vec,
        generate_linegraph_data_Jac(method_list[i], datatype, resolution)[[1]])
  sd_vec = append(sd_vec,
              generate_linegraph_data_Jac(method_list[i], datatype, resolution)[[2]])
  
  number_vec = append(number_vec,
            generate_linegraph_data_Jac(method_list[i], datatype, resolution)[[3]])
  
  method_vec = append(method_vec,
                      generate_linegraph_data_Jac(method_list[i], datatype, resolution)[[4]])

 }
output_df = data.frame(mean_vec, sd_vec, number_vec, method_vec)
return(output_df)
}


linegraph_SparSim_coarse_df_Jac = generate_linegraph_df_Jac("SparSim", "coarse")
linegraph_SparSim_fine_df_Jac = generate_linegraph_df_Jac("SparSim", "fine")
linegraph_PBMC_coarse_df_Jac = generate_linegraph_df_Jac("PBMC", "coarse")
linegraph_PBMC_fine_df_Jac = generate_linegraph_df_Jac("PBMC", "fine")




## SparSim clustering performance evaluation


generate_linegraph_data_SparSim_clustering = function(Method,  resolution){

meansil_vec =c()
meanpurity_vec= c()
meanari_vec = c()
meannmi_vec = c()
meanknn_vec = c()
number_vec = c()
method_vec = c()
for (j in 1:10 ) {
  number_set = seq(300,3000,300)
  sil_list = c()
  ari_list = c()
  nmi_list = c()
  purity_list = c()
  knn_list = c()

for (i in 1:24) {
  
  if(resolution == "fine"){
  load(file =paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_fine_", 
                    i, ".RData"))
  workdata=get(paste0("SparSim_fine_", i))
  gene.name<-rownames(workdata)  
  
  if(Method == "Mcadet"){
    gene = get(paste0(Method, "_", "SparSim", "_fine_3000_", i))[1:number_set[j]]
  }else{
    gene = get(paste0(Method, "_", "SparSim", "_fine_3000_", i))[1:number_set[j]]
  }}else{
    load(file =paste0("C:/Users/13541/Desktop/Thesis/SPARSim/SparSim_Mcadet_data/SparSim_", 
                      i, ".RData"))
    workdata=get(paste0("SparSim_", i))
    gene.name<-rownames(workdata)  
    
    if(Method == "Mcadet"){
      gene = get(paste0(Method, "_", "SparSim", "_3000_", i))[1:number_set[j]]
    }else{
      gene = get(paste0(Method, "_", "SparSim", "_3000_", i))[1:number_set[j]]
    }
  }
  
    cell.info<-c(rep("Type 1",  60),
                 rep("Type 2",  180),
                 rep("Type 3",  240),
                 rep("Type 4",  360),
                 rep("Type 5",  420),
                 rep("Type 6",  480),
                 rep("Type 7",  600),
                 rep("Type 8",  960),
                 rep("Type 9",  1200),
                 rep("Type 10", 1500))

  
  eva_value = eva(gene, cell.info, gene.name, workdata,
                  n_components= 15, 
                  kmeans.centers= 10, 
                  nstart = 10, 
                  k = 30)
  sil_list = append(sil_list,  eva_value[1])
  purity_list = append(purity_list, eva_value[2])
  ari_list = append(ari_list, eva_value[3])
  nmi_list = append(nmi_list, eva_value[4])
  knn_list = append(knn_list, eva_value[5])
  
  if(resolution == "fine"){
  rm(list = paste0("SparSim_fine_", i))
  rm(workdata)
  gc()}else{
    rm(list = paste0("SparSim_", i))
    rm(workdata)
    gc()
  }
  print(c(j,i))
}

  meansil_vec = append(meansil_vec, mean(sil_list))
  meanpurity_vec = append(meanpurity_vec, mean(purity_list))
  meanari_vec = append(meanari_vec, mean(ari_list))
  meannmi_vec = append(meannmi_vec, mean(nmi_list))
  meanknn_vec = append(meanknn_vec, mean(knn_list))
  
  
  number_vec = append(number_vec, number_set[j])
  
  if(Method == "Mcadet"){
    method_vec = append(method_vec, "Mcadet")
  }else{
    method_vec = append(method_vec, Method)
  }
  
  
}

obj = list(meansil_vec, meanpurity_vec, meanari_vec,
           meannmi_vec, meanknn_vec, number_vec, method_vec)
return(obj)

}



## PBMC clustering performance evaluation

generate_linegraph_data_PBMC_clustering = function(Method,  resolution, kmeancenter){
  
  meansil_vec =c()
  meanpurity_vec= c()
  meanari_vec = c()
  meannmi_vec = c()
  meanknn_vec = c()
  number_vec = c()
  method_vec = c()
  for (j in 1:10 ) {
    number_set = seq(300,3000,300)
    sil_list = c()
    ari_list = c()
    nmi_list = c()
    purity_list = c()
    knn_list = c()
    
    for (p in 1:8) {
      for (q in 1:3) {
        donor.set<-c(1,2,3,4,5,6,7,8)
        time.set<-c(0,3,7)
        
        if(resolution == "fine"){
          load(file =paste0("C:/Users/13541/Desktop/Thesis/pbmc/data_fine/pbmc_fine_", 
                            donor.set[p], "_", time.set[q], ".RData"))
          load(file =paste0("C:/Users/13541/Desktop/Thesis/pbmc/data_fine/pbmc_fine_clabel_", 
                            donor.set[p], "_", time.set[q], ".RData"))
          
          workdata=get(paste0("pbmc_fine_", donor.set[p], "_", time.set[q]))
          gene.name<-rownames(workdata)  
          
          if(Method == "Mcadet"){
            gene = get(paste0(Method, "_", "pbmc", "_fine_3000_", donor.set[p], "_", time.set[q]))[1:number_set[j]]
          }else{
            gene = get(paste0(Method, "_", "PBMC", "_fine_3000_", donor.set[p], "_", time.set[q]))[1:number_set[j]]
          }}else{
            load(file =paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmcdata_", 
                              donor.set[p], "_", time.set[q], ".RData"))
            load(file =paste0("C:/Users/13541/Desktop/Thesis/pbmc/data/pbmccell1_", 
                              donor.set[p], "_", time.set[q], ".RData"))
            workdata=get(paste0("pbmcdata_", donor.set[p], "_", time.set[q]))
            gene.name<-rownames(workdata)  
            
            if(Method == "Mcadet"){
              gene = get(paste0(Method, "_", "pbmc", "_3000_", donor.set[p], "_", time.set[q]))[1:number_set[j]]
            }else{
              gene = get(paste0(Method, "_", "PBMC", "_3000_", donor.set[p], "_", time.set[q]))[1:number_set[j]]
            }
          }
        
        
        
        if(resolution == "fine"){
          cell.info <- get(paste0( "pbmc_fine_clabel_",donor.set[p], "_", time.set[q]))
        }else{
          cell.info <-get(paste0( "pbmccell1_",donor.set[p], "_", time.set[q]))
        }
        
        
        
        eva_value = eva(gene, cell.info, gene.name, workdata,
                        n_components= 15, 
                        kmeans.centers= kmeancenter, 
                        nstart = 10, 
                        k = 30)
        
        sil_list = append(sil_list,  eva_value[1])
        purity_list = append(purity_list, eva_value[2])
        ari_list = append(ari_list, eva_value[3])
        nmi_list = append(nmi_list, eva_value[4])
        knn_list = append(knn_list, eva_value[5])
        
        if(resolution == "fine"){
          rm(list = paste0("pbmc_fine_", donor.set[p], "_", time.set[q]))
          rm(workdata)
          gc()}else{
            rm(list = paste0("pbmcdata_",donor.set[p], "_", time.set[q]))
            rm(workdata)
            gc()
          }
        print(c(j,p,q))
      }
    }
    
    meansil_vec = append(meansil_vec, mean(sil_list))
    meanpurity_vec = append(meanpurity_vec, mean(purity_list))
    meanari_vec = append(meanari_vec, mean(ari_list))
    meannmi_vec = append(meannmi_vec, mean(nmi_list))
    meanknn_vec = append(meanknn_vec, mean(knn_list))
    
    
    number_vec = append(number_vec, number_set[j])
    
    if(Method == "Mcadet"){
      method_vec = append(method_vec, "Mcadet")
    }else{
      method_vec = append(method_vec, Method)
    }
    
    
  }
  
  obj = list(meansil_vec, meanpurity_vec, meanari_vec,
             meannmi_vec, meanknn_vec, number_vec, method_vec)
  return(obj)
  
}

################# generate final dataframe for evaluation results


generate_linegraph_df_clustering <- function(resolution, datatype, 
                                                     kmeancenter){
  
  meansil_vec =c()
  meanpurity_vec= c()
  meanari_vec = c()
  meannmi_vec = c()
  meanknn_vec = c()
  number_vec = c()
  method_vec = c()
  
  for (s in 1:7) {
    method_list = c("Mcadet", "scryFS", "NBDrop", "M3Drop",
                    "SeuratDisp", "SeuratVst", "random")
    if(datatype == "SparSim"){
    result = generate_linegraph_data_SparSim_clustering(method_list[s], resolution)
    }else{
      result = generate_linegraph_data_PBMC_clustering(method_list[s], resolution, kmeancenter)
      
    }
    meansil_vec = append(meansil_vec, result[[1]])
                     
    meanpurity_vec = append(meanpurity_vec, result[[2]])
    
    meanari_vec = append(meanari_vec, result[[3]])
    
    meannmi_vec = append(meannmi_vec, result[[4]])
    
    meanknn_vec = append(meanknn_vec, result[[5]])
    
    number_vec = append(number_vec, result[[6]])
    
    method_vec = append(method_vec, result[[7]])
    print(method_list[s])
  }
  output_df = data.frame(meansil_vec, meanpurity_vec, meanari_vec,
                         meannmi_vec, meanknn_vec, number_vec, method_vec)
  return(output_df)
}


linegraph_SparSim_coarse_df_clustering <- 
  generate_linegraph_df_clustering("coarse", "SparSim", 10)


linegraph_SparSim_fine_df_clustering <- 
  generate_linegraph_df_clustering("fine", "SparSim", 10)


linegraph_PBMC_coarse_df_clustering <-  
  generate_linegraph_df_clustering("coarse", "PBMC", 8)


linegraph_PBMC_fine_df_clustering <-  
  generate_linegraph_df_clustering("fine", "PBMC", 11)



### save and load the data

save(linegraph_SparSim_coarse_df_clustering, file = 
"C:/Users/13541/Desktop/Thesis/Results/Line_graph_trend/linegraph_SparSim_coarse_df_clustering.RData")
save(linegraph_SparSim_fine_df_clustering, file = 
       "C:/Users/13541/Desktop/Thesis/Results/Line_graph_trend/linegraph_SparSim_fine_df_clustering.RData")
save(linegraph_PBMC_coarse_df_clustering, file = 
       "C:/Users/13541/Desktop/Thesis/Results/Line_graph_trend/linegraph_PBMC_coarse_df_clustering.RData")
save(linegraph_PBMC_fine_df_clustering, file = 
       "C:/Users/13541/Desktop/Thesis/Results/Line_graph_trend/linegraph_PBMC_fine_df_clustering.RData")
save(linegraph_SparSim_coarse_df_Jac, file = 
       "C:/Users/13541/Desktop/Thesis/Results/Line_graph_trend/linegraph_SparSim_coarse_df_Jac.RData")
save(linegraph_SparSim_fine_df_Jac, file = 
       "C:/Users/13541/Desktop/Thesis/Results/Line_graph_trend/linegraph_SparSim_fine_df_Jac.RData")
save(linegraph_PBMC_coarse_df_Jac, file = 
       "C:/Users/13541/Desktop/Thesis/Results/Line_graph_trend/linegraph_PBMC_coarse_df_Jac.RData")
save(linegraph_PBMC_fine_df_Jac, file = 
       "C:/Users/13541/Desktop/Thesis/Results/Line_graph_trend/linegraph_PBMC_fine_df_Jac.RData")



load(file ="C:/Users/css22/Desktop/Thesis1/Results/Line_graph_trend/linegraph_SparSim_coarse_df_clustering.RData")
load(file ="C:/Users/css22/Desktop/Thesis1/Results/Line_graph_trend/linegraph_SparSim_fine_df_clustering.RData")
load(file ="C:/Users/css22/Desktop/Thesis1/Results/Line_graph_trend/linegraph_PBMC_coarse_df_clustering.RData")
load(file ="C:/Users/css22/Desktop/Thesis1/Results/Line_graph_trend/linegraph_PBMC_fine_df_clustering.RData")
load(file ="C:/Users/css22/Desktop/Thesis1/Results/Line_graph_trend/linegraph_SparSim_coarse_df_Jac.RData")
load(file ="C:/Users/css22/Desktop/Thesis1/Results/Line_graph_trend/linegraph_SparSim_fine_df_Jac.RData")
load(file ="C:/Users/css22/Desktop/Thesis1/Results/Line_graph_trend/linegraph_PBMC_coarse_df_Jac.RData")
load(file ="C:/Users/css22/Desktop/Thesis1/Results/Line_graph_trend/linegraph_PBMC_fine_df_Jac.RData")



linegraph_PBMC_coarse_df_Jac


## function of generating figures


library(ggplot2)

generate_linegraph_figure <- function(DF, evatype,
                                      plot_title, y_title, show_legend, label){
  
  if(evatype == "Jaccard"){
  DF$method_vec[DF$method_vec == "scryFS"] = "Scry"
  DF$method_vec[DF$method_vec == "SeuratDisp"] = "Seurat Disp"
  DF$method_vec[DF$method_vec == "SeuratVst"] = "Seurat Vst"
  
  DF$method_vec = factor(DF$method_vec, levels = c("Mcadet", "Scry",
                                                   "Seurat Disp", "Seurat Vst",
                                                   "M3Drop", "NBDrop", "random"))
  }else{
  
  
  DF = cbind(apply(DF[,1:5],1,mean), DF[,6:7])
  colnames(DF)[1] = "mean_vec"
  
  DF$method_vec[DF$method_vec == "scryFS"] = "Scry"
  DF$method_vec[DF$method_vec == "SeuratDisp"] = "Seurat Disp"
  DF$method_vec[DF$method_vec == "SeuratVst"] = "Seurat Vst"
  DF$method_vec = factor(DF$method_vec, levels = c("Mcadet", "Scry",
                                                   "Seurat Disp", "Seurat Vst",
                                                   "M3Drop", "NBDrop", "random"))
    
  }
  
  
  
  
  
  
  
  if(show_legend){ 
    out_figure = ggplot(data =  DF, 
                        mapping = aes(x = number_vec, y = mean_vec, 
                                      colour = method_vec)) + geom_line(linewidth = 1.2)+
      geom_point(size = 3)+
      labs(x = "", y = y_title, 
           title = plot_title )+
      theme_bw() + scale_colour_discrete(name = "Method", 
                                         labels = c("Mcadet", "Scry",
                                                    "Seurat Disp", "Seurat Vst",
                                                    "M3Drop", "NBDrop", "Random"))+
      scale_x_continuous(limits=c(200, 3000), breaks=seq(200, 3000, 200))+
      theme( panel.grid = element_blank(),  
             panel.grid.major = element_line(color = "grey", linetype = "dashed"),  
             panel.grid.minor = element_line(color = "grey", linetype = "dashed"),
             axis.title = element_text(face ="bold", size = 16),
             plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
             legend.position = "bottom",
             legend.text = element_text(size = 10, face = "bold"),
             legend.title = element_blank(), 
             legend.background = element_rect(color = "black", 
                                              linewidth = 1.5,
                                              fill = "lightblue"),
             axis.text.x = element_text(face = "bold", size = 12, color = "black"), 
             axis.text.y = element_text(face = "bold", size = 12),
             panel.border = element_rect(color = "black", fill = NA, size = 1.5))+
      annotate("text", x = -Inf, y = Inf, label = label, hjust = -1, vjust = 1.5, size = 8, fontface = "bold")
    
    
    
    }else{
      out_figure = ggplot(data =  DF, 
                          mapping = aes(x = number_vec, y = mean_vec, 
                                        colour = method_vec)) + geom_line(linewidth = 1.2)+
        geom_point(size = 3)+
        labs(x = "", y = y_title, 
             title = plot_title )+
        theme_bw() + scale_colour_discrete(name = "Method", 
                                           labels = c("Mcadet", "Scry",
                                                      "Seurat Disp", "Seurat Vst",
                                                      "M3Drop", "NBDrop", "Random"))+
        scale_x_continuous(limits=c(200, 3000), breaks=seq(200, 3000, 200))+
        theme( panel.grid = element_blank(),  
               panel.grid.major = element_line(color = "grey", linetype = "dashed"),  
               panel.grid.minor = element_line(color = "grey", linetype = "dashed"),
               axis.title = element_text(face ="bold", size = 16),
               plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
               legend.position = "None",
               legend.text = element_text(size = 10, face = "bold"),
               legend.title = element_blank(), 
               legend.background = element_rect(color = "black", 
                                                linewidth = 1.5,
                                                fill = "lightblue"),
               axis.text.x = element_text(face = "bold", size = 12, color = "black"), 
               axis.text.y = element_text(face = "bold", size = 12),
               panel.border = element_rect(color = "black", fill = NA, size = 1.5)) +
        annotate("text", x = -Inf, y = Inf, label = label, hjust = -1, vjust = 1.5, size = 8, fontface = "bold")
              
               
             }
  
  return(out_figure)
  
}





line_PBMC_coarse_Jaccard <- generate_linegraph_figure(DF = linegraph_PBMC_coarse_df_Jac, evatype = "Jaccard",
plot_title = "PBMC coarse-resolution datasets",  y_title = "Mean Jaccard similarity", show_legend = F, label = "A") +
  geom_vline(aes(xintercept = 1800), colour = "#990000", linetype = "dashed", linewidth = 1.1)


line_PBMC_fine_clustering <- generate_linegraph_figure(DF = linegraph_PBMC_fine_df_clustering, evatype = "clustering",
                                                      plot_title = "PBMC fine-resolution datasets",  
                                                      y_title = "", 
                                                    show_legend = F, label = "B") +
  geom_vline(aes(xintercept = 3000), colour = "#990000", linetype = "dashed", linewidth = 1.1)






line_simulated_coarse_clustering  <- generate_linegraph_figure(DF = linegraph_SparSim_coarse_df_clustering , 
                                                           evatype = "clustering ",
                                                      plot_title = "Simulated coarse-resolution datasets",  
                                                      y_title = "Mean averaged clustering metrics", show_legend = T, 
                                                      label = "C") +
  geom_vline(aes(xintercept = 1800), colour = "#990000", linetype = "dashed", linewidth = 1.1)


line_simulated_fine_clustering <- generate_linegraph_figure(DF = linegraph_SparSim_fine_df_clustering, evatype = "clustering",
                                                    plot_title = "Simulated fine-resolution datasets",  y_title = "", 
                                                    show_legend = F, label = "D") +
  geom_vline(aes(xintercept = 2100), colour = "#990000", linetype = "dashed", linewidth = 1.1)







(line_PBMC_coarse_clustering + line_PBMC_fine_clustering) /
  (line_simulated_coarse_clustering + line_simulated_fine_clustering )




