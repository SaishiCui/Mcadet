
length_output=function(method, datatype, resolution){
  number_list = c()
  true_number_list = c()
if(datatype == "SparSim"){  
for (i in 1:24) {
if(resolution == "coarse"){
  if(method == "Mcadet"){
  number_list= append(number_list, length(get(paste0(method, "_",
                                              datatype, "_", i))$gene))}else{
                                                
  number_list= append(number_list, 
                    length(get(paste0(method, "_",datatype, "_", i))))                                             
                                              }

}else{
  if(method == "Mcadet"){
  number_list= append(number_list, length(get(paste0(method, "_",
                                              datatype, "_fine_", i))$gene))}else{
                                                
  number_list= append(number_list, 
                    length(get(paste0(method, "_",datatype, "_fine_", i))))                                             
  }
 }
  true_number_list = append(true_number_list, length(true_SparSim_vector))
}
  
}else{
  for (p in 1:8) {
    for (k in 1:3) {
      
      donor.set<-c(1,2,3,4,5,6,7,8)
      time.set<-c(0,3,7)
      
    if(resolution == "coarse"){
      
      if(method == "Mcadet"){
      number_list= append(number_list, length(get(paste0(method, "_",
      "pbmc", "_", donor.set[p], "_", time.set[k]))$gene))}else{
      number_list= append(number_list, length(get(paste0(method, "_",
      datatype, "_", donor.set[p], "_", time.set[k]))))        
      }
    }else{
      if(method == "Mcadet"){
      number_list= append(number_list, length(get(paste0(method, "_",
      "pbmc", "_fine_", donor.set[p], "_", time.set[k]))$gene))}else{
      number_list= append(number_list, length(get(paste0(method, "_",
      datatype, "_fine_", donor.set[p], "_", time.set[k]))))}
    }
    }
    
    
    if(resolution == "coarse"){
    true_number_list = append(true_number_list, length(get(paste0("DEgene_",
                                        donor.set[p], "_", time.set[k]))))}else{
      
    true_number_list = append(true_number_list,  length(get(paste0("DEgene_fine_",
                                                   donor.set[p], "_", time.set[k]))))
      
    }
  }
 
}
 
  return(list(quantile(number_list), mean(number_list), sd(number_list), 
              round(mean(true_number_list)) ))
  
  
}
  




method_list = c("Mcadet", "NBDrop", "M3Drop",
                "scryFS", "BreFS","SeuratVst", "SeuratMvp")
datatype_list = c("SparSim", "PBMC")
resolution_list = c("coarse", "fine")
first_q = c()
second_q = c()
third_q = c()
mean = c()
sd = c()
true_number = c()
method_vec = c()
datatype_vec = c()
resolution_vec = c()
for (j in 1:2) {
  for (k in 1:2) {
    for (i in 1:7) {
      first_q = append(first_q ,length_output(method_list[i], datatype_list[j], resolution_list[k])[[1]][2])
      second_q = append(second_q, length_output(method_list[i], datatype_list[j], resolution_list[k])[[1]][3])
      third_q = append(third_q, length_output(method_list[i], datatype_list[j], resolution_list[k])[[1]][4])
      mean = append(mean,round(length_output(method_list[i], datatype_list[j], resolution_list[k])[[2]]))
      sd = append(sd, round(length_output(method_list[i], datatype_list[j], resolution_list[k])[[3]],2))
      true_number = append(true_number, length_output(method_list[i], datatype_list[j], resolution_list[k])[[4]])
      method_vec = append(method_vec, method_list[i])
      datatype_vec = append(datatype_vec, datatype_list[j])
      resolution_vec = append(resolution_vec, resolution_list[k])
    }
    
  }
  
}

out_df <- data.frame("first quantile" = first_q,
           "median" = second_q,
           "third quantile" = third_q,
           "mean" = mean,
           "sd" = sd,
          "true number"=  true_number,
           "Method" = method_vec,
           "Data" = datatype_vec,
           "resolution" = resolution_vec)





