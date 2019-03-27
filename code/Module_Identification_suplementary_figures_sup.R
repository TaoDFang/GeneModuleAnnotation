jaccard_index=c()
jaccard_module=c()
for(i in 1:length(tao_module_names)){
  selected_pathwayIDs=Tao_results[Tao_results$module==tao_module_names[i],"pathway_id"]
  selected_pathwayIDs_vec=1:length(selected_pathwayIDs)
  print(i)
  if(length(selected_pathwayIDs)>1){
    modulei_jaccard_index=sapply(1:length(selected_pathwayIDs),function(j){
      #print(j)
      if(j < length(selected_pathwayIDs)){
        compared_pathwayIDs=selected_pathwayIDs[selected_pathwayIDs_vec[-(1:j)]]
        selected_genes=all_genes[gene_pathway_matrix[,selected_pathwayIDs[j]]==1] ##!!! use common genes or not????
        return(max(sapply(compared_pathwayIDs, function(x){
          compared_genes=all_genes[gene_pathway_matrix[,x]==1]
          jaccard=length(intersect(selected_genes,compared_genes))/length(union(selected_genes,compared_genes))
          if(jaccard==1){
            print(selected_pathwayIDs[j])
            print(x)
          }
          return(jaccard)
        }
        )))
      }
    })
    if(max(unlist(modulei_jaccard_index))==1){
      jaccard_module=c(jaccard_module,tao_module_names[i])
    }
    jaccard_index=c(jaccard_index,median(unlist(modulei_jaccard_index)))
  }
}
jaccard_index=unlist(jaccard_index)
saveRDS(jaccard_index,file="tao_jaccard_index_median_0.5.rds")
saveRDS(jaccard_module,file="tao_jaccard_module_0.5.rds")

daniel_jaccard_index=c()
daniel_jaccard_module=c()
Daniel_network_GO_results=Daniel_network_results[grepl("GO",Daniel_network_results$termId),]
daniel_GO_module_names=unique(Daniel_network_GO_results$module)

for(i in 1:length(daniel_GO_module_names)){
  selected_pathwayIDs=Daniel_network_GO_results[Daniel_network_GO_results$module==daniel_GO_module_names[i],"termId"]
  selected_pathwayIDs=selected_pathwayIDs[selected_pathwayIDs %in% raw_go_list_names ]
  selected_pathwayIDs_vec=1:length(selected_pathwayIDs)
  print(i)
  if(length(selected_pathwayIDs)>1){
    modulei_jaccard_index=sapply(1:length(selected_pathwayIDs),function(j){
      if(j< length(selected_pathwayIDs)){
        compared_pathwayIDs=selected_pathwayIDs[selected_pathwayIDs_vec[-(1:j)]]
        #selected_genes=all_genes[gene_pathway_matrix[,selected_pathwayIDs[j]]==1]
        selected_genes=raw_go_list[[selected_pathwayIDs[j]]]
        return(max(sapply(compared_pathwayIDs, function(x){
          #compared_genes=all_genes[gene_pathway_matrix[,x]==1]
          compared_genes=raw_go_list[[x]]
          jaccard=length(intersect(selected_genes,compared_genes))/length(union(selected_genes,compared_genes))
          #if(jaccard==1){
          #  print(selected_pathwayIDs[j])
          #  print(x)
          #}
          return(jaccard)
        }
        )))
      }
    })
    if(max(unlist(modulei_jaccard_index))==1){
      daniel_jaccard_module=c(daniel_jaccard_module,daniel_GO_module_names[i])
    }
    daniel_jaccard_index=c(daniel_jaccard_index,median(unlist(modulei_jaccard_index)))
  }
}

daniel_jaccard_index=unlist(daniel_jaccard_index)
#saveRDS(daniel_jaccard_index,file="daniel_jaccard_index.rds")
#saveRDS(daniel_jaccard_index,file="daniel_jaccard_index_median.rds")
#saveRDS(daniel_jaccard_module,file = "daniel_jaccard_module.rds")
