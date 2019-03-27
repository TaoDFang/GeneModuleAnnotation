library(ggplot2)
library(reshape2)
library(Matrix)
library(ribiosMath)
library(plotly)
library(cogena)
library(Matrix)
library(glmnet)
library(ribiosAnnotation)
library(igraph)
library(data.table)
library(RColorBrewer)
library(ribiosPlot)
library(R2HTML)
#setwd("/media/sf_sharded_ubuntu/Disease_module_identification_DREAM_change")

# set the workplace, mainly to use some preprocessed data
setwd("~/Documents/DreamChallengeModuleAnnotation/HPC_Disease_module_identification_DREAM_change/")

## set path to R scripts
script_path="~/Documents/DreamChallengeModuleAnnotation/GeneModuleAnnotation/code/"

## set path to laterest results
result_path="~/Documents/DreamChallengeModuleAnnotation/Results/"


source(paste(script_path,"Module_Identification_funcs.R",sep = ""))

# get root nodes from hierarchical tree
go_ontology <- readRDS("go_ontology.rds")
go_roots=c("GO:0003674","GO:0005575","GO:0008150")
go_sub_roots=sapply(go_roots, function(x){
  as_ids(neighbors(go_ontology, x, mode = "out"))
})
go_sub_roots=unlist(go_sub_roots)

reactome_ontology=readRDS("human_reactome_ontology.rds")
reactome_roots=which(sapply(sapply(V(reactome_ontology), function(x) neighbors(reactome_ontology,x, mode="in")), length) == 0)
reactome_roots=names(reactome_roots)

go_ontology_names=V(go_ontology)$pathway_names
names(go_ontology_names)=as_ids(V(go_ontology))
reactome_ontology_names=V(reactome_ontology)$pathway_names
names(reactome_ontology_names)=as_ids(V(reactome_ontology))

######################################

module_genesets=scan("dream_consensus_modules.gmt",what = "",sep = "\n")
module_genesets=strsplit(module_genesets,split = "\t")

module_names=sapply(module_genesets, function(x){x[1]})
network_names=lapply(module_names, function(x){
  strsplit(x,split = "_")[[1]][1]
})
unique_network_names=unique(network_names)

module_genesets=lapply(module_genesets,function(x){
  x_genes=x[-c(1,2)]
  return(x_genes)
})
names(module_genesets)=module_names


#######################

dream_consensus_modules_results=data.frame()

#for(net_index in 1:length(unique_network_names)){
for(net_index in 1:1){  
  print(net_index)
  Daniel_results_name="dream_consensus_modules.functional_enrichment.txt"
  Daniel_results=read.csv(Daniel_results_name,header = TRUE,sep = "\t",stringsAsFactors = FALSE)
  Daniel_network_results=Daniel_results[Daniel_results$network==unique_network_names[net_index],]
  Daniel_network_simple_results=Daniel_network_results[,c("pathwayDb","network","module","term","termId","P.noncentral","P.noncentral.fdr" )]
  
  network_genesets=module_genesets[network_names %in% unique_network_names[net_index]]
  network_all_genes=vector()
  for (i in 1:length(network_genesets)){
    network_all_genes=union(network_all_genes,network_genesets[[i]])
  }
  
  load("gene_pathway_matrix.rda")     # 20227 17599
  all_genes=rownames(gene_pathway_matrix)
  all_pathways=colnames(gene_pathway_matrix)

  
  for(module_index in 1:length(network_genesets)){
    module_genes=network_genesets[[module_index]]
    module_name=names(network_genesets[module_index])
    Daniel_module_results=Daniel_network_simple_results[Daniel_network_results$module==module_name,]
    
    module_labels=rep(0,length(all_genes))          #len:20244
    names(module_labels)=all_genes
    module_common_genes=intersect(all_genes,module_genes) 
    if(length(module_common_genes)>1){
      module_labels[module_common_genes]=1
      cvfit=glmnet(gene_pathway_matrix,module_labels,lambda = 0.007956622,alpha = 0.5)   #s5;0.007956622, dream[[1]]:
      coef=coef(cvfit, s = "lambda.min")
      non0index=coef@i[-1]   #remove intercept
      non0coef=coef@x[-1]
      selected_index=non0index[which(non0coef>0)]
      selected_pathways=all_pathways[selected_index]
      selected_coef=non0coef[which(non0coef>0)]
      names(selected_coef)=selected_pathways
      
      if(length(selected_pathways)>0){
        selected_pathway_names=from_id2name(selected_pathways )
        names(selected_coef)=selected_pathway_names
        fisher_exact_test_results=fisher_exact_test(selected_pathways,module_common_genes,gene_pathway_matrix )
        selected_pathways_fisher_pvalue=fisher_exact_test_results$selected_pathways_fisher_pvalue
        selected_pathways_num_genes=fisher_exact_test_results$selected_pathways_num_genes
        
        module_results=data.frame(module=rep(module_name,length(selected_pathways)),
                                  pathway_id=selected_pathways,
                                  pathway_name=selected_pathway_names,
                                  fish_Pvalue=format(selected_pathways_fisher_pvalue,scientific=TRUE,digits=3),
                                  reg_coeff=format(selected_coef,scientific=TRUE,digits=3),
                                  GenesInPathway=selected_pathways_num_genes,
                                  Go_Reactome_root_id=unlist(find_root_ids(selected_pathways )),
                                  Go_Reactome_root_names=unlist(find_root_names(selected_pathways )),
                                  step2root=unlist(get_steps(selected_pathways,go_ontology,reactome_ontology,go_roots,reactome_roots)),
                                  row.names = NULL)
        dream_consensus_modules_results=rbind(dream_consensus_modules_results,module_results)
        
      }
    }
    
  
    
  }
  
  
}

write.csv(dream_consensus_modules_results, file = paste(result_path,"dream_consensus_modules_results_alpha0.5.csv",sep = ""),
          row.names=FALSE, col.names = TRUE)
