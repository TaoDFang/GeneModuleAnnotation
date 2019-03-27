

############regression functions
regression_selected_pathways <-function(gene_input,lambda=0.007956622,alpha=0.5){
  setwd("/pstore/home/fangt3/Disease_module_identification_DREAM_change")
  data_path="/pstore/home/fangt3/Disease_module_identification_DREAM_change/"
  go_ontology <- readRDS(paste(data_path,"go_ontology.rds",sep = ""))
  reactome_ontology=readRDS(paste(data_path,"human_reactome_ontology.rds",sep = ""))
  go_ontology_names=V(go_ontology)$pathway_names
  names(go_ontology_names)=as_ids(V(go_ontology))
  reactome_ontology_names=V(reactome_ontology)$pathway_names
  names(reactome_ontology_names)=as_ids(V(reactome_ontology))
  
  load(paste(data_path,"gene_pathway_matrix.rda",sep = ""))     # 20227 17599
  all_genes=rownames(gene_pathway_matrix)
  all_pathways=colnames(gene_pathway_matrix)
  
  module_labels=rep(0,length(all_genes))          #len:20244F
  names(module_labels)=all_genes
  module_common_genes=intersect(all_genes,gene_input) 
   
  if(length(module_common_genes)>1){                     # should set a lower thereshold for num of module common genes, more than 50%
    module_labels[module_common_genes]=1
    cvfit=glmnet(gene_pathway_matrix,module_labels,lambda = lambda,alpha =alpha)   #
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
      
      pathway2gene_p=pathway2gene_plot_new(gene_pathway_matrix,selected_coef,selected_pathways,selected_pathway_names,selected_pathways_num_genes,module_common_genes)
      
      new_order=order(selected_coef,decreasing = TRUE)
      return(list(selected_pathways=selected_pathways[new_order],
                  selected_pathway_names=selected_pathway_names[new_order],
                  selected_pathways_fisher_pvalue=selected_pathways_fisher_pvalue[new_order],
                  pathway2gene_p=pathway2gene_p))
    }else{
      return(NULL)}
  }
  
}
