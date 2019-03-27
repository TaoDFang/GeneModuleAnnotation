library(ggplot2)
library(reshape2)
library(Matrix)
library(ribiosMath)
library(plotly)
library(cogena)
library(Matrix)
library(glmnet)

#setwd("/media/sf_sharded_ubuntu/Disease_module_identification_DREAM_change")
setwd("/pstore/home/fangt3/Disease_module_identification_DREAM_change")
source("Module_Identification_funcs.R")

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

for(net_index in 1:length(unique_network_names)){
  print(net_index)
  network_genesets=module_genesets[network_names %in% unique_network_names[net_index]]
  network_all_genes=vector()
  for (i in 1:length(network_genesets)){
    network_all_genes=union(network_all_genes,network_genesets[[i]])
  }
  
  load("gene_pathway_matrix.rda")     # 20227 17599
  all_genes=rownames(gene_pathway_matrix)
  all_pathways=colnames(gene_pathway_matrix)
  
  common_genes=intersect(network_all_genes,all_genes)        # 16637
  gene_pathway_matrix=gene_pathway_matrix[common_genes,]
  all_genes=rownames(gene_pathway_matrix)
  all_pathways=colnames(gene_pathway_matrix)
  
  for(module_index in 1:length(network_genesets)){
    module_genes=network_genesets[[module_index]]
    module_name=names(network_genesets[module_index])
    
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
        
        results=matrix(nrow = length(selected_pathways),ncol = 7)
        rownames(results)=selected_pathways
        colnames(results)=c("pathway_name","fish_Pvalue","reg_coeff","GenesInPathway","Go_Reactome_root_id","Go_Reactome_root_names","step2root")
        results[,"step2root"]=unlist(get_steps(selected_pathways,go_ontology,reactome_ontology,go_roots,reactome_roots))
        results[,"Go_Reactome_root_id"]=unlist(find_root_ids(selected_pathways ))
        results[,"Go_Reactome_root_names"]=unlist(find_root_names(selected_pathways ))
        results[,"pathway_name"]=selected_pathway_names
        results[,"reg_coeff"]=format(selected_coef,scientific=TRUE,digits=3)
        results[,"fish_Pvalue"]=format(selected_pathways_fisher_pvalue,scientific=TRUE,digits=3)
        results[,"GenesInPathway"]=selected_pathways_num_genes
        
        
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
        #####################plot################################
        pathway_colors=colorRampPalette(brewer.pal(9, "Set1"))(length(selected_pathways))
        names(pathway_colors)=selected_pathway_names
        pathway2coeff_p=pathway2coeff_plot(selected_pathway_names,selected_coef,pathway_colors)
        log_p_value=-log(selected_pathways_fisher_pvalue)
        pathway2fisher_p=pathway2fisher_plot(selected_pathway_names,log_p_value,pathway_colors)
        pathway2gene_p=pathway2gene_plot(gene_pathway_matrix,selected_pathways,selected_pathway_names,selected_pathways_num_genes,module_common_genes)
        
        if(length(results[grepl("HSA",results[,"Go_Reactome_root_id"]),"Go_Reactome_root_names"])>0)
        {
          plot_reactome_p=plot_reactome_statistics(results)
          ggsave(paste("HTML_results/",module_name,"_reactome_statistics.PNG",sep=""),  width =10 ,height =1, 
                 plot = plot_reactome_p)
        }
        
        if(length(results[grepl("GO",results[,"Go_Reactome_root_id"]),"Go_Reactome_root_names"])>0)
        {
          plot_go_p=plot_go_statistics(results)
          ggsave(paste("HTML_results/",module_name,"_go_statistics.PNG",sep=""),  width =10 ,height =2.5, 
                 plot = plot_go_p)
        }
        
   
        ggsave(paste("HTML_results/",module_name,"_pathway2coef.PNG",sep=""), width =20 ,height =3.5, 
               plot = pathway2coeff_p)
        ggsave(paste("HTML_results/",module_name,"_pathway2fisher.PNG",sep=""), width =20 ,height =3.5, 
               plot = pathway2fisher_p)
        ggsave(paste("HTML_results/",module_name,"_pathway2gene.PNG",sep=""),  width =20 ,height =3.5, 
               plot = pathway2gene_p)
        
        ggsave(paste("HTML_results/",module_name,"_pathway2coef.PNG",sep=""), width =20 ,height =3.5, 
               plot = pathway2coeff_p)
        ggsave(paste("HTML_results/",module_name,"_pathway2fisher.PNG",sep=""), width =20 ,height =3.5, 
               plot = pathway2fisher_p)
        
        
        
        
        patterns=c("signaling","pathway")
        patterns=paste(patterns,collapse ="|")
        pattern_pathways=selected_pathway_names[grepl(patterns,selected_pathway_names,ignore.case = TRUE)]
        HTMLStart(outdir="./HTML_results/", file=module_name,
                  extension="html", echo=FALSE,HTMLframe = FALSE,CSSFile = "R2HTML_tao.css")
        HTML.title(module_name, HR=1,CSSFile = "R2HTML_tao.css")
        HTML.title("data_frame", HR=2,CSSFile = "R2HTML_tao.css")
        HTML(results,CSSFile = "R2HTML_tao.css")
        HTML.title("Pathways containing 'signaling ' or other words in their names", HR=2,CSSFile = "R2HTML_tao.css")
        HTML(paste(pattern_pathways,sep = "<br/>"))
        
        if(length(results[grepl("GO",results[,"Go_Reactome_root_id"]),"Go_Reactome_root_names"])>0)
        {
          HTML.title("go_statistics plot", HR=2,CSSclass="")
          HTMLInsertGraph(paste(module_name,"_go_statistics.PNG",sep=""),file =paste("./HTML_results/",module_name,".html",sep = "") ,
                          WidthHTML =1200 )
        }
        
        if(length(results[grepl("HSA",results[,"Go_Reactome_root_id"]),"Go_Reactome_root_names"])>0)
        {
          HTML.title("reactome_statistics plot", HR=2,CSSclass="")
          HTMLInsertGraph(paste(module_name,"_reactome_statistics.PNG",sep=""),file =paste("./HTML_results/",module_name,".html",sep = "") ,
                          WidthHTML =1200 )
        }
        
        
        HTML.title("pathway2fisher plot", HR=2,CSSclass="")
        HTMLInsertGraph(paste(module_name,"_pathway2fisher.PNG",sep=""),file =paste("./HTML_results/",module_name,".html",sep = ""),
                        WidthHTML =1200 )
        # pathway2fisher_p
        # HTMLplot(Width = 1200, Height = 350,CSSclass="") 
        HTML.title("pathway2coeff plot", HR=2,CSSclass="")
        HTMLInsertGraph(paste(module_name,"_pathway2coef.PNG",sep=""),file =paste("./HTML_results/",module_name,".html",sep = ""),
                        WidthHTML =1200 )
        # pathway2coeff_p
        # HTMLplot(Width = 1200, Height = 350,CSSclass="")
        HTML.title("pathway2gene plot", HR=2,CSSclass="")
        HTMLInsertGraph(paste(module_name,"_pathway2gene.PNG",sep=""),file =paste("./HTML_results/",module_name,".html",sep = "") ,
                        WidthHTML =1200 )
        
        # pathway2gene_p
        # HTMLplot(Width = 1200, Height = 350,CSSclass="")
        HTMLStop()
        
        #####################plot################################
        
      }
    }
    
  }
  
  
}

write.csv(dream_consensus_modules_results, file = "dream_consensus_modules_results.csv",
          row.names=FALSE, col.names = TRUE)
