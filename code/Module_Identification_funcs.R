# GO:0008150 bp root node, GO:0003674 mf root node, GO:0005575 cc root node
get_steps=function(selected_pathways,go_ontology,reactome_ontology,go_roots,reactome_roots){
  steps=lapply(selected_pathways,function(x){
    if(grepl("GO",x)){
      min(distances(go_ontology,v=go_roots,to=x))
    }else{
      min(distances(reactome_ontology,v=reactome_roots,to=x))
    }
  })
  return(steps)
}

find_root_ids=function(selected_pathways){
  find_roots=lapply(selected_pathways, function(x){
    if(grepl("GO",x)){
      parent_nodes=ego(go_ontology,order = gorder(go_ontology),x,mode = "in")
      find_roots=intersect(go_sub_roots,names(unlist(parent_nodes)))
    }else{
      parent_nodes=ego(reactome_ontology,order = gorder(reactome_ontology),x,mode = "in")
      find_roots= intersect(reactome_roots,names(unlist(parent_nodes)))
    }
    if(length(find_roots)>1){
      find_roots=paste0(unlist(find_roots),collapse = "#")
    }
    return(find_roots)
  })
  return(find_roots)
  
}

find_root_names=function(selected_pathways){
  find_roots=lapply(selected_pathways, function(x){
    if(grepl("GO",x)){
      parent_nodes=ego(go_ontology,order = gorder(go_ontology),x,mode = "in")
      find_roots=intersect(go_sub_roots,names(unlist(parent_nodes)))
      find_roots=go_ontology_names[find_roots]
    }else{
      parent_nodes=ego(reactome_ontology,order = gorder(reactome_ontology),x,mode = "in")
      find_roots= intersect(reactome_roots,names(unlist(parent_nodes)))
      find_roots=reactome_ontology_names[find_roots]
    }
    if(length(find_roots)>1){
      find_roots=paste0(unlist(find_roots),collapse  = "#")
    }
    return(find_roots)
  })
  
  return(find_roots)
}


find_sources=function(selected_pathways){
  sources=lapply(selected_pathways, function(x){
    if(grepl("GO",x)){
      source="GO"
    }else{
      source="Reactome"
    }
    
    return(source)
  })
  
  return(sources)
}

find_pathway_links<-function(selected_pathways,selected_pathway_names){
  links=lapply(1:length(selected_pathways), function(i){
    #print(i)
    if(grepl("GO",selected_pathways[i])){
      #http://amigo.geneontology.org/amigo/medial_search?q=GO%3A0005520
      url=paste("http://amigo.geneontology.org/amigo/medial_search?q=GO%3A",strsplit(selected_pathways[i],split=":")[[1]][2],sep = "")
      links=paste("<a href=\"",url,"\">",selected_pathway_names[i],"</a></p>",sep = "")
    }else{
      #https://reactome.org/content/query?q=R-HSA-1810476&species=Homo+sapiens&species=Entries+without+species&cluster=true
      url=paste("https://reactome.org/content/query?q=",selected_pathways[i],"&species=Homo+sapiens&species=Entries+without+species&cluster=true",sep = "")
      links=paste("<a href=\"",url,"\">",selected_pathway_names[i],"</a></p>",sep = "")
    }
    
    return(links)
  })
  return(links)
}

find_classes=function(selected_pathways,go_ontology,reactome_ontology,reactome_roots,reactome_ontology_names){
  find_roots=lapply(selected_pathways, function(x){
    if(grepl("GO",x)){
      ## GO:0008150 bp root node, GO:0003674 mf root node, GO:0005575 cc root node
      go_roots=c("GO:0003674","GO:0005575","GO:0008150")
      parent_nodes=ego(go_ontology,order = gorder(go_ontology),x,mode = "in")
      find_roots=intersect(go_roots,names(unlist(parent_nodes)))
      find_roots=go_ontology_names[find_roots]
    }else{
      parent_nodes=ego(reactome_ontology,order = gorder(reactome_ontology),x,mode = "in")
      find_roots= intersect(reactome_roots,names(unlist(parent_nodes)))
      find_roots=reactome_ontology_names[find_roots]
    }
    if(length(find_roots)>1){
      find_roots=paste0(unlist(find_roots),collapse  = "#")
    }
    return(find_roots)
  })
  
  return(find_roots)
}


from_id2name=function(selected_pathways){
  names=lapply(selected_pathways, function(x){
    if(grepl("GO",x)){
      unname(go_ontology_names[x])
    }else{
      unname(reactome_ontology_names[x])
    }
  })
  return(unname(unlist(names)))
}



full_fisher_exact_test=function(input_genes){
  load("gene_pathway_matrix.rda")     # 17176,6493
  all_genes=rownames(gene_pathway_matrix)
  all_pathways=colnames(gene_pathway_matrix)
  common_genes=intersect(all_genes,input_genes)
  fisher_test=fisher_exact_test(all_pathways,common_genes,gene_pathway_matrix)
  #all.equal(all_pathways,names(fisher_test2$selected_pathways_fisher_pvalue))
  fisher_test_orderPathways=all_pathways[order(fisher_test$selected_pathways_fisher_pvalue)]
  fisher_test_orderPathways_names=from_id2name(fisher_test_orderPathways)
  return(list(fisher_test_orderPathways=fisher_test_orderPathways,
              fisher_test_orderPathways_names=fisher_test_orderPathways_names,
              fisher_test_orderPvalue=fisher_test$selected_pathways_fisher_pvalue[order(fisher_test$selected_pathways_fisher_pvalue)]))
}



fisher_exact_test=function(selected_pathways,module1_common_genes,gene_pathway_matrix){ 
  selected_pathways_fisher_pvalue=vector()
  selected_pathways_num_genes=vector()
  for(index in 1:length(selected_pathways)){
    fisher_pathway=selected_pathways[index]
    #          genes_in_common      genes_in_reference
    #pathway          a                   c
    #no_pathway       b                   d
    a=sum(gene_pathway_matrix[module1_common_genes,fisher_pathway])
    b=length(module1_common_genes)-a
    c=sum(gene_pathway_matrix[,fisher_pathway])-a
    d=length(all_genes)-a-c-b
    contigency_table=matrix(c(a,b,c,d),nrow = 2)
    fisher_result=fisher.test(contigency_table)
    selected_pathways_fisher_pvalue[fisher_pathway]=fisher_result[['p.value']]
    selected_pathways_num_genes[fisher_pathway]=sum(gene_pathway_matrix[,fisher_pathway])
  }
  return(list(selected_pathways_fisher_pvalue=selected_pathways_fisher_pvalue,selected_pathways_num_genes=selected_pathways_num_genes))
}


# pathway2coeff_plot=function(selected_pathway_names,selected_coef,pathway_colors){
#   pathway2coeff_df=data.frame(pathway=selected_pathway_names,coeff=selected_coef,names=selected_pathway_names)
#   pathway2coeff_p<-ggplot(data=pathway2coeff_df, aes(x=pathway, y=coeff,fill=names)) +
#     geom_bar(stat="identity")+
#     #geom_text(aes(label=selected_coef), vjust=0, hjust=0,size=1)+
#     coord_flip()+
#     scale_x_discrete(limits=selected_pathway_names[order(selected_coef)])+                # !!!!!!!!!
#     scale_fill_manual(values=pathway_colors[selected_pathway_names[order(selected_coef)]])+
#     ylab("regression coefficients")+
#     theme(axis.text.y = element_text(size=20),
#           axis.text.x = element_text(size=20),
#           legend.position="none",
#           axis.title.x = element_text(size=20),
#           axis.title.y = element_text(size=20))
# }

pathway2coeff_plot=function(selected_pathway_names,selected_coef,pathway_colors){
  pathway2coeff_df=data.frame(pathway=selected_pathway_names,coeff=selected_coef,names=selected_pathway_names)
  pathway2coeff_p<-ggplot(data=pathway2coeff_df, aes(x=pathway, y=coeff,fill=names)) +
    geom_bar(stat="identity")+
    #geom_text(aes(label=selected_coef), vjust=0, hjust=0,size=1)+
    coord_flip()+
    scale_x_discrete(limits=selected_pathway_names[order(selected_coef)])+                # !!!!!!!!!
    scale_fill_manual(values=pathway_colors[selected_pathway_names[order(selected_coef)]])+
    ylab("regression coefficients")+
    theme(legend.position="none")
}


pathway2fisher_plot=function(selected_pathway_names,log_p_value,pathway_colors){
  pathway2fisher_df=data.frame(pathway=selected_pathway_names,coeff=log_p_value,names=selected_pathway_names)
  pathway2fisher_p<-ggplot(data=pathway2fisher_df, aes(x=pathway, y=coeff,fill=names)) +
    geom_bar(stat="identity")+
    #geom_text(aes(label=log_p_value), vjust=0, hjust=0,size=1)+
    coord_flip()+
    scale_x_discrete(limits=selected_pathway_names[order(log_p_value)])+
    scale_fill_manual(values=pathway_colors[selected_pathway_names[order(log_p_value)]])+
    ylab("negative log fisher exact test pvalues")+
    theme(legend.position="none")
}


pathway2gene_plot_new=function(gene_pathway_matrix,selected_coef,selected_pathways,selected_pathway_names,selected_pathways_num_genes,s5_module_genes){
  pathway2gene=matrix(0,nrow = length(selected_pathway_names),ncol = length(s5_module_genes))
  mode(pathway2gene) <- "integer"
  rownames(pathway2gene)=selected_pathway_names
  colnames(pathway2gene)=s5_module_genes
  for (j in 1:length(selected_pathway_names)) {
    module_gene_in_pathway=intersect(all_genes[which(gene_pathway_matrix[,selected_pathways[j]]==1)],s5_module_genes)
    pathway2gene[selected_pathway_names[j],module_gene_in_pathway]=1
  }
  
  if(length(selected_coef)>1){
    # pathway2gene=pathway2gene[order(as.numeric(selected_coef),decreasing = TRUE),]
    pathway2gene <- as.matrix(pathway2gene)
    rtmatColSum <- colSums(pathway2gene)
    pathway2gene <- pathway2gene[, order(rtmatColSum, decreasing = TRUE)]
    ## myheat <- function(x,...) biosHeatmap(x, Rowv=FALSE, Colv=FALSE, dendro="none", col="royalbluered", labRow=NA, labCol=NA,...)
    p2gRowOrd <- cascadeOrder(pathway2gene)
    p2gColOrd <- cascadeOrder(t(pathway2gene))
    pathway2gene_reorder=pathway2gene[p2gRowOrd, ]    #module_index=9
    #pathway2gene_reorder=pathway2gene[p2gRowOrd, p2gColOrd]    #module_index=9
    #write.csv(pathway2gene,file = "cascadeOrder_example.csv",row.names = TRUE,col.names = TRUE)
  }else{
    pathway2gene_reorder=pathway2gene
    p2gRowOrd=1
  }
  
  
  # selected_pathways_with_gene_num=vector()
  # selected_pathway_names_newOrder=selected_pathway_names[order(as.numeric(selected_coef),decreasing = TRUE)]
  # selected_pathways_num_genes_newOrder=selected_pathways_num_genes[order(as.numeric(selected_coef),decreasing = TRUE)]
  selected_pathways_with_gene_num=vector()
  selected_pathway_names_newOrder=rownames(pathway2gene_reorder)
  selected_pathways_num_genes_newOrder=selected_pathways_num_genes[p2gRowOrd]
  
  for (i in 1:length(selected_pathway_names)) {
    selected_pathways_with_gene_num[i]=paste(selected_pathway_names_newOrder[i],"(",selected_pathways_num_genes_newOrder[i],")",sep = "")
  }
  rownames(pathway2gene_reorder)=selected_pathways_with_gene_num
  
  if(nrow(pathway2gene_reorder)>1){
    melted_pathway2gene=melt(pathway2gene_reorder[rev(rownames(pathway2gene_reorder)),])
    melted_pathway2gene$value=factor(melted_pathway2gene$value)
  }else{
    melted_pathway2gene=data.frame(Var1=rep(rownames(pathway2gene_reorder),ncol(pathway2gene_reorder)),Var2=colnames(pathway2gene_reorder),value=pathway2gene_reorder[,])
    melted_pathway2gene$value=factor(melted_pathway2gene$value)
  }
  
  #melted_pathway2gene=melt(pathway2gene)
  pathway2gene_p=ggplot(data = melted_pathway2gene, aes(x=Var2, y=Var1, fill=value)) +
    geom_tile(color="black") +
    ggtitle("Module genes in pathways")+
    xlab("genes in the module")+
    ylab("selected pathways")+
    scale_fill_manual(values=c("white","blue"))+
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=0.5),
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
  
  return(pathway2gene_p)
}


pathway2gene_plot=function(gene_pathway_matrix,selected_pathways,selected_pathway_names,selected_pathways_num_genes,s5_module_genes){
  pathway2gene=matrix(0,nrow = length(selected_pathway_names),ncol = length(s5_module_genes))
  rownames(pathway2gene)=selected_pathway_names
  colnames(pathway2gene)=s5_module_genes
  for (j in 1:length(selected_pathway_names)) {
    module_gene_in_pathway=intersect(all_genes[which(gene_pathway_matrix[,selected_pathways[j]]==1)],s5_module_genes)
    pathway2gene[selected_pathway_names[j],module_gene_in_pathway]=1
  }
  
  selected_pathways_with_gene_num=vector()
  for (i in 1:length(selected_pathway_names)) {
    selected_pathways_with_gene_num[i]=paste(selected_pathway_names[i],"(",selected_pathways_num_genes[i],")",sep = "")
  }
  rownames(pathway2gene)=selected_pathways_with_gene_num
  
  melted_pathway2gene=melt(pathway2gene)
  
  pathway2gene_p=ggplot(data = melted_pathway2gene, aes(x=Var2, y=Var1, fill=value)) +
    geom_tile(color='white') +
    xlab("genes in the module")+
    ylab("selected pathways")+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  
  return(pathway2gene_p)
}

pathway2gene_plot_ribios=function(gene_pathway_matrix,selected_pathways,selected_pathway_names,selected_pathways_num_genes,s5_module_genes){
  pathway2gene=matrix(0,nrow = length(selected_pathway_names),ncol = length(s5_module_genes))
  rownames(pathway2gene)=selected_pathway_names
  colnames(pathway2gene)=s5_module_genes
  for (j in 1:length(selected_pathway_names)) {
    module_gene_in_pathway=intersect(all_genes[which(gene_pathway_matrix[,selected_pathways[j]]==1)],s5_module_genes)
    pathway2gene[selected_pathway_names[j],module_gene_in_pathway]=1
  }
  
  jpeg("GO_bp_module_heatmap.jpg",quality = 100,pointsize = 30,width = 2000,height = 1500)
  rtmat <- as.matrix(pathway2gene)
  rtmatColSum <- colSums(rtmat)
  rtmat <- rtmat[, order(rtmatColSum, decreasing = TRUE)]
  library(ribiosPlot)
  rtmatCascade <- cascadeOrder(rtmat)
  biosHeatmap(rtmat[rtmatCascade,], Colv=FALSE, Rowv=FALSE,
              col="blackyellow",
              xlab="genes in the module",
              ylab="selected pathways",
              cexCol=1.25,
              trace = "both",vline = FALSE)
  
  dev.off()
  
  return(pathway2gene_p)
}

########### node pathway barplot ################
plot_reactome_statistics=function(results){
  plot_reactome=results[grepl("HSA",results[,"Go_Reactome_root_id"]),"Go_Reactome_root_names"]
  plot_reactome=unlist(sapply(plot_reactome, function(x){
    strsplit(x,split = "#")
  }))
  unique_plot_reactome=unique(plot_reactome)
  unique_plot_reactome_count=sapply(unique_plot_reactome,function(x){
    sum(plot_reactome==x)
  } )
  
  plot_reactome_df=data.frame(pathways=unique_plot_reactome,pathway_count=unique_plot_reactome_count,names=unique_plot_reactome)
  plot_reactome_p<-ggplot(data=plot_reactome_df, aes(x=pathways, y=pathway_count,fill=names)) +
    geom_bar(stat="identity")+
    coord_flip()+
    scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(length(unique_plot_reactome)))+
    theme(legend.position="none")
  
  return(plot_reactome_p)
}

#################################################
plot_go_statistics=function(results){
  plot_go=results[grepl("GO",results[,"Go_Reactome_root_id"]),"Go_Reactome_root_names"]
  plot_go=unlist(sapply(plot_go, function(x){
    strsplit(x,split = "#")
  }))
  unique_plot_go=unique(plot_go)
  unique_plot_go_count=sapply(unique_plot_go,function(x){
    sum(plot_go==x)
  } )
  
  plot_go_df=data.frame(pathways=unique_plot_go,pathway_count=unique_plot_go_count,names=unique_plot_go)
  plot_go_p<-ggplot(data=plot_go_df, aes(x=pathways, y=pathway_count,fill=names)) +
    geom_bar(stat="identity")+
    coord_flip()+
    scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(length(unique_plot_go)))+
    theme(legend.position="none")
  
  return(plot_go_p)
}



get_subgraphs=function(results){
  module_subgraphs=list()
  for (i in 1:nrow(results)){
    print(i)
    if(grepl("GO",rownames(results)[i])){
      module_subgraphs[[rownames(results)[i]]]=induced.subgraph(go_ontology,names(unlist(ego(go_ontology,order = length(V(go_ontology)),rownames(results)[i],mode = "in"))))
    }else{
      module_subgraphs[[rownames(results)[i]]]=induced.subgraph(reactome_ontology,names(unlist(ego(reactome_ontology,order = length(V(reactome_ontology)),rownames(results)[i],mode = "in"))))
    }
    
  }
  return(module_subgraphs)
}