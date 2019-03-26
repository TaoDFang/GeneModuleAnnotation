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
  
  
  
  ## check here again !!!!!!!!!!!!!!!!!!!!!!!very important 
  # common_genes=intersect(network_all_genes,all_genes)        # 16637
  # gene_pathway_matrix=gene_pathway_matrix[common_genes,]
  # all_genes=rownames(gene_pathway_matrix)
  # all_pathways=colnames(gene_pathway_matrix)
  
  
  
  
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
        
        #https://www.ncbi.nlm.nih.gov/gene?term=(als2cr8[gene])%20AND%20(Homo%20sapiens[orgn])%20AND%20alive[prop]%20NOT%20newentry[gene]&sort=weight})
        #gene_url=paste("https://www.ncbi.nlm.nih.gov/gene?term=(",module_genes_metaInfo[,"Module Gene Symbol"][2],"[gene])%20AND%20(Homo%20sapiens[orgn])%20AND%20alive[prop]%20NOT%20newentry[gene]&sort=weight",sep = "")
        #gene_link=paste("<a href=\"",gene_url,"\">",module_genes_metaInfo[,"Module Gene Symbol"][2],"</a></p>",sep = "")
        #https://www.ncbi.nlm.nih.gov/gene?term=(als2cr8[gene])%20AND%20(Homo%20sapiens[orgn])%20AND%20alive[prop]%20NOT%20newentry[gene]&sort=weight
        #https://www.ncbi.nlm.nih.gov/gene?term=(cacul1[gene])%20AND%20(Homo%20sapiens[orgn])%20AND%20alive[prop]%20NOT%20newentry[gene]&sort=weight
        #################################################################
        
        most_specific_annotations=matrix(nrow = length(selected_pathways),ncol = 5)
        #rownames(most_specific_annotations)=selected_pathways
        colnames(most_specific_annotations)=c("Coeff","P-value","Source","Term","Class")
        most_specific_annotations[,"Coeff"]=format(selected_coef,scientific=TRUE,digits=3)
        most_specific_annotations[,"P-value"]=format(selected_pathways_fisher_pvalue,scientific=TRUE,digits=3)
        most_specific_annotations[,"Source"]=unlist(find_sources(selected_pathways))
        most_specific_annotations[,"Term"]=unlist(find_pathway_links(selected_pathways ,selected_pathway_names))
        
        most_specific_annotations[,"Class"]= unlist(find_classes(selected_pathways,go_ontology,reactome_ontology,reactome_roots,reactome_ontology_names))
        if(nrow(most_specific_annotations)>1){
          most_specific_annotations=most_specific_annotations[order(as.numeric(most_specific_annotations[,"Coeff"]),decreasing = TRUE),]
        }
        #<p>This text contains <sup>superscript</sup> text.</p>
        most_specific_annotations_newColnames=sapply(1:ncol(most_specific_annotations), function(j){paste(colnames(most_specific_annotations)[j],".<sup>",j,"</sup>",sep = "")})
        colnames(most_specific_annotations)=most_specific_annotations_newColnames  
  ###################################################################
       
       # pathway2gene_p=pathway2gene_plot_new(gene_pathway_matrix,selected_coef,selected_pathways,selected_pathway_names,selected_pathways_num_genes,module_common_genes)
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
        
        
        pathway2gene_p=pathway2gene_plot_new(gene_pathway_matrix,selected_coef,selected_pathways,selected_pathway_names,selected_pathways_num_genes,module_common_genes)
        
       
        ggsave(paste("HTML_results/",module_name,"_pathway2gene.PNG",sep=""),  width =20 ,height =3.5, 
               plot = pathway2gene_p)
        
      }
    }
    
    
    
    ############################### HTML plot#################################
    
    #####################plot################################
    module_metaInfo=matrix(nrow = 4,ncol = 1)
    rownames(module_metaInfo)=c("Assigned name","Network","Module ID","Module size")
    module_metaInfo["Network",]=strsplit(module_name,split="_")[[1]][1]
    module_metaInfo["Module ID",]=module_name
    module_metaInfo["Module size",]=paste(length(module_genes),"genes",sep = " ")
    ###############################################################
    rabiosGeneInfo=annotateAnyIDs(module_genes)
    module_genes_metaInfo=rabiosGeneInfo[module_genes,c("GeneID","Input","GeneName")]
    rownames(module_genes_metaInfo)=module_genes
    colnames(module_genes_metaInfo)=c("Gene ID","Gene Symbol","Gene Name")
    module_genes_metaInfo=module_genes_metaInfo[order(module_genes_metaInfo[,"Gene Symbol"]),]
    #https://www.ncbi.nlm.nih.gov/gene/?term=79800
    gene_urls=sapply(module_genes_metaInfo[,"Gene ID"], function(x){paste("https://www.ncbi.nlm.nih.gov/gene/?term=",x,sep = "")})
    gene_links=sapply(1:length(gene_urls), function(i){paste("<a href=\"",gene_urls[i],"\">",module_genes_metaInfo[,"Gene ID"][i],"</a></p>",sep = "")})
    #gene_urls=sapply(module_genes_metaInfo[,"Module Gene Symbol"], function(x){paste("https://www.ncbi.nlm.nih.gov/gene?term=(",x,"[gene])%20AND%20(Homo%20sapiens[orgn])%20AND%20alive[prop]%20NOT%20newentry[gene]&sort=weight",sep = "")})
    #gene_links=sapply(1:length(gene_urls), function(i){paste("<a href=\"",gene_urls[i],"\">",module_genes_metaInfo[,"Module Gene Symbol"][i],"</a></p>",sep = "")})
    module_genes_metaInfo[,"Gene ID"]=gene_links  
    
    HTMLStart(outdir="./HTML_results/", file=module_name,
              extension="html", echo=FALSE,HTMLframe = FALSE,CSSFile = "R2HTML_tao.css")
    #[PPI-STRING_Consensus_mod1](./htmls/PPI-STRING_Consensus_mod1.html)
    # #<a href="./htmls/PPI-STRING_Consensus_mod347.html">PPI-STRING_Consensus_mod347</a><br />
    back_to_main='<a href="../index.html">Back to main page</a>'
    HTML(back_to_main,CSSFile = "R2HTML_tao.css")
    HTML.title("DREAM Module Identification Challenge – Consensus modules", HR=1,CSSFile = "R2HTML_tao.css")
    
    HTML.title(module_metaInfo["Module ID",], HR=2,CSSFile = "R2HTML_tao.css")
    HTML(module_metaInfo,CSSFile = "R2HTML_tao.css")
    
    HTML.title("Module genes", HR=2,CSSFile = "R2HTML_tao.css")
    HTML.title("This module comprises the following genes:", HR=3,CSSFile = "R2HTML_tao.css")
    HTML(module_genes_metaInfo,CSSFile = "R2HTML_tao.css",row.names = FALSE)
    
    HTML.title("Functional annotation", HR=2,CSSFile = "R2HTML_tao.css")
    HTML("<p>Modules were tested for enrichment in functional and pathway annotations using two complementary approaches:
         <br>    1. To select a small number of specific / non-redundant annotations for each module, a regression-based approach was used;
         <br>    2. To obtain the complete set of enriched annotations, an extension of Fisher’s exact test that takes annotation bias into account was employed (Wallenius’ non-central hypergeometric distribution).
         </p>",CSSFile = "R2HTML_tao.css")
    
    if(length(selected_pathways)>0){
      HTML.title("Most specific annotations for this module", HR=3,CSSFile = "R2HTML_tao.css")
      HTML(most_specific_annotations,CSSFile = "R2HTML_tao.css",row.names=FALSE)
      #<p>This text contains <sup>superscript</sup> text.</p>
      HTML("<p><sup>1</sup>Regression coefficient</p>",CSSFile = "R2HTML_tao.css")
      HTML("<p><sup>2</sup>Fisher’s exact test nominal P-value</p>",CSSFile = "R2HTML_tao.css")
      HTML("<p><sup>3</sup>Annotation source (Reactome, GO biological process (BP), molecular function (MF) and cellular component (CC))</p>",CSSFile = "R2HTML_tao.css")
      HTML("<p><sup>4</sup>GO category or Reactome pathway</p>",CSSFile = "R2HTML_tao.css")
      HTML("<p><sup>5</sup>High-level branch of annotation tree</p>",CSSFile = "R2HTML_tao.css")
      
      
      HTML.title("Gene membership", HR=3,CSSFile = "R2HTML_tao.css")
      HTMLInsertGraph(paste(module_name,"_pathway2gene.PNG",sep=""),file =paste("./HTML_results/",module_name,".html",sep = "") ,
                      WidthHTML =1200 )
    }else{
      
      HTML("<p>No specific annotations were found for this module</p>",CSSFile = "R2HTML_tao.css")
    }
    
    
    
    HTML.title("All enriched annotations", HR=3,CSSFile = "R2HTML_tao.css")
    
    if(nrow(Daniel_module_results)){
      ############################################################
      Daniel_module_GO_results=Daniel_module_results[grepl("GO",Daniel_module_results$termId),]
      if(nrow(Daniel_module_GO_results)){
        Daniel_module_GO_results=Daniel_module_GO_results[Daniel_module_GO_results$P.noncentral<0.05,]
        Daniel_module_GO_results=Daniel_module_GO_results[order(Daniel_module_GO_results$P.noncentral),]
        K=20
        if(nrow(Daniel_module_GO_results)>K){
          Daniel_module_GO_results=Daniel_module_GO_results[1:20,]
        }
        
        Daniel_module_GO_links=lapply(1:length(Daniel_module_GO_results$termId), function(i){
          url=paste("http://amigo.geneontology.org/amigo/medial_search?q=GO%3A",strsplit(Daniel_module_GO_results$termId[i],split=":")[[1]][2],sep = "")
          links=paste("<a href=\"",url,"\">",Daniel_module_GO_results$term[i],"</a></p>",sep = "")
        })
        Daniel_module_GO_results$term=unlist(Daniel_module_GO_links)
        Daniel_module_GO_results_2html=Daniel_module_GO_results[,c("P.noncentral","P.noncentral.fdr","term" )]
        Daniel_module_GO_results_2html=format(Daniel_module_GO_results_2html,scientific=TRUE,digits=3)
        colnames(Daniel_module_GO_results_2html)=c(paste("P-value","<sup>1<sup>",sep = ""),paste("FDR","<sup>2<sup>",sep = ""),"Term")
        ##########################################################
        HTML.title("Gene Ontology", HR=4,CSSFile = "R2HTML_tao.css")
        HTML(Daniel_module_GO_results_2html,CSSFile = "R2HTML_tao.css",row.names=FALSE)
        HTML("<p><sup>1</sup>1Nominal enrichment p-value (Wallenius’ noncentral hypergeometric distribution)
             </p>",CSSFile = "R2HTML_tao.css")
        HTML("<p><sup>2</sup>FDR corrected p-value (Benjamini-Hochberg)
             ",CSSFile = "R2HTML_tao.css")
      }
      
      Daniel_module_REACTOME_results=Daniel_module_results[grepl("REACTOME",Daniel_module_results$term),]
      if(nrow(Daniel_module_REACTOME_results)){
        Daniel_module_REACTOME_results=Daniel_module_REACTOME_results[Daniel_module_REACTOME_results$P.noncentral<0.05,]
        Daniel_module_REACTOME_results=Daniel_module_REACTOME_results[order(Daniel_module_REACTOME_results$P.noncentral),]
        K1=20
        if(nrow(Daniel_module_REACTOME_results)>K1){
          Daniel_module_REACTOME_results=Daniel_module_REACTOME_results[1:20,]
        }
        
        Daniel_module_REACTOME_results_2html=Daniel_module_REACTOME_results[,c("P.noncentral","P.noncentral.fdr","term" )]
        reactome_terms=gsub("REACTOME_","",Daniel_module_REACTOME_results_2html$term)
        reactome_terms=gsub("_"," ",reactome_terms)
        Daniel_module_REACTOME_results_2html$term=reactome_terms
        Daniel_module_REACTOME_results_2html=format(Daniel_module_REACTOME_results_2html,scientific=TRUE,digits=3)
        colnames(Daniel_module_REACTOME_results_2html)=c(paste("P-value","<sup>1<sup>",sep = ""),paste("FDR","<sup>2<sup>",sep = ""),"Term")
        ##########################################################
        HTML.title("Reactome", HR=4,CSSFile = "R2HTML_tao.css")
        HTML(Daniel_module_REACTOME_results_2html,CSSFile = "R2HTML_tao.css",row.names=FALSE)
        HTML("<p><sup>1</sup>1Nominal enrichment p-value (Wallenius’ noncentral hypergeometric distribution)
</p>",CSSFile = "R2HTML_tao.css")
        HTML("<p><sup>2</sup>FDR corrected p-value (Benjamini-Hochberg)
",CSSFile = "R2HTML_tao.css")
      }
      
      
      Daniel_module_MP_results=Daniel_module_results[grepl("MP",Daniel_module_results$termId),]
      if(nrow(Daniel_module_MP_results)){
        Daniel_module_MP_results=Daniel_module_MP_results[Daniel_module_MP_results$P.noncentral<0.05,]
        Daniel_module_MP_results=Daniel_module_MP_results[order(Daniel_module_MP_results$P.noncentral),]
        K=20
        if(nrow(Daniel_module_MP_results)>K){
          Daniel_module_MP_results=Daniel_module_MP_results[1:20,]
        }
        
        #http://www.informatics.jax.org/vocab/mp_ontology/MP:0005124
        Daniel_module_MP_links=lapply(1:length(Daniel_module_MP_results$termId), function(i){
          url=paste("http://www.informatics.jax.org/vocab/mp_ontology/",Daniel_module_MP_results$termId[i],sep = "")
          links=paste("<a href=\"",url,"\">",Daniel_module_MP_results$term[i],"</a></p>",sep = "")
        })
        Daniel_module_MP_results$term=unlist(Daniel_module_MP_links)
        Daniel_module_MP_results_2html=Daniel_module_MP_results[,c("P.noncentral","P.noncentral.fdr","term" )]
        Daniel_module_MP_results_2html=format(Daniel_module_MP_results_2html,scientific=TRUE,digits=3)
        colnames(Daniel_module_MP_results_2html)=c(paste("P-value","<sup>1<sup>",sep = ""),paste("FDR","<sup>2<sup>",sep = ""),"Term")
        ##########################################################
        HTML.title("Mouse mutant phenotypes", HR=4,CSSFile = "R2HTML_tao.css")
        HTML(Daniel_module_MP_results_2html,CSSFile = "R2HTML_tao.css",row.names=FALSE)
        HTML("<p><sup>1</sup>1Nominal enrichment p-value (Wallenius’ noncentral hypergeometric distribution)
                 </p>",CSSFile = "R2HTML_tao.css")
        HTML("<p><sup>2</sup>FDR corrected p-value (Benjamini-Hochberg)
                 ",CSSFile = "R2HTML_tao.css")
      }
      
      
      }else{
        HTML("<p>No enriched annotations found for this module</p>",CSSFile = "R2HTML_tao.css")
      }
    
    HTMLStop()
    #####################plot################################
    
  }
  
  
}

write.csv(dream_consensus_modules_results, file = "dream_consensus_modules_results_alpha0.5.csv",
          row.names=FALSE, col.names = TRUE)
