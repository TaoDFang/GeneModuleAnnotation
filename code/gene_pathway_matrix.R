library(ggplot2)
library(reshape2)
library(Matrix)
library(ribiosMath)
library(plotly)

### using data from go_preprocess_tree.R and reactome_preprocess.R

#setwd("/media/sf_sharded_ubuntu/Disease_module_identification_DREAM_change")
setwd("/pstore/home/fangt3/Disease_module_identification_DREAM_change")

go_genesets_raw=scan("tao_go.gmt",what = "",sep = "\n")
go_genesets=strsplit(go_genesets_raw,split = "\t")


reactome_genesets_raw=scan("ReactomePathways.gmt",what = "",sep = "\n")
reactome_genesets=strsplit(reactome_genesets_raw,split = "\t")

go_all_genes=vector()
for (i in 1:length(go_genesets)){
  go_all_genes=union(go_all_genes,go_genesets[[i]][-c(1,2)])
}
save(go_all_genes,file = "go_all_genes.rda")

load("go_all_genes.rda")  #16146 unique genes

reactome_all_genes=vector()
for (i in 1:length(reactome_genesets)){
  reactome_all_genes=union(reactome_all_genes,reactome_genesets[[i]][-c(1,2,3)])
}
save(reactome_all_genes,file = "reactome_all_genes.rda")

load("reactome_all_genes.rda")  #10669 unique genes

go_all_pathways_id=vector()
for (i in 1:length(go_genesets)){
  go_all_pathways_id=union(go_all_pathways_id,go_genesets[[i]][1])
}
save(go_all_pathways_id,file = "go_all_pathways_id.rda")

load("go_all_pathways_id.rda")    #4464 pathways

reactome_all_pathways_id=vector()
for (i in 1:length(reactome_genesets)){
  reactome_all_pathways_id=union(reactome_all_pathways_id,reactome_genesets[[i]][2])
}
save(reactome_all_pathways_id,file = "reactome_all_pathways_id.rda")

load("reactome_all_pathways_id.rda")    #2029 pathways

library(Matrix)
library(glmnet)

go_gene_pathway_matrix=Matrix(0,nrow = length(go_all_genes),ncol = length(go_all_pathways_id))
rownames(go_gene_pathway_matrix)=go_all_genes
colnames(go_gene_pathway_matrix)=go_all_pathways_id
for(i in 1:length(go_all_pathways_id)){
  pathway_name=go_all_pathways_id[i]
  pathway_genes=go_genesets[[i]][-c(1,2)]
  go_gene_pathway_matrix[pathway_genes, pathway_name]=1
  print(i)
}
save(go_gene_pathway_matrix,file="go_gene_pathway_matrix.rda")

load("go_gene_pathway_matrix.rda")   #16146 * 4464


reactome_gene_pathway_matrix=Matrix(0,nrow = length(reactome_all_genes),ncol = length(reactome_all_pathways_id))
rownames(reactome_gene_pathway_matrix)=reactome_all_genes
colnames(reactome_gene_pathway_matrix)=reactome_all_pathways_id
for(i in 1:length(reactome_all_pathways_id)){
  pathway_name=reactome_all_pathways_id[i]
  pathway_genes=reactome_genesets[[i]][-c(1,2,3)]
  reactome_gene_pathway_matrix[pathway_genes, pathway_name]=1
  print(i)
}
save(reactome_gene_pathway_matrix,file="reactome_gene_pathway_matrix.rda")

load("reactome_gene_pathway_matrix.rda")    # 10669 * 2029



all_genes=union(go_all_genes,reactome_all_genes)   #17176
#write.csv(all_genes,file="ModuleIdentification_all_genes.csv")
all_pathways=union(go_all_pathways_id,reactome_all_pathways_id)  #6493


gene_pathway_matrix=Matrix(0,nrow = length(all_genes),ncol = length(all_pathways))
rownames(gene_pathway_matrix)=all_genes
colnames(gene_pathway_matrix)=all_pathways
for(i in 1:length(all_pathways)){
  pathway_name=all_pathways[i]
  if(pathway_name %in% go_all_pathways_id){
    pathway_genes=rownames(go_gene_pathway_matrix)[go_gene_pathway_matrix[,pathway_name]==1]
  }else{
    pathway_genes=rownames(reactome_gene_pathway_matrix)[reactome_gene_pathway_matrix[,pathway_name]==1]
  }
  gene_pathway_matrix[pathway_genes,pathway_name]=1
  print(i)
}
save(gene_pathway_matrix,file="gene_pathway_matrix.rda")














load('go_gene_pathway_matrix_old.rda')
load("gene_pathway_matrix.rda")     # 17176 6493
all_genes=rownames(gene_pathway_matrix)
all_pathways=colnames(gene_pathway_matrix)

common_genes=intersect(genes,all_genes)        # 14904
gene_pathway_matrix=gene_pathway_matrix[common_genes,]
all_genes=rownames(gene_pathway_matrix)
all_pathways=colnames(gene_pathway_matrix)


# from go_preprocess_tree.R in virtural machine
library(data.table)
#library(ribiosAnnotation)
library(ontologyIndex)
library(RamiGO)
library(ontoCAT)
#library(Ramigo)

#setwd("/media/sf_sharded_ubuntu/Disease_module_identification_DREAM_change")
setwd("/pstore/home/fangt3/Disease_module_identification_DREAM_change")

go_gaf_filename="goa_human.gaf"
go_gaf=fread(go_gaf_filename,header = FALSE,skip = 23,fill=TRUE, na.strings="NA",stringsAsFactors = FALSE)
go_gaf_simple=go_gaf[,c(3,5)]
go_gaf_simple=go_gaf_simple[which(go_gaf_simple[,1] !=""),]
go_gaf_simple=as.data.frame(go_gaf_simple)

#overlap_genes=intersect(fdata$GeneSymbol,go_gaf_simple[,1])   #11163/11515

all_go_ids=unique(go_gaf_simple[,2])
go_list=vector("list",length = length(all_go_ids))
names(go_list)=all_go_ids

for (i in 1:length(all_go_ids)){
  go_list[[i]]=go_gaf_simple[go_gaf_simple[,2]==all_go_ids[i],1]
}

saveRDS(go_list, "raw_go_list.rds")

