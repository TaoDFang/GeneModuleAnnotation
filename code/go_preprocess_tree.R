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

test=go_list
rmGO_index=vector()
for (i in 1:length(all_go_ids)){
  if (length(go_list[[i]])<10 || length(go_list[[i]])>200){
    rmGO_index=c(rmGO_index,i)
  }
}

go_list=go_list[-(rmGO_index)]

saveRDS(go_list, "go_list.rds")
go_list <- readRDS("go_list.rds")
go_list_names=names(go_list)


go_ontology_filename="go-basic.obo"
go_ontology_obo=get_OBO(go_ontology_filename,propagate_relationships = "is_a",extract_tags = "minimal")
for (j in 1:length(go_list)){
  go_name=go_ontology_obo$name[go_ontology_obo$id==go_list_names[j]]
  go_list[[j]]=c(go_name,go_list[[j]])
}


saveRDS(go_list, "go_list_gmt.rds")
go_list_gmt <- readRDS("go_list_gmt.rds")


for (i in seq_along(go_list_gmt) ){ 
  cat(names(go_list_gmt[i]), go_list_gmt[[i]], file="tao_go.gmt", append=TRUE, sep = "\t")
  cat("\n", append=TRUE, file="tao_go.gmt")
}



go_ontology_filename="go-basic.obo"
go_ontology_obo=get_OBO(go_ontology_filename,propagate_relationships = "is_a",extract_tags = "minimal")
go_ontology_obo_children=go_ontology_obo$children
go_ontology_obo_list=list()
count=1
for (i in 1:length(go_ontology_obo_children)) {
  if(length(go_ontology_obo_children[[i]])>0 && grepl("GO",names(go_ontology_obo_children[i]))){
    for (j in 1:length(go_ontology_obo_children[[i]])) {
      go_ontology_obo_list[[count]]=list(names(go_ontology_obo_children[i]),go_ontology_obo_children[[i]][j])
      count=count+1
    }
  }
}
go_ontology_obo_frame=matrix(nrow = length(go_ontology_obo_list),ncol = 2)
for (i in 1:length(go_ontology_obo_list)) {
  go_ontology_obo_frame[i,]=unlist(go_ontology_obo_list[[i]])
}

go_ontology=as.matrix(go_ontology_obo_frame)
go_ontology=graph_from_edgelist(go_ontology,directed = TRUE)

go_ontology_obo_pathway_names=go_ontology_obo$name
names(go_ontology_obo_pathway_names)=go_ontology_obo$id

V(go_ontology)$pathway_names=go_ontology_obo_pathway_names[match(as_ids(V(go_ontology)),names(go_ontology_obo_pathway_names))]
#V(human_reactome_ontology)$pathway_names=reactome_pathwaysName[ match(as_ids(V(human_reactome_ontology)),reactome_pathwaysName[,1]),2]

root_nodes=which(sapply(sapply(V(go_ontology), function(x) neighbors(go_ontology,x, mode="in")), length) == 0)
go_ontology_subgraphs=list()
for(i in 1:length(root_nodes)){
  a=dfs(go_ontology, names(root_nodes)[i], neimode = "out",
        unreachable = FALSE, order = TRUE, order.out = FALSE, father = TRUE,
        dist = TRUE)
  order_node_index=!is.na(a$order)
  order_node=a$order[order_node_index]
  go_ontology_subgraphs[[names(root_nodes)[i]]]=induced.subgraph(go_ontology,order_node)
}

saveRDS(go_ontology, "go_ontology.rds")
saveRDS(go_ontology_subgraphs, "go_ontology_subgraphs.rds")








library(cogena)
go_tao_genesets=gmt2list("tao_go.gmt")

# test ontoCAT package
go_ontology_filename="/media/sf_sharded_ubuntu/Disease_module_identification_DREAM_change/go-basic.obo"
go_ontology=getOntology(go_ontology_filename)
getAllTerms(go_ontology)
getTermNameById(go_ontology,"GO:2000813")
getTermChildrenById(go_ontology,"GO:2000813")
getTermParentsById(go_ontology,"GO:2000813")

# GO:0008150 bp root node, GO:0003674 mf root node, GO:0005575 cc root node
getTermChildrenById(go_ontology,"GO:0003674")
isRootById(go_ontology,"GO:0008150")
showHierarchyDownToTermById(go_ontology,"GO:2000813")
