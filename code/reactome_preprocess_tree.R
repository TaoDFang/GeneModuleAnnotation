library(igraph)
library()

setwd("/pstore/home/fangt3/Disease_module_identification_DREAM_change")

reactome_relations_name="ReactomePathwaysRelation.txt"
reactome_relations=read.csv(reactome_relations_name,header = FALSE,sep = "\t")

# check the pathway overlapping between hierarchy tree and gmt file

human_reactome_relations=reactome_relations[grep("HSA",reactome_relations[,1]),]
human_reactome_relations=as.matrix(human_reactome_relations)
human_reactome_ontology=graph_from_edgelist(human_reactome_relations,directed = TRUE)

# get pathwaynames from reactome pathway id
reactome_pathwaysName=read.csv("ReactomePathways.txt",header = FALSE,sep = "\t",stringsAsFactors = FALSE)
V(human_reactome_ontology)$pathway_names=reactome_pathwaysName[ match(as_ids(V(human_reactome_ontology)),reactome_pathwaysName[,1]),2]

#plot(reactome_ontology)
#human_reactome_components=components(human_reactome_ontology, mode ="weak")
#count_components(reactome_ontology, mode ="weak")

#check root nodes
root_nodes=which(sapply(sapply(V(human_reactome_ontology), function(x) neighbors(human_reactome_ontology,x, mode="in")), length) == 0)

human_reactome_subgraphs=list()

for(i in 1:length(root_nodes)){
  a=dfs(human_reactome_ontology, names(root_nodes)[i], neimode = "out",
        unreachable = FALSE, order = TRUE, order.out = FALSE, father = TRUE,
        dist = TRUE)
  order_node_index=!is.na(a$order)
  order_node=a$order[order_node_index]
  human_reactome_subgraphs[[names(root_nodes)[i]]]=induced.subgraph(human_reactome_ontology,order_node)
}


saveRDS(human_reactome_ontology, "human_reactome_ontology.rds")
saveRDS(human_reactome_subgraphs, "human_reactome_subgraphs.rds")
#pdf(file="test.pdf")
#plot(human_reactome_subgraphs[[1]])
#dev.off()
