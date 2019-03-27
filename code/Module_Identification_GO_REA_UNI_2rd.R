library(ggplot2)
library(reshape2)
library(Matrix)
library(ribiosMath)
library(plotly)
library(Matrix)
library(glmnet)
library(R2HTML)
library(igraph)
library(RColorBrewer)


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

# read module information
raw_modules=scan("DSDSpectral.consensus_top50.1_ppi_anonym_v2.txt",what = "",sep = "\n")
modules_list=strsplit(raw_modules,"[[:space:]]+")
load("gene2module.rda")
genes=names(gene2module)     #17397 unique genes


load("gene_pathway_matrix.rda")     # 17176 6494
all_genes=rownames(gene_pathway_matrix)
all_pathways=colnames(gene_pathway_matrix)

common_genes=intersect(genes,all_genes)        # 14904
gene_pathway_matrix=gene_pathway_matrix[common_genes,]
all_genes=rownames(gene_pathway_matrix)
all_pathways=colnames(gene_pathway_matrix)


#s5_module_genes=dream_genesets[[1]]
s5_module_genes=c("TRPC4AP","CDC37","TNIP1","IKBKB","NKIRAS2","NFKBIA","TIMM50","RELB","TNFAIP3","NFKBIB","HSPA1A","NFKBIE",
                  "SPAG9","NFKB2","ERLIN1","REL","TNIP2","TUBB6","MAP3K8")

module1_label=rep(0,length(all_genes))          #len:20244
names(module1_label)=all_genes
#module1_common_genes=intersect(all_genes,modules_list[[1]][c(-1,-2)])       # len: 33
module1_common_genes=intersect(all_genes,s5_module_genes) 
#module1_common_genes=intersect(all_genes,genecodis_module[[1]]) 
### !!! need to check how many common genes are in this specific modules
### !!! neeed to integrate more data 
print("number of genes in the specific module: \n")
print(length(s5_module_genes))
print("number of common genes: \n")
print(length(module1_common_genes))
module1_label[module1_common_genes]=1
#fit=glmnet(gene_pathway_matrix,module1_label,alpha = 1)
cvfit=glmnet(gene_pathway_matrix,module1_label,lambda = 0.007956622,alpha = 0.5)   #s5;0.007956622, dream[[1]]:
#cvfit=glmnet(gene_pathway_matrix,module1_label,family = "binomial",alpha = 0.5)
#cvfit=cv.glmnet(gene_pathway_matrix,module1_label,alpha = 1)
#cvfit=cv.glmnet(gene_pathway_matrix,module1_label,alpha = 0.5)
#cvfit=cv.glmnet(gene_pathway_matrix,module1_label,alpha = 0)
#cvfit=cv.glmnet(gene_pathway_matrix,module1_label,family = "binomial",alpha = 0.5)   # binomial model works so bad so may be it is not a good idea to just label 0 an 1
#save(cvfit,file='cvfit.rda')
#load('cvfit.rda')
plot(cvfit,xvar="lamda")
cvfit$lambda.min
coef=coef(cvfit, s = "lambda.min")
non0index=coef@i[-1]   #remove intercept
non0coef=coef@x[-1]
selected_index=non0index[which(non0coef>0)]
selected_pathways=all_pathways[selected_index]
selected_coef=non0coef[which(non0coef>0)]

results=matrix(nrow = length(selected_pathways),ncol = 7)
rownames(results)=selected_pathways
colnames(results)=c("pathway_name","fish_Pvalue","reg_coeff","GenesInPathway","Go_Reactome_root_id","Go_Reactome_root_names","step2root")


results[,"step2root"]=unlist(get_steps(selected_pathways,go_ontology,reactome_ontology,go_roots,reactome_roots))
results[,"Go_Reactome_root_id"]=unlist(find_root_ids(selected_pathways ))
results[,"Go_Reactome_root_names"]=unlist(find_root_names(selected_pathways ))


selected_pathway_names=from_id2name(selected_pathways )

names(selected_coef)=selected_pathway_names


fisher_exact_test_results=fisher_exact_test(selected_pathways )
selected_pathways_fisher_pvalue=fisher_exact_test_results$selected_pathways_fisher_pvalue
selected_pathways_num_genes=fisher_exact_test_results$selected_pathways_num_genes




results[,"pathway_name"]=selected_pathway_names
results[,"reg_coeff"]=format(selected_coef,scientific=TRUE,digits=3)
results[,"fish_Pvalue"]=format(selected_pathways_fisher_pvalue,scientific=TRUE,digits=3)


#color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
#pathway_colors=sample(color,length(selected_pathways))
#names(pathway_colors)=selected_pathway_names

pathway_colors=brewer.pal(length(selected_pathways),"Set3")
names(pathway_colors)=selected_pathway_names


pathway2coeff_p=pathway2coeff_plot(selected_pathway_names,selected_coef,pathway_colors)
#ggsave("GO_REA_UNI_test_pathway2coef.pdf", 
#       plot = pathway2coeff_p)



log_p_value=-log(selected_pathways_fisher_pvalue)
pathway2fisher_p=pathway2fisher_plot(selected_pathway_names,log_p_value,pathway_colors)
#ggsave("GO_REA_UNI_test_pathway2fisher.pdf", 
#       plot = pathway2fisher_p)



results[,"GenesInPathway"]=selected_pathways_num_genes
pathway2gene_p=pathway2gene_plot(gene_pathway_matrix,selected_pathways,selected_pathway_names,selected_pathways_num_genes,s5_module_genes)
#ggsave("GO_REA_UNI_test_pathway2gene.pdf", 
#       plot = pathway2gene_p)



module_subgraphs=get_subgraphs(results )
for (i in 1:length(module_subgraphs)) {
  jpeg(paste("./HTML_results/moduletest_",i,".PNG",sep = ""))
  plot(module_subgraphs[[i]])
  dev.off()
  
}


patterns=c("signaling","pathways")
patterns=paste(patterns,collapse ="|")
pattern_pathways=selected_pathway_names[grepl(patterns,selected_pathway_names)]

HTMLStart(outdir="./HTML_results/", file="moduletest_html",
          extension="html", echo=FALSE,HTMLframe = FALSE,CSSFile = "R2HTML_tao.css")
HTML.title("Module in the paper", HR=1,CSSFile = "R2HTML_tao.css")

#HTMLhr()

HTML.title("data_frame", HR=2,CSSFile = "R2HTML_tao.css")
HTML(results,CSSFile = "R2HTML_tao.css")

#HTMLhr()

HTML.title("Pathways containing 'signaling ' or other words in their names", HR=2,CSSFile = "R2HTML_tao.css")
HTML(paste(pattern_pathways,sep = "<br/>"))


HTML.title("pathway2fisher plot", HR=2,CSSclass="")
pathway2fisher_p
HTMLplot(Width = 1200, Height = 350,CSSclass="") 

HTML.title("pathway2coeff plot", HR=2,CSSclass="")
pathway2coeff_p
HTMLplot(Width = 1200, Height = 350,CSSclass="")

HTML.title("pathway2gene plot", HR=2,CSSclass="")
pathway2gene_p
HTMLplot(Width = 1200, Height = 350,CSSclass="")

# for (i in 1:length(module_subgraphs)) {
#   HTMLInsertGraph(paste("moduletest_",i,".PNG",sep = ""),
#                   file ="./HTML_results/moduletest_html.html" )
# }


HTMLStop()



## heatmap of reactome  for all modules

write.csv(dream_consensus_modules_results, file = "dream_consensus_modules_results.csv",
          row.names=FALSE, col.names = TRUE)


DC_results=read.csv(file = "dream_consensus_modules_results_old.csv", header = TRUE,sep = ",",stringsAsFactors = FALSE)

DC_reactome_results=DC_results[grepl("HSA",DC_results[,"Go_Reactome_root_id"]),]

DC_reactome_str=DC_reactome_results[grepl("STRING",DC_reactome_results[,"module"]),]

unique_modules=unique(DC_reactome_str[,"module"])


reactome_ontology=readRDS("human_reactome_ontology.rds")
reactome_roots=which(sapply(sapply(V(reactome_ontology), function(x) neighbors(reactome_ontology,x, mode="in")), length) == 0)
reactome_roots_pathways=V(reactome_ontology)$pathway_names[reactome_roots]


reactome_heatMap_m=matrix(0,nrow = length(unique_modules),ncol = length(reactome_roots_pathways))
reactome_heatMap=data.frame(reactome_heatMap_m,row.names = unique_modules)
colnames(reactome_heatMap)=reactome_roots_pathways

for (i in 1:length(unique_modules)) {
  #i=143
  module_name=unique_modules[i]
  module_pathways=DC_reactome_str[DC_reactome_str[,"module"]==module_name,"Go_Reactome_root_names"]
  tmp=sapply(module_pathways,function(x){
    if(grepl("#",x)){
      x=strsplit(x,split="#")
    }
    return(x)
  })
  tmp=unlist(tmp)
  unique_tmp=unique(tmp)
  for(j in 1:length(unique_tmp)){
    print(sum(tmp==unique_tmp[j]) )
    reactome_heatMap[module_name,unique_tmp[j]]=sum(tmp==unique_tmp[j]) 
  }
}


jpeg("reactome_module_heatmap.jpg",quality = 100,pointsize = 30,width = 1500,height = 1500)
#heatmap.2(as.matrix(reactome_heatMap),col=redgreen,
#          vline = TRUE,hline = TRUE,margins = c(12,8),
#         labRow = FALSE,srtCol=45)


rtmat <- as.matrix(reactome_heatMap)
rtmatColSum <- colSums(rtmat)
rtmat <- rtmat[, order(rtmatColSum, decreasing = TRUE)]
library(ribiosPlot)
rtmatCascade <- cascadeOrder(rtmat)
biosHeatmap(rtmat[rtmatCascade,], Colv=FALSE, Rowv=FALSE,
            col="royalbluered", labRow="", scale="row",
            xlab="Reactome pathway cateogies",
            cexCol=1.25, ylab="PPI-STRING-Modules")

dev.off()
#blackyellow
#heatmap.2(sel, col=redgreen(75), scale="row", ColSideColors=col,
#          key=TRUE, symkey=FALSE, density.info="none",cexRow=1,cexCol=1,margins=c(12,8),trace="none",srtCol=45)
#graphics.off()
#image(as.matrix(reactome_heatMap[1,]))
          


## heatmap of MGI Go terms   for all modules
MGI_GO_bin=read.csv(file = "MGI_GO_slim/map2MGIslim.txt",header = TRUE,sep = "\t",stringsAsFactors = FALSE)
MGI_GO_bp=unique(MGI_GO_bin[MGI_GO_bin$aspect=="P","GOSlim_bin"])
MGI_GO_mf=unique(MGI_GO_bin[MGI_GO_bin$aspect=="F","GOSlim_bin"])
MGI_GO_cc=unique(MGI_GO_bin[MGI_GO_bin$aspect=="C","GOSlim_bin"])

GO_categories=c(MGI_GO_bp,MGI_GO_cc,MGI_GO_mf,"unknow")

DC_results=read.csv(file = "dream_consensus_modules_results.csv", header = TRUE,sep = ",",stringsAsFactors = FALSE)

DC_GO_results=DC_results[grepl("GO",DC_results[,"Go_Reactome_root_id"]),]

DC_GO_str=DC_GO_results[grepl("STRING",DC_GO_results[,"module"]),]

go_ontology <- readRDS("go_ontology.rds")
go_ontology_names=V(go_ontology)$pathway_names
names(go_ontology_names)=as_ids(V(go_ontology))

#test1=unique(MGI_GO_bin$GO_id)
#test2=unique(DC_GO_str$pathway_id)
#test3=intersect(unique(MGI_GO_bin$GO_id),unique(DC_GO_str$pathway_id))
#test4=setdiff(unique(DC_GO_str$pathway_id),unique(MGI_GO_bin$GO_id))
#test4=setdiff(unique(DC_GO_results$pathway_id),unique(MGI_GO_bin$GO_id))
#test5=go_ontology_names[test4]

unique_modules=unique(DC_GO_str[,"module"])



GO_heatMap_m=matrix(0,nrow = length(unique_modules),ncol = length(GO_categories))
GO_heatMap=data.frame(GO_heatMap_m,row.names = unique_modules)
colnames(GO_heatMap)=GO_categories

for (i in 1:length(unique_modules)) {
  #i=143
  print(i)
  module_name=unique_modules[i]
  module_pathways=DC_GO_str[DC_GO_str[,"module"]==module_name,"Go_Reactome_root_id"]
  tmp=sapply(module_pathways,function(x){
    if(grepl("#",x)){
      x=strsplit(x,split="#")
    }
    return(x)
  })
  tmp=unlist(tmp)
  tmp_parents=MGI_GO_bin[match(tmp,MGI_GO_bin$GO_id),"GOSlim_bin"]
  tmp_parents=replace(tmp_parents,is.na(tmp_parents),"unknow")
  unique_tmp_parents=unique(tmp_parents)
  for(j in 1:length(unique_tmp_parents)){
    #print( sum(tmp_parents==unique_tmp_parents[j]))
    GO_heatMap[module_name,unique_tmp_parents[j]]=sum(sum(tmp_parents==unique_tmp_parents[j])) 
  }
}


jpeg("MGI_GO_module_heatmap.jpg",quality = 100,pointsize = 30,width = 1500,height = 1500)
#heatmap.2(as.matrix(reactome_heatMap),col=redgreen,
#          vline = TRUE,hline = TRUE,margins = c(12,8),
#         labRow = FALSE,srtCol=45)


rtmat <- as.matrix(GO_heatMap)
rtmatColSum <- colSums(rtmat)
rtmat <- rtmat[, order(rtmatColSum, decreasing = TRUE)]
library(ribiosPlot)
rtmatCascade <- cascadeOrder(rtmat)
biosHeatmap(rtmat[rtmatCascade,], Colv=FALSE, Rowv=FALSE,
            col="royalbluered", labRow="", scale="row",
            xlab="MGI GO  pathway cateogies",
            cexCol=1.25, ylab="PPI-STRING-Modules")

dev.off()
#blackyellow
#heatmap.2(sel, col=redgreen(75), scale="row", ColSideColors=col,
#          key=TRUE, symkey=FALSE, density.info="none",cexRow=1,cexCol=1,margins=c(12,8),trace="none",srtCol=45)
#graphics.off()
#image(as.matrix(reactome_heatMap[1,]))          




## heatmap of Go   for all modules

DC_results=read.csv(file = "dream_consensus_modules_results_old.csv", header = TRUE,sep = ",",stringsAsFactors = FALSE)

DC_GO_results=DC_results[grepl("GO",DC_results[,"Go_Reactome_root_id"]),]

DC_GO_str=DC_GO_results[grepl("STRING",DC_GO_results[,"module"]),]

unique_modules=unique(DC_GO_str[,"module"])

go_ontology <- readRDS("go_ontology.rds")
go_roots=c("GO:0003674","GO:0005575","GO:0008150")
go_sub_roots=sapply(go_roots, function(x){
  as_ids(neighbors(go_ontology, x, mode = "out"))
})
go_sub_roots=unlist(go_sub_roots)
go_ontology_names=V(go_ontology)$pathway_names
names(go_ontology_names)=as_ids(V(go_ontology))
go_sub_roots_pathways=go_ontology_names[go_sub_roots]


GO_heatMap_m=matrix(0,nrow = length(unique_modules),ncol = length(go_sub_roots_pathways))
GO_heatMap=data.frame(GO_heatMap_m,row.names = unique_modules)
colnames(GO_heatMap)=go_sub_roots_pathways

for (i in 1:length(unique_modules)) {
  #i=143
  module_name=unique_modules[i]
  module_pathways=DC_GO_str[DC_GO_str[,"module"]==module_name,"Go_Reactome_root_names"]
  tmp=sapply(module_pathways,function(x){
    if(grepl("#",x)){
      x=strsplit(x,split="#")
    }
    return(x)
  })
  tmp=unlist(tmp)
  unique_tmp=unique(tmp)
  for(j in 1:length(unique_tmp)){
    print(sum(tmp==unique_tmp[j]) )
    GO_heatMap[module_name,unique_tmp[j]]=sum(tmp==unique_tmp[j]) 
  }
}


jpeg("GO_module_heatmap.jpg",quality = 100,pointsize = 30,width = 2000,height = 1500)
rtmat <- as.matrix(GO_heatMap)
rtmatColSum <- colSums(rtmat)
rtmat <- rtmat[, order(rtmatColSum, decreasing = TRUE)]
library(ribiosPlot)
rtmatCascade <- cascadeOrder(rtmat)
biosHeatmap(rtmat[rtmatCascade,], Colv=FALSE, Rowv=FALSE,
            col="royalbluered", labRow="", scale="row",
            xlab="GO pathway cateogies",
            cexCol=1.25, ylab="PPI-STRING-Modules")

dev.off()

GO_bp=go_ontology_names[sapply(go_roots[3], function(x){
  as_ids(neighbors(go_ontology, x, mode = "out"))
})]
GO_bp_heatmap=GO_heatMap[,GO_bp]
jpeg("GO_bp_module_heatmap.jpg",quality = 100,pointsize = 30,width = 2000,height = 1500)
rtmat <- as.matrix(GO_bp_heatMap)
rtmatColSum <- colSums(rtmat)
rtmat <- rtmat[, order(rtmatColSum, decreasing = TRUE)]
library(ribiosPlot)
rtmatCascade <- cascadeOrder(rtmat)
biosHeatmap(rtmat[rtmatCascade,], Colv=FALSE, Rowv=FALSE,
            col="royalbluered", labRow="", scale="row",
            xlab="GO bp pathway cateogies",
            cexCol=1.25, ylab="PPI-STRING-Modules")

dev.off()
