library(ggplot2)
library(ribiosPlot)
library(ribiosPlot)
library(RColorBrewer)
setwd("~/Disease_module_identification_DREAM_change/")

Daniel_results_name="dream_consensus_modules.functional_enrichment.txt"
Daniel_results=read.csv(Daniel_results_name,header = TRUE,sep = "\t",stringsAsFactors = FALSE)
Daniel_network_results=Daniel_results[Daniel_results$network=="PPI-STRING",]
GO_Rec_index=c(grep("GO",Daniel_network_results$termId),grep("REACTOME",Daniel_network_results$term))
Daniel_network_results=Daniel_network_results[GO_Rec_index,]

Tao_results= read.csv(file = "dream_consensus_modules_results_alpha0.5.csv",header = TRUE,stringsAsFactors = FALSE) # dream_consensus_modules_results_old.csv dream_consensus_modules_results_alpha0.5.csv
Tao_reactome_results=Tao_results[grep("HSA",Tao_results$pathway_id),]


###############################Reactome classes associated with modules (analogous to Fig. 3f)
reactome_classes=unique(Tao_reactome_results$Go_Reactome_root_names)
reactome_classes=unique(unlist(sapply(reactome_classes, function(x){strsplit(x,split = "#")[[1]]})))
reactome_colors=colorRampPalette(brewer.pal(9, "Set1"))(length(reactome_classes))
names(reactome_colors)=reactome_classes
reactome_class2moduleCount=data.frame(reactome_classes=reactome_classes,module_num=rep(0,length(reactome_classes)),names=reactome_classes)
rownames(reactome_class2moduleCount)=reactome_classes
for(i in 1:length(reactome_classes)){
  reactome_class2moduleCount[reactome_classes[i],"module_num"]=sum(grepl(reactome_classes[i],Tao_reactome_results$Go_Reactome_root_names,fixed = TRUE))
}

disease_index=grep("Disease",rownames(reactome_class2moduleCount))
reactome_class2moduleCount_noDisease=reactome_class2moduleCount[-disease_index,]
reactome_classes_noDisease=reactome_classes[-disease_index]
reactome_colors_noDisease=reactome_colors[-disease_index]
reactome_class2moduleCount_ggplot=ggplot(data = reactome_class2moduleCount_noDisease,aes(x=reactome_classes,y=module_num,fill=names))+
  geom_bar(stat = "identity")+
  coord_flip()+
  scale_x_discrete(limits=reactome_classes_noDisease[order(reactome_class2moduleCount_noDisease$module_num)])+                # !!!!!!!!!
  scale_fill_manual(values=reactome_colors_noDisease[reactome_classes[order(reactome_class2moduleCount_noDisease$module_num)]])+
  ylab("Module count")+
  xlab("Reactome classes")+
  theme(legend.position="none")


jpeg("reactome_class2moduleCount_ggplot0.5.jpg",quality = 100,pointsize = 30,width = 1000,height = 500)
reactome_class2moduleCount_ggplot
dev.off()


########################
#Distribution of (a) #significant terms (FDR corrected P-value from noncentral hypergeometric) and
#(b) #selected terms (regression). E.g., two boxplots or histograms.

tao_module_names=unique(Tao_results$module)
hist_tao=data.frame(module=tao_module_names,pathway_num=rep(0,length(tao_module_names)),row.names = tao_module_names)  
for(i in 1:length(tao_module_names)){
 hist_tao[tao_module_names[i],"pathway_num"]=sum(Tao_results$module==tao_module_names[i]) 
}


daniel_module_names=unique(Daniel_network_results$module)
hist_daniel=data.frame(module=daniel_module_names,pathway_num=rep(0,length(daniel_module_names)),row.names=daniel_module_names)
for (i in 1:length(daniel_module_names)) {
  hist_daniel[daniel_module_names[i],"pathway_num"]=sum(Daniel_network_results$module==daniel_module_names[i])
}

ggplot(hist_tao, aes(x=pathway_num)) + geom_histogram()
ggplot(hist_daniel,aes(x=pathway_num))+geom_histogram()

combined_hist_frame=data.frame(name=c(rep('tao',nrow(hist_tao)),rep('daniel',nrow(hist_daniel))),
                               pathway_num=c(hist_tao$pathway_num,hist_daniel$pathway_num))

jpeg("Hist_plot0.5.jpg",quality = 100,pointsize = 30,width = 1000,height = 500)
ggplot(combined_hist_frame,aes(x=pathway_num))+
  geom_histogram(
    bins=30,
    colour="black", fill="white"
  )+
  facet_wrap(~name,ncol=1)+
  theme_bw()+
  ggtitle("Distribution of num for selected terms ")
dev.off()


##########################Scatterplot: #significant terms (Fisher) vs. #selected terms for each module
common_module_names=intersect(tao_module_names,daniel_module_names)
scatter_plot=data.frame(module=common_module_names,row.names = common_module_names,
                        tao_pathway_num=hist_tao[common_module_names,]$pathway_num,
                        daniel_pathway_num=hist_daniel[common_module_names,]$pathway_num
                        )


jpeg("scatter_plot0.5.jpg",quality = 100,pointsize = 30,width = 1000,height = 500)
ggplot(scatter_plot, aes(x=daniel_pathway_num, y=tao_pathway_num)) + geom_point(shape=1) +
  theme_bw()
 # geom_abline(intercept = 0)    #show sparsity
dev.off()

cor(scatter_plot$tao_pathway_num,scatter_plot$daniel_pathway_num)



###################################### 
#overlap GO selecetd terms between two methods
Daniel_network_GO_results=Daniel_network_results[grepl("GO",Daniel_network_results$termId),]
daniel_GO_module_names=unique(Daniel_network_GO_results$module)
Tao_go_results=Tao_results[grep("GO",Tao_results$pathway_id),]
tao_go_module_names=unique(Tao_go_results$module)
common_module_names=intersect(tao_go_module_names,daniel_GO_module_names)

tao_go_id=unique(Tao_go_results$pathway_id)  #1949
daniel_go_id=unique(Daniel_network_GO_results$termId) #10838
total_overlap_go_id=intersect(tao_go_id,daniel_go_id)  #1872

overlap_go_terms_frame=data.frame(matrix(0,nrow = length(common_module_names),ncol = 3),row.names = common_module_names)
colnames(overlap_go_terms_frame)=c('tao_go_pathway_num','daniel_go_pathway_num','overlap_go_pathway_num')

for(i in 1:length(common_module_names)){
  tao_selected_pathwayIDs=Tao_go_results[Tao_go_results$module==common_module_names[i],"pathway_id"]
  overlap_go_terms_frame[common_module_names[i],"tao_go_pathway_num"]=length(tao_selected_pathwayIDs)
  daniel_selected_pathwayIDs=Daniel_network_GO_results[Daniel_network_GO_results$module==common_module_names[i],"termId"]
  overlap_go_terms_frame[common_module_names[i],"daniel_go_pathway_num"]=length(daniel_selected_pathwayIDs)
  overlap_go_terms_frame[common_module_names[i],"overlap_go_pathway_num"]=length(intersect(tao_selected_pathwayIDs,daniel_selected_pathwayIDs))
}



jpeg("overlap_go_terms_scatter_plot0.5.jpg",quality = 100,pointsize = 30,width = 1000,height = 500)
ggplot(overlap_go_terms_frame, aes(x=tao_go_pathway_num, y=overlap_go_pathway_num)) + geom_point(shape=1) +
  theme_bw()
# geom_abline(intercept = 0)    #show sparsity
dev.off()

##################################
#Show that the selected terms are less redundant than the ones from Fisher: 
#E.g., for each significant term, compute max Jaccard index across all other significant terms of that module. 
#Show distribution of these max Jaccard indexes for all significant terms and modules. 
#Do the same for the selected terms.
load("gene_pathway_matrix.rda")     # 
all_genes=rownames(gene_pathway_matrix)
all_pathways=colnames(gene_pathway_matrix)

raw_go_list=readRDS("raw_go_list.rds")
raw_go_list_names=names(raw_go_list)


jaccard_index=c()
jaccard_module=c()
for(i in 1:length(common_module_names)){
  selected_pathwayIDs=Tao_go_results[Tao_go_results$module==common_module_names[i],"pathway_id"]
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
     jaccard_module=c(jaccard_module,common_module_names[i])
   }
   jaccard_index=c(jaccard_index,unlist(modulei_jaccard_index))
 }
}
jaccard_index=unlist(jaccard_index)
#saveRDS(jaccard_index,file="tao_jaccard_index_median_0.5.rds")
saveRDS(jaccard_module,file="tao_jaccard_module_0.5.rds")
saveRDS(jaccard_index,file="tao_jaccard_index_0.5.rds")

jaccard_index=readRDS("tao_jaccard_index_0.5.rds")    

jaccard_index_frame=data.frame(jaccard=jaccard_index)

ggplot(jaccard_index_frame, aes(x=jaccard)) + 
  geom_histogram(      # Histogram with density instead of count on y-axis
    bins=30,
    colour="black", fill="white") +
  #facet_wrap()
  theme_bw()+
  ggtitle("jaccard index distributioin(Tao's result)")
#geom_density(alpha=.2, fill="#FF6666") 



daniel_jaccard_index=c()
daniel_jaccard_module=c()

#Daniel_network_GO_results_sort=Daniel_network_GO_results[order(Daniel_network_GO_results$P.noncentral),]

for(i in 1:length(common_module_names)){
  #temp_daniel_results=Daniel_network_GO_results_sort[Daniel_network_GO_results_sort$module==common_module_names[i],]
  selected_pathwayIDs=Daniel_network_GO_results[Daniel_network_GO_results$module==common_module_names[i],"termId"]
  selected_pathwayIDs=selected_pathwayIDs[selected_pathwayIDs %in% raw_go_list_names ] # !!!!!!!problem here ???
  tao_selected_pathwayIDs=Tao_go_results[Tao_go_results$module==common_module_names[i],"pathway_id"]
  selected_pathwayIDs=selected_pathwayIDs[1:length(tao_selected_pathwayIDs)]
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
      daniel_jaccard_module=c(daniel_jaccard_module,common_module_names[i])
    }
    daniel_jaccard_index=c(daniel_jaccard_index,unlist(modulei_jaccard_index))
  }
}

daniel_jaccard_index=unlist(daniel_jaccard_index)
saveRDS(daniel_jaccard_index,file="daniel_jaccard_index.rds")
#saveRDS(daniel_jaccard_index,file="daniel_jaccard_index_median.rds")
saveRDS(daniel_jaccard_module,file = "daniel_jaccard_module.rds")

daniel_jaccard_index=readRDS("daniel_jaccard_index.rds")    

daniel_jaccard_index_frame=data.frame(jaccard=daniel_jaccard_index)
ggplot(daniel_jaccard_index_frame, aes(x=jaccard)) + 
  geom_histogram(      # Histogram with density instead of count on y-axis
    bins=30,
    colour="black", fill="white") +
  #facet_wrap()
  theme_bw()+
  ggtitle("Jaccard index distributioin(Daniel's result)")
#geom_density(alpha=.2, fill="#FF6666") 


combined_jaccard_index_frame=data.frame(name=c(rep("tao",length(jaccard_index)),rep("daniel",length(daniel_jaccard_index))),
                                        jaccard=c(jaccard_index,daniel_jaccard_index))
ggplot(combined_jaccard_index_frame, aes(x=jaccard)) + 
  geom_histogram(      # Histogram with density instead of count on y-axis
    bins=30,
    colour="black", fill="white") +
  facet_wrap(~name,ncol = 1)+
  theme_bw()+
  ggtitle("Jaccard index distributioin(Daniel's result)")
  
tao_multiply_jaccard_index=rep(jaccard_index,round(length(daniel_jaccard_index)/length(jaccard_index)))
multiply_combined_jaccard_index_frame=data.frame(name=c(rep("tao",length(tao_multiply_jaccard_index)),rep("daniel",length(daniel_jaccard_index))),
                                        jaccard=c(tao_multiply_jaccard_index,daniel_jaccard_index))
ggplot(multiply_combined_jaccard_index_frame, aes(x=jaccard)) + 
  geom_histogram(      # Histogram with density instead of count on y-axis
    bins=50,
    colour="black", fill="white") +
  facet_wrap(~name,ncol = 1)+
  theme_bw()+
  ggtitle("Jaccard index distributioin(Daniel's result)")



#For each module associated with a Reactome class (columns), 
#count of module genes of a given UniProt class (in rows) (analogous to Fig. 3g)

# 
# library(ribiosAnnotation)
# all_genes_ribios=annotateAnyIDs(all_genes)
# gene_pathway_matrix_all_genes=data.frame(gene_symbol=all_genes)
# gene_pathway_matrix_all_geneIDs=data.frame(gene_id=all_genes_ribios$GeneID)
# write.table(gene_pathway_matrix_all_genes,file = "gene_pathway_matrix_all_genes.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
# write.table(gene_pathway_matrix_all_geneIDs,file = "gene_pathway_matrix_all_geneIDs.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
# 

consensus_genesets=scan("dream_consensus_modules.gmt",what = "",sep = "\n")
consensus_genesets=strsplit(consensus_genesets,split = "\t")
consensus_names=sapply(consensus_genesets, function(x){x[1]})
consensus_genesets=lapply(consensus_genesets,function(x){
  x_genes=x[-c(1,2)]
  return(x_genes)
})
names(consensus_genesets)=consensus_names

tao_reactome_module_names=unique(Tao_reactome_results$module)
tao_reactome_genesets=consensus_genesets[tao_reactome_module_names]
tao_reactome_all_genes=vector()
for (i in 1:length(tao_reactome_genesets)){
  tao_reactome_all_genes=union(tao_reactome_all_genes,tao_reactome_genesets[[i]])
}


uniprot_frame=read.csv(file="uniprot_frame.csv",header = TRUE,row.names = 1)   #19982 * 721
uniprot_all_genes=rownames(uniprot_frame)
uniprot_reactome_common_genes=intersect(uniprot_all_genes,tao_reactome_all_genes)
uniprot_frame=uniprot_frame[uniprot_reactome_common_genes,]
#uniprot_all_genes=rownames(uniprot_frame)
#uniprot_all_keywords=colnames(uniprot_frame)


#
test=hist(colSums(uniprot_frame),200)
sum(colSums(uniprot_frame)>100 )#& colSums(uniprot_frame)<150)
sum(rowSums(uniprot_frame)>50)

uniprot_frame_filted=uniprot_frame[,colSums(uniprot_frame)>50]
uniprot_all_genes=rownames(uniprot_frame_filted)
uniprot_all_keywords=colnames(uniprot_frame_filted)
uniprot_reactome_common_genes=intersect(uniprot_all_genes,tao_reactome_all_genes)

## cluster uniprot keywords
#raw=unname(uniprot_frame)
uniprot_frame_d <- dist(t(uniprot_frame_filted), method = "manhattan") # distance matrix
uniprot_frame_fit <- hclust(uniprot_frame_d, method="ward.D2") 
plot(uniprot_frame_fit,ann = FALSE)
clusterCut <- cutree(uniprot_frame_fit,h = 1000)   #k=40

clusterCut_groups=unique(clusterCut)
clusterCut_name_groups=lapply(clusterCut_groups, function(i){colnames(uniprot_frame_filted)[clusterCut==clusterCut_groups[i]]})
saveRDS(clusterCut_name_groups,"clusterCut_name_groups.rds")

uniprot_class_tao_defined=clusterCut_name_groups
names(uniprot_class_tao_defined)=c("other","Transport","Secreted_Signal","Transcription","G.protein.coupled.receptor","Nucleotide.binding",
     "Phosphoprotein","Ubl.conjugation_Isopeptide.bond","Developmental.protein_Differentiation","Polymorphism",
     "Transmembrane" ,"Transit.peptide_Mitochondrion","Complete/Reference.proteome","Nucleus",
     "Glycoprotein","Zinc","Hydrolase","Transferase","Membrane","Disease.mutation","Endoplasmic.reticulum",
     "Coiled.coil","Alternative.splicing","Metal.binding","Repeat","X3D.structure","Direct.protein.sequencing",
     "Cytoplasm","Acetylation","Disulfide.bond","Cell.membrane")

#library(ape)
#hcd <- as.dendrogram(uniprot_frame_fit)
#plot(as.phylo(uniprot_frame_fit), type = "fan",show.tip.label = FALSE)
#uniprot_frame_reorder=uniprot_frame[,fit$order]



uniprot_hist=data.frame(matrix(0,nrow = length(names(uniprot_class_tao_defined)),ncol = length(tao_reactome_module_names)),
                        row.names =names(uniprot_class_tao_defined))
colnames(uniprot_hist)=tao_reactome_module_names
for(i in 1:length(tao_reactome_module_names)){
  temp_genes=tao_reactome_genesets[[tao_reactome_module_names[i]]]
  temp_common_genes=intersect(temp_genes,uniprot_reactome_common_genes)
  for (j in 1:length(names(uniprot_class_tao_defined))) {
    temp_uniprot_classes=uniprot_class_tao_defined[[names(uniprot_class_tao_defined[j])]]
    if(length(temp_uniprot_classes)>1){
      uniprot_hist[names(uniprot_class_tao_defined)[j],tao_reactome_module_names[i]]=sum(colSums(uniprot_frame_filted[temp_common_genes,temp_uniprot_classes]) )
    }else{
      uniprot_hist[names(uniprot_class_tao_defined)[j],tao_reactome_module_names[i]]=sum(uniprot_frame_filted[temp_common_genes,temp_uniprot_classes] )
    }
  }
}

#test=uniprot_hist[apply(uniprot_hist, 1, max)>30,]
test=uniprot_hist

library(pheatmap)

reactome_classes=unique(Tao_reactome_results$Go_Reactome_root_names)
reactome_classes=unique(unlist(sapply(reactome_classes, function(x){strsplit(x,split = "#")[[1]]})))

module2reactomeClasses=data.frame(row.names = c("module","reactome_class"))

for(i in 1:length(tao_reactome_module_names)){
  #i=198
  temp_classes=unique(Tao_reactome_results$Go_Reactome_root_names[Tao_reactome_results$module==tao_reactome_module_names[i]])
  temp_classes=unique(unlist(sapply(temp_classes, function(x){strsplit(x,split = "#")[[1]]})))
  for (j in 1:length(temp_classes)) {
    module2reactomeClasses=rbind(module2reactomeClasses,data.frame(module=tao_reactome_module_names[i],reactome_class=temp_classes[j]))
  }
}

test1=NULL
test1_classes=c()
test1_class_num=c()
for (i in 1:length(reactome_classes)) {
    temp_modules=module2reactomeClasses$module[module2reactomeClasses$reactome_class==reactome_classes[i]] 
    if(length(temp_modules)>1){
      temp_frame=test[,temp_modules]
    }else{
      temp_frame=data.frame(matrix(test[,temp_modules],ncol = 1))
    }
    test1_classes=c(test1_classes,rep(reactome_classes[i],ncol(temp_frame)))
    test1_class_num=c(test1_class_num,ncol(temp_frame))
    colnames(temp_frame)=paste(reactome_classes[i],1:ncol(temp_frame),sep="_")
    if(i>1){
      if(ncol(temp_frame)>1){
        d <- dist(t(temp_frame), method = "euclidean") # distance matrix
        fit <- hclust(d, method="centroid") 
        test1=cbind(test1,temp_frame[,fit$order])
      }else{
        test1=cbind(test1,temp_frame)
      }
      
    }else{
      test1=temp_frame
    }
}

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
reactome_class_colors=sample(color,length(reactome_classes))
names(reactome_class_colors)=reactome_classes
hist_colors=unlist(sapply(1:length(reactome_classes), function(i){rep(reactome_class_colors[i],test1_class_num[i])}))
names(hist_colors)=colnames(test1)
annotation=data.frame(reactome_classes=test1_classes)
rownames(annotation)=colnames(test1)

anno_colors=list(reactome_classes=reactome_class_colors)

jpeg("uniprotClass_ReactomeModule_heatmap1.0.jpg",quality = 100,pointsize = 30,width = 2000,height = 1500)
pheatmap(as.matrix(test1),cluster_cols = FALSE,cluster_rows = TRUE,show_colnames  = FALSE,annotation = annotation,annotation_colors = anno_colors)
dev.off()



jpeg("uniprotClass_ReactomeModule_biosHeatmap1.0.jpg",quality = 100,pointsize = 30,width = 2000,height = 1500)
pheatmap(as.matrix(test1[cascadeOrder(test1),]),cluster_cols = FALSE,cluster_rows = FALSE,show_colnames  = FALSE,annotation = annotation,annotation_colors = anno_colors)
dev.off()
