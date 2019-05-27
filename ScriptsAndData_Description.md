The raw descriptions can be found in "How to use TG_GATEs R code.docx" file

# Data sources and link
“ReactomePathways.gmt”  is downloaded from : https://reactome.org/download-data/Specialized data formats/Reactome Pathways Gene Set
“ReactomePathwaysRelation.txt” and “ReactomePathways.txt” are downloaded from https://reactome.org/download-data/Pathways

“goa_human.gaf” is download from  http://geneontology.org/page/download-go-annotations
Goa/gaf formart: ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/README

“go-basic.obo” is downloaded from http://www.geneontology.org/page/download-ontology



# Description of code

go_preprocess_tree.R is used for go data preprocess and tree function test

reactome_preprocess_tree.R is used for reactome data preprocess and tree function test

Gene_pathway_matrix.R is used to construct gene pathway matrix and get raw go gene list(raw_go_list.rda)

“Module_Identification_GO_REA_PPI_STRING_test.R” is used to test model on DreamChallenge PPT-STRING_Consensus networks,
which will generate html fils whose format are specified in https://drive.google.com/drive/folders/1tuFt83ot_P2kTb9vPkH1bbL0Oo08rBdZ and a txt file contains all annotation results

Module_Identification_GO_REA_PPI_STRING_test_onlyTxtfile.R is used if you only want get some results in txt file format

Module_Identification_GO_REA_PPI_STRING_test_onlyTxtfile_forDanielsPaper.R is used to get txt file format for Daniel's paper


Module_Identification_suplementary_figures.R is used to draw supplementary figures as described in :https://drive.google.com/drive/folders/1tuFt83ot_P2kTb9vPkH1bbL0Oo08rBdZ   or  https://docs.google.com/document/d/1w-RcyMFWbx3fictyHthqdfu2sDgf0pvdnb0H-Z_5iNY/edit#
Module_Identification_suplementary_figures_sup.R  is used to save some old code for Module_Identification_suplementary_figures.R script

regression_selected_pathways_funcs.R  is  used to write a function to use the model


“Module_Identification_funcs.R”  is used to store all the supplementary  functions for module identification project



DreamChallenge_website and DreamChallenge folder in local computer is used to build a website to use our results on DreamChallenge data. DreamChallenge_website is original folder and DreamChallenge is a cleaned folders uploaded to Github.








4. “Module_Identification_GO_REA_UNI_2rd.R” is used to test new model. I.e. map to hierarchical trees and draw heatmap of modules for reactome terms

5. . “Module_Identification_funcs.R”  is used to store all the supplementary  functions for module identification project
6. “Module_Identification_GO_REA_UNI_all_2rd.R ” is used to get module identification for all modules

7. “Module_Identification_GO_REA_PPI_STRING_test.R” is used to test new model after 08.14 meeting
8. “Module_Identification_GO_REA_PPI_STRING.R” is used to get module identification for all modules
9. Module_Identification_suplementary_figures.R is used to draw supplementary figures as described in :https://docs.google.com/document/d/1w-RcyMFWbx3fictyHthqdfu2sDgf0pvdnb0H-Z_5iNY/edit#
Module_Identification_suplementary_figures_sup.R  is used to save some old code for Module_Identification_suplementary_figures.R script

Interpret_TGGATEs.R  is  used to interpret TG-GATEs model based on the regression based pathway selection methods
regression_selected_pathways_funcs.R  is  used to write a function to use the model
last_interpret.Rmd
