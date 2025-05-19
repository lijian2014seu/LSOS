##### description: Gene set annotation (categorization) and importance ranking

##### 01Start:  Gene set annotation (categorization) #####
go_enrich_res <- read.table( file = "WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_GOenrich_res.txt", sep="\t")
kegg_enrich_res <- read.table( file = "WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_KEGGenrich_res.txt", sep="\t")

geneset_all_df_PPI <- rbind(go_enrich_res, kegg_enrich_res)
#Bone-related
bone_enrich_res <- unique(geneset_all_df_PPI[c(grep(" bone", geneset_all_df_PPI$Description, ignore.case = TRUE), grep("^bone ", geneset_all_df_PPI$Description, ignore.case = TRUE)), ])
dim(bone_enrich_res)

bone_enrich_res <- unique(rbind(bone_enrich_res, geneset_all_df_PPI[c(grep(" stone", geneset_all_df_PPI$Description, ignore.case = TRUE), grep("^stone ", geneset_all_df_PPI$Description, ignore.case = TRUE)), ]))
dim(bone_enrich_res)

bone_enrich_res <- unique(rbind(bone_enrich_res, geneset_all_df_PPI[c(grep(" osteoclast", geneset_all_df_PPI$Description, ignore.case = TRUE), grep("^osteoclast ", geneset_all_df_PPI$Description, ignore.case = TRUE)), ]))
dim(bone_enrich_res)

bone_enrich_res <- unique(rbind(bone_enrich_res, geneset_all_df_PPI[c(grep(" collagen", geneset_all_df_PPI$Description, ignore.case = TRUE), grep("^collagen ", geneset_all_df_PPI$Description, ignore.case = TRUE)), ]))
dim(bone_enrich_res)
bone_enrich_res <- unique(rbind(bone_enrich_res, geneset_all_df_PPI[c(grep("chondroitin sulfate", geneset_all_df_PPI$Description, ignore.case = TRUE), grep("^collagen ", geneset_all_df_PPI$Description, ignore.case = TRUE)), ]))
dim(bone_enrich_res)

#Muscle-related 
muscle_enrich_res <- unique(geneset_all_df_PPI[c(grep(" muscle", geneset_all_df_PPI$Description, ignore.case = TRUE), grep("^muscle ", geneset_all_df_PPI$Description, ignore.case = TRUE)), ])
dim(muscle_enrich_res)
muscle_enrich_res$Description

#Cardiovascular
cardiovascular_enrich_res <- unique(geneset_all_df_PPI[grep("cardia|vascular|endocardium|cardio|cardiac|VEGF", geneset_all_df_PPI$Description, ignore.case = TRUE), ])
dim(cardiovascular_enrich_res)
cardiovascular_enrich_res$Description


#Coagulation-related
Coagulation_enrich_res <- unique(geneset_all_df_PPI[grep("Coagulation|hemostasis|platelet|thrombosis", geneset_all_df_PPI$Description, ignore.case = TRUE), ])
dim(Coagulation_enrich_res)
Coagulation_enrich_res$Description

#Mitochondria-related
mitochon_enrich_res <- unique(geneset_all_df_PPI[grep("mitochon|NADPH", geneset_all_df_PPI$Description, ignore.case = TRUE), ])
dim(mitochon_enrich_res)
mitochon_enrich_res <- unique(rbind(mitochon_enrich_res,geneset_all_df_PPI[grep("NADPH", geneset_all_df_PPI$Description, ignore.case = TRUE), ]))
dim(mitochon_enrich_res)
mitochon_enrich_res$Description

#RNA
ncRNAs_enrich_res <- unique(geneset_all_df_PPI[grep("RNA ", geneset_all_df_PPI$Description, ignore.case = TRUE), ])
dim(ncRNAs_enrich_res)
ncRNAs_enrich_res$Description

#neurological (neurons, axons, telencephalon, synapses, neural precursors, astrocytes, beta amyloid)
neuron_enrich_res <- unique(geneset_all_df_PPI[grep("neuron|axonogenesis|telencephalon|neural|nerve|astrocyte|amyloid-beta", geneset_all_df_PPI$Description, ignore.case = TRUE), ])
dim(neuron_enrich_res)
neuron_enrich_res$Description

#Immunity-related
Immune_enrich_res <- unique(geneset_all_df_PPI[grep("neutrophil|mast cell|leukocyte|immun|interleukin|T cell|lymphocyte|B cell|cytokine|phagocytosis|natural killer cell|T-helper 2|mononuclear|antigen|Th1|Chemokine|C-type lectin receptor|Hematopoietic", geneset_all_df_PPI$Description, ignore.case = TRUE), ])
dim(Immune_enrich_res)
Immune_enrich_res$Description

#inflammatory-related
inflam_enrich_res <- unique(geneset_all_df_PPI[grep("inflammatory|cytokine|chemokine|NF-kappaB|TNF|Toll|MAPK|JAK", geneset_all_df_PPI$Description, ignore.case = TRUE), ])
dim(inflam_enrich_res)
inflam_enrich_res$Description


#Metabolism-related
catabolic_enrich_res <- unique(geneset_all_df_PPI[grep("catabolic|metabolic", geneset_all_df_PPI$Description, ignore.case = TRUE), ])
dim(catabolic_enrich_res)
catabolic_enrich_res <- unique(rbind(catabolic_enrich_res,geneset_all_df_PPI[grep("metabolic", geneset_all_df_PPI$Description, ignore.case = TRUE), ]))
dim(catabolic_enrich_res)
catabolic_enrich_res$Description

#Stress response, reactive oxygen species
stress_enrich_res <- unique(geneset_all_df_PPI[grep("reactive oxygen", geneset_all_df_PPI$Description, ignore.case = TRUE), ])
dim(stress_enrich_res)
stress_enrich_res$Description

#Pathogens,Viral,Bacterial
pathogens_enrich_res <- unique(geneset_all_df_PPI[grep("Pathogen|Vir(al|us)|Bacteria|infection", geneset_all_df_PPI$Description, ignore.case = TRUE), ])
dim(pathogens_enrich_res)
pathogens_enrich_res$Description

#Basic cell function and structure (death, apoptosis, differentiation; cytoskeleton: actin, cytoskeleton, microtubules)
cell_enrich_res <- unique(geneset_all_df_PPI[grep("cell death|differentiation|necro|apopto|cytoskeleton|actin|microtubule", geneset_all_df_PPI$Description, ignore.case = TRUE), ])
dim(cell_enrich_res)
cell_enrich_res$Description
##### 01End:  Gene set annotation (categorization) #####


##### 02Start: Ranking of geneset importance  #####
enrich_res_list <- ls()[grep("_enrich_res",ls())]
enrich_res_geneset <- lapply(enrich_res_list, function(x){x <- eval(parse(text=x)); genesets <- x$Description; return(genesets)})
names(enrich_res_geneset) <- enrich_res_list

sig_geneset_interact_scores_df <- read.table("WGCNA_PPI_analysis/result/gene_set_interact/sig_geneset_interact_scores.txt", sep="\t")
sub_sig_geneset_interact_scores_df <- sig_geneset_interact_scores_df[(sig_geneset_interact_scores_df$V1%in%unlist(enrich_res_geneset) & sig_geneset_interact_scores_df$V2%in%unlist(enrich_res_geneset)),]

geneset_importance_list <- lapply(enrich_res_list, function(x){
       genesets <- unlist(enrich_res_geneset[x]);
	   sub_sig_geneset_interact_scores_df1 <- sub_sig_geneset_interact_scores_df[(!sub_sig_geneset_interact_scores_df$V1%in%genesets | !sub_sig_geneset_interact_scores_df$V2%in%genesets),];       
	   geneset_importance_df <- data.frame(geneset=genesets);
       geneset_importance_df$importance <- unlist(lapply(geneset_importance_df$geneset, function(y){sum(sub_sig_geneset_interact_scores_df1$V3[sub_sig_geneset_interact_scores_df1$V1==y|sub_sig_geneset_interact_scores_df1$V2==y])}));
	   geneset_importance_df <- geneset_importance_df[geneset_importance_df$importance!=0,];
       return(geneset_importance_df);
	   })
names(geneset_importance_list) <- enrich_res_list
geneset_importance_df <- do.call(rbind, geneset_importance_list)
geneset_importance_df$type <- stringr::str_split_i(rownames(geneset_importance_df),"_",1)
ggplot(geneset_importance_df, aes(x=type, y=importance)) + geom_boxplot() + geom_signif(comparisons = list(c("inflam", "cell")),
                y_position = max(geneset_importance_df$importance) + 1,
                tip_length = 0.02)

write.table(sub_sig_geneset_interact_scores_df, , file="WGCNA_PPI_analysis/result/gene_set_interact//sub_sig_geneset_interact_scores_df.txt", row.names = F, sep="\t", quote = F)				
write.table(geneset_importance_df, file="WGCNA_PPI_analysis/result/gene_set_interact//geneset_importance_type.txt", row.names = F, sep="\t", quote = F)
##### 02End: Ranking of geneset importance  #####