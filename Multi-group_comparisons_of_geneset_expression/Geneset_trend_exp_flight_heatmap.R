# Description: Gene set expression change trends in M90, M180-1 and M180-2 (Extended Data Fig. 9, Extended Data Fig. 10 and Supplementary Fig. 7)

# Load packages
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(stringr)

# Load genesets enriched by PPI nodes 
go_enrich_res <- read.table( file = "WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_GOenrich_res.txt", sep="\t")
kegg_enrich_res <- read.table( file = "WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_KEGGenrich_res.txt", sep="\t")
geneset_all_df_PPI <- rbind(go_enrich_res, kegg_enrich_res)

# Load count data of 90-day experiment (M90),nobatch_expression_profile.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 90-day spaceflight_M90, RNA-seq)
combat_mat_01 <- read.table("90-day spaceflight_M90/RNA-seq/DE_analysis/nobatch_expression_profile.txt", header=T, row.names=1, as.is=T, sep=" ")
SampleSheet_01 <- read.csv("90-day spaceflight_M90/RNA-seq/DE_analysis/SampleSheet.csv")
combat_mat_01_log <- log2(combat_mat_01+1)
# Load count data of 180-day experiment (M180-1),nobatch_expression_profile.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 180-day spaceflight_180-1, RNA-seq)
combat_mat_02 <- read.table("180-day spaceflight_M180-1/RNA-seq/DE_analysis/nobatch_expression_profile.txt", header=T, row.names=1, as.is=T, sep=" ")
SampleSheet_02 <- read.csv("180-day spaceflight_M180-1/RNA-seq/DE_analysis/SampleSheet.csv")
combat_mat_02_log <- log2(combat_mat_02+1)


pattern_gene_list_01 <- readRDS(file="Differential methylation analysis/result/90-day spaceflight_M90/pattern_gene_list.Rds")
pattern_gene_list_02 <- readRDS(file="Differential methylation analysis/result/180-day spaceflight_M180-1/pattern_gene_list.Rds")
DEG_list_trends_01 <- readRDS("90-day spaceflight_M90/RNA-seq/DE_analysis/DEG_list_trends.Rds")
DEG_list_trends_02 <- readRDS("180-day spaceflight_M180-1/RNA-seq/DE_analysis/DEG_list_trends.Rds")


DEG_T1_to_T3_01 <- read.table( file = "90-day spaceflight_M90/RNA-seq/DE_analysis/T1_to_T3_DEG_table.txt", header = T, sep = "\t", row.names = 1)
DEG_T1_to_T3_01[unlist(DEG_list_trends_01), "change_process"] <-   unlist(lapply(names(DEG_list_trends_01), function(x){rep(x, length(DEG_list_trends_01[[x]]))}))
DEG_T1_to_T3_01$change_process[is.na(DEG_T1_to_T3_01$change_process)] <- "stable_stable"
DEG_T1_to_T3_01$baseline <- "Unchanged"
DEG_T1_to_T3_01$baseline[intersect(grep("_down",DEG_T1_to_T3_01$change_process), which(DEG_T1_to_T3_01$logFC<0))] <- "DEG_over"
DEG_T1_to_T3_01$baseline[intersect(grep("_up",DEG_T1_to_T3_01$change_process), which(DEG_T1_to_T3_01$logFC>0))] <- "DEG_over"
DEG_T1_to_T3_01$baseline[DEG_T1_to_T3_01$change_process%in%c("up_stable","down_stable")] <- "DEG_unrecovered"
DEG_T1_to_T3_01$baseline[DEG_T1_to_T3_01$baseline=="Unchanged"&DEG_T1_to_T3_01$change_process!="stable_stable"] <- "Other DEGs"
table(DEG_T1_to_T3_01$baseline)


DEG_T1_to_T3_02 <- read.table( file = "180-day spaceflight_M180-1/RNA-seq/DE_analysis/T1_to_T3_DEG_table.txt", header = T, sep = "\t", row.names = 1)
DEG_T1_to_T3_02[unlist(DEG_list_trends_02), "change_process"] <-   unlist(lapply(names(DEG_list_trends_02), function(x){rep(x, length(DEG_list_trends_02[[x]]))}))
DEG_T1_to_T3_02$change_process[is.na(DEG_T1_to_T3_02$change_process)] <- "stable_stable"
DEG_T1_to_T3_02$baseline <- "Unchanged"
DEG_T1_to_T3_02$baseline[intersect(grep("_down",DEG_T1_to_T3_02$change_process), which(DEG_T1_to_T3_02$logFC<0))] <- "DEG_over"
DEG_T1_to_T3_02$baseline[intersect(grep("_up",DEG_T1_to_T3_02$change_process), which(DEG_T1_to_T3_02$logFC>0))] <- "DEG_over"
DEG_T1_to_T3_02$baseline[DEG_T1_to_T3_02$change_process%in%c("up_stable","down_stable")] <- "DEG_unrecovered"
DEG_T1_to_T3_02$baseline[DEG_T1_to_T3_02$baseline=="Unchanged"&DEG_T1_to_T3_02$change_process!="stable_stable"] <- "Other DEGs"
table(DEG_T1_to_T3_02$baseline)


write.table(DEG_T1_to_T3_01, "90-day spaceflight_M90/RNA-seq/DE_analysis/T1_to_T3_DEG_table_anno.txt",  sep = "\t", quote=F)
write.table(DEG_T1_to_T3_02, "180-day spaceflight_M180-1/RNA-seq/DE_analysis/T1_to_T3_DEG_table_anno.txt",  sep = "\t", quote=F)


#term <- "bone resorption"
outdir <- "Multi-group comparisons of geneset expression/result/"
library(pheatmap)
for(term in geneset_all_df_PPI$Description){
   interested_geneset <- geneset_all_df_PPI[geneset_all_df_PPI$Description==term, "geneID"]
   interested_geneset <- unique(unlist(strsplit(interested_geneset,"/")))
   term <- str_split_i(term, " - ", 1)
   word_count <- str_count(term, "\\w+")

   if(word_count > 10){
      first_ten_words <- str_extract_all(term, "\\w+")[[1]][1:10]
      trimmed_text <- paste(first_ten_words, collapse = " ")
	  term <- trimmed_text
	}
   sample_order <- order(SampleSheet_01$Sample_Group,SampleSheet_01$Sample_Well)
   sub_exp_mat <- combat_mat_01_log[interested_geneset,sample_order]
  anno_col <- SampleSheet_01[,c("Sample_Group", "Sample_Well")]
  rownames(anno_col) <-  paste0("X90_",SampleSheet_01$Sample_Name)
  colnames(anno_col) <- c("Time", "Subject")
  anno_col$Subject <- paste0(anno_col$Subject,"_90")
  exp_trend <- mapply(function(x){names(DEG_list_trends_01[mapply(function(y){x%in%y}, DEG_list_trends_01)])}, interested_geneset)
  exp_trend <- exp_trend[mapply(length,exp_trend)!=0]
  meth_trend <- mapply(function(x){names(pattern_gene_list_01[mapply(function(y){x%in%y}, pattern_gene_list_01)])}, interested_geneset)
  meth_trend <- meth_trend[mapply(length,meth_trend)!=0]
  anno_row <- data.frame(matrix(ncol = 2, nrow = nrow(sub_exp_mat)))
  rownames(anno_row) <- rownames(sub_exp_mat)
  colnames(anno_row) <- c("meth_trend", "exp_trend")
  anno_row$rebound <- DEG_T1_to_T3_01[rownames(anno_row),"baseline"]
  anno_row$meth_trend[match(names(meth_trend),rownames(anno_row))] <- as.character(meth_trend)
  anno_row$meth_trend[is.na(anno_row$meth_trend)] <- "stable_stable"
  anno_row$exp_trend <- "stable_stable"
  anno_row$exp_trend[match(names(exp_trend),rownames(anno_row))] <- as.character(exp_trend)
  outfile <- paste0("exp_flight01_",term,".pdf")
  my_palette <- c(colorRampPalette(c("#536BE1", "white"))(100),colorRampPalette(c("white", "#FE4C5B"))(100))
  my_color <- brewer.pal(9, "Pastel2")
  gaps_col <- c(3,6)
  num_gene1 <- nrow(sub_exp_mat)
  if(num_gene1<100){
         pheatmap::pheatmap(sub_exp_mat, annotation_col = anno_col, annotation_row = anno_row,cluster_cols = F,cluster_rows = F, scale = "row", color = my_palette, filename = paste0(outdir,outfile), gaps_col = gaps_col, annotation_colors = list(
           Time  = c("T1" = my_color[1], "T2" = my_color[2], "T3" = my_color[3]),
           Subject = c("H01_90" = my_color[4], "H02_90" = my_color[5], "H03_90" = my_color[6]),
		   meth_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7]),
		   exp_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7]),
		   rebound = c("DEG_unrecovered" = my_color[1], "Other DEGs" = my_color[9], "DEG_over" = my_color[4], "Unchanged" = my_color[3])))
  }else{
  pheatmap::pheatmap(sub_exp_mat, annotation_col = anno_col, annotation_row = anno_row,cluster_cols = F,cluster_rows = F, scale = "row", color = my_palette, filename = paste0(outdir,outfile), gaps_col = gaps_col, width =  1*log(num_gene1), height =  2.5*log2(num_gene1), annotation_colors = list(
           Time = c("T1" = my_color[1], "T2" = my_color[2], "T3" = my_color[3]),
           Subject = c("H01_90" = my_color[4], "H02_90" = my_color[5], "H03_90" = my_color[6]),
		   meth_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7]),
		   exp_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7]),
		   rebound = c("DEG_unrecovered" = my_color[1], "Other DEGs" = my_color[9], "DEG_over" = my_color[4], "Unchanged" = my_color[3])))
    }


  sample_order <- order(SampleSheet_02$Sample_Group,SampleSheet_02$Sample_Well)
  sub_exp_mat <- combat_mat_02_log[interested_geneset,sample_order]
  anno_col <- SampleSheet_02[,c("Sample_Group", "Sample_Well")]
  rownames(anno_col) <- paste0("X180_",SampleSheet_02$Sample_Name)
  colnames(anno_col) <- c("Time", "Subject")
  anno_col$Subject <- paste0(anno_col$Subject,"_180-1")
  exp_trend <- mapply(function(x){names(DEG_list_trends_02[mapply(function(y){x%in%y}, DEG_list_trends_02)])}, interested_geneset)
  exp_trend <- exp_trend[mapply(length,exp_trend)!=0]
  meth_trend <- mapply(function(x){names(pattern_gene_list_02[mapply(function(y){x%in%y}, pattern_gene_list_02)])}, interested_geneset)
  meth_trend <- meth_trend[mapply(length,meth_trend)!=0]
  anno_row <- data.frame(matrix(ncol = 2, nrow = nrow(sub_exp_mat)))
  rownames(anno_row) <- rownames(sub_exp_mat)
  colnames(anno_row) <- c("meth_trend", "exp_trend")
  anno_row$rebound <- DEG_T1_to_T3_02[rownames(anno_row),"baseline"]
  anno_row$meth_trend[match(names(meth_trend),rownames(anno_row))] <- as.character(meth_trend)
  anno_row$meth_trend[is.na(anno_row$meth_trend)] <- "stable_stable"
  anno_row$exp_trend <- "stable_stable"
  anno_row$exp_trend[match(names(exp_trend),rownames(anno_row))] <- as.character(exp_trend)
  outfile <- paste0("exp_flight02_",term,".pdf")
  my_color <- brewer.pal(9, "Pastel2")
  gaps_col <- c(3,6)
  num_gene1 <- nrow(sub_exp_mat)
  if(num_gene1<100){
         pheatmap::pheatmap(sub_exp_mat, annotation_col = anno_col, annotation_row = anno_row,cluster_cols = F,cluster_rows = F, scale = "row", color = my_palette, filename = paste0(outdir,outfile), gaps_col = gaps_col, annotation_colors = list(
           Time = c("T1" = my_color[1], "T2" = my_color[2], "T3" = my_color[3]),
           Subject = c("H01_180-1" = my_color[4], "H02_180-1" = my_color[5], "H03_180-1" = my_color[6]),
		   meth_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7]),
		   exp_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7]),
		   rebound = c("DEG_unrecovered" = my_color[1], "Other DEGs" = my_color[9], "DEG_over" = my_color[4], "Unchanged" = my_color[3])))
  }else{                      
  pheatmap::pheatmap(sub_exp_mat, annotation_col = anno_col, annotation_row = anno_row,cluster_cols = F,cluster_rows = F, scale = "row", color = my_palette, filename = paste0(outdir,outfile), gaps_col = gaps_col, width =  1*log(num_gene1), height =  2.5*log2(num_gene1), annotation_colors = list(
           Time = c("T1" = my_color[1], "T2" = my_color[2], "T3" = my_color[3]),
           Subject = c("H01_180-1" = my_color[4], "H02_180-1" = my_color[5], "H03_180-1" = my_color[6]),
		   meth_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7]),
		   exp_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7]),
		   rebound = c("DEG_unrecovered" = my_color[1], "Other DEGs" = my_color[9], "DEG_over" = my_color[4], "Unchanged" = my_color[3])))
    }

}


# Load count data of 180-day experiment (M180-2),nobatch_expression_profile.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 180-day spaceflight_180-2, RNA-seq)
combat_mat_03 <- read.table("180-day spaceflight_M180-2/RNA-seq/DE_analysis/nobatch_expression_profile.txt", header=T, row.names=1, as.is=T, sep=" ")
SampleSheet_03 <- read.csv("180-day spaceflight_M180-2/RNA-seq/DE_analysis/SampleSheet.csv")
SampleSheet_03 <- SampleSheet_03[!grepl("Torbit", SampleSheet_03$Sample_Group), ]
combat_mat_03 <- combat_mat_03[,SampleSheet_03$Sample_Name]
combat_mat_03_log <- log2(combat_mat_03+1)

pattern_gene_list_03 <- readRDS(file="Differential methylation analysis/result/180-day spaceflight_M180-1/pattern_gene_list.Rds")
DEG_list_trends_03 <- readRDS("180-day spaceflight_M180-2/RNA-seq/DE_analysis/DEG_list_trends.Rds")

DEG_T1_to_T3_03 <- read.table( file = "180-day spaceflight_M180-2/RNA-seq/DE_analysis/T1_to_T3_DEG_table.txt", header = T, sep = "\t", row.names = 1)
DEG_T1_to_T3_03[unlist(DEG_list_trends_03), "change_process"] <-   unlist(lapply(names(DEG_list_trends_03), function(x){rep(x, length(DEG_list_trends_03[[x]]))}))
DEG_T1_to_T3_03$change_process[is.na(DEG_T1_to_T3_03$change_process)] <- "stable_stable"
DEG_T1_to_T3_03$baseline <- "Unchanged"
DEG_T1_to_T3_03$baseline[intersect(grep("_down",DEG_T1_to_T3_03$change_process), which(DEG_T1_to_T3_03$logFC<0))] <- "DEG_over"
DEG_T1_to_T3_03$baseline[intersect(grep("_up",DEG_T1_to_T3_03$change_process), which(DEG_T1_to_T3_03$logFC>0))] <- "DEG_over"
DEG_T1_to_T3_03$baseline[DEG_T1_to_T3_03$change_process%in%c("up_stable","down_stable")] <- "DEG_unrecovered"
DEG_T1_to_T3_03$baseline[DEG_T1_to_T3_03$baseline=="Unchanged"&DEG_T1_to_T3_03$change_process!="stable_stable"] <- "Other DEGs"
table(DEG_T1_to_T3_03$baseline)


write.table(DEG_T1_to_T3_03, "180-day spaceflight_M180-2/RNA-seq/DE_analysis/T1_to_T3_DEG_table_anno.txt",  sep = "\t", quote=F)

library(RColorBrewer)
library(stringr)
#term <- "bone resorption"
outdir <- "geneset/"
library(pheatmap)
for(term in geneset_all_df_PPI$Description){
   interested_geneset <- geneset_all_df_PPI[geneset_all_df_PPI$Description==term, "geneID"]
   interested_geneset <- unique(unlist(strsplit(interested_geneset,"/")))
   term <- str_split_i(term, " - ", 1)
   word_count <- str_count(term, "\\w+")

   if(word_count > 10){
      first_ten_words <- str_extract_all(term, "\\w+")[[1]][1:10]
      trimmed_text <- paste(first_ten_words, collapse = " ")
	  term <- trimmed_text
	}

  sample_order <- order(SampleSheet_03$Sample_Group,SampleSheet_03$Sample_Well)
  sub_exp_mat <- combat_mat_03_log[interested_geneset,sample_order]
  anno_col <- SampleSheet_03[,c("Sample_Group", "Sample_Well")]
  rownames(anno_col) <- SampleSheet_03$Sample_Name
  colnames(anno_col) <- c("Time", "Subject")
  anno_col$Subject <- paste0(anno_col$Subject,"_180-2")
  exp_trend <- mapply(function(x){names(DEG_list_trends_03[mapply(function(y){x%in%y}, DEG_list_trends_03)])}, interested_geneset)
  exp_trend <- exp_trend[mapply(length,exp_trend)!=0]
  meth_trend <- mapply(function(x){names(pattern_gene_list_03[mapply(function(y){x%in%y}, pattern_gene_list_03)])}, interested_geneset)
  meth_trend <- meth_trend[mapply(length,meth_trend)!=0]
  anno_row <- data.frame(matrix(ncol = 2, nrow = nrow(sub_exp_mat)))
  rownames(anno_row) <- rownames(sub_exp_mat)
  colnames(anno_row) <- c("meth_trend", "exp_trend")
  anno_row$rebound <- DEG_T1_to_T3_03[rownames(anno_row),"baseline"]
  anno_row$meth_trend[match(names(meth_trend),rownames(anno_row))] <- as.character(meth_trend)
  anno_row$meth_trend[is.na(anno_row$meth_trend)] <- "stable_stable"
  anno_row$exp_trend <- "stable_stable"
  anno_row$exp_trend[match(names(exp_trend),rownames(anno_row))] <- as.character(exp_trend)
  outfile <- paste0( "180-2/exp_flight03_",term,".pdf")
  my_palette <- c(colorRampPalette(c("#536BE1", "white"))(100),colorRampPalette(c("white", "#FE4C5B"))(100))
  my_color <- brewer.pal(9, "Pastel2")
  gaps_col <- c(3,6)
  num_gene1 <- nrow(sub_exp_mat)
  if(num_gene1<100){
         pheatmap::pheatmap(sub_exp_mat, annotation_col = anno_col, annotation_row = anno_row,cluster_cols = F,cluster_rows = F, scale = "row", color = my_palette, filename = paste0(outdir,outfile), gaps_col = gaps_col, annotation_colors = list(
           Time = c("T1" = my_color[1], "T2" = my_color[2], "T3" = my_color[3]),
           Subject = c("H01_180-2" = my_color[4], "H02_180-2" = my_color[5], "H03_180-2" = my_color[6]),
		   meth_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7]),
		   exp_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7]),
		   rebound = c("DEG_unrecovered" = my_color[1], "Other DEGs" = my_color[9], "DEG_over" = my_color[4], "Unchanged" = my_color[3])))
  }else{                      
  pheatmap::pheatmap(sub_exp_mat, annotation_col = anno_col, annotation_row = anno_row,cluster_cols = F,cluster_rows = F, scale = "row", color = my_palette, filename = paste0(outdir,outfile), gaps_col = gaps_col, width =  1*log(num_gene1), height =  2.5*log2(num_gene1), annotation_colors = list(
           Time = c("T1" = my_color[1], "T2" = my_color[2], "T3" = my_color[3]),
           Subject = c("H01_180-2" = my_color[4], "H03_180-2" = my_color[5], "H03_180-2" = my_color[6]),
		   meth_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7]),
		   exp_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7]),
		   rebound = c("DEG_unrecovered" = my_color[1], "Other DEGs" = my_color[9], "DEG_over" = my_color[4], "Unchanged" = my_color[3])))
    }

}
