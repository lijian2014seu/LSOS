##### description:
# Correlate gene expression log‑fold changes with methylation changes in specified gene regions:
#  - Load DMP and DEG results for two missions and multiple timepoint comparisons
#  - Annotate DMP probes with enhancer overlaps and promoter/genebody features
#  - For each comparison and gene group (all genes, DEGs, and baseline subgroups), compute
#    median probe Δβ per gene and correlate against gene expression logFC via Pearson’s r
#  - Summarize correlation statistics and plot bar charts of r and P‑values

# Working directory need to be adjusted according to your path.
setwd("/home/lqwang/Program")
source("./main/utils/annotate_enhancer_hg38.R")

# Command line args.
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-f", "--flag"), type="character", help="flag"),
  make_option(c("-e", "--exp1"), type="character", help="exp1"),
  make_option(c("-E", "--exp2"), type="character", help="exp2"),
  make_option(c("-r", "--gene_region"), type="character", help="gene_region")
)
parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
opts <- parse_args(parser)
flag <- opts$flag
exp1 <- opts$exp1
exp2 <- opts$exp2
gene_region <- opts$gene_region

library(ggplot2)
library(openxlsx)
library(rtracklayer)
library(GenomicRanges)

DMP_result_path = paste0("./DMP_result", flag)
Corr_result_path = paste0("./Corr_result", flag)
DEG_result_path = "DEG_result"

##### 01 Load DMP and DEG tables #####
df_DMP_M90 = readRDS(file.path(DMP_result_path, exp1, paste0("DMP_res_", exp1, ".Rds")))
df_DMP_M180_1 = readRDS(file.path(DMP_result_path, exp2, paste0("DMP_res_", exp2, ".Rds")))

# extract each time‐comparison
df_DMP_M90_T1_to_T2 <- df_DMP_M90[["T1 to T2"]]
df_DMP_M90_T2_to_T3 <- df_DMP_M90[["T2 to T3"]]
df_DMP_M90_T1_to_T3 <- df_DMP_M90[["T1 to T3"]]
df_DMP_M180_1_T1_to_T2 <- df_DMP_M180_1[["T1 to T2"]]
df_DMP_M180_1_T2_to_T3 <- df_DMP_M180_1[["T2 to T3"]]
df_DMP_M180_1_T1_to_T3 <- df_DMP_M180_1[["T1 to T3"]]

# Load DEG information
DEG_list_trends_M90 <- readRDS( file = file.path(DEG_result_path, "M90/DEG_list_trends.Rds" ))
DEG_T1_to_T3_M90 <- read.table( file = file.path(DEG_result_path, "M90/T1_to_T3_DEG_table.txt" ))
DEG_T1_to_T2_M90 <- read.table( file = file.path(DEG_result_path, "M90/T1_to_T2_DEG_table.txt" ))
DEG_T2_to_T3_M90 <- read.table( file = file.path(DEG_result_path, "M90/T2_to_T3_DEG_table.txt" ))
DEG_T1_to_T3_M90[ unlist(DEG_list_trends_M90), "change_process" ] <- unlist(lapply(names(DEG_list_trends_M90), 
                                                                                   function(x){rep(x, length(DEG_list_trends_M90[[x]]))}))
DEG_T1_to_T3_M90$change_process[ is.na(DEG_T1_to_T3_M90$change_process) ] <- "stable_stable"

DEG_list_trends_M180_1 <- readRDS( file = file.path(DEG_result_path, "M180-1/DEG_list_trends.Rds" ))
DEG_T1_to_T3_M180_1 <- read.table( file = file.path(DEG_result_path, "M180-1/T1_to_T3_DEG_table.txt" ))
DEG_T1_to_T2_M180_1 <- read.table( file = file.path(DEG_result_path, "M180-1/T1_to_T2_DEG_table.txt" ))
DEG_T2_to_T3_M180_1 <- read.table( file = file.path(DEG_result_path, "M180-1/T2_to_T3_DEG_table.txt" ))
DEG_T1_to_T3_M180_1[ unlist(DEG_list_trends_M180_1), "change_process" ] <- unlist(lapply(names(DEG_list_trends_M180_1), 
                                                                                         function(x){rep(x, length(DEG_list_trends_M180_1[[x]]))}))
DEG_T1_to_T3_M180_1$change_process[ is.na(DEG_T1_to_T3_M180_1$change_process) ] <- "stable_stable"


##### 02 Classification of DEGs #####
DEG_T1_to_T3_M90$baseline <- "Unchanged"
DEG_T1_to_T3_M90$baseline[ intersect(grep("_down", DEG_T1_to_T3_M90$change_process), which(DEG_T1_to_T3_M90$logFC < 0)) ] <- "DEG_over"
DEG_T1_to_T3_M90$baseline[ intersect(grep("_up", DEG_T1_to_T3_M90$change_process), which(DEG_T1_to_T3_M90$logFC > 0)) ] <- "DEG_over"
DEG_T1_to_T3_M90$baseline[ DEG_T1_to_T3_M90$change_process %in% c("up_stable", "down_stable") ] <- "DEG_unrecovered"
DEG_T1_to_T3_M90$baseline[ DEG_T1_to_T3_M90$baseline == "Unchanged" & DEG_T1_to_T3_M90$change_process != "stable_stable" ] <- "Other DEGs"
table(DEG_T1_to_T3_M90$baseline)

DEG_T1_to_T3_M180_1$baseline <- "Unchanged"
DEG_T1_to_T3_M180_1$baseline[ intersect(grep("_down", DEG_T1_to_T3_M180_1$change_process), which(DEG_T1_to_T3_M180_1$logFC < 0)) ] <- "DEG_over"
DEG_T1_to_T3_M180_1$baseline[ intersect(grep("_up", DEG_T1_to_T3_M180_1$change_process), which(DEG_T1_to_T3_M180_1$logFC > 0)) ] <- "DEG_over"
DEG_T1_to_T3_M180_1$baseline[ DEG_T1_to_T3_M180_1$change_process %in% c("up_stable", "down_stable") ] <- "DEG_unrecovered"
DEG_T1_to_T3_M180_1$baseline[ DEG_T1_to_T3_M180_1$baseline == "Unchanged" & DEG_T1_to_T3_M180_1$change_process != "stable_stable" ] <- "Other DEGs"
table(DEG_T1_to_T3_M180_1$baseline)

promoter_region = c("TSS1500", "TSS200", "5'UTR", "1stExon")

res_df <- c()

# Define the DEG and DMP data frame pairs and time comparisons
comparisons <- list(
  list(DEG = "DEG_T1_to_T2_M90", DMP = "df_DMP_M90_T1_to_T2", M = "M90", time = "T2 vs T1"),
  list(DEG = "DEG_T1_to_T2_M180_1", DMP = "df_DMP_M180_1_T1_to_T2", M = "M180-1", time = "T2 vs T1"),
  list(DEG = "DEG_T2_to_T3_M90", DMP = "df_DMP_M90_T2_to_T3", M = "M90", time = "T3 vs T2"),
  list(DEG = "DEG_T2_to_T3_M180_1", DMP = "df_DMP_M180_1_T2_to_T3", M = "M180-1", time = "T3 vs T2"),
  list(DEG = "DEG_T1_to_T3_M90", DMP = "df_DMP_M90_T1_to_T3", M = "M90", time = "T3 vs T1"),
  list(DEG = "DEG_T1_to_T3_M180_1", DMP = "df_DMP_M180_1_T1_to_T3", M = "M180-1", time = "T3 vs T1")
)

baseline_types <- c("Unchanged", "DEG_over", "DEG_unrecovered", "Other DEGs")

# Loop through each comparison
for (comp in comparisons) {
  M_type <- comp$M; time_point <- comp$time
  DEG_df <- get(comp$DEG); DMP_df <- get(comp$DMP)

  DMP_df <- annotate_enhancer_for_probes(DMP_df)

  DEG_T1_to_T3_df <- get(paste0("DEG_T1_to_T3_", ifelse(M_type == "M90", "M90", "M180_1")))

  # Correlation for all genes
  all_promoter_probes <- rownames(DMP_df)[DMP_df$gene %in% rownames(DEG_df) & DMP_df$feature %in% promoter_region]
  
  all_enhancer_probes <- rownames(DMP_df)[DMP_df$gene %in% rownames(DEG_df) & DMP_df$in_enhancer == TRUE]
  
  all_genebody_probes <- rownames(DMP_df)[DMP_df$gene %in% rownames(DEG_df) & DMP_df$feature %in% c("Body")]
  

  if (gene_region == "promoter") {
    sub_probe_annotation <- DMP_df[all_promoter_probes, c("gene", "feature", "deltaBeta")]
  } else if (gene_region == "enhancer") {
    sub_probe_annotation <- DMP_df[all_enhancer_probes, c("gene", "feature", "deltaBeta")]
  } else if (gene_region == "genebody") {
    sub_probe_annotation <- DMP_df[all_genebody_probes, c("gene", "feature", "deltaBeta")]
  } else {
    stop("!!!")
  }
  
  gene_promoter_deltaBeta_all <- tapply(sub_probe_annotation$deltaBeta, as.character(sub_probe_annotation$gene), median)
  genes_all <- intersect(names(gene_promoter_deltaBeta_all), rownames(DEG_df))
  
  if (length(genes_all) > 1) {
    gene_promoter_deltaBeta_all <- gene_promoter_deltaBeta_all[genes_all]
    gene_expression_logFC_all <- DEG_df[genes_all, "logFC"]
    df_all <- data.frame(gene_expression_logFC = gene_expression_logFC_all, gene_promoter_deltaBeta = gene_promoter_deltaBeta_all)
    res_all <- cor.test(df_all$gene_expression_logFC, df_all$gene_promoter_deltaBeta)
    res_df <- rbind(res_df, c(res_all$statistic, res_all$estimate, res_all$p.value, res_all$parameter, res_all$conf.int))
    rownames(res_df)[nrow(res_df)] <- paste("All genes in", M_type, time_point, sep = " ")
  } else {
    cat(paste("Not enough genes for correlation (All genes) in", M_type, time_point, "\n", sep = " "))
    stop()
  }


  # Correlation for DEGs (baseline != "Unchanged")
  de_genes <- rownames(DEG_T1_to_T3_df)[DEG_T1_to_T3_df$baseline != "Unchanged"]
  promoter_probes_de <- rownames(DMP_df)[DMP_df$gene %in% rownames(DEG_df) & DMP_df$feature %in% promoter_region]

  enhancer_probes_de <- rownames(DMP_df)[DMP_df$gene %in% rownames(DEG_df) & DMP_df$in_enhancer == TRUE]

  if (gene_region == "promoter") {
    sub_probe_annotation_de <- DMP_df[all_promoter_probes, c("gene", "feature", "deltaBeta")]
  } else if (gene_region == "enhancer") {
    sub_probe_annotation_de <- DMP_df[all_enhancer_probes, c("gene", "feature", "deltaBeta")]
  } else if (gene_region == "genebody") {
    sub_probe_annotation_de <- DMP_df[all_genebody_probes, c("gene", "feature", "deltaBeta")]
  } else {
    stop("!!!")
  }

  
  gene_promoter_deltaBeta_de <- tapply(sub_probe_annotation_de$deltaBeta, as.character(sub_probe_annotation_de$gene), median)
  genes_de <- intersect(names(gene_promoter_deltaBeta_de), rownames(DEG_df))
  genes_de <- intersect(genes_de, de_genes)

  if (length(genes_de) > 1) {
    gene_promoter_deltaBeta_de <- gene_promoter_deltaBeta_de[genes_de]
    gene_expression_logFC_de <- DEG_df[genes_de, "logFC"]
    df_de <- data.frame(gene_expression_logFC = gene_expression_logFC_de, gene_promoter_deltaBeta = gene_promoter_deltaBeta_de)
    res_de <- cor.test(df_de$gene_expression_logFC, df_de$gene_promoter_deltaBeta)
    res_df <- rbind(res_df, c(res_de$statistic, res_de$estimate, res_de$p.value, res_de$parameter, res_de$conf.int))
    rownames(res_df)[nrow(res_df)] <- paste("DEGs in", M_type, time_point, sep = " ")
  } else {
    cat(paste("Not enough genes for correlation (DEGs) in", M_type, time_point, "\n"))
  }

  # Loop through baseline types
  for (baseline in baseline_types) {
    genes0 <- rownames(DEG_T1_to_T3_df)[DEG_T1_to_T3_df$baseline == baseline]
    promoter_probes_baseline <- rownames(DMP_df)[DMP_df$gene %in% rownames(DEG_df) & DMP_df$feature %in% promoter_region]
    
    enhancer_probes_baseline <- rownames(DMP_df)[DMP_df$gene %in% rownames(DEG_df) & DMP_df$in_enhancer == TRUE]

    if (gene_region == "promoter") {
      sub_probe_annotation_baseline <- DMP_df[all_promoter_probes, c("gene", "feature", "deltaBeta")]
    } else if (gene_region == "enhancer") {
      sub_probe_annotation_baseline <- DMP_df[all_enhancer_probes, c("gene", "feature", "deltaBeta")]
    } else if (gene_region == "genebody") {
      sub_probe_annotation_baseline <- DMP_df[all_genebody_probes, c("gene", "feature", "deltaBeta")]
    } else {
      stop("!!!")
    }
    
    gene_promoter_deltaBeta_baseline <- tapply(sub_probe_annotation_baseline$deltaBeta, as.character(sub_probe_annotation_baseline$gene), median)
    genes_baseline <- intersect(names(gene_promoter_deltaBeta_baseline), rownames(DEG_df))
    genes_baseline <- intersect(genes_baseline, genes0)

    if (length(genes_baseline) > 1) {
      gene_promoter_deltaBeta_baseline <- gene_promoter_deltaBeta_baseline[genes_baseline]
      gene_expression_logFC_baseline <- DEG_df[genes_baseline, "logFC"]
      df_baseline <- data.frame(gene_expression_logFC = gene_expression_logFC_baseline, gene_promoter_deltaBeta = gene_promoter_deltaBeta_baseline)
      res_baseline <- cor.test(df_baseline$gene_expression_logFC, df_baseline$gene_promoter_deltaBeta)
      res_df <- rbind(res_df, c(res_baseline$statistic, res_baseline$estimate, res_baseline$p.value, res_baseline$parameter, res_baseline$conf.int))
      rownames(res_df)[nrow(res_df)] <- paste(gsub("_", " ", baseline), "genes in", M_type, time_point, sep = " ")
    } else {
      cat(paste("Not enough genes for correlation (", baseline, ") in", M_type, time_point, "\n"))
    }
  }
}

colnames(res_df) <- c("statistic", "coefficient", "P-value", "parameter", "conf.int_bottom ", "conf.int_top")
res_df1 <- format(res_df, scientific = TRUE, digits = 3)
res_df1 <- as.data.frame(res_df1)
res_df1$coefficient <- as.numeric(res_df1$coefficient)
res_df1$`P-value` <- as.numeric(res_df1$`P-value`)
res_df1 <- as.data.frame(res_df1)
res_df1$gene_group <- rownames(res_df1)

# Plotting bar charts to show the results of correlation analyses
show_cols_T2_vs_T1 = c("All genes in M180-1 T2 vs T1", "All genes in M90 T2 vs T1",
                       "DEGs in M90 T2 vs T1", "DEGs in M180-1 T2 vs T1", 
                       "Unchanged genes in M180-1 T2 vs T1", "Unchanged genes in M90 T2 vs T1")

p_T2_vs_T1 <- ggplot( res_df1[show_cols_T2_vs_T1, ], aes(x = gsub(" T2 vs T1", "", gene_group), y = coefficient, fill = `P-value`) ) +
    geom_bar(stat = "identity") + geom_text(aes( label = paste0("r = ", coefficient) ), vjust = 1, size = 3) +
    theme( axis.text.x = element_text(angle = 45, hjust = 1) ) +
	geom_text( aes(label = paste0("P = ", `P-value`)), vjust = 2, size = 3 ) +
	labs(title = paste0("Correlation of gene expression and ", gene_region, " methylation in adaptation"), x = "gene group")


show_cols_T3_vs_T2 = c("DEGs in M180-1 T3 vs T2", "DEGs in M90 T3 vs T2",
                       "DEG unrecovered genes in M180-1 T3 vs T2", "DEG unrecovered genes in M90 T3 vs T2", 
                       "Other DEGs genes in M180-1 T3 vs T2", "Other DEGs genes in M90 T3 vs T2", 
                       "DEG over genes in M180-1 T3 vs T2", "DEG over genes in M90 T3 vs T2")
p_T3_vs_T2 <- ggplot(res_df1[show_cols_T3_vs_T2, ], aes(x = gsub(" T3 vs T2", "", gene_group), y = coefficient, fill = `P-value`)) +
    geom_bar(stat = "identity") + geom_text(aes(label = paste0("r = ", coefficient)), vjust = 1, size = 3) +
    theme( axis.text.x = element_text(angle = 45, hjust = 1) ) +
	geom_text( aes(label = paste0("P = ", `P-value`)), vjust = 2, size = 3 )+
	labs(title = paste0("Correlation of gene expression and ", gene_region, " methylation in recovery"), x="gene group")


show_cols_T3_vs_T1 = c("DEGs in M180-1 T3 vs T1", "DEGs in M90 T3 vs T1",
                       "DEG unrecovered genes in M180-1 T3 vs T1", "DEG unrecovered genes in M90 T3 vs T1", 
                       "Other DEGs genes in M180-1 T3 vs T1", "Other DEGs genes in M90 T3 vs T1", 
                       "DEG over genes in M180-1 T3 vs T1", "DEG over genes in M90 T3 vs T1")
p_T3_vs_T1 <- ggplot(res_df1[show_cols_T3_vs_T1, ], aes(x = gsub(" T3 vs T1", "", gene_group), y = coefficient, fill = `P-value`)) +
    geom_bar(stat = "identity") + geom_text(aes(label = paste0("r = ", coefficient)), vjust = 1, size = 3) +
    theme( axis.text.x = element_text(angle = 45, hjust = 1) ) +
	geom_text( aes(label = paste0("P = ", `P-value`)), vjust = 2, size = 3 )+
	labs(title = paste0("Correlation of gene expression and ", gene_region, " methylation in recovery"), x="gene group")

ggsave(filename = paste0("Corr_Plot_adaptation_", gene_region, ".pdf"), plot = p_T2_vs_T1, path = Corr_result_path, width = 10, height = 8, device = "pdf")
ggsave(filename = paste0("Corr_Plot_recovery_", gene_region, ".pdf"), plot = p_T3_vs_T2, path = Corr_result_path, width = 10, height = 8, device = "pdf")
ggsave(filename = paste0("Corr_Plot_recovery_T3vsT1", gene_region, ".pdf"), plot = p_T3_vs_T1, path = Corr_result_path, width = 10, height = 8, device = "pdf")
