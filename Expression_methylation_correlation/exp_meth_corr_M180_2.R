##### description:
# Correlate gene expression log‑fold changes with methylation changes in specified gene regions for a single experiment:
#  - Load DMP (differentially methylated positions) and DEG (differentially expressed genes) results for timepoint comparisons T1→T2 and T2→T3 for M180‑2  
#  - Classify DEGs into “DEG_over”, “DEG_unrecovered”, “Other DEGs”, and “Unchanged” based on expression trends and logFC  
#  - Annotate probes with promoter, enhancer, or gene body context via `annotate_enhancer_for_probes()`  
#  - For all genes and each DEG subgroup, compute median Δβ per gene and perform Pearson correlation against gene expression logFC  
#  - Summarize correlation statistics and generate bar plots of correlation coefficients and P‑values :contentReference[oaicite:0]{index=0}

# Working directory need to be adjusted according to your path.
setwd("/home/lqwang/Program")
source("./main/utils/annotate_enhancer_hg38.R")

# Command line args.
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-f", "--flag"), type="character", help="flag"),
  make_option(c("-e", "--exp"), type="character", help="exp"),
  make_option(c("-r", "--gene_region"), type="character", help="gene_region")
)
parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
opts <- parse_args(parser)
flag <- opts$flag
exp <- opts$exp
gene_region <- opts$gene_region

library(ggplot2)
library(openxlsx)
library(rtracklayer)
library(GenomicRanges)

DMP_result_path = paste0("./DMP_result", flag)
Corr_result_path = paste0("./Corr_result", flag)
DEG_result_path = "DEG_result"

##### 01 Load DMP and DEG tables #####
df_DMP_M180_2 = readRDS(file.path(DMP_result_path, exp, paste0("DMP_res_", exp, ".Rds")))

df_DMP_M180_2_T1_to_T2 <- df_DMP_M180_2[["T1 to T2"]]
df_DMP_M180_2_T2_to_T3 <- df_DMP_M180_2[["T2 to T3"]]

# Load DEG information
DEG_list_trends_M180_2 <- readRDS( file = file.path(DEG_result_path, exp, "DEG_list_trends.Rds" ))
DEG_T1_to_T3_M180_2 <- read.table( file = file.path(DEG_result_path, exp, "T1_to_T3_DEG_table.txt" ))
DEG_T1_to_T2_M180_2 <- read.table( file = file.path(DEG_result_path, exp, "T1_to_T2_DEG_table.txt" ))
DEG_T2_to_T3_M180_2 <- read.table( file = file.path(DEG_result_path, exp, "T2_to_T3_DEG_table.txt" ))
DEG_T1_to_T3_M180_2[ unlist(DEG_list_trends_M180_2), "change_process" ] <- unlist(lapply(names(DEG_list_trends_M180_2), 
                                                                                         function(x){rep(x, length(DEG_list_trends_M180_2[[x]]))}))
DEG_T1_to_T3_M180_2$change_process[ is.na(DEG_T1_to_T3_M180_2$change_process) ] <- "stable_stable"


##### 02 Classification of DEGs #####
DEG_T1_to_T3_M180_2$baseline <- "Unchanged"
DEG_T1_to_T3_M180_2$baseline[ intersect(grep("_down", DEG_T1_to_T3_M180_2$change_process), which(DEG_T1_to_T3_M180_2$logFC < 0)) ] <- "DEG_over"
DEG_T1_to_T3_M180_2$baseline[ intersect(grep("_up", DEG_T1_to_T3_M180_2$change_process), which(DEG_T1_to_T3_M180_2$logFC > 0)) ] <- "DEG_over"
DEG_T1_to_T3_M180_2$baseline[ DEG_T1_to_T3_M180_2$change_process %in% c("up_stable", "down_stable") ] <- "DEG_unrecovered"
DEG_T1_to_T3_M180_2$baseline[ DEG_T1_to_T3_M180_2$baseline == "Unchanged" & DEG_T1_to_T3_M180_2$change_process != "stable_stable" ] <- "Other DEGs"
table(DEG_T1_to_T3_M180_2$baseline)

promoter_region = c("TSS1500", "TSS200", "5'UTR", "1stExon")

res_df <- c()

# Define the DEG and DMP data frame pairs and time comparisons
comparisons <- list(
  list(DEG = "DEG_T1_to_T2_M180_2", DMP = "df_DMP_M180_2_T1_to_T2", M = "M180-2", time = "T2 vs T1"),
  list(DEG = "DEG_T2_to_T3_M180_2", DMP = "df_DMP_M180_2_T2_to_T3", M = "M180-2", time = "T3 vs T2")
)

baseline_types <- c("Unchanged", "DEG_over", "DEG_unrecovered", "Other DEGs")

# Loop through each comparison
for (comp in comparisons) {
  M_type <- comp$M; time_point <- comp$time
  DEG_df <- get(comp$DEG); DMP_df <- get(comp$DMP)

  DMP_df <- annotate_enhancer_for_probes(DMP_df)

  DEG_T1_to_T3_df <- get(paste0("DEG_T1_to_T3_", ifelse(exp == "M180-2", "M180_2", NULL)))

  # Correlation for all genes
  all_promoter_probes <- rownames(DMP_df)[DMP_df$gene %in% rownames(DEG_df) & DMP_df$feature %in% promoter_region]
  
  all_enhancer_probes <- rownames(DMP_df)[DMP_df$gene %in% rownames(DEG_df) & DMP_df$in_enhancer == TRUE]
  
  if (gene_region == "promoter") {
    sub_probe_annotation <- DMP_df[all_promoter_probes, c("gene", "feature", "deltaBeta")]
  } else if (gene_region == "enhancer") {
    sub_probe_annotation <- DMP_df[all_enhancer_probes, c("gene", "feature", "deltaBeta")]
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
    sub_probe_annotation_de <- DMP_df[promoter_probes_de, c("gene", "feature", "deltaBeta")]
  } else if (gene_region == "enhancer") {
    sub_probe_annotation_de <- DMP_df[enhancer_probes_de, c("gene", "feature", "deltaBeta")]
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
    rownames(res_df)[nrow(res_df)] <- paste("DEGs in", M_type, time_point)
  } else {
    cat(paste("Not enough genes for correlation (DEGs) in", M_type, time_point, "\n"))
  }

  # Loop through baseline types
  for (baseline in baseline_types) {
    genes0 <- rownames(DEG_T1_to_T3_df)[DEG_T1_to_T3_df$baseline == baseline]
    promoter_probes_baseline <- rownames(DMP_df)[DMP_df$gene %in% rownames(DEG_df) & DMP_df$feature %in% promoter_region]
    
    enhancer_probes_baseline <- rownames(DMP_df)[DMP_df$gene %in% rownames(DEG_df) & DMP_df$in_enhancer == TRUE]

    if (gene_region == "promoter") {
      sub_probe_annotation_baseline <- DMP_df[promoter_probes_baseline, c("gene", "feature", "deltaBeta")]
    } else if (gene_region == "enhancer") {
      sub_probe_annotation_baseline <- DMP_df[enhancer_probes_baseline, c("gene", "feature", "deltaBeta")]
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
      rownames(res_df)[nrow(res_df)] <- paste(gsub("_", " ", baseline), "genes in", M_type, time_point)
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
show_cols_T2_vs_T1 = c("All genes in M180-2 T2 vs T1",
                       "DEGs in M180-2 T2 vs T1", 
                       "Unchanged genes in M180-2 T2 vs T1")

p_T2_vs_T1 <- ggplot( res_df1[show_cols_T2_vs_T1, ], aes(x = gsub(" T2 vs T1", "", gene_group), y = coefficient, fill = `P-value`) ) +
    geom_bar(stat = "identity") + geom_text(aes( label = paste0("r = ", coefficient) ), vjust = 1, size = 3) +
    theme( axis.text.x = element_text(angle = 45, hjust = 1) ) +
	geom_text( aes(label = paste0("P = ", `P-value`)), vjust = 2, size = 3 ) +
	labs(title = paste0("Correlation of gene expression and ", gene_region, " methylation in adaptation"), x = "gene group")

ggsave(filename = paste0("Corr_Plot_adaptation ", gene_region, "_M180_2.pdf"), plot = p_T2_vs_T1, path = Corr_result_path, width = 10, height = 8, device = "pdf")
