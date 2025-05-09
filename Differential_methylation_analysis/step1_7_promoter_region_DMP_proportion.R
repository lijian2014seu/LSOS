### description: Count up/down DMPs by genomic feature (gene region) and export pie charts; compile feature counts into tables, compute promoter‑region percentages, and plot a bar chart of promoter DMP proportions across experiments. ###

# Working directory need to be adjusted according to your path.
setwd("/home/lqwang/Program")
source("./main/utils/annotate_enhancer_hg38.R")

# Command line args.
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-f", "--flag"), type="character", help="flag")
)
parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)
flag <- opts$flag

message("step1_7_promoter_region_DMP_proportion.R", flag)

library(ggplot2)

DMP_result_path = paste0("./DMP_result", flag)

missions = list(list(name = "M90", 
                     beta_path = paste0("./90day_spaceflight_M90/DNA_methylation/beta/nobatch_beta", flag, ".Rds"), 
                     sample_sheet_path = "./90day_spaceflight_M90/DNA_methylation/rawdata/SampleSheet.csv"),
                list(name = "M180-1",
                     beta_path = paste0("./180day_spaceflight_M180_1/DNA_methylation/beta/nobatch_beta", flag, ".Rds"),
                     sample_sheet_path = "./180day_spaceflight_M180_1/DNA_methylation/rawdata/SampleSheet.csv")
                )
comparisons = list(c("T1", "T2"), c("T2", "T3"))
gene_region_count_list = list()
gene_region = c("TSS1500", "TSS200", "3'UTR", "5'UTR", "1stExon", "Body", "IGR")
enhancer_proportion_list = list()

for (exp in missions) {
  for (comp in comparisons) {
    df_DMP_T1_to_T2 <- readRDS(file.path(DMP_result_path, 
                                         exp$name, 
                                         paste0("DMP_res_", exp$name, ".Rds")))[[paste(comp[[1]], "to", comp[[2]], sep = " ")]]
    df_DMP_T1_to_T2 <- annotate_enhancer_for_probes(df_DMP_T1_to_T2)

    

    up_DMP_T1_to_T2 <- df_DMP_T1_to_T2[df_DMP_T1_to_T2$change_deltaBeta0.05 %in% c("up"), ]
    up_enhancer_num <- table(up_DMP_T1_to_T2$in_enhancer)[["TRUE"]]
    up_num <- nrow(up_DMP_T1_to_T2)
    enhancer_proportion_list[[paste("up", comp[[2]], "vs", comp[[1]], exp$name, "enhancer_proportion", sep = "_")]] = up_enhancer_num / up_num

    up_T1_to_T2_DMP_geneRegion <- sapply(gene_region, function(x) {
      temp <- up_DMP_T1_to_T2[grepl(x, up_DMP_T1_to_T2$feature), ]
      return (dim(temp)[1])
    })
    up_T1_to_T2_geneRegion_count <- up_T1_to_T2_DMP_geneRegion[up_T1_to_T2_DMP_geneRegion != 0]
    gene_region_count_list[[paste("up", comp[[2]], "vs", comp[[1]], exp$name, "geneRegion_count", sep = "_")]] = up_T1_to_T2_geneRegion_count



    down_DMP_T1_to_T2 <- df_DMP_T1_to_T2[df_DMP_T1_to_T2$change_deltaBeta0.05 %in% c("down"), ]
    down_enhancer_num <- table(down_DMP_T1_to_T2$in_enhancer)[["TRUE"]]
    down_num <- nrow(down_DMP_T1_to_T2)
    enhancer_proportion_list[[paste("down", comp[[2]], "vs", comp[[1]], exp$name, "enhancer_proportion", sep = "_")]] = down_enhancer_num / down_num

    down_T1_to_T2_DMP_geneRegion <- sapply(gene_region, function (x) {
      temp <- down_DMP_T1_to_T2[grepl(x, down_DMP_T1_to_T2$feature), ]
      return (dim(temp)[1])
    })
    down_T1_to_T2_geneRegion_count <- down_T1_to_T2_DMP_geneRegion[down_T1_to_T2_DMP_geneRegion != 0]
    gene_region_count_list[[paste("down", comp[[2]], "vs", comp[[1]], exp$name, "geneRegion_count", sep = "_")]] = down_T1_to_T2_geneRegion_count
  }
}



rowname_df_geneRegion = c("down_T2_vs_T1_M90_geneRegion_count", "up_T2_vs_T1_M90_geneRegion_count",
               "down_T2_vs_T1_M180-1_geneRegion_count", "up_T2_vs_T1_M180-1_geneRegion_count",
               "down_T3_vs_T2_M90_geneRegion_count", "up_T3_vs_T2_M90_geneRegion_count",
               "down_T3_vs_T2_M180-1_geneRegion_count", "up_T3_vs_T2_M180-1_geneRegion_count")
geneRegion_count_df <- rbind(gene_region_count_list[[rowname_df_geneRegion[1]]],
                             gene_region_count_list[[rowname_df_geneRegion[2]]], 
                             gene_region_count_list[[rowname_df_geneRegion[3]]],
                             gene_region_count_list[[rowname_df_geneRegion[4]]],                             
                             gene_region_count_list[[rowname_df_geneRegion[5]]], 
                             gene_region_count_list[[rowname_df_geneRegion[6]]],
                             gene_region_count_list[[rowname_df_geneRegion[7]]],                              
                             gene_region_count_list[[rowname_df_geneRegion[8]]])
rownames(geneRegion_count_df) = rowname_df_geneRegion
write.table(geneRegion_count_df, file = file.path(DMP_result_path, "geneRegion_count_df_deltaBeta0.05.txt"), sep="\t", quote = FALSE)

geneRegion_count_percent <- apply(geneRegion_count_df, 1, function(x) {x / sum(x)})
promoter_DMP_percent_df <- data.frame(percent = apply(geneRegion_count_percent , 2, function (x) {sum(x[c("TSS1500", "TSS200", "1stExon", "5'UTR")])})) 
promoter_DMP_percent_df$change_type <- rep(c("down", "up"), 4)
promoter_DMP_percent_df$Experiment_time <- gsub("^(down_|up_)|_geneRegion_count$", "", rownames(promoter_DMP_percent_df))
time_levels <- unique(promoter_DMP_percent_df$Experiment_time)
promoter_DMP_percent_df$Experiment_time <- factor(
  promoter_DMP_percent_df$Experiment_time,
  levels = time_levels
)
p <- ggplot(promoter_DMP_percent_df, aes(x = Experiment_time, y = percent, fill = change_type )) + 
     geom_bar(stat = "identity", position = "dodge") + 
     scale_fill_manual(values = c(down = "#95BDE0", up = "#C1548E")) + 
     ylab("Percentage of promoter region probes (%)") + 
     geom_text(aes(label = round(percent, 3)), vjust = -0.7, position = position_dodge(0.9))
ggsave(filename = "bar_promoter_region.pdf", plot = p, path = DMP_result_path, width = 8, height = 8, device = "pdf")
cat(">> [DONE] Plot saved to", file.path(DMP_result_path, "bar_promoter_region.pdf"), "\n")



rowname_df_enhancer = c("down_T2_vs_T1_M90_enhancer_proportion", "up_T2_vs_T1_M90_enhancer_proportion",
               "down_T2_vs_T1_M180-1_enhancer_proportion", "up_T2_vs_T1_M180-1_enhancer_proportion",
               "down_T3_vs_T2_M90_enhancer_proportion", "up_T3_vs_T2_M90_enhancer_proportion",
               "down_T3_vs_T2_M180-1_enhancer_proportion", "up_T3_vs_T2_M180-1_enhancer_proportion")
enhancer_proportion_df <- rbind(enhancer_proportion_list[[rowname_df_enhancer[1]]],
                                enhancer_proportion_list[[rowname_df_enhancer[2]]], 
                                enhancer_proportion_list[[rowname_df_enhancer[3]]],
                                enhancer_proportion_list[[rowname_df_enhancer[4]]],                             
                                enhancer_proportion_list[[rowname_df_enhancer[5]]], 
                                enhancer_proportion_list[[rowname_df_enhancer[6]]],
                                enhancer_proportion_list[[rowname_df_enhancer[7]]],                              
                                enhancer_proportion_list[[rowname_df_enhancer[8]]]
                              )
rownames(enhancer_proportion_df) = rowname_df_enhancer
write.table(enhancer_proportion_df, file = file.path(DMP_result_path, "enhancer_proportion_df_deltaBeta0.05.txt"), sep="\t", quote = FALSE)

enhancer_DMP_percent_df <- data.frame(percent = enhancer_proportion_df)
enhancer_DMP_percent_df$change_type <- rep(c("down", "up"), 4)
enhancer_DMP_percent_df$Experiment_time <- gsub("^(down_|up_)|_enhancer_count$", "", rownames(enhancer_DMP_percent_df))
time_levels <- unique(enhancer_DMP_percent_df$Experiment_time)
enhancer_DMP_percent_df$Experiment_time <- factor(
  enhancer_DMP_percent_df$Experiment_time,
  levels = time_levels
)
p <- ggplot(enhancer_DMP_percent_df, aes(x = Experiment_time, y = percent, fill = change_type )) + 
     geom_bar(stat = "identity", position = "dodge") + 
     scale_fill_manual(values = c(down = "#95BDE0", up = "#C1548E")) + 
     ylab("Percentage of enhancer region probes (%)") + 
     geom_text(aes(label = round(percent, 3)), vjust = -0.7, position = position_dodge(0.9))
ggsave(filename = "bar_enhancer_region.pdf", plot = p, path = DMP_result_path, width = 8, height = 8, device = "pdf")
cat(">> [DONE] Plot saved to", file.path(DMP_result_path, "bar_enhancer_region.pdf"), "\n")
