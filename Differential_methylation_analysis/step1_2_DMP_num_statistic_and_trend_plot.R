### description: Statistics on the number of DMPs and presentation of the changing trend of DMPs in T2 vs T1 ###

# Working directory need to be adjusted according to your path.
setwd("/home/lqwang/Program")

# Command line args.
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-f", "--flag"), type="character", help="flag")
)
parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)
flag <- opts$flag

library(RColorBrewer)
library(ggplot2)
library(ggalluvial)

message("step1_2_DMP_num_statistic_and_trend_plot.R", flag)

if (flag %in% c("_3_3", "_6", "_1_1_1_1_1_1")) {
  missions <- c("M13", "M15", "M33", "M90", "M180-1", "M180-2")
} else if (flag %in% c("_2_1_rev", "_2_1")) {
  missions <- c("M90", "M180-1", "M180-2")
} else {
  message("NOT PRE-DEFINED FLAG! STOP RUNNING!")
  quit()
}

DMP_result_path = paste0("./DMP_result", flag)

# Define mission information.
dataset_info <- list(
  list(name = "M13", T1 = "T1", T2 = "T2", T3 = "T3", M113 = "M113", D10pro = "D10_pro", days = 13),
  list(name = "M15", T1 = "T1", T2 = "T2", T3 = "T3", M113 = "M113", D10pro = "D10_pro",days = 15),
  list(name = "M33", T1 = "T1", T2 = "T2", T3 = "T3", M113 = "M113", D10pro = "D10_pro", days = 33),
  list(name = "M90", T1 = "T1", T2 = "T2", T3 = "T3", M113 = "M113", D10pro = "D10_pro", days = 90),
  list(name = "M180-1", T1 = "T1", T2 = "T2", T3 = "T3", M113 = "M113", D10pro = "D10_pro", days = 180),
  list(name = "M180-2", T1 = "T1", T2 = "T2", T3 = "T3", M113 = "M113", D10pro = "D10_pro", days = 180)
)

dataset_info <- Filter(function(x) x$name %in% missions, dataset_info)


## 01Start: Statistics on the number of hypermethylated and hypomethylated DMPs ##
message("\n[01] Starting DMP count statistics...")

result_df <- data.frame()
for (exp in dataset_info) {
  message(paste("\nProcessing dataset:", exp$name))

  # Read DMP results.
  df_DMP_list <- readRDS(file.path(DMP_result_path, exp$name, paste0("DMP_res_", exp$name, ".Rds")))
  df_DMP_T1_to_T2 <- df_DMP_list[[paste0(exp$T1, " to ", exp$T2)]]
  df_DMP_T2_to_T3 <- df_DMP_list[[paste0(exp$T2, " to ", exp$T3)]]
  df_DMP_T1_to_T3 <- df_DMP_list[[paste0(exp$T1, " to ", exp$T3)]]
  df_DMP_T1_to_M113 <- df_DMP_list[[paste0(exp$T1, " to ", exp$M113)]]
  df_DMP_T1_to_D10pro <- df_DMP_list[[paste0(exp$T1, " to ", exp$D10pro)]]

  # Classify DMPs by methylation changing direction
  DMP_list_deltaBeta0.05 <- list(up_T1_to_T2 = rownames(df_DMP_T1_to_T2[df_DMP_T1_to_T2$change_deltaBeta0.05 == "up",]),
                                down_T1_to_T2 = rownames(df_DMP_T1_to_T2[df_DMP_T1_to_T2$change_deltaBeta0.05 == "down",]),
					                      stable_T1_to_T2 = rownames(df_DMP_T1_to_T2[df_DMP_T1_to_T2$change_deltaBeta0.05 == "stable",]),

                                up_T2_to_T3 = rownames(df_DMP_T2_to_T3[df_DMP_T2_to_T3$change_deltaBeta0.05 == "up",]),
                                down_T2_to_T3 = rownames(df_DMP_T2_to_T3[df_DMP_T2_to_T3$change_deltaBeta0.05 == "down",]),
					                      stable_T2_to_T3 = rownames(df_DMP_T2_to_T3[df_DMP_T2_to_T3$change_deltaBeta0.05 == "stable",]),

                                up_T1_to_T3 = rownames(df_DMP_T1_to_T3[df_DMP_T1_to_T3$change_deltaBeta0.05 == "up",]),
                                down_T1_to_T3 = rownames(df_DMP_T1_to_T3[df_DMP_T1_to_T3$change_deltaBeta0.05 == "down",]),
					                      stable_T1_to_T3 = rownames(df_DMP_T1_to_T3[df_DMP_T1_to_T3$change_deltaBeta0.05 == "stable",]),
                                
                                up_T1_to_M113 = rownames(df_DMP_T1_to_M113[df_DMP_T1_to_M113$change_deltaBeta0.05 == "up",]),
                                down_T1_to_M113 = rownames(df_DMP_T1_to_M113[df_DMP_T1_to_M113$change_deltaBeta0.05 == "down",]),
					                      stable_T1_to_M113 = rownames(df_DMP_T1_to_M113[df_DMP_T1_to_M113$change_deltaBeta0.05 == "stable",]),

                                up_T1_to_D10pro = rownames(df_DMP_T1_to_D10pro[df_DMP_T1_to_D10pro$change_deltaBeta0.05 == "up",]),
                                down_T1_to_D10pro = rownames(df_DMP_T1_to_D10pro[df_DMP_T1_to_D10pro$change_deltaBeta0.05 == "down",]),
					                      stable_T1_to_D10pro = rownames(df_DMP_T1_to_D10pro[df_DMP_T1_to_D10pro$change_deltaBeta0.05 == "stable",])
                                )

  # Save classification results.
  saveRDS(DMP_list_deltaBeta0.05, file.path(DMP_result_path, exp$name, "DMP_list_deltaBeta0.05.Rds"))
  lapply(DMP_list_deltaBeta0.05, length)
  message(paste("Saved DMP classification for", exp$name))

  # Aggregate DMP counts for each mission.
  counts <- sapply(DMP_list_deltaBeta0.05, length)
  row_df <- data.frame(Dataset = exp$name, Days = exp$days, t(counts), check.names = FALSE)
  result_df <- rbind(result_df, row_df)
}

# Export summary statistics.
write.csv(result_df, file = file.path(DMP_result_path, "DMP_counts_summary.csv"), row.names = FALSE, quote = FALSE)
message("\n[01] DMP count summary saved to DMP_counts_summary.csv")
## 01End: Statistics on the number of hypermethylated and hypomethylated DMPs ##


## 02Start: Count the number of DMPs and categorise them according to the change process in the T2_vs_T1 and T3_vs_T2 ##
message("\n[02] Analyzing DMP temporal trends...")

for (exp in dataset_info) {
  message(paste("\nProcessing trend analysis for:", exp$name))

  # Define possible state transitions.
  state <- c("up", "stable", "down")
  state_groups <- unique(t(combn(c(1: 3, 1: 3), 2)))
  trends <- t(apply(state_groups, 1, function(x) {return(c(state[x]))}))
  trends <- as.data.frame(trends)
  colnames(trends) <- c("T1_to_T2","T2_to_T3")

  # Load DMP classification data.
  DMP_list_deltaBeta0.05 <- readRDS(file.path(DMP_result_path, exp$name, "DMP_list_deltaBeta0.05.Rds"))

  # Identify trend-specific probes.
  DMP_list_trends <- apply(trends, 1, function(x) {
    group1 <- paste0(x[1], "_", colnames(trends)[1])
    group2 <- paste0(x[2], "_", colnames(trends)[2])
    trend_probes <- intersect(DMP_list_deltaBeta0.05[[group1]], DMP_list_deltaBeta0.05[[group2]])
    if (length(trend_probes)) {return (trend_probes)} else {list()}
  })
  names(DMP_list_trends) <- paste0(trends$T1_to_T2, "_", trends$T2_to_T3)

  # Remove stable-stable trends.
  DMP_list_trends <- DMP_list_trends[-grep("stable_stable", names(DMP_list_trends))]
  saveRDS(DMP_list_trends, file = file.path(DMP_result_path, exp$name, "DMP_list_trend.Rds"))
  message(paste("Saved trend analysis for", exp$name))

  # Generate trend statistics.
  DMP_list_trend <- readRDS(file.path(DMP_result_path, exp$name, "DMP_list_trend.Rds"))
  DMP_num <- mapply(length, DMP_list_trend)
  DMP_num_df <- data.frame(trend = names(DMP_num), number = DMP_num)
  DMP_num_df$T1_to_T2 <- stringr::str_split_i(DMP_num_df$trend, "_", 1)
  DMP_num_df$T2_to_T3 <- stringr::str_split_i(DMP_num_df$trend, "_", 2)
  write.csv(DMP_num_df, file = file.path(DMP_result_path, exp$name, "DMP_number_statistic.csv"), row.names = F, quote = FALSE)

  # Prepare data for visualization.
  DMP_num_df <- DMP_num_df[DMP_num_df$T1_to_T2 != "stable", ]
  DMP_num_df$T1_to_T2 <- factor(DMP_num_df$T1_to_T2, levels = c("up", "stable", "down"))
  DMP_num_df$T2_to_T3 <- factor(DMP_num_df$T2_to_T3, levels = c("up", "stable", "down"))

  # Creating an alluvium map.
  stratum_fill_color = colors <- c(
    if(any(DMP_num_df$T1_to_T2 == "up"    & DMP_num_df$number != 0)) "#73B9CE",
    if(any(DMP_num_df$T1_to_T2 == "down"  & DMP_num_df$number != 0)) "#C1548E",
    if(any(DMP_num_df$T2_to_T3 == "up"    & DMP_num_df$number != 0)) "#73B9CE",
    if(any(DMP_num_df$T2_to_T3 == "stable"& DMP_num_df$number != 0)) "darkgray",
    if(any(DMP_num_df$T2_to_T3 == "down"  & DMP_num_df$number != 0)) "#C1548E"
  )	
  p_alluvium <- ggplot(data = DMP_num_df, aes(axis1 = T1_to_T2, axis2 = T2_to_T3, y = number)) + 
    geom_alluvium(aes(fill = T1_to_T2),
                  width = 0.1,
                  aes.bind = "flows") +
    geom_stratum(fill = stratum_fill_color,
                width=0.1) +
    scale_fill_manual(breaks = c("up", "stable", "down"),
                      values = c("up" = "#C1548E", "stable" = "darkgray", "down" = "#73B9CE"), 
                      labels = c("up", "stable", "down")) +
    scale_color_manual(breaks = c("up", "stable", "down"), 
                      values = c("#C1548E", "darkgray", "#73B9CE"), 
                      labels = c("up", "stable", "down")) +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) + 
    scale_x_continuous(breaks = c(1, 2), 
                      labels = c("T1_to_T2", "T2_to_T3"),
                      expand = expansion(mult = c(0, 0))) + 
    coord_cartesian(clip="off") + 
    theme_classic()

  ggsave(filename = "Alluvium_Plot.pdf", plot = p_alluvium, path = file.path(DMP_result_path, exp$name), width = 8, height = 8, device = "pdf")
  message(paste("Alluvium plot saved for", exp$name))
}
message("\n[02] Trend analysis completed for all datasets")


## 03Start: Trend plot for differentially methylated probes ##
message("\n[03] Generating methylation trend plots...")

# Trend_matrix: Col: gene; Rownames: time by order; Value: mean of methlation beta value of the genes.
trend_plot <- function(trend_matrix, my_col, title) {
  time <- colnames(trend_matrix)
  trend_matrix <- t(apply(trend_matrix, 1, function(x) scale(x, center = F)))
  t <- 1: length(time)
  data <- as.data.frame(t(trend_matrix))
  data$t <- t
  data_long <- reshape::melt(data, id.vars = "t")
  x <- 1: length(time)
  quantile_value <- apply(trend_matrix, 2, function(x) {
    temp <- round(quantile(x, seq(0, 1, 0.01)), 3)
    return (temp)
  })
  pdata.list <- list()
  for (k in 1: 50) {
    pdata <- data.frame(x, lower = quantile_value[k, ], upper = quantile_value[101 - k, ])
    pdata.list[[k]] <- pdata
  }
  myPalette <- colorRampPalette(my_col)(45)
  plot.trend <- ggplot() + 
    geom_rect(aes(xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf), fill = "white", alpha = 0.4) + 
    geom_rect(aes(xmin = 1, xmax = length(time), ymin = -Inf, ymax = Inf), fill = "#DEF6F3", alpha = 0.4) + 
    geom_rect(aes(xmin = length(time), xmax = Inf, ymin = -Inf, ymax = Inf), fill = "white", alpha = 0.4) + 
    theme(panel.grid = element_blank(), panel.border = element_rect(fill = NA, color = "black", size = 1))
  for(k in 6: 50) {
    pdata <- pdata.list[[k]]
    plot.trend <- plot.trend + geom_ribbon(data = pdata, aes(ymin = lower, ymax = upper, x = x), fill = myPalette[k - 5], alpha = 1)
  }
  plot.trend <- plot.trend + 
    theme_bw() + 
    xlab("Time") + 
    ylab("deltaBeta(vs T1)") + 
    labs(title = title) + 
    theme(panel.grid = element_blank()) + 
    scale_x_continuous(breaks = seq(1, length(time), 1), labels = c("T1", "T2", "T3"))
  return (plot.trend)
}

for (exp in dataset_info) {
  message(paste("\nGenerating trend plots for:", exp$name))

  # Load DMP data.
  df_DMP_list <- readRDS(file.path(DMP_result_path, exp$name, paste0("DMP_res_", exp$name, ".Rds")))
  df_DMP_T1_to_T2 <- df_DMP_list[[paste0(exp$T1, " to ", exp$T2)]]
  df_DMP_T2_to_T3 <- df_DMP_list[[paste0(exp$T2, " to ", exp$T3)]]
  df_DMP_T1_to_T3 <- df_DMP_list[[paste0(exp$T1, " to ", exp$T3)]]

  # Prepare deltaBeta data.
  probes <- rownames(df_DMP_T1_to_T2)
  deltaBeta_df <- data.frame(probe_name = probes, T1_value = 0, T2_value = df_DMP_T1_to_T2$deltaBeta,
                            T3_value = df_DMP_T1_to_T3[probes, "deltaBeta"], change_type = df_DMP_T1_to_T2$change_deltaBeta0.05)
  write.csv(deltaBeta_df, file = file.path(DMP_result_path, exp$name, "Changing_trend_of_DMPs.csv"), row.names = F)

  # Generate trend plots.
  deltaBeta_df_up <- deltaBeta_df[deltaBeta_df$change_type == "up", 2: 4]
  my_col <- c("#F0D8E5", "#C1548E")
  up_plot <- trend_plot(deltaBeta_df_up, my_col, title = "Hyper-methylated probes in T2 vs T1")
  deltaBeta_df_down <- deltaBeta_df[deltaBeta_df$change_type == "down", 2: 4]
  my_col <- c("#C6E4F5", "#6BACD1")
  down_plot <- trend_plot(deltaBeta_df_down, my_col, title = "Hypo-methylated probes in T2 vs T1")

  # Combine and save plots.
  combined_plot <- ggpubr::ggarrange(up_plot, down_plot, ncol = 1, nrow = 2)
  ggsave(filename = "DMP_trend_plots.pdf", plot = combined_plot, 
          path = file.path(DMP_result_path, exp$name), width = 10, height = 8, device = "pdf")
  message(paste("Trend plots saved for", exp$name))
}

message("\n[03] Trend visualization completed")

## 03End: Trend plot for differentially methylated probes ##


## 04Start: Methylation changing trend of genes ##

# Count the overall methylation changes in the promoter region of a gene, only if all probes are consistently hypermethylated or hypomethylated will the gene be labelled up or down accordingly.
message("\n[04] Analyzing gene-level methylation trends...")

for (exp in dataset_info) {
  message(paste("\nProcessing gene trends for:", exp$name))

  # Generate all possible state transition combinations.
  state <- c("up", "stable", "down")
  state_groups <- unique(t(combn(c(1: 3, 1: 3), 2)))
  trends <- t(apply(state_groups, 1, function(x) {return(c(state[x]))}))
  trends <- as.data.frame(trends)
  colnames(trends) <- c("T1_to_T2", "T2_to_T3")
  trends <- trends[!apply(trends, 1, function(x) {x[1] == x[2] & x[2] == "stable"}), ]

  # Load DMP classification data.
  message("Loading DMP classification results...")
  DMP_list_deltaBeta0.05 <- readRDS(file.path(DMP_result_path, exp$name, "DMP_list_deltaBeta0.05.Rds"))
  df_DMP_list <- readRDS(file.path(DMP_result_path, exp$name, paste0("DMP_res_", exp$name, ".Rds")))

  # Identify probes with consistent trends.
  DMP_trends_list <- apply(trends, 1, function(x) {
    group1 <- paste0(x[1], "_", colnames(trends)[1])
    group2 <- paste0(x[2], "_", colnames(trends)[2])
    trend_genes <- intersect(DMP_list_deltaBeta0.05[[group1]], DMP_list_deltaBeta0.05[[group2]])
    if (length(trend_genes)) {return (trend_genes)} else {list()}
  })
  names(DMP_trends_list) <- paste0(trends$T1_to_T2, "_", trends$T2_to_T3)

  # Map probes to genes in promoter regions.
  message("Mapping probes to promoter-associated genes...")
  df_DMP_T1_to_T2 <- df_DMP_list[[paste0(exp$T1, " to ", exp$T2)]]
  DMG_list_deltaBeta0.05 <- DMP_list_deltaBeta0.05
  for (i in 1: length(DMP_list_deltaBeta0.05)) {
    probe_CpG_info <- df_DMP_T1_to_T2[DMP_list_deltaBeta0.05[[i]], ]
    promoter_region <- c("TSS1500", "TSS200", "5'UTR", "1stExon")
    DM_genes <- as.character(probe_CpG_info[probe_CpG_info$feature %in% promoter_region, "gene"])
    DM_genes <- DM_genes[DM_genes != ""]
    DM_genes <- unique(na.omit(DM_genes))
    DMG_list_deltaBeta0.05[[i]] <- DM_genes
  }

  # Resolve conflicting gene annotations.
  message("Resolving conflicting gene annotations...")
  DMG_list_deltaBeta0.05$stable_T1_to_T2 <- DMG_list_deltaBeta0.05$stable_T1_to_T2[!DMG_list_deltaBeta0.05$stable_T1_to_T2 %in% c(DMG_list_deltaBeta0.05$up_T1_to_T2, DMG_list_deltaBeta0.05$down_T1_to_T2)]
  common_genes <- intersect(DMG_list_deltaBeta0.05$up_T1_to_T2, DMG_list_deltaBeta0.05$down_T1_to_T2)
  DMG_list_deltaBeta0.05$up_T1_to_T2 <- DMG_list_deltaBeta0.05$up_T1_to_T2[!DMG_list_deltaBeta0.05$up_T1_to_T2 %in% common_genes]
  DMG_list_deltaBeta0.05$down_T1_to_T2 <- DMG_list_deltaBeta0.05$down_T1_to_T2[!DMG_list_deltaBeta0.05$down_T1_to_T2 %in% common_genes]
  DMG_list_deltaBeta0.05$stable_T2_to_T3 <- DMG_list_deltaBeta0.05$stable_T2_to_T3[!DMG_list_deltaBeta0.05$stable_T2_to_T3 %in% c(DMG_list_deltaBeta0.05$up_T2_to_T3, DMG_list_deltaBeta0.05$down_T2_to_T3)]
  common_genes <- intersect(DMG_list_deltaBeta0.05$up_T2_to_T3, DMG_list_deltaBeta0.05$down_T2_to_T3)
  DMG_list_deltaBeta0.05$up_T2_to_T3 <- DMG_list_deltaBeta0.05$up_T2_to_T3[!DMG_list_deltaBeta0.05$up_T2_to_T3 %in% common_genes]
  DMG_list_deltaBeta0.05$down_T2_to_T3 <- DMG_list_deltaBeta0.05$down_T2_to_T3[!DMG_list_deltaBeta0.05$down_T2_to_T3 %in% common_genes]

  mapply(length, DMG_list_deltaBeta0.05)

  # Identify trend-specific gene groups.
  message("Identifying trend-specific gene groups...")
  patterns <- names(DMP_trends_list)
  pattern_gene_list <- lapply(patterns, function(x) {
    x <- unlist(strsplit(x, "_"))
    group1 <- paste0(x[1], "_", "T1_to_T2")
    group2 <- paste0(x[2], "_", "T2_to_T3")
    trend_genes <- intersect(DMG_list_deltaBeta0.05[[group1]], DMG_list_deltaBeta0.05[[group2]])
    if (length(trend_genes)) {return (trend_genes)} else {list()}
  })
  names(pattern_gene_list) <- patterns
  mapply(length, pattern_gene_list)

  # Visualize gene overlaps.
  message("Generating Venn diagram...")
  p_venn <- ggvenn::ggvenn(pattern_gene_list, c("up_stable", "up_down", "down_up", "down_stable"))
  ggsave(filename = "venn.pdf", plot = p_venn, 
          path = file.path(DMP_result_path, exp$name), width = 12, height = 12, device = "pdf")
  
  message("Saving pattern gene list...")
  saveRDS(pattern_gene_list, file=file.path(DMP_result_path, exp$name, "pattern_gene_list.Rds"))
  message(paste("Gene trend analysis saved for", exp$name))
}

message("\n[04] Gene-level analysis completed")

## 04End: Methylation changing trend of genes ##