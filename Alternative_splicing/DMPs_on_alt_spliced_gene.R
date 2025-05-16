### description:
# Analyze hyper‑ and hypo‑methylated DMPs overlapping alternatively spliced (AS) genes:
#  - Load DMP tables for M90 and M180‑1 across timepoint comparisons
#  - Read AS gene lists and linked probe IDs from external files
#  - For each AS gene, extract DMP annotations (gene, feature, methylation change)
#  - Perform hypergeometric tests for enrichment of “Body” region probes
#  - Summarize results in Excel and generate bar‑ and pie‑charts of DMP counts by region

suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-f", "--flag"), type="character", help="flag"),
  make_option(c("-d", "--direction_flag"), type="character", help="direction_flag")
)
parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)
flag <- opts$flag
direction_flag <- opts$direction_flag

# Working directory need to be adjusted according to your path.
setwd("/home/lqwang/Program")

library(openxlsx)
library(ggplot2)

DMP_result_path = paste0("./DMP_result", flag)
spling_diff_gene_list_dir = "./main/Alternative_splicing/result/"

AS_result_path = paste0("./AS_result", flag)
if (!dir.exists(AS_result_path)) {
    dir.create(AS_result_path, recursive = TRUE, showWarnings = FALSE)
    message("Created directory:", AS_result_path)
}

exp1 = "M90"
exp2 = "M180-1"

if (direction_flag == "up_down") {
  direction = c("up", "down")
} else if (direction_flag == "up") {
  direction = c("up")
} else if (direction_flag == "down") {
  direction = c("down")
} else {
  stop("!!!")
}

#Load DMP information
df_DMP_list_M90 = readRDS(file.path(DMP_result_path, exp1, paste0("DMP_res_", exp1, ".Rds")))
df_DMP_list_M180_1 = readRDS(file.path(DMP_result_path, exp2, paste0("DMP_res_", exp2, ".Rds")))

df_DMP_M90_T1_to_T2 <- df_DMP_list_M90[["T1 to T2"]]
df_DMP_M90_T2_to_T3 <- df_DMP_list_M90[["T2 to T3"]]
df_DMP_M90_T1_to_T3 <- df_DMP_list_M90[["T1 to T3"]]
df_DMP_M180_1_T1_to_T2 <- df_DMP_list_M180_1[["T1 to T2"]]
df_DMP_M180_1_T2_to_T3 <- df_DMP_list_M180_1[["T2 to T3"]]
df_DMP_M180_1_T1_to_T3 <- df_DMP_list_M180_1[["T1 to T3"]]


compositions = c("T1_to_T2_M90", "T2_to_T3_M90", "T1_to_T2_M180_1", "T2_to_T3_M180_1", "T1_to_T3_M180_1")
splice_diff_path_list = c("T1_to_T2_M90_splice_diff.xlsx", "T2_to_T3_M90_splice_diff.xlsx", "T1_to_T2_M180_1_splice_diff.xlsx", 
                          "T2_to_T3_M180_1_splice_diff.xlsx", "T1_to_T3_M180_1_splice_diff.xlsx")
splice_diff_probes_path_list = c("T1_to_T2_M90_splice_diff_probes.txt", "T2_to_T3_M90_splice_diff_probes.txt", 
                                 "T1_to_T2_M180_1_splice_diff_probes.txt", "T2_to_T3_M180_1_splice_diff_probes.txt",
                                 "T1_to_T3_M180_1_splice_diff_probes.txt")
df_DMP_get_list = c("df_DMP_M90_T1_to_T2", "df_DMP_M90_T2_to_T3", "df_DMP_M180_1_T1_to_T2", "df_DMP_M180_1_T2_to_T3", "df_DMP_M180_1_T1_to_T3")

splicing_gene_meth_df_list = list()

splicing_gene_meth_list_list = list()

p_body_list = list()

for (idx in c(1, 2, 3, 4, 5)) {
  splice_diff_path = splice_diff_path_list[idx]
  splice_diff_probes_path = splice_diff_probes_path_list[idx]
  df_DMP = get(df_DMP_get_list[idx])

  spling_diff_gene <- readxl::read_xlsx(file.path(spling_diff_gene_list_dir, splice_diff_path))
  spling_diff_gene <- spling_diff_gene[, 1: 16]
  splicing_gene_meth <- read.table(file.path(spling_diff_gene_list_dir, splice_diff_probes_path))
  spling_diff_gene$meth_probes <- splicing_gene_meth$V3
  splicing_gene_meth_list <- apply(spling_diff_gene, 1, function(x) {
    gene <- x[[3]]
    probes <- unlist(strsplit(x[[17]], ","))[-1]

    probe_anno <- df_DMP[probes, ]

    probe_anno$gene <- as.character(probe_anno$gene)
    probe_gene <- names(which.max(table(probe_anno$gene)))
    if (sum( probe_anno$gene == gene ) > 0) {
      probe_anno <- probe_anno[probe_anno$gene == gene, ]
      print(dim(probe_anno))
    } else {
      print(c(gene, probe_gene))
      probe_anno$gene[probe_anno$gene == probe_gene] <- gene
      probe_anno <- probe_anno[probe_anno$gene == gene, ]
      print(dim(probe_anno))
    }
    print(table(probe_anno$feature[ which(probe_anno$change_deltaBeta0.05 %in% direction) ]))
    probe_anno$test_id <- x[[1]]
    probe_anno$gene_id <- x[[2]]
    return(probe_anno)
  })
  lapply(splicing_gene_meth_list, function(probe_anno){table(probe_anno$feature[ which(probe_anno$change_deltaBeta0.05 %in% direction) ])})
  splicing_gene_meth_df <- do.call(rbind, splicing_gene_meth_list)

  splicing_gene_meth_list_list[[compositions[idx]]] = splicing_gene_meth_list
  splicing_gene_meth_df_list[[compositions[idx]]] = splicing_gene_meth_df
  
  # Hypergeometric test for the significance of methylation changes occurring in probes on alternatively spliced genes.
  N = length(df_DMP$change_deltaBeta0.05[df_DMP$change_deltaBeta0.05 %in% direction])
  M = table(df_DMP$feature[df_DMP$change_deltaBeta0.05 %in% direction])["Body"]

  n = sum(splicing_gene_meth_df$change_deltaBeta0.05 %in% direction)
  k = sum(splicing_gene_meth_df$change_deltaBeta0.05 %in% direction & splicing_gene_meth_df$feature == "Body")
  k/n

  q = k
  k = n
  n = N-M
  m = M

  p_body <- phyper( q-1, m, n, k, lower.tail= FALSE )
  print(p_body)
  if (p_body < 1e-3) {
    label <- format(p_body, scientific = TRUE, digits = 2)
  } else {
    label <- formatC(p_body, format = "f", digits = 3)
  }
  p_body_list[[compositions[idx]]] = label
}


# Saving the result tables.
wb <- createWorkbook()

addWorksheet(wb, "AS_gene_Probe_M90_T1_to_T2")
addWorksheet(wb, "AS_gene_Probe_M90_T2_to_T3")
addWorksheet(wb, "AS_gene_Probe_M180_1_T1_to_T2")
addWorksheet(wb, "AS_gene_Probe_M180_1_T2_to_T3")
addWorksheet(wb, "AS_gene_Probe_M180_1_T1_to_T3")

writeData(wb, "AS_gene_Probe_M90_T1_to_T2", splicing_gene_meth_df_list[["T1_to_T2_M90"]], rowNames = T)
writeData(wb, "AS_gene_Probe_M90_T2_to_T3", splicing_gene_meth_df_list[["T2_to_T3_M90"]], rowNames = T)
writeData(wb, "AS_gene_Probe_M180_1_T1_to_T2", splicing_gene_meth_df_list[["T1_to_T2_M180_1"]], rowNames = T)
writeData(wb, "AS_gene_Probe_M180_1_T2_to_T3", splicing_gene_meth_df_list[["T2_to_T3_M180_1"]], rowNames = T)
writeData(wb, "AS_gene_Probe_M180_1_T1_to_T3", splicing_gene_meth_df_list[["T1_to_T3_M180_1"]], rowNames = T)

saveWorkbook(wb, file.path(AS_result_path, "splicing_gene_meth_df.xlsx"), overwrite = TRUE)

#Barplot of AS gene number 
splicing_diff_num <- c(length(splicing_gene_meth_list_list[["T1_to_T2_M90"]]), length(splicing_gene_meth_list_list[["T2_to_T3_M90"]]), 0,
                       length(splicing_gene_meth_list_list[["T1_to_T2_M180_1"]]), length(splicing_gene_meth_list_list[["T2_to_T3_M180_1"]]),
                       length(splicing_gene_meth_list_list[["T1_to_T3_M180_1"]]))

exp_names <- c("M90_T2_vs_T1", "M90_T3_vs_T2", "M90_T3_vs_T1",
               "M180_1_T1_to_T2", "M180_1_T2_to_T3", "M180_1_T1_to_T3")

splicing_diff_num_df <- data.frame(test_name = factor(exp_names, levels = exp_names), splicing_diff_num = splicing_diff_num, 
                                   Experiments = c(rep("M90", 3), rep("M180-1", 3)))

p_bar = ggplot(data = splicing_diff_num_df, mapping = aes(x = test_name, y = splicing_diff_num, fill = Experiments)) + 
       geom_bar(stat = 'identity') + 
	     geom_text(aes(label = splicing_diff_num), vjust = -0.5) + 
       scale_fill_manual(values = c("#A71462", "#8695C2"))

ggsave(filename = "Bar_Plot_AS.pdf", plot = p_bar, path = AS_result_path, width = 10, height = 6, device = "pdf")


#Statistics on the number of sites with methylation changes
tM90_T1_to_T2 <- table(splicing_gene_meth_df_list[["T1_to_T2_M90"]]$feature[splicing_gene_meth_df_list[["T1_to_T2_M90"]]$change_deltaBeta0.05 %in% direction])
tM90_T2_to_T3 <- table(splicing_gene_meth_df_list[["T2_to_T3_M90"]]$feature[splicing_gene_meth_df_list[["T2_to_T3_M90"]]$change_deltaBeta0.05 %in% direction])

tM180_1_T1_to_T2 <- table(splicing_gene_meth_df_list[["T1_to_T2_M180_1"]]$feature[splicing_gene_meth_df_list[["T1_to_T2_M180_1"]]$change_deltaBeta0.05 %in% direction])
tM180_1_T2_to_T3 <- table(splicing_gene_meth_df_list[["T2_to_T3_M180_1"]]$feature[splicing_gene_meth_df_list[["T2_to_T3_M180_1"]]$change_deltaBeta0.05 %in% direction])
tM180_1_T1_to_T3 <- table(splicing_gene_meth_df_list[["T1_to_T3_M180_1"]]$feature[splicing_gene_meth_df_list[["T1_to_T3_M180_1"]]$change_deltaBeta0.05 %in% direction])

splicing_meth_num_mat <- rbind(tM90_T1_to_T2, tM90_T2_to_T3, tM180_1_T1_to_T2, tM180_1_T2_to_T3, tM180_1_T1_to_T3)
splicing_meth_num_df <- reshape::melt(splicing_meth_num_mat)
p_values <- c(p_body_list[["T1_to_T2_M90"]], p_body_list[["T2_to_T3_M90"]], p_body_list[["T1_to_T2_M180_1"]], 
                    p_body_list[["T2_to_T3_M180_1"]], p_body_list[["T1_to_T3_M180_1"]])
# p_values[3: 4] <- c("8.16e-06", "1.27e-07")
ggplot(data = splicing_meth_num_df, mapping = aes(x = X1, y = value, fill = X2)) + geom_bar(stat = 'identity', position = 'stack') +
  annotate(geom = "text", x = 1, y = 42, label = p_values[1]) +
  annotate(geom = "text", x = 2, y = 55, label = p_values[2]) +
  annotate(geom = "text", x = 3, y = 90, label = p_values[3]) +
  annotate(geom = "text", x = 4, y = 90, label = p_values[4]) +
  annotate(geom = "text", x = 4, y = 90, label = p_values[5])

#Drawing pie charts of the distribution of methylation-altered probes in functional regions of genes
splicing_meth_num_df <- splicing_meth_num_df[ order(splicing_meth_num_df$X1, splicing_meth_num_df$X2), ]
splicing_meth_num_df$percentage <- unlist( tapply( splicing_meth_num_df$value, factor(splicing_meth_num_df$X1), function (x) {x / sum(x)} ) ) * 100
splicing_meth_num_df$X1 <- gsub("tM90_", "90-days flight ", splicing_meth_num_df$X1)
splicing_meth_num_df$X1 <- gsub("tM180_1_", "180-days flight ", splicing_meth_num_df$X1)
plot <- vector("list", length = length(unique(splicing_meth_num_df$X1)))
names(plot) <- unique(splicing_meth_num_df$X1)
t <- tapply(splicing_meth_num_df$value, factor(splicing_meth_num_df$X1), function(x) { sum(x) })
names(p_values) <- unique(splicing_meth_num_df$X1)
colors <- c("#9E77AA", "#C2B3A2", "#6FB492", "#8694BF", "#D176A1", "#6EA1C4","#B25084", "#D0D1DC")
for(i in unique(splicing_meth_num_df$X1)){
  data <- splicing_meth_num_df[splicing_meth_num_df$X1 == i, ]
  plot[[i]] <- ggplot(data, aes(x = "", y = percentage, fill = as.factor(X2))) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = colors) +
    coord_polar(theta = "y") +
    labs(fill = "Gene Regions") +
    xlab("percentage") +
    ylab(paste0("N=", t[i], ", P(in Body)=", p_values[i])) +
    ggtitle(i)
}
p_pie = ggpubr::ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], nrow = 2, ncol=2, common.legend=T, legend="right")
ggsave(filename = paste0("Pie_Plot_AS_", paste(direction, collapse = "_"), ".pdf"), plot = p_pie, path = AS_result_path, width = 10, height = 6, device = "pdf")
