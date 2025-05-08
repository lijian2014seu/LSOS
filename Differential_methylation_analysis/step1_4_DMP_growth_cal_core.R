### description: Identify genome regions with positive correlation between DMP density and inflight duration in T2 vs T1 ###

# Working directory need to be adjusted according to your path.
setwd("/home/lqwang/Program")

# Command line args.
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-f", "--flag"), type="character", help="flag"),
  make_option(c("-g", "--group1"), type="character", help="group1"),
  make_option(c("-G", "--group2"), type="character", help="group2"),
  make_option(c("-w", "--window_size"), type="numeric", help="window_size"),
  make_option(c("-d", "--direction"), type="character", help="direction")
)
parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)
flag <- opts$flag
window_size <- opts$window_size
direction <- opts$direction
group1 <- opts$group1
group2 <- opts$group2

message("step2_3_2_DMP_growth_cal_core.R", flag, window_size, direction, group1, group2)

DMP_growth_res_path = "./DMP_growth_result"

density_all_missions_gw = readRDS(file.path(DMP_growth_res_path, 
                                        paste("density_across_missions_genomewide", group1, group2, 
                                              window_size, direction, ".Rds", sep = "_")))

density_gw_M13 = density_all_missions_gw[["M13"]]
density_gw_M15 = density_all_missions_gw[["M15"]]
density_gw_M33 = density_all_missions_gw[["M33"]]
density_gw_M90 = density_all_missions_gw[["M90"]]
density_gw_M180_1 = density_all_missions_gw[["M180-1"]]
density_gw_M180_2 = density_all_missions_gw[["M180-2"]]

density_all = data.frame(density_gw_M13$chr, density_gw_M13$start, density_gw_M13$end, 
                density_gw_M13$value, 
                density_gw_M15$value, 
                density_gw_M33$value,
                density_gw_M90$value, 
                density_gw_M180_1$value, 
                density_gw_M180_2$value)
colnames(density_all) <- c("chr_hg38", "start_hg38", "end_hg38", "M13", "M15", "M33", "M90", "M180-1", "M180-2")
rownames(density_all) <- rownames(density_gw_M13)


clean_dataframe <- function(df) {
  rows_with_zero <- apply(df[, c("M13", "M15", "M33", "M90", "M180-1", "M180-2")], 1, function (row) any(row == 0))
  df <- df[!rows_with_zero, ]
  
  # Positve correlation between DMP density and inflight duration. In other words, we want to identify windows (genome region) like this: DMP density monotonically rises as mission duration increases.
  condition1 <- apply(df[, c("M13", "M15", "M33", "M90", "M180-1")], 1, function(row) all(diff(row) > 0))
  condition2 <- apply(df[, c("M13", "M15", "M33", "M90", "M180-2")], 1, function(row) all(diff(row) > 0))

  df <- df[condition1 | condition2, ]
  return(df)
}

growthcore <- clean_dataframe(density_all)

saveRDS(growthcore, file.path(DMP_growth_res_path, paste("growthcore", group1, group2, window_size, direction, ".Rds", sep = "_")))
message("Done!")