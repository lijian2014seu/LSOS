### description: Calculate DMP density on sliding windows across the genome in T2 vs T1 ###

# Working directory need to be adjusted according to your path.
setwd("/home/lqwang/Program")
source("./main/utils/Calculate_Genomic_Density.R")

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

library(circlize)
library(data.table)

cytoband_hg38_path = "./main/Differential_methylation_analysis/genome_conver/cytoBand.txt"

DMP_result_path = paste0("./DMP_result", flag)
bg_density_450k_path = paste0("./main/Differential_methylation_analysis/genome_conver/", "450k_density_", window_size, ".Rds")
bg_density_epic_path = paste0("./main/Differential_methylation_analysis/genome_conver/", "epic_density_", window_size, ".Rds")

DMP_growth_res_path = "./DMP_growth_result"
if (!dir.exists(DMP_growth_res_path)) {
    dir.create(DMP_growth_res_path, recursive = TRUE, showWarnings = FALSE)
    message("Created directory:", DMP_growth_res_path)
}

comp = c(group1, group2)
comp_name = paste0(comp[[1]], "_", comp[[2]])

# Read background probe info of 450K and EPIC. A dataframe records the number of CpG probes in a fixed length genome bp window, sliding along the genome. 
bg_density_450k <- readRDS(bg_density_450k_path)
bg_density_epic <- readRDS(bg_density_epic_path)

track_height = 0.07

# Define mission information.
dataset_paths <- list(
  list(name = "M13", T1 = "T1", T1 = "T2", T1 = "T3", color = "black"),
  list(name = "M15", T1 = "T1", T1 = "T2", T1 = "T3", color = "black"),
  list(name = "M33", T1 = "T1", T1 = "T2", T1 = "T3", color = "black"), 
  list(name = "M90", T1 = "T1", T2 = "T2", T3 = "T3", color = "black"),
  list(name = "M180-1", T1 = "T1", T2 = "T2", T3 = "T3", color = "black"),
  list(name = "M180-2", T1 = "T1", T2 = "T2", T3 = "T3", color = "black")
)

# Initialize the genome information of 22 chromosomes and their bp range.
circos.initializeWithIdeogram(
  cytoband = cytoband_hg38_path,
  plotType = c("axis", "labels"),
  ideogram.height = 0.05
)

# Calculate DMP counts in each window for all the missions, divided by background CpG probe counts, resulting in DMP density per window.
result_list = list()
for (exp in dataset_paths) {
  message("Processing ", exp$name, "...")
  if ( exp$name %in% c("M13", "M15", "M33") ) {
    bg_density = bg_density_450k
  } else {
    bg_density = bg_density_epic
  }

  DMP1 <- readRDS(file.path(DMP_result_path, exp$name, paste0("DMP_res_", exp$name, ".Rds")))
  DMP2 <- DMP1[[paste0(comp[[1]], " to ", comp[[2]])]]
  DMP2 <- DMP2[DMP2$change_deltaBeta0.05 == direction, ]
  DMP2 <- DMP2[!(DMP2$CHR_hg38 %in% c("chrX", "chrY")), ]
  DMP3 <- data.frame(seqnames = DMP2$CHR_hg38, start = DMP2$MAPINFO_hg38,
                      end = DMP2$MAPINFO_hg38, row.names = rownames(DMP2)
                    )
  DMP4 <- Calculate_Genomic_Density(data = DMP3, window.size = window_size, overlap = TRUE, count_by = "number")[[1]]
  
  if (!all(rownames(DMP4) %in% rownames(bg_density))) {
    print(rownames(DMP4))
    print(rownames(bg_density))
    stop("inconsistance between bg_density and DMP_density!")
  }

  DMP4$value <- DMP4$value / (bg_density[row.names(DMP4), ]$value + 1)
  result_list[[exp$name]] = DMP4
}

saveRDS(result_list, file.path(DMP_growth_res_path, 
                               paste("density_across_missions_genomewide", group1, group2, window_size, direction, ".Rds", sep = "_")))
