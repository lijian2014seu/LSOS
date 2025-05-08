### description: Plot global changes in DNA methylation across 6 missions in T2 vs T1 ###

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
  make_option(c("-l", "--up_limit"), type="numeric", help="up_limit"),
  make_option(c("-d", "--direction"), type="character", help="direction")
)
parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)
flag <- opts$flag
window_size <- opts$window_size
up_limit <- opts$up_limit
direction <- opts$direction
group1 <- opts$group1
group2 <- opts$group2


if (flag %in% c("_3_3", "_6", "_1_1_1_1_1_1")) {
  missions <- c("M13", "M15", "M33", "M90", "M180-1", "M180-2")
} else if (flag %in% c("_2_1_rev", "_2_1")) {
  missions <- c("M90", "M180-1", "M180-2")
} else {
  message("NOT PRE-DEFINED FLAG! STOP RUNNING!")
  quit()
}

library(circlize)
library(data.table)

track_height = 0.1

cytoband_hg38_path = "./main/Differential_methylation_analysis/genome_conver/cytoBand.txt"
probe_mapping = readRDS("./main/Differential_methylation_analysis/genome_conver/probe_annotation_hg19tohg38.Rds")

DMP_result_path = paste0("./DMP_result", flag)
DMP_growth_res_path = "./DMP_growth_result"

bg_density_450k_path <- paste0("./main/Differential_methylation_analysis/genome_conver/", "450k_density_", window_size, ".Rds")
bg_density_epic_path <- paste0("./main/Differential_methylation_analysis/genome_conver/", "epic_density_", window_size, ".Rds")

density_all_missions_gw <- readRDS(file.path(DMP_growth_res_path, 
                                        paste("density_across_missions_genomewide", group1, group2, 
                                              window_size, direction, ".Rds", sep = "_")))
growthcore <- readRDS(file.path(DMP_growth_res_path, paste("growthcore", group1, group2, window_size, direction, ".Rds", sep = "_")))

comp = c(group1, group2)


# Define mission information.
dataset_paths <- list(
  list(name = "M180-1", T1 = "T1", T2 = "T2", T3 = "T3", color = "#C1548E"),
  list(name = "M180-2", T1 = "T1", T2 = "T2", T3 = "T3", color = "#B1689B"),
  list(name = "M90", T1 = "T1", T2 = "T2", T3 = "T3", color = "#A27CA8"),
  list(name = "M33", T1 = "T1", T1 = "T2", T1 = "T3", color = "#9291B5"),
  list(name = "M15", T1 = "T1", T1 = "T2", T1 = "T3", color = "#82A5C1"),
  list(name = "M13", T1 = "T1", T1 = "T2", T1 = "T3", color = "#73B9CE")
)
dataset_paths <- Filter(function(x) x$name %in% missions, dataset_paths)


Plot_Genomic_Density <- 
  function (data, ymax, ylim.force = FALSE, window.size = NULL, overlap = TRUE, 
            count_by = "number", col = ifelse(area, "grey", "#562e2e"), 
            lwd = par("lwd"), lty = par("lty"), type = "l", 
            area = TRUE, area.baseline = NULL, baseline = 0, border = NA, 
            ...) {
    circos.genomicTrackPlotRegion(data, ylim = c(0,ymax), numeric.column = 4,
        panel.fun = function(region, value, ...) {
            circos.genomicLines(region, value, col = col, ytop = 0.8, ybottom = 0.2,
                lwd = lwd, lty = lty, type = type, border = border, 
                area = area, baseline = baseline)
        }, ...)
  }


comp_name <- paste0(comp[[1]], "_", comp[[2]])
message("Processing ", comp_name, "...")
pdf(file.path(DMP_growth_res_path, paste0(comp_name, "_", direction, "_DMP_growth_", window_size, "_", up_limit, ".pdf")), width = 8, height = 8)
par(mar = c(1, 1, 2, 1), cex.main = 1.5)
circos.par(gap.after = 0.3, start.degree = 90, track.margin = c(0, 0), clock.wise = TRUE)
circos.initializeWithIdeogram(
  cytoband = cytoband_hg38_path,
  plotType = c("axis", "labels"),
  ideogram.height = 0.05
)
bg_density_450k <- readRDS(bg_density_450k_path)
bg_density_epic <- readRDS(bg_density_epic_path)


circos.genomicTrackPlotRegion(
  growthcore, ylim = c(0.1, 0.2),
  panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value, ytop = 0.5, ybottom = 0, 
                        col = "#C1548E",
                        border = NA)
  },
  track.height = 0.05,
  bg.border = NA,
  bg.col = "grey"
)

for (exp in dataset_paths) {
  message("Processing ", exp$name, "...")
  
  DMP1 <- density_all_missions_gw[[exp$name]]

  Plot_Genomic_Density(data = DMP1, ymax = up_limit, col = exp$color, bg.col = "#F0F0F0", bg.border = NA, 
                                window.size = window_size, count_by = "number", track.height = track_height, 
                                track.margin = c(0, 0))
}

circos.clear()
title(main = paste0("DMP_distribution_", comp_name, "_", direction), line = -1)
dev.off()