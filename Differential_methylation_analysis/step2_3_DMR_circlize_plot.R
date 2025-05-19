### description:
# Generate circos plots of DMR distributions for each mission and comparison:
#  - Load DMR results (hg38 coordinates) per mission & comparison
#  - Exclude sex chromosomes, split hyper- vs hypo-methylated regions
#  - Plot genome-wide DMR distribution via circlize_plot utility ###

# Working directory need to be adjusted according to your path.
setwd("/home/lqwang/Program")
source("./main/utils/circlize_plot.R")

# Command line args.
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-f", "--flag"), type="character", help="flag"),
  make_option(c("-g", "--maxGap450k"), type="character", help="maxGmaxGap450kap"),
  make_option(c("-G", "--maxGapEPIC"), type="character", help="maxGapEPIC")
)
parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)
flag <- opts$flag
maxGap450k <- opts$maxGap450k
maxGapEPIC <- opts$maxGapEPIC

message("step2_3_DMR_circlize_plot.R", flag, maxGap450k, maxGapEPIC)

library(circlize)


DMR_result_path = paste0("./DMR_result", flag)
comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3"))
cytoband_hg38_path = "./main/Differential_methylation_analysis/genome_conver/cytoBand.txt"

dataset_paths <- list(
  list( name = "M13" ),
  list( name = "M15" ),
  list( name = "M33" ), 
  list( name = "M90" ),
  list( name = "M180-1" ),
  list( name = "M180-2" )
)

# Loop over missions and comparisons.
for (exp in dataset_paths) {
    if (exp$name %in% c("M13", "M15", "M33")) {
      gap = maxGap450k
    } else if (exp$name %in% c("M90", "M180-1", "M180-2")) {
      gap = maxGapEPIC
    } else {
      stop("NOT VALID arraytype!")
    }
    for (comp in comparisons) {
        comp_name <- paste0(comp[1], "_to_", comp[2])

        # Safely load DMR RDS; skip if absent
        DMR <- tryCatch(
              {
                readRDS(file.path(DMR_result_path, exp$name, paste0("DMR_", comp_name, gap, ".Rds")))
              }, error = function(e) 
              {
                message(e$message)
                message(file.path(DMR_result_path, exp$name, paste0("DMR_", comp_name, gap, ".Rds")))
                message("Skip.....")
                return (NULL)
              }
            )
        if (is.null(DMR)) { next }
        DMR <- DMR[, c("chr_hg38", "start_hg38", "end_hg38", "value")]
        names(DMR)[1] <- "chr"
        DMR <- DMR[DMR$chr != "chrX" & DMR$chr != "chrY", ]

        DMR_hyper <- DMR[DMR$value > 0, ]
        message(exp$name, " ", comp_name, " hyper DMRs: ", nrow(DMR_hyper))
        DMR_hypo <- DMR[DMR$value < 0, ]
        message(exp$name, " ", comp_name, " hypo DMRs: ", nrow(DMR_hypo))

        # Generate and save circos plot
        plot_circlize(DMR, 
                      file.path(DMR_result_path, exp$name, paste0("genome_distribution_", comp[1], "_to_", comp[2], ".pdf")), 
                      paste0("DMRs_", exp$name, "_", comp[2], "vs", comp[1]), 
                      cytoband_hg38_path)
        message(exp$name, " ", comp_name, " Done!")
    }
}
