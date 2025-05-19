### description:
# Compare DMR and DMP results for a given mission and group pair to check consistance between these results and gurantee correctness:
#  - Load DMP and DMR data
#  - Match CpG probes within each DMR region
#  - Summarize per‑DMR deltaBeta: total, hyper‑ and hypo‑methylated counts
#  - Print a summary line for each DMR ###

# Working directory need to be adjusted according to your path.
setwd("/home/lqwang/Program")

# Command line args.
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-m", "--mission"), type="character", help="mission"),
  make_option(c("-t", "--group1"), type="character", help="group1"),
  make_option(c("-T", "--group2"), type="character", help="group2"),
  make_option(c("-f", "--flag"), type="character", help="flag"),
  make_option(c("-g", "--maxGap450k"), type="numeric", help="maxGap450k"),
  make_option(c("-G", "--maxGapEPIC"), type="numeric", help="maxGapEPIC")
)
parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)
mission <- opts$mission
group1 <- opts$group1
group2 <- opts$group2
flag <- opts$flag
maxGap450k <- opts$maxGap450k
maxGapEPIC <- opts$maxGapEPIC

message("step2_2_check_DMR_DMP_result_match.R", mission, group1, group2, flag, maxGap450k, maxGapEPIC)

library(dplyr)
library(tibble)
library(purrr)

if (mission %in% c("M13", "M15", "M33")) {
    gap = maxGap450k
} else if (mission %in% c("M90", "M180-1", "M180-2")) {
    gap = maxGapEPIC
} else {
    stop("NOT VALID GAP!")
}

DMP_result_path = paste0("./DMP_result", flag)
DMR_result_path = paste0("./DMR_result", flag)

comp = c(group1, group2)

paste0(comp[[1]], "_to_", comp[[2]])
paste0(comp[[1]], " to ", comp[[2]])

# Read and subset DMP table for the comparison
DMP_all = readRDS(file.path(DMP_result_path, mission, paste0("DMP_res_", mission, ".Rds")))
DMP = DMP_all[[paste0(comp[[1]], " to ", comp[[2]])]]
DMR = readRDS(file.path(DMR_result_path, mission, paste0("DMR_", comp[[1]], "_to_", comp[[2]], gap, ".Rds")))

DMP = DMP[, c("deltaBeta", "change_deltaBeta0.05", "CHR", "MAPINFO")]
DMP$CHR = paste0("chr", DMP$CHR)
DMR = DMR[, c("value", "seqnames", "start", "end")]
colnames(DMR) = c("value", "CHR", "start", "end")

DMR_list = list()
for (DMR_idx in rownames(DMR)) {
   DMR_list[[DMR_idx]] = c()
}

# Annotate DMR by mapping overlapping DMP probes
DMR_list <- DMR %>%
  tibble::rownames_to_column("DMR_idx") %>%
  rowwise() %>%
  mutate(probes = list(
    DMP %>%
      filter( CHR == CHR, MAPINFO >= start, MAPINFO <= end ) %>%
      select(deltaBeta, change_deltaBeta0.05)
    )
  ) %>%
  ungroup() %>%
  select(DMR_idx, value, probes) %>%
  mutate(
    deltaBeta_sum = map_dbl(probes, ~ sum(.x$deltaBeta, na.rm = TRUE)),
    deltaBeta_up = map_int(probes, ~ sum(.x$deltaBeta >  0, na.rm = TRUE)),
    deltaBeta_down = map_int(probes, ~ sum(.x$deltaBeta <  0, na.rm = TRUE))
)

for (i in seq_len(nrow(DMR_list))) {
  cat(
    "DMR_idx:",        DMR_list$DMR_idx[i],        ", ",
    "value:",          DMR_list$value[i],          ", ",
    "deltaBeta_sum:",  DMR_list$deltaBeta_sum[i],  ", ",
    "deltaBeta_up:",   DMR_list$deltaBeta_up[i],   ", ",
    "deltaBeta_down:", DMR_list$deltaBeta_down[i], "\n"
  )
}
