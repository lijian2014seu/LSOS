### description: Perform DMR analysis across specified missions:
#  - Load normalized beta data and sample metadata
#  - Identify DMRs via Bumphunter with array‑specific maxGap
#  - Annotate and filter DMRs (≥3 CpGs), convert to hg38 coordinates
#  - Save per‑comparison results and compile summary workbook ###

# Working directory need to be adjusted according to your path.
setwd("/home/lqwang/Program")

# Command line args.
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-f", "--flag"), type="character", help="flag"),
  make_option(c("-g", "--maxGap450k"), type="numeric", help="maxGap450k"),
  make_option(c("-G", "--maxGapEPIC"), type="numeric", help="maxGapEPIC")
)
parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)
flag <- opts$flag
maxGap450k <- opts$maxGap450k
maxGapEPIC <- opts$maxGapEPIC

message("step2_1_DMR_analysis.R", flag, maxGap450k, maxGapEPIC)

if (flag %in% c("_3_3", "_6", "_1_1_1_1_1_1")) {
  missions <- c("M13", "M15", "M33", "M90", "M180-1", "M180-2")
} else if (flag %in% c("_2_1_rev", "_2_1")) {
  missions <- c("M90", "M180-1", "M180-2")
} else {
  message("NOT PRE-DEFINED FLAG! STOP RUNNING!")
  quit()
}

library(ChAMP)
library(stringr)
library(dplyr)
library(tibble)
library(openxlsx)
library(bumphunter)
set.seed(123)

# Define output paths.
DMP_result_path = paste0("./DMP_result", flag)
DMR_result_path = paste0("./DMR_result", flag)
xlsx_file = file.path(DMR_result_path, paste0("DMR_summary", flag, ".xlsx"))

if (!dir.exists(DMR_result_path)) {
    dir.create(DMR_result_path, recursive = TRUE, showWarnings = FALSE)
    message("Created directory:", DMR_result_path)
}

# Define mission information.
dataset_paths <- list(
  list(
    name = "M13",
    data_path = "./13day_spaceflight_M13/DNA_methylation/",
    comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3")),
    data_type = "450K"
  ),
  list(
    name = "M15",
    data_path = "./15day_spaceflight_M15/DNA_methylation/",
    comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3")),
    data_type = "450K"
  ),
  list(
    name = "M33",
    data_path = "./33day_spaceflight_M33/DNA_methylation/",
    comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3")),
    data_type = "450K"
  ), 
  list(
    name = "M90",
    data_path = "./90day_spaceflight_M90/DNA_methylation/",
    comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3")),
    data_type = "EPIC"
  ),
  list(
    name = "M180-1",
    data_path = "./180day_spaceflight_M180_1/DNA_methylation/",
    comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3")),
    data_type = "EPIC"
  ),
  list(
    name = "M180-2",
    data_path = "./180day_spaceflight_M180_2/DNA_methylation/",
    comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3")),
    data_type = "EPIC"
  )
)

dataset_paths <- Filter(function(x) x$name %in% missions, dataset_paths)

probe_annotation_hg19tohg38 <- readRDS("./main/Differential_methylation_analysis/genome_conver/probe_annotation_hg19tohg38.Rds")

wb <- createWorkbook()

# Main processing loop.
for (exp in dataset_paths) {
  message("\n", strrep("-", 40))
  message("Processing experiment: ", exp$name)
  message("Comparison: ", paste(exp$comparison, collapse = " vs "))

  # Load preprocessed beta matrix.
  message("\n[1/6] Loading normalized beta matrix...")
  myCombat <- readRDS(paste0(exp$data_path, "beta/nobatch_beta", flag, ".Rds"))

  # Read sample metadata.
  SampleSheet <- read.csv(paste0(exp$data_path, "rawdata/SampleSheet.csv"))
  table(SampleSheet$Sample_Group)

  # Filter samples for comparison groups.
  message("[2/6] Filtering comparison groups...")
  df_DMP_list <- readRDS(file.path(DMP_result_path, exp$name, paste0("DMP_res_", exp$name, ".Rds")))
  

  for (comp in exp$comparisons) {
    if (exp$data_type == "450K") {
      maxGap = maxGap450k
    } else if (exp$data_type == "EPIC") {
      maxGap = maxGapEPIC
    } else {
      stop("NOT VALID arraytype!")
    }
    message("Processing comparison between ", comp[1], " and ", comp[2])
    df_DMP <- df_DMP_list[[paste0(comp[1], " to ", comp[2])]]
    pheno <- factor(SampleSheet$Sample_Group[SampleSheet$Sample_Group %in% comp], levels=comp)

    # Run DMR detection.
    message("[3/6] Identifying DMRs using Bumphunter...") 
    DMR <- tryCatch(  
      champ.DMR(beta = as.matrix(myCombat[, SampleSheet$Sample_Group %in% comp]),
                                pheno = pheno,
                                method = "Bumphunter",
                                compare.group = NULL,
                                arraytype = exp$data_type,
                                maxGap = maxGap,
                                cores = 16),
      error = function(e) {
        message("DMR detection failed: ", e$message)
        return(NULL)
      }
    )
    if (is.null(DMR)) {
      message("Aborting DMR analysis for this experiment.")
      next
    }

    if ("value" %in% colnames(DMR$BumphunterDMR)) {
      up_DMR <- sum(DMR$BumphunterDMR$value > 0, na.rm = TRUE)
      down_DMR <- sum(DMR$BumphunterDMR$value < 0, na.rm = TRUE)
      message("\n", exp$name, "_", comp[1], "_", comp[2], " Hypermethylated DMRs: ", up_DMR)
      message(exp$name, "_", comp[1], "_", comp[2], " Hypomethylated DMRs: ", down_DMR, "\n")
    }

    # Save DMR results.
    message("[5/6] Saving output files...")
    if (!dir.exists(file.path(DMR_result_path, exp$name))) {
      dir.create(file.path(DMR_result_path, exp$name), recursive = TRUE, showWarnings = FALSE)
      message("Created directory:", file.path(DMR_result_path, exp$name), "\n")
    }

    # Map CpG sites to DMRs.
    message("[6/6] Annotating DMR-associated CpGs...")
    index <- apply(DMR$BumphunterDMR, 1, function(x) which(df_DMP$CHR ==
                                                    gsub("chr", "", x[1]) & df_DMP$MAPINFO >= as.numeric(x[2]) &
                                                    df_DMP$MAPINFO <= as.numeric(x[3])))
    
    # Filter DMRs with ≥3 CpG sites.
    DMR_CpG_list <- lapply(index, function(x) df_DMP[x, ])
    index1 <- which(mapply(nrow, DMR_CpG_list) >= 3)

    # Create final annotation table.
    DMR_CpG_table <- do.call(rbind, DMR_CpG_list[index1])
    if (is.null(DMR_CpG_table)) {next}
    DMR_CpG_table $ DMRindex <- str_split_fixed(rownames(DMR_CpG_table), "[.]", 2)[, 1]
    rownames(DMR_CpG_table) <- str_split_fixed(rownames(DMR_CpG_table), "[.]", 2)[, 2]

    dmr_hg38_map <- DMR_CpG_table %>%
      group_by(DMRindex) %>%
      summarise(
        chr_hg38 = first(CHR_hg38),
        start_hg38 = min(MAPINFO_hg38),
        end_hg38 = max(MAPINFO_hg38),
        n_CpG_hg38 = n(),
        .groups = "drop"
      ) %>%
      filter(n_CpG_hg38 >= 3)

    DMR$BumphunterDMR <- DMR$BumphunterDMR %>%
      rownames_to_column("DMRindex") %>% 
      inner_join(
        dmr_hg38_map %>% select(DMRindex, chr_hg38, start_hg38, end_hg38),
        by = "DMRindex"
      ) %>% 
      column_to_rownames("DMRindex")

    saveRDS(DMR$BumphunterDMR, file.path(DMR_result_path, exp$name, 
                                                  paste0("DMR_", comp[1], "_to_", comp[2], "_", maxGap,".Rds")))

    sheet_name <- paste0("DMR_", exp$name, "_", comp[2], "_vs_", comp[1])
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet = sheet_name, x = DMR$BumphunterDMR, rowNames = TRUE)
  }                                  
  message("Analysis completed for ", exp$name, "\n", strrep("-", 40))
}

message("Writing DMR result to xlsx...")
saveWorkbook(wb, xlsx_file, overwrite = TRUE)
message("DMR summary xslx finished: ", xlsx_file)