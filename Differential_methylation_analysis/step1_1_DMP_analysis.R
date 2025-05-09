### description: Differential methylation analysis for 90-day space flight (M90 / 01) and  180-day space flight (M180-1 / 02) ###

# Working directory need to be adjusted according to your path.
setwd("/home/lqwang/Program")
source("./main/utils/consist_judge_fun.R")

# Command line args.
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-f", "--flag"), type="character", help="flag")
)
parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)
flag <- opts$flag

message("step1_1_DMP_analysis.R", flag)

library(ChAMP)
library(ggvenn)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(tibble)
library(stringr)
library(data.table)


if (flag %in% c("_3_3", "_6", "_1_1_1_1_1_1")) {
  missions <- c("M13", "M15", "M33", "M90", "M180-1", "M180-2")
} else if (flag %in% c("_2_1_rev", "_2_1")) {
  missions <- c("M90", "M180-1", "M180-2")
} else {
  message("NOT PRE-DEFINED FLAG! STOP RUNNING!")
  quit()
}

DMP_result_path = paste0("./DMP_result", flag)

# Define mission information and group comparisons need to calculate DMPs.
dataset_paths <- list(
  list(
    name = "M13",
    data_path = "./13day_spaceflight_M13/DNA_methylation/",
    comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3"), c("T1", "D10_pro")),
    data_type = "450K",
    name1 = "M13"
  ),
  list(
    name = "M15",
    data_path = "./15day_spaceflight_M15/DNA_methylation/",
    comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3"), c("T1", "D10_pro")),
    data_type = "450K",
    name1 = "M15"
  ),
  list(
    name = "M33",
    data_path = "./33day_spaceflight_M33/DNA_methylation/",
    comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3"), c("T1", "D10_pro")),
    data_type = "450K",
    name1 = "M33"
  ),
  list(
    name = "M90",
    data_path = "./90day_spaceflight_M90/DNA_methylation/",
    comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3")),
    data_type = "EPIC",
    name1 = "M90"
  ),
  list(
    name = "M180-1",
    data_path = "./180day_spaceflight_M180_1/DNA_methylation/",
    comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3")),
    data_type = "EPIC",
    name1 = "M180_1"
  ),
  list(
    name = "M180-2",
    data_path = "./180day_spaceflight_M180_2/DNA_methylation/",
    comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3"), c("T1", "M113")),
    data_type = "EPIC",
    name1 = "M180_2"
  )
)

all_results <- list()

# Load chr and pos mapping file from hg19 to hg38.
probe_annotation_hg19tohg38 <- readRDS("./main/Differential_methylation_analysis/genome_conver/probe_annotation_hg19tohg38.Rds")

dataset_paths <- Filter(function(x) x$name %in% missions, dataset_paths)

# Calculate DMPs for each mission.
for(exp in dataset_paths){

  message("\n===== processing mission ", exp$name, "=====\n", "\n")
  
  # Load nobatch data and samplesheet.
  beta <- readRDS(paste0(exp$data_path, "beta/nobatch_beta", flag, ".Rds"))
  pd <- read.csv(paste0(exp$data_path, "rawdata/SampleSheet.csv"))
  rownames(pd) <- pd$Sample_Name

  # Check consistance between beta matrix and samplesheet.
  if(!all(colnames(beta) == paste0(rownames(pd), "_", exp$name1))){
    stop(paste(exp$name, "Sample Information Mismatch Error!"))
  }
  
  # Champ.DMP().
  dmp_res <- tryCatch(
    champ.DMP(
      beta = beta,
      pheno = pd$Sample_Group,
      compare.group = NULL,
      adjPVal = 1,
      arraytype = exp$data_type),
    error = function(e) {
      message("DMP detection failed: ", e$message)
      return(NULL)
    }
  )
  if (is.null(dmp_res)) {next}

  # Filter probes with fluctuation over 0.05
  r <- apply(beta, 1, function(x) {max(x) - min(x)})
  filter_beta <- beta[r > 0.05,]
  filter_beta <- as.data.frame(filter_beta)

  # Calculate DMP for each comparison, and write all the outputs to a list named exp_result with item name of "group1 to group2".
  exp_results <- list()
  for(comp in exp$comparisons){
    comp_name <- paste(comp[1], "to", comp[2], sep = " ")
    df <- dmp_res[[paste0(comp[1], "_to_", comp[2])]]
    
    # Filter probes with a consistant change direction among all the samples.
    consist_judge_fun_comp1_to_comp2 <- consist_judge_fun(comp[1], comp[2], filter_beta, pd)

    # Define significant DMPs with a criteria of p value less than 0.05 and avg deltaBeta more than 0.05.
    up_probes <- intersect(names(consist_judge_fun_comp1_to_comp2)[consist_judge_fun_comp1_to_comp2 == "Smaller"], 
                                  rownames(df)[df$P.Value < 0.05])
    down_probes <- intersect(names(consist_judge_fun_comp1_to_comp2)[consist_judge_fun_comp1_to_comp2 == "Greater"], 
                                  rownames(df)[df$P.Value < 0.05])
    df$change <- "stable"
    df[up_probes, "change"] <- "up"
    df[down_probes, "change"] <- "down"
    table(df$change)

    for (cutoff in c(0.05, 0.1)) {
      df$logFC <- df$deltaBeta
      new_col <- paste0("change_deltaBeta", cutoff)
      df[, new_col] <- df$change
      deltaBeta_vec <- df[df$change != "stable", "deltaBeta"]
      change_vec <- ifelse (deltaBeta_vec > cutoff, "up", ifelse (deltaBeta_vec < (-cutoff), "down", "stable"))
      df[df$change != "stable", new_col] <- change_vec
      table(df[, new_col])
    }
    
    # Convert genome position annotation from hg19 to hg38.
    df <- df %>%
    rownames_to_column("ProbeID") %>%
    left_join(probe_annotation_hg19tohg38 %>%
                rownames_to_column("ProbeID") %>%
                select(ProbeID, CHR_hg38, MAPINFO_hg38),
              by = "ProbeID") %>%
    column_to_rownames("ProbeID")
    message(exp$name, "---", comp_name, " mapping done!", "\n")

    exp_results[[comp_name]] <- df
    message("---", comp_name, "analysis completed ---\n")
    table(df$change)
  }

  # Make result directory for each mission.
  if (!dir.exists(file.path(DMP_result_path, exp$name))) {
    dir.create(file.path(DMP_result_path, exp$name), recursive = TRUE, showWarnings = FALSE)
    message("Created directory:", file.path(DMP_result_path, exp$name), "\n")
  }

  all_results[[exp$name]] <- exp_results
  message("\n Now saving RDS... \n")
  saveRDS(exp_results, file.path(DMP_result_path, exp$name, paste0("DMP_res_", exp$name, ".Rds")))
  message("Finish saving ", exp$name, " RDS to ", file.path(DMP_result_path, exp$name, paste0("DMP_res_", exp$name, ".Rds")), "\n")
}

# Saving DMP analysis result for all probes to xlsx file.
wb <- createWorkbook()
for(exp_name in names(all_results)){
  for(comp_name in names(all_results[[exp_name]])){
    comp <- strsplit(comp_name, " to ", fixed = TRUE)[[1]]
    sheet_name <- paste(comp[2], "vs", comp[1], "DMP", exp_name, sep = "_")
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, all_results[[exp_name]][[comp_name]], rowNames = TRUE)
    message("Finish writing sheet ", sheet_name, "\n")
  }
}

message("Now writing .xlsx file to ", DMP_result_path, "/Results of differential methylation analysis of DNA methylation probes.xlsx, which is quite slow...", "\n")
saveWorkbook(wb, file.path(DMP_result_path, "Results of differential methylation analysis of DNA methylation probes.xlsx"), overwrite = TRUE)
message("Finish writing .xlsx file to ", DMP_result_path, "/Results of differential methylation analysis of DNA methylation probes.xlsx", "\n")
