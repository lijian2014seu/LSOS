### description:
# Perform comprehensive DMR intersection and functional enrichment across missions:
#  - Load and preprocess DMR (“up”/“down”) sets (hg38) for each mission and comparison
#  - Identify overlaps, containments, and exclusives via GRanges
#  - Map overlapping/exclusive DMRs to DMP-derived genes
#  - Run GO enrichment and collate results into summary workbook and circos figures ###

# Working directory need to be adjusted according to your path.
setwd("/home/lqwang/Program")
source("./main/utils/enrich_function.R")
source("./main/utils/circlize_plot.R")
source("./main/utils/DMP_genome_region_function_annotation.R")

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

message("step2_4_DMR_intersection.R", flag, maxGap450k, maxGapEPIC)

library(dplyr); library(GenomicRanges); library(openxlsx); library(readxl); library(purrr)

cytoband_hg38_path = "./main/Differential_methylation_analysis/genome_conver/cytoBand.txt"
DMP_result_path = paste0("./DMP_result", flag)
DMR_result_path = paste0("./DMR_result", flag)

# Function grange_containment: compute within/contains/equal relationships
grange_containment <- function (query, exp_1, subject, exp_2) {
  empty_df <- data.frame(
    chr_hg38    = character(0),
    start_hg38  = integer(0),
    end_hg38    = integer(0),
    width       = integer(0),
    strand      = character(0),
    value       = numeric(0),
    DMR_index   = character(0),
    group       = character(0),
    stringsAsFactors = FALSE
  )
  rownames(empty_df) <- character(0)

  if (nrow(query) == 0 || nrow(subject) == 0) {
    return(list(DMR1_which_be_contained_DMR2 = empty_df, 
                DMR2_which_be_contained_DMR1 = empty_df,
                DMR1_which_contains_DMR2 = empty_df,
                DMR2_which_contains_DMR1 = empty_df,
                DMR1_equal_with_DMR2 = empty_df))
  }
  makeGR <- function (df) {
    makeGRangesFromDataFrame(df, 
                              keep.extra.columns = TRUE,
                              seqnames.field = "chr_hg38",
                              start.field = "start_hg38", 
                              end.field = "end_hg38")
  }
  formatting <- function (df) {
    df <- df %>%
      group_by(seqnames, start, end) %>%
      summarise(        
        width = mean(width),
        strand = paste(strand, sep = "+"),
        value = mean(value),
        DMR_index = paste(DMR_index, sep = "+"),
        .groups = "drop"
      ) %>% 
      as.data.frame() %>%
      rename(., c("chr_hg38" = "seqnames", 
                  "start_hg38" = "start", 
                  "end_hg38" = "end")) 
    rownames(df) = paste(df$chr_hg38, df$start_hg38, df$end_hg38, sep="_")
    df
  }
  DMR1 <- makeGR(query)
  DMR2 <- makeGR(subject)
  
  hits <- findOverlaps(DMR1, DMR2, type = "within")
  if (length(hits) != 0) {
    DMR1_which_be_contained_DMR2 <- formatting(as.data.frame(DMR1[unique(queryHits(hits))]))
    DMR1_which_be_contained_DMR2$group <- paste(exp_2, "contain", exp_1, sep = "")
    DMR2_which_contains_DMR1 <- formatting(as.data.frame(DMR2[unique(subjectHits(hits))]))
    DMR2_which_contains_DMR1$group <- paste(exp_2, "contain", exp_1, sep = "")
  } else {
    DMR1_which_be_contained_DMR2 <- empty_df
    DMR2_which_contains_DMR1 <- empty_df
  }

  hits_rev <- findOverlaps(DMR2, DMR1, type = "within")
  if (length(hits_rev) != 0) {
    DMR2_which_be_contained_DMR1 <- formatting(as.data.frame(DMR2[unique(queryHits(hits_rev))]))
    DMR2_which_be_contained_DMR1$group <- paste(exp_1, "contain", exp_2, sep = "")
    DMR1_which_contains_DMR2 <- formatting(as.data.frame(DMR1[unique(subjectHits(hits_rev))]))
    DMR1_which_contains_DMR2$group <- paste(exp_1, "contain", exp_2, sep = "")
  } else {
    DMR2_which_be_contained_DMR1 <- empty_df
    DMR1_which_contains_DMR2 <- empty_df
  }

  hits_equal <- findOverlaps(DMR1, DMR2, type = "equal")
  if (length(hits_equal) != 0) {
    DMR1_equal_with_DMR2 <- formatting(as.data.frame(DMR1[unique(queryHits(hits_equal))]))
    DMR1_equal_with_DMR2$group <- paste(exp_1, "contain", exp_2, sep = "")
  } else {
    DMR1_equal_with_DMR2 <- empty_df
  }

  return(list(DMR1_which_be_contained_DMR2 = DMR1_which_be_contained_DMR2, 
              DMR2_which_be_contained_DMR1 = DMR2_which_be_contained_DMR1,
              DMR1_which_contains_DMR2 = DMR1_which_contains_DMR2,
              DMR2_which_contains_DMR1 = DMR2_which_contains_DMR1,
              DMR1_equal_with_DMR2 = DMR1_equal_with_DMR2))
}

# Read DMR data, split into hyper/hypo (exclude chrX/Y)
read_DMR_data <- function (file_path, mission_name) {
  DMR_raw = readRDS(file_path)[, c("chr_hg38", "start_hg38", "end_hg38", "width", "strand", "value")]
  DMR_raw$DMR_index = paste(mission_name, rownames(DMR_raw), sep = "_")
  DMR_raw$width = DMR_raw$end_hg38 - DMR_raw$start_hg38
  rownames(DMR_raw) = paste(DMR_raw$chr_hg38, DMR_raw$start_hg38, DMR_raw$end_hg38, sep="_")
  DMR_up = DMR_raw[DMR_raw$value > 0 & !(DMR_raw$chr_hg38 %in% c("chrX", "chrY")), ]
  DMR_down = DMR_raw[DMR_raw$value < 0 & !(DMR_raw$chr_hg38 %in% c("chrX", "chrY")), ]
  return ( list(DMR_up = DMR_up, DMR_down = DMR_down) )
}

# Retrieve unique genes of overlapping DMPs within DMR regions for experiments
get_genes_for_DMP_in_DMR <- function (DMP_result_path, DMR, exp1, exp2, group1, group2) {
  up_overlap_M90_DMPs <- filter_DMP_in_a_particular_genome_region(
                                DMP_result_path = file.path(DMP_result_path, exp1, paste0("DMP_res_", exp1, ".Rds")),
                                DMR_list = DMR,
                                g1 = group1, g2 = group2, condition = c("up", "down")
                          )
  up_overlap_M180_1_DMPs <- filter_DMP_in_a_particular_genome_region(
                                DMP_result_path = file.path(DMP_result_path, exp2, paste0("DMP_res_", exp2, ".Rds")),
                                DMR_list = DMR,
                                g1 = group1, g2 = group2, condition = c("up", "down")
                          )
  up_overlap_M90_M180_1_genes <- unique(c(up_overlap_M90_DMPs$gene, up_overlap_M180_1_DMPs$gene))
  as.character(up_overlap_M90_M180_1_genes)
}

get_genes_for_DMP_in_DMR_single <- function (DMP_result_path, DMR, exp, group1, group2) {
  up_overlap_M90_DMPs <- filter_DMP_in_a_particular_genome_region(
                                DMP_result_path = file.path(DMP_result_path, exp, paste0("DMP_res_", exp, ".Rds")),
                                DMR_list = DMR,
                                g1 = group1, g2 = group2, condition = c("up", "down")
                          )
  as.character(unique(up_overlap_M90_DMPs$gene))
}

# Define the data structure of DMR_intersection_number sheet
count_pair <- function(name1, all_up1, all_down1, ct_up2, ct_down2, ex_up1, ex_down1, ex_gene_up1, ex_gene_down1, ex_fun_up1, ex_fun_down1,
                       name2, all_up2, all_down2, ct_up1, ct_down1, ex_up2, ex_down2, ex_gene_up2, ex_gene_down2, ex_fun_up2, ex_fun_down2,
                       op_up, op_down, op_gene_up, op_gene_down, op_fun_up, op_fun_down, eq_all, group1, group2) {
  safe_nrow <- function(x) {
    nr <- tryCatch(nrow(x), error = function(e) NULL)
    if (is.null(nr) || length(nr) != 1 || is.na(nr)) {
      return(0L)
    }
    return(as.integer(nr))
  }
  data.frame(
    Experiment_group = paste(name1, "vs", name2, sep=""),
    Time_group = paste(group2, "vs", group1, sep=""),
    Number_of_hypermethylated_DMR_in_first_experiment = safe_nrow(all_up1),
    Number_of_hypomethylated_DMR_in_first_experiment = safe_nrow(all_down1),
    Number_of_hypermethylated_DMR_in_second_experiment = safe_nrow(all_up2),
    Number_of_hypomethylated_DMR_in_second_experiment = safe_nrow(all_down2),
    Number_of_hypermethylated_DMRs_in_the_second_experiment_that_included_DMRs_in_the_first_experiment = safe_nrow(ct_up2),
    Number_of_hypomethylated_DMRs_in_the_second_experiment_that_included_DMRs_in_the_first_experiment = safe_nrow(ct_down2),
    Number_of_hypermethylated_DMRs_in_the_first_experiment_that_included_DMRs_in_the_second_experiment = safe_nrow(ct_up1),
    Number_of_hypomethylated_DMRs_in_the_first_experiment_that_included_DMRs_in_the_second_experiment  = safe_nrow(ct_down1),
    Number_of_DMRs_completely_overlapping = safe_nrow(eq_all),
    Number_of_hypermethylated_DMRs_in_the_first_experiment_exclusively = safe_nrow(ex_up1),
    Number_of_hypomethylated_DMRs_in_the_first_experiment_exclusively = safe_nrow(ex_down1),
    Number_of_DMRs_in_the_first_experiment_exclusively = safe_nrow(ex_up1) + safe_nrow(ex_down1),
    Number_of_hypermethylated_DMRs_in_the_second_experiment_exclusively = safe_nrow(ex_up2),
    Number_of_hypomethylated_DMRs_in_the_second_experiment_exclusively = safe_nrow(ex_down2),
    Number_of_DMRs_in_the_second_experiment_exclusively = safe_nrow(ex_up2) + safe_nrow(ex_down2),
    Number_of_overlapped_hypermethylated_DMRs = safe_nrow(op_up),
    Number_of_overlapped_hypomethylated_DMRs = safe_nrow(op_down),
    Number_of_overlapped_DMRs = safe_nrow(op_up) + safe_nrow(op_down),
    Number_of_hypermethylated_genes_on_overlapping_DMRs = safe_nrow(op_gene_up),
    Number_of_hypomethylated_genes_on_overlapping_DMRs = safe_nrow(op_gene_down),
    up_gene_fun_overlapping = op_fun_up,
    down_gene_fun_overlapping = op_fun_down,
    Number_of_hypermethylated_genes_on_DMRs_exclusively_in_the_first_experiment = safe_nrow(ex_gene_up1),
    Number_of_hypomethylated_genes_on_DMRs_exclusively_in_the_first_experiment = safe_nrow(ex_gene_down1),
    up_gene_fun_exclusively_in_the_first_experiment = ex_fun_up1,
    down_gene_fun_exclusively_in_the_first_experiment = ex_fun_down1,
    Number_of_hypermethylated_genes_on_DMRs_exclusively_in_the_second_experiment = safe_nrow(ex_gene_up2),
    Number_of_hypomethylated_genes_on_DMRs_exclusively_in_the_second_experiment = safe_nrow(ex_gene_down2),
    up_gene_fun_exclusively_in_the_second_experiment = ex_fun_up2,
    down_gene_fun_exclusively_in_the_second_experiment = ex_fun_down2,
    stringsAsFactors = FALSE
  )
}


comparisons = list(c("T1", "T2"), c("T2", "T3"), c("T1", "T3"))
exp_comparisons = list(c("M90", "M180-1"), c("M90", "M180-2"), c("M180-1", "M180-2"))
res_list = list()
wb <- createWorkbook()

# Main loop: for each comparison and mission pair.
for (comp in comparisons) {
  group1 = comp[[1]]; group2 = comp[[2]]

  # Load DMRs for M90, M180-1, M180-2
  DMR_M90_real = read_DMR_data(file.path(DMR_result_path, "M90", paste0("DMR_", group1, "_to_", group2, maxGapEPIC, ".Rds")), "M90")
  DMR_M180_1_real = read_DMR_data(file.path(DMR_result_path, "M180-1", paste0("DMR_", group1, "_to_", group2, maxGapEPIC, ".Rds")), "M180-1")
  DMR_M180_2_real = read_DMR_data(file.path(DMR_result_path, "M180-2", paste0("DMR_", group1, "_to_", group2, maxGapEPIC, ".Rds")), "M180-2")

  # Identify CpGs in DMR regions and write to workbook
  M90_CpGs = filter_DMP_in_a_particular_genome_region(
                              DMP_result_path = file.path(DMP_result_path, "M90","DMP_res_M90.Rds"),
                              DMR_list = rbind(DMR_M90_real[["DMR_up"]], DMR_M90_real[["DMR_down"]]),
                              g1 = group1, g2 = group2, condition = c("up", "down", "stable")
                        )

  sheet_name_M90 = paste("CpG_on_DMR_in_M90", group2, "vs", group1, sep = "_")
  addWorksheet(wb, sheet_name_M90)
  writeData(wb, sheet = sheet_name_M90, x = M90_CpGs, rowNames = TRUE)

  # GO enrichment for up/down gene sets
  M90_hyper_enrich = go_enrich( gene_symbols = as.character(unique(M90_CpGs[M90_CpGs$change_deltaBeta0.05 == "up", ]$gene)), gene_type = "DMR_hypermethylated gene" )
  M90_hypo_enrich = go_enrich( gene_symbols = as.character(unique(M90_CpGs[M90_CpGs$change_deltaBeta0.05 == "down", ]$gene)), gene_type = "DMR_hypomethylated gene" )
  sheet_name_M90 = paste("DMR_M90", group2, "vs", group1, "gene_fun", sep = "_")
  addWorksheet(wb, sheet_name_M90)
  writeData(wb, sheet = sheet_name_M90, x = rbind(M90_hyper_enrich, M90_hypo_enrich), rowNames = FALSE)

  # Repeat for M180-1 and M180-2
  M180_1_CpGs = filter_DMP_in_a_particular_genome_region(
                          DMP_result_path = file.path(DMP_result_path, "M180-1", "DMP_res_M180-1.Rds"),
                          DMR_list = rbind(DMR_M180_1_real[["DMR_up"]], DMR_M180_1_real[["DMR_down"]]),
                          g1 = group1, g2 = group2, condition = c("up", "down", "stable")
                    )
  sheet_name_M180_1_cpg = paste("CpG_on_DMR_in_M180-1", group2, "vs", group1, sep = "_")
  addWorksheet(wb, sheet_name_M180_1_cpg)
  writeData(wb, sheet = sheet_name_M180_1_cpg, x = M180_1_CpGs, rowNames = TRUE)

  M180_1_hyper_enrich = go_enrich( gene_symbols = as.character(unique(M180_1_CpGs[M180_1_CpGs$change_deltaBeta0.05 == "up", ]$gene)), gene_type = "DMR_hypermethylated gene" )
  M180_1_hypo_enrich = go_enrich( gene_symbols = as.character(unique(M180_1_CpGs[M180_1_CpGs$change_deltaBeta0.05 == "down", ]$gene)), gene_type = "DMR_hypomethylated gene" )
  sheet_name_M180_1_enrich = paste("DMR_M180-1", group2, "vs", group1, "gene_fun", sep = "_")
  addWorksheet(wb, sheet_name_M180_1_enrich)
  writeData(wb, sheet = sheet_name_M180_1_enrich, x = rbind(M180_1_hyper_enrich, M180_1_hypo_enrich), rowNames = FALSE)


  M180_2_CpGs = filter_DMP_in_a_particular_genome_region(
                          DMP_result_path = file.path(DMP_result_path, "M180-2", "DMP_res_M180-2.Rds"),
                          DMR_list = rbind(DMR_M180_2_real[["DMR_up"]], DMR_M180_2_real[["DMR_down"]]),
                          g1 = group1, g2 = group2, condition = c("up", "down", "stable")
                    )
  sheet_name_M180_2_cpg = paste("CpG_on_DMR_in_M180-2", group2, "vs", group1, sep = "_")
  addWorksheet(wb, sheet_name_M180_2_cpg)
  writeData(wb, sheet = sheet_name_M180_2_cpg, x = M180_2_CpGs, rowNames = TRUE)

  M180_2_hyper_enrich = go_enrich( gene_symbols = as.character(unique(M180_2_CpGs[M180_2_CpGs$change_deltaBeta0.05 == "up", ]$gene)), gene_type = "DMR_hypermethylated gene" )
  M180_2_hypo_enrich = go_enrich( gene_symbols = as.character(unique(M180_2_CpGs[M180_2_CpGs$change_deltaBeta0.05 == "down", ]$gene)), gene_type = "DMR_hypomethylated gene" )
  sheet_name_M180_2_enrich = paste("DMR_M180-2", group2, "vs", group1, "gene_fun", sep = "_")
  addWorksheet(wb, sheet_name_M180_2_enrich)
  writeData(wb, sheet = sheet_name_M180_2_enrich, x = rbind(M180_2_hyper_enrich, M180_2_hypo_enrich), rowNames = FALSE)

  # Pairwise intersection analyses between missions
  message("Processing comparison ", group2, " vs ", group1)
  for (exp_comp in exp_comparisons) {
    DMR_M90 = get(paste("DMR", gsub("-", "_", exp_comp[[1]]), "real", sep = "_"))
    DMR_M90_up = DMR_M90[["DMR_up"]]; DMR_M90_down = DMR_M90[["DMR_down"]]
    DMR_M180_1 = get(paste("DMR", gsub("-", "_", exp_comp[[2]]), "real", sep = "_"))
    DMR_M180_1_up = DMR_M180_1[["DMR_up"]]; DMR_M180_1_down = DMR_M180_1[["DMR_down"]]
    message("Processing", exp_comp[[1]], exp_comp[[2]], "overlap")

    # Compute containment, exclusives, overlaps for up/down sets
    tmp = grange_containment(DMR_M90_up, exp_comp[[1]], DMR_M180_1_up, exp_comp[[2]])
    up_DMR_of_M90_which_be_contained_M180_1 = tmp[["DMR1_which_be_contained_DMR2"]]; up_DMR_of_M180_1_which_be_contained_M90 = tmp[["DMR2_which_be_contained_DMR1"]]; 
    up_DMR_of_M90_which_contains_M180_1 = tmp[["DMR1_which_contains_DMR2"]]; up_DMR_of_M180_1_which_contains_M90 = tmp[["DMR2_which_contains_DMR1"]]
    up_M90_equal_with_M180_1 = tmp[["DMR1_equal_with_DMR2"]]
    tmp = grange_containment(DMR_M90_down, exp_comp[[1]], DMR_M180_1_down, exp_comp[[2]])
    down_DMR_of_M90_which_be_contained_M180_1 = tmp[["DMR1_which_be_contained_DMR2"]]; down_DMR_of_M180_1_which_be_contained_M90 = tmp[["DMR2_which_be_contained_DMR1"]]; 
    down_DMR_of_M90_which_contains_M180_1 = tmp[["DMR1_which_contains_DMR2"]]; down_DMR_of_M180_1_which_contains_M90 = tmp[["DMR2_which_contains_DMR1"]]
    down_M90_equal_with_M180_1 = tmp[["DMR1_equal_with_DMR2"]]

    
    sheet_name_M90_contain_M180_1 = paste0("ct_DMR_", exp_comp[[1]], ">", exp_comp[[2]], "_", group2, "vs", group1)
    addWorksheet(wb, sheet_name_M90_contain_M180_1)
    writeData(wb, sheet = sheet_name_M90_contain_M180_1, x = rbind(up_DMR_of_M90_which_contains_M180_1, down_DMR_of_M90_which_contains_M180_1), rowNames = FALSE)

    sheet_name_M180_1_contain_M90 = paste0("ct_DMR_", exp_comp[[2]], ">", exp_comp[[1]], "_", group2, "vs", group1)
    addWorksheet(wb, sheet_name_M180_1_contain_M90)
    writeData(wb, sheet = sheet_name_M180_1_contain_M90, x = rbind(up_DMR_of_M180_1_which_contains_M90, down_DMR_of_M180_1_which_contains_M90), rowNames = FALSE)
    

    up_M90_exclusive_M180_1 = DMR_M90_up[!(rownames(DMR_M90_up) %in% c(rownames(up_DMR_of_M90_which_be_contained_M180_1), rownames(up_DMR_of_M90_which_contains_M180_1))), ]
    down_M90_exclusive_M180_1 = DMR_M90_down[!(rownames(DMR_M90_down) %in% c(rownames(down_DMR_of_M90_which_be_contained_M180_1), rownames(down_DMR_of_M90_which_contains_M180_1))), ]
    up_M180_1_exclusive_M90 = DMR_M180_1_up[!(rownames(DMR_M180_1_up) %in% c(rownames(up_DMR_of_M180_1_which_be_contained_M90), rownames(up_DMR_of_M180_1_which_contains_M90))), ]
    down_M180_1_exclusive_M90 = DMR_M180_1_down[!(rownames(DMR_M180_1_down) %in% c(rownames(down_DMR_of_M180_1_which_be_contained_M90), rownames(down_DMR_of_M180_1_which_contains_M90))), ]


    sheet_name_M90_exclusive_M180_1 = paste0("un_DMR_", exp_comp[[1]], "(X", exp_comp[[2]], ")_", group2, "vs", group1)
    addWorksheet(wb, sheet_name_M90_exclusive_M180_1)
    writeData(wb, sheet = sheet_name_M90_exclusive_M180_1, x = rbind(up_M90_exclusive_M180_1, down_M90_exclusive_M180_1), rowNames = FALSE)

    sheet_name_M180_1_exclusive_M90 = paste0("un_DMR_", exp_comp[[2]], "(X", exp_comp[[1]], ")_", group2, "vs", group1)
    addWorksheet(wb, sheet_name_M180_1_exclusive_M90)
    writeData(wb, sheet = sheet_name_M180_1_exclusive_M90, x = rbind(up_M180_1_exclusive_M90, down_M180_1_exclusive_M90), rowNames = FALSE)


    up_overlap_M90_M180_1 = bind_rows(up_DMR_of_M90_which_be_contained_M180_1, up_DMR_of_M180_1_which_be_contained_M90) %>%
      { if (nrow(.) == 0) { data.frame( chr_hg38 = character(0), start_hg38 = numeric(0), end_hg38 = numeric(0), value = numeric(0), p.value = numeric(0) ) } else . } %>%
      group_by(chr_hg38, start_hg38, end_hg38) %>%
      summarise(across(where(is.numeric), mean),
      .groups = "drop")
    down_overlap_M90_M180_1 = bind_rows(down_DMR_of_M90_which_be_contained_M180_1, down_DMR_of_M180_1_which_be_contained_M90) %>%
      { if (nrow(.) == 0) { data.frame( chr_hg38 = character(0), start_hg38 = numeric(0), end_hg38 = numeric(0), value = numeric(0), p.value = numeric(0) ) } else . } %>%
      group_by(chr_hg38, start_hg38, end_hg38) %>%
      summarise(across(where(is.numeric), mean),
      .groups = "drop")

    # Identify genes and functions in DMR
    up_overlap_M90_M180_1_genes = get_genes_for_DMP_in_DMR(DMP_result_path, up_overlap_M90_M180_1, exp_comp[[1]], exp_comp[[2]], group1, group2)
    down_overlap_M90_M180_1_genes = get_genes_for_DMP_in_DMR(DMP_result_path, down_overlap_M90_M180_1, exp_comp[[1]], exp_comp[[2]], group1, group2)
    up_M90_exclusive_M180_1_genes = get_genes_for_DMP_in_DMR_single(DMP_result_path, up_M90_exclusive_M180_1, exp_comp[[1]], group1, group2)
    down_M90_exclusive_M180_1_genes = get_genes_for_DMP_in_DMR_single(DMP_result_path, down_M90_exclusive_M180_1, exp_comp[[1]], group1, group2)
    up_M180_1_exclusive_M90_genes = get_genes_for_DMP_in_DMR_single(DMP_result_path, up_M180_1_exclusive_M90, exp_comp[[2]], group1, group2)
    down_M180_1_exclusive_M90_genes = get_genes_for_DMP_in_DMR_single(DMP_result_path, down_M180_1_exclusive_M90, exp_comp[[2]], group1, group2)

    print(up_overlap_M90_M180_1_genes)
    print(class(up_overlap_M90_M180_1_genes))
    print(nrow(up_overlap_M90_M180_1_genes))
    print(down_overlap_M90_M180_1_genes)
    print(nrow(down_overlap_M90_M180_1_genes))


    up_overlap_M90_M180_1_funs_org = go_enrich( gene_symbols = up_overlap_M90_M180_1_genes, gene_type = "overlap" )
    up_overlap_M90_M180_1_funs = paste(up_overlap_M90_M180_1_funs_org$Description, collapse = "; ")

    down_overlap_M90_M180_1_funs_org = go_enrich( gene_symbols = down_overlap_M90_M180_1_genes, gene_type = "overlap" )
    down_overlap_M90_M180_1_funs = paste(down_overlap_M90_M180_1_funs_org$Description, collapse = "; ")
    
    up_M90_exclusive_M180_1_funs_org = go_enrich( gene_symbols = up_M90_exclusive_M180_1_genes, gene_type = "exclusive" )
    up_M90_exclusive_M180_1_funs = paste(up_M90_exclusive_M180_1_funs_org$Description, collapse = "; ")
    
    down_M90_exclusive_M180_1_funs_org = go_enrich( gene_symbols = down_M90_exclusive_M180_1_genes, gene_type = "exclusive" )
    down_M90_exclusive_M180_1_funs = paste(down_M90_exclusive_M180_1_funs_org$Description, collapse = "; ")
    
    up_M180_1_exclusive_M90_funs_org = go_enrich( gene_symbols = up_M180_1_exclusive_M90_genes, gene_type = "exclusive" )
    up_M180_1_exclusive_M90_funs = paste(up_M180_1_exclusive_M90_funs_org$Description, collapse = "; ")
    
    down_M180_1_exclusive_M90_funs_org = go_enrich( gene_symbols = down_M180_1_exclusive_M90_genes, gene_type = "exclusive" )
    down_M180_1_exclusive_M90_funs = paste(down_M180_1_exclusive_M90_funs_org$Description, collapse = "; ")


    # plot circle overlap/exclusive/... DMR figures
    plot_circlize(DMR = rbind(up_overlap_M90_M180_1, down_overlap_M90_M180_1)[, c("chr_hg38", "start_hg38", "end_hg38", "value")], 
                  file_path = file.path(DMR_result_path, paste(exp_comp[[1]], exp_comp[[2]], "overlap", group1, "to", group2, ".pdf", sep = "_")),
                  title_name = paste0(exp_comp[[1]], " & ", exp_comp[[2]], " Overlap DMR"),
                  cytoband_hg38_path = cytoband_hg38_path)

    plot_circlize(DMR = rbind(up_M90_exclusive_M180_1, down_M90_exclusive_M180_1)[, c("chr_hg38", "start_hg38", "end_hg38", "value")], 
                  file_path = file.path(DMR_result_path, paste(exp_comp[[1]], exp_comp[[2]], "exclusive", group1, "to", group2, ".pdf", sep = "_")),
                  title_name = paste0(exp_comp[[1]], "-exclusive DMR_", group2, "vs", group1),
                  cytoband_hg38_path = cytoband_hg38_path)

    plot_circlize(DMR = rbind(up_M180_1_exclusive_M90, down_M180_1_exclusive_M90)[, c("chr_hg38", "start_hg38", "end_hg38", "value")], 
                  file_path = file.path(DMR_result_path, paste(exp_comp[[2]], exp_comp[[1]], "exclusive", group1, "to", group2, ".pdf", sep = "_")),
                  title_name = paste0(exp_comp[[2]], "-exclusive DMR_", group2, "vs", group1),
                  cytoband_hg38_path = cytoband_hg38_path)

    tbl_M90_M180_1 <- count_pair(name1 = exp_comp[[1]], all_up1 = DMR_M90_up, all_down1 = DMR_M90_down, 
                                  ct_up2 = up_DMR_of_M180_1_which_contains_M90, ct_down2 = down_DMR_of_M180_1_which_contains_M90,
                                  ex_up1 = up_M90_exclusive_M180_1, ex_down1 = down_M90_exclusive_M180_1,
                                  ex_gene_up1 = up_M90_exclusive_M180_1_genes, ex_gene_down1 = down_M90_exclusive_M180_1_genes,
                                  ex_fun_up1 = up_M90_exclusive_M180_1_funs, ex_fun_down1 = down_M90_exclusive_M180_1_funs,
                                  name2 = exp_comp[[2]], all_up2 = DMR_M180_1_up, all_down2 = DMR_M180_1_down, 
                                  ct_up1 = up_DMR_of_M90_which_contains_M180_1, ct_down1 = down_DMR_of_M90_which_contains_M180_1,
                                  ex_up2 = up_M180_1_exclusive_M90, ex_down2 = down_M180_1_exclusive_M90,
                                  ex_gene_up2 = up_M180_1_exclusive_M90_genes, ex_gene_down2 = down_M180_1_exclusive_M90_genes,
                                  ex_fun_up2 = up_M180_1_exclusive_M90_funs, ex_fun_down2 = down_M180_1_exclusive_M90_funs, 
                                  eq_all = rbind(up_M90_equal_with_M180_1, down_M90_equal_with_M180_1),
                                  op_up = up_overlap_M90_M180_1, op_down = down_overlap_M90_M180_1, 
                                  op_gene_up = up_overlap_M90_M180_1_genes, op_gene_down = down_overlap_M90_M180_1_genes,
                                  op_fun_up = up_overlap_M90_M180_1_funs, op_fun_down = down_overlap_M90_M180_1_funs,
                                  group1 = group1, group2 = group2)

    res_list[[length(res_list) + 1]] <- tbl_M90_M180_1
  }
}

combined_res <- do.call(rbind, res_list)
addWorksheet(wb, "DMR_intersection_number")
writeData(wb, sheet = "DMR_intersection_number", x = combined_res, rowNames = FALSE)


DMR_xlsx_path <- file.path(table_summary_path, paste0("DMR_summary", flag, ".xlsx"))
sheets_pre <- excel_sheets(DMR_xlsx_path)
wb_pre <- loadWorkbook(DMR_xlsx_path)


current_sheets <- names(wb)

idx_cpg <- grep("^CpG", current_sheets)
idx_fun <- grep("gene_fun$", current_sheets)
other_idx <- setdiff(seq_along(current_sheets), unique(c(idx_cpg, idx_fun)))
sheets_cur <- c(c("DMR_intersection_number"), current_sheets[idx_cpg], current_sheets[idx_fun], current_sheets[other_idx])

sheets_all <- c(c("DMR_intersection_number"), sheets_pre, current_sheets[idx_cpg], current_sheets[idx_fun], current_sheets[other_idx])

merged_list <- list()
for (sh in sheets_all) {
  if (sh %in% sheets_pre) {
    data <- readWorkbook(wb_pre, sheet = sh)
  } else {
    data <- readWorkbook(wb, sheet = sh)
  }
  merged_list[[sh]] <- data
}

wb_out <- createWorkbook()
for (sh in names(merged_list)) {
  addWorksheet(wb_out, sh)
  writeData(wb_out, sh, merged_list[[sh]])
}

message("Writing DMR result to xlsx, which is really slow...")
saveWorkbook(wb_out, file = file.path(DMR_result_path, "Table S2. Results of quantitative and functional enrichment analysis of DMRs.xlsx"), overwrite = TRUE)
message("DMR summary xslx finished: ", file.path(DMR_result_path, "Table S2. Results of quantitative and functional enrichment analysis of DMRs.xlsx"))
