### description: Batch removal and dimensionality reduction visualization for DNA methylation data. ###

# This script offers options for different batch removal strategy. In the paper, we combine M90 and M180-1 to remove batch across mission and M180-2 as separate validation for the different sampling time of T2, referring to _2_1. For Fig .1a, we filtered common probes for all the six missions and remove batch effect for all to ensure a fair compare of DMPs across different experiments, referring to _6.

# Working directory need to be adjusted according to your path.
setwd("/home/lqwang/Program")

# Command line args.
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-f", "--flag"), type="character", help="flag"),
  make_option(c("-p", "--perplexity"), type="numeric", help="flag"),
  make_option(c("-s", "--seed"), type="numeric", help="flag")
)
parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)
flag <- opts$flag
perplexity <- opts$perplexity
seeds <- opts$seeds


message("step2_DNA_meth_batch_remove.R", flag, perplexity, seeds)

library(Rtsne)
library(sva)
library(ChAMP)
library(ggplot2)
library(ggpubr)


# Set the mission groups of batch removal.
if (flag ==  "_3_3") {
  mission_combinations <- list(c("M13", "M15", "M33"), c("M90", "M180_1", "M180_2"))
} else if (flag == "_6") {
  mission_combinations <- list(c("M13", "M15", "M33", "M90", "M180_1", "M180_2"))
} else if (flag == "_2_1") {
  mission_combinations <- list(c("M90", "M180_1"), c("M180_2"))
} else if (flag == "_2_1_rev") {
  mission_combinations <- list(c("M90", "M180_2"), c("M180_1"))
} else if (flag == "_1_1_1_1_1_1") {
  mission_combinations <- list(c("M13"), c("M15"), c("M33"), c("M90"), c("M180_1"), c("M180_2"))
} else {
  message("NOT PRE-DEFINED FLAG! STOP RUNNING!")
  quit()
}

vis_save_path = file.path("./remove_batch_vis", flag)
if (!dir.exists(vis_save_path)) {
  dir.create(vis_save_path, recursive = TRUE, showWarnings = FALSE)
  message("Created directory:", vis_save_path)
}

# Set data directory.
dataset_paths <- list(
  M13 = list(
    norm_path = "./13day_spaceflight_M13/DNA_methylation/beta/Norm_beta.Rds",
    SampleSheet_dir = "./13day_spaceflight_M13/DNA_methylation/rawdata/SampleSheet.csv"
  ),
  M15 = list(
    norm_path = "./15day_spaceflight_M15/DNA_methylation/beta/Norm_beta.Rds",
    SampleSheet_dir = "./15day_spaceflight_M15/DNA_methylation/rawdata/SampleSheet.csv"
  ),
  M33 = list(
    norm_path = "./33day_spaceflight_M33/DNA_methylation/beta/Norm_beta.Rds",
    SampleSheet_dir = "./33day_spaceflight_M33/DNA_methylation/rawdata/SampleSheet.csv"
  ),
  M90 = list(
    norm_path = "./90day_spaceflight_M90/DNA_methylation/beta/Norm_beta.Rds",
    SampleSheet_dir = "./90day_spaceflight_M90/DNA_methylation/rawdata/SampleSheet.csv"
  ),
  M180_1 = list(
    norm_path = "./180day_spaceflight_M180_1/DNA_methylation/beta/Norm_beta.Rds",
    SampleSheet_dir = "./180day_spaceflight_M180_1/DNA_methylation/rawdata/SampleSheet.csv"
  ),
  M180_2 = list(
    norm_path = "./180day_spaceflight_M180_2/DNA_methylation/beta/Norm_beta.Rds",
    SampleSheet_dir = "./180day_spaceflight_M180_2/DNA_methylation/rawdata/SampleSheet.csv"
  )
)



# For each group/combination, execute...
for (combination in mission_combinations) {
  mission_data <- list()
  sample_sheets <- list()

  # Read DNA methylation beta matrix and sample sheet of each mission together into a key-value list.
  for (mission_code in combination) {
    mission_path <- dataset_paths[[mission_code]]
    mission_data[[mission_code]] <- readRDS(mission_path$norm_path)
    message("read beta matrix of ", mission_code, " from ", mission_path$norm_path, " successfully!")
    sample_sheets[[mission_code]] <- read.csv(mission_path$SampleSheet_dir)
    message("read sample sheet of ", mission_code, " from ", mission_path$SampleSheet_dir, " successfully!")
  }

  common_probes <- Reduce(intersect, lapply(mission_data, rownames))

  # Merge all beta matrixes into one big matrix.
  beta_mat_all <- do.call(cbind, lapply(combination, function(code) {
    mission_beta <- mission_data[[code]][common_probes, ]
    colnames(mission_beta) <- paste(colnames(mission_beta), code, sep = "_")
    return(mission_beta)
  }))

  # Merge all sample sheets into one big sheet.
  SampleSheet_all <- do.call(rbind, lapply(combination, function(code) {
    ss <- sample_sheets[[code]]
    ss$Sample_Name <- paste0(ss$Sample_Name, "_", code)
    ss$Subject <- paste0(ss$Sample_Well, "_", code)
    ss$Experiment <- code
    ss
  }))

  message(
  "Successfully merged:\n",
  "1. Beta matrix: ", nrow(beta_mat_all), " probes x ", ncol(beta_mat_all), " samples\n",
  "   Example column names: ", paste(head(colnames(beta_mat_all)), collapse = ", "), if(ncol(beta_mat_all)>3) "...\n" else "\n",
  "2. SampleSheet: ", nrow(SampleSheet_all), " rows x ", ncol(SampleSheet_all), " columns\n",
  "   Example sample names: ", paste(head(SampleSheet_all$Sample_Name), collapse = ", "), if(nrow(SampleSheet_all)>3) "..." else ""
  )

  # Record the number of subjects and time points.
  n_subjects <- length(unique(SampleSheet_all$Subject))
  n_times <- length(unique(SampleSheet_all$Sample_Group))

  # PART ONE：visualization before batch removal.
  set.seed(seeds)
  tsne_out <- Rtsne(t(beta_mat_all), pca=FALSE, perplexity=perplexity, theta=0.0)
  tsnes <- as.data.frame(tsne_out$Y)
  colnames(tsnes) <- c("tSNE1", "tSNE2")

  colors_subject <- head(c("#C1548E", "#B7A637", "#37BA5F", "#90BAFF",
                            "#F590E8", "#73B9CE", "#FF7F50", "#FF6B6B",
                            "#4B0082", "#FFA500", "#00CED1", "#808080",
                            "#6A5ACD", "#E6B566", "#8B0000", "#2E8B57",
                            "#BC8F8F"), n_subjects)

  colors_time <- head(c("#C1548E", "#37BA5F", "#90BAFF","#B7A637","#F590E8"), n_times)

  tsnes$Subject <- SampleSheet_all$Subject
  tsnes$Time <- SampleSheet_all$Sample_Group
  tsnes$Experiment <- SampleSheet_all$Experiment
  tsnes$Sample_Name <- SampleSheet_all$Sample_Name

  write.csv(tsnes, file = file.path(vis_save_path, paste0("tsne_pre_combat_", paste(combination, collapse="_"), ".csv")), row.names = FALSE)

  p1 <- ggplot(tsnes, aes(tSNE1, tSNE2)) + 
        geom_point(aes(color=Subject), size=3) +
        scale_color_manual(values=colors_subject) +
        theme(panel.background = element_rect(fill = "white", color = "black"))+  
        theme(panel.grid.major = element_line(color = "gray"), panel.grid.minor = element_line(color = "gray"))

  p2 <- ggplot(tsnes, aes(tSNE1, tSNE2)) + 
        geom_point(aes(color=Time, shape=Experiment), size=3) +
        scale_color_manual(values=colors_time) +
        theme(panel.background = element_rect(fill = "white", color = "black"))+  
        theme(panel.grid.major = element_line(color = "gray"), panel.grid.minor = element_line(color = "gray"))

  pre_combat_plot <- ggpubr::ggarrange(p1, p2, nrow=1)
  ggsave(file.path(vis_save_path, paste0("tsne_pre_combat_", paste(combination, collapse="_"), ".pdf")), pre_combat_plot, width=14, height=6)
  message("tsne figure before batch removal saved to ", file.path(vis_save_path, paste0("tsne_pre_combat_", paste(combination, collapse="_"), ".pdf")))

  # PART2: batch removal.

  # Inter-mission: remove individual difference among subjects in one mission.
  combat_data <- list()
  for (code in combination) {
    pd <- sample_sheets[[code]]
    message("Start removing individual difference for mission ", code)
    combat_data[[code]] <- champ.runCombat(
                            beta=as.data.frame(mission_data[[code]]),
                            pd=pd,
                            batchname=c("Sample_Well")
                          )
    message("Successfully removed individual difference for mission ", code)
  }
  combat_combined <- do.call(cbind, lapply(combination, function(code) combat_data[[code]][common_probes, ]))
  colnames(combat_combined) <- SampleSheet_all$Sample_Name

  # Cross-mission: remove experimental batch effect across all the missions.
  message("Start removing batch effect of experiment for mission combination ", paste(combination, collapse=", "))
  if (length(combination) > 1) {  
    combat_combined_logit <- logit2(combat_combined)
    combat_combined_corrected <- ComBat(
                                  dat=combat_combined_logit,
                                  batch=SampleSheet_all$Experiment,
                                  par.prior=TRUE
                                )
                                
    final_beta <- ilogit2(combat_combined_corrected)
    message("Successfully removed batch effect of experiment for mission combination ", paste(combination, collapse=", "))

  } else {
    final_beta <- combat_combined
    message("Only one experiment, skip exp batch effect stage")
  }

  # PART THREE：visualization after batch removal.
  set.seed(seeds)
  tsne_out_post <- Rtsne(t(final_beta), pca=FALSE, perplexity=perplexity, theta=0.0)
  tsnes_post <- as.data.frame(tsne_out_post$Y)
  colnames(tsnes_post) <- c("tSNE1", "tSNE2")

  tsnes_post$Subject <- SampleSheet_all$Subject
  tsnes_post$Time <- SampleSheet_all$Sample_Group
  tsnes_post$Experiment <- SampleSheet_all$Experiment
  tsnes_post$Sample_Name <- SampleSheet_all$Sample_Name
  write.csv(tsnes_post, file = file.path(vis_save_path, paste0("tsne_post_combat_", paste(combination, collapse="_"), ".csv")), row.names = FALSE)

  p1_post <- ggplot(tsnes_post, aes(tSNE1, tSNE2)) + 
              geom_point(aes(color=Subject), size=3) +
              scale_color_manual(values=colors_subject) +
              theme(panel.background = element_rect(fill = "white", color = "black"))+  
              theme(panel.grid.major = element_line(color = "gray"), panel.grid.minor = element_line(color = "gray"))

  p2_post <- ggplot(tsnes_post, aes(tSNE1, tSNE2)) + 
              geom_point(aes(color=Time, shape=Experiment), size=3) +
              scale_color_manual(values=colors_time) +
              theme(panel.background = element_rect(fill = "white", color = "black"))+  
              theme(panel.grid.major = element_line(color = "gray"), panel.grid.minor = element_line(color = "gray"))

  post_combat_plot <- ggpubr::ggarrange(p1_post, p2_post, nrow=1)
  ggsave(file.path(vis_save_path, paste0("tsne_post_combat_", paste(combination, collapse="_"), ".pdf")), 
        post_combat_plot, width=14, height=6)
  message("tsne figure after batch removal saved to ", file.path(vis_save_path, paste0("tsne_post_combat_", paste(combination, collapse="_"), ".pdf")))

  # Save plot.pdf and beta_nobatch.Rds.
  final_plot <- ggpubr::ggarrange(pre_combat_plot, post_combat_plot, ncol=1)
  ggsave(file.path(vis_save_path, paste0("combined_tsne_post_combat_", paste(combination, collapse="_"), ".pdf")), 
        final_plot, width=14, height=10)

  start_col <- 1
  for (code in combination) {
    n_samples <- ncol(mission_data[[code]])
    end_col <- start_col + n_samples - 1
    day_number <- gsub("^M(\\d+).*", "\\1", code)
    mission_dir <- file.path(
      paste0(day_number, "day_spaceflight_", code),
      "DNA_methylation",
      "beta"
    )
    saveRDS(
      object = final_beta[, start_col:end_col],
      file = file.path(mission_dir, paste0("nobatch_beta", flag, ".Rds"))
    )
    message("RDS of ", code, " nobatch beta saved to ", mission_dir, paste0("/nobatch_beta", flag, ".Rds"))
    start_col <- end_col + 1
  }
}
