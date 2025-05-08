### Processing of raw DNA methylation data from Illumina 450K and EPIC platform. ###

# Working directory need to be adjusted according to your path.
setwd("/home/lqwang/Program")

library(ChAMP)

# Set input and output paths.
dataset_paths <- list(
  list(
    input_dir = "./13day_spaceflight_M13/DNA_methylation/rawdata",
    beta_path = "./13day_spaceflight_M13/DNA_methylation/beta/beta.Rds",
    norm_path = "./13day_spaceflight_M13/DNA_methylation/beta/Norm_beta.Rds",
    data_type = "450K"
  ),
  list(
    input_dir = "./15day_spaceflight_M15/DNA_methylation/rawdata",
    beta_path = "./15day_spaceflight_M15/DNA_methylation/beta/beta.Rds",
    norm_path = "./15day_spaceflight_M15/DNA_methylation/beta/Norm_beta.Rds",
    data_type = "450K"
  ),
  list(
    input_dir = "./33day_spaceflight_M33/DNA_methylation/rawdata",
    beta_path = "./33day_spaceflight_M33/DNA_methylation/beta/beta.Rds",
    norm_path = "./33day_spaceflight_M33/DNA_methylation/beta/Norm_beta.Rds",
    data_type = "450K"
  ),
  list(
    input_dir = "./90day_spaceflight_M90/DNA_methylation/rawdata",
    beta_path = "./90day_spaceflight_M90/DNA_methylation/beta/beta.Rds",
    norm_path = "./90day_spaceflight_M90/DNA_methylation/beta/Norm_beta.Rds",
    data_type = "EPIC"
  ),
  list(
    input_dir = "./180day_spaceflight_M180_1/DNA_methylation/rawdata",
    beta_path = "./180day_spaceflight_M180_1/DNA_methylation/beta/beta.Rds",
    norm_path = "./180day_spaceflight_M180_1/DNA_methylation/beta/Norm_beta.Rds",
    data_type = "EPIC"
  ),
  list(
    input_dir = "./180day_spaceflight_M180_2/DNA_methylation/rawdata",
    beta_path = "./180day_spaceflight_M180_2/DNA_methylation/beta/beta.Rds",
    norm_path = "./180day_spaceflight_M180_2/DNA_methylation/beta/Norm_beta.Rds",
    data_type = "EPIC"
  )
)

# Load raw DNA methylation data
for (path in dataset_paths) {
    if (!dir.exists(dirname(path$norm_path))) {
        dir.create(dirname(path$norm_path), recursive = TRUE, showWarnings = FALSE)
        message("Created directory:", dirname(path$norm_path))
    }
    myload <- champ.load(
        directory = path$input_dir, methValue = "B", 
        filterXY = FALSE, filterDetP = TRUE, detPcut = 0.01, 
        filterBeads = TRUE, beadCutoff = 0.05, filterNoCG = FALSE, 
        filterSNPs = TRUE, filterMultiHit = TRUE, arraytype = path$data_type
        )
    message(path$input_dir, " loaded!")
    saveRDS(myload, file = path$beta_path)
    message("Myload has been saved to ", path$beta_path)
    message()
}


# Normalization of DNA methylation data
for (path in dataset_paths) {
    myload <- readRDS(path$beta_path)
    myNorm <- champ.norm(
        beta = myload$beta, rgSet = myload$rgSet, mset = myload$mset,
        resultsDir = "./CHAMP_Normalization/", method = "BMIQ",
        plotBMIQ = TRUE, arraytype = path$data_type, cores = 8
        )
    message(path$beta_path, " has been normalized!")
    saveRDS(myNorm, file = path$norm_path)
    message("MyNorm has been saved to ", path$norm_path)
    message()
}
