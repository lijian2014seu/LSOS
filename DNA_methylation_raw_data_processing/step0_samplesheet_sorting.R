### Sort samplesheet by row, to make sure sample_well and sample_group in the order of H01, H02, H03 and T1, T2, T3. ###

# Working directory need to be adjusted according to your path.
setwd("/home/lqwang/Program")

# Rawdata including .idat files and samplesheets should be downloaded from online database, you can refer to the paper for link. Files should be organized as the format below for execution.
file_paths <- c(
  "./13day_spaceflight_M13/DNA_methylation/rawdata/SampleSheet.csv",
  "./15day_spaceflight_M15/DNA_methylation/rawdata/SampleSheet.csv",
  "./33day_spaceflight_M33/DNA_methylation/rawdata/SampleSheet.csv",
  "./90day_spaceflight_M90/DNA_methylation/rawdata/SampleSheet.csv",
  "./180day_spaceflight_M180_1/DNA_methylation/rawdata/SampleSheet.csv",
  "./180day_spaceflight_M180_2/DNA_methylation/rawdata/SampleSheet.csv"
)

# Rewrite samplesheet in the right order.
for (file_path in file_paths) {
  tryCatch({
    if (!file.exists(file_path)) {
      warning(paste("File not exits:", file_path), immediate. = TRUE)
      next
    }
    data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
    if (!"Sample" %in% colnames(data)) {
      warning(paste("File", file_path, "lack column Sample, skipped"), immediate. = TRUE)
      next
    }
    
    sorted_data <- data[order(data$Sample), ]
    
    write.csv(sorted_data, file = file_path, row.names = FALSE, na = "", quote = F)
    message(paste("Successfully processed file:", file_path))
  }, error = function(e) {
    warning(paste("Error occurs when processing file", file_path, ":", e$message), immediate. = TRUE)
  })
}

message("done")