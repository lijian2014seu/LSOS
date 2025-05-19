##### description: Gene set methylation change trends in M90 and M180-1
setwd("/home/lqwang/Program")
#load packages
library(readxl)
library(RColorBrewer)
library(ggplot2)

source("./main/utils/consist_judge_fun.R")
flag = "_2_1"
DMP_result_path = paste0("./DMP_result", flag)
Multi_com_res_path = "./MultiCompare_result"
interested_terms = c("Osteoclast differentiation", "blood coagulation", "NF-kappa B signaling pathway", "Platelet activation")
#Load genesets enriched by PPI nodes 
go_enrich_res <- read.table( file = "./main/WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_GOenrich_res.txt", sep = "\t")
kegg_enrich_res <- read.table( file = "./main/WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_KEGGenrich_res.txt", sep = "\t")
geneset_all_df_PPI <- rbind(go_enrich_res, kegg_enrich_res)

#trend_matrix: col:gene; rownames:time by order; value:mean of methlation beta value of the genes.
trend_plot <- function (trend_matrix) {
  time <- colnames(trend_matrix)
  trend_matrix <- t(apply(trend_matrix, 1, function (x) scale(x, center = F)))

  t <- 1: length(time)
  data <- as.data.frame(t(trend_matrix))
  data$t <- t
  data_long <- reshape::melt(data, id.vars = "t")
  x <- 1: length(time)
  quantile_value <- apply(trend_matrix, 2, function (x) {
    temp <- round(quantile(x, seq(0, 1, 0.01)), 3)
    return(temp)
  })
  pdata.list <- list()
  for(k in 1: 50){
    pdata <- data.frame(x, lower = quantile_value[k, ], upper = quantile_value[102 - k, ])
    pdata.list[[k]] <- pdata
  }
  my_col <- c("#C6E4F5", "#2BAAF3")
  myPalette <- colorRampPalette(my_col)(45)
  plot.trend<-ggplot() +
    geom_rect(aes(xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf), fill = "#FFFFCC", alpha = 0.4) +
    geom_rect(aes(xmin = 1, xmax = length(time), ymin = -Inf, ymax = Inf), fill = "#DEF6F3", alpha = 0.4) +
    geom_rect(aes(xmin = length(time), xmax = Inf, ymin = -Inf, ymax = Inf), fill = "#FFFFCC", alpha = 0.4) +
    theme(panel.grid = element_blank(), panel.border = element_rect(fill = NA, color = "black", size = 1))

  for(k in 6: 50) {
    pdata <- pdata.list[[k]]
    plot.trend <- plot.trend + geom_ribbon(data = pdata, aes(ymin = lower, ymax = upper, x = x), fill = myPalette[k - 5], alpha = 1)
  }
  plot.trend <- plot.trend +
    theme_bw() +
    xlab("") +
    theme(panel.grid = element_blank()) + 
    scale_x_continuous(breaks = seq(1, length(time), 1), labels = time)
  return (plot.trend)
}


missions = list(list(name = "M90", 
                     beta_path = paste0("./90day_spaceflight_M90/DNA_methylation/beta/nobatch_beta", flag, ".Rds"), 
                     sample_sheet_path = "./90day_spaceflight_M90/DNA_methylation/rawdata/SampleSheet.csv"),
                list(name = "M180-1",
                     beta_path = paste0("./180day_spaceflight_M180_1/DNA_methylation/beta/nobatch_beta", flag, ".Rds"),
                     sample_sheet_path = "./180day_spaceflight_M180_1/DNA_methylation/rawdata/SampleSheet.csv")
                )
for (exp in missions) {
  df_DMP_T1_to_T2 <- readRDS(file.path(DMP_result_path, exp$name, paste0("DMP_res_", exp$name, ".Rds")))[["T1 to T2"]]
  myCombat <- readRDS(exp$beta_path)
  SampleSheet <- read.csv(exp$sample_sheet_path)
  plot_list_flight <- c()
  for (term in interested_terms) {
    interested_geneset <- geneset_all_df_PPI[geneset_all_df_PPI$Description == term, "geneID"]
    interested_geneset <- unique(unlist(strsplit(interested_geneset, "/")))
    
    promoter_region <- c("TSS1500", "TSS200", "5'UTR", "1stExon")
    interested_geneset_annotation <- df_DMP_T1_to_T2[df_DMP_T1_to_T2$gene%in%interested_geneset & df_DMP_T1_to_T2$feature %in% promoter_region, ]
    interested_geneset_meth_probe <- as.data.frame(myCombat[rownames(interested_geneset_annotation), ])
    interested_geneset_meth_probe$gene <- interested_geneset_annotation$gene

    interested_geneset_meth_mat <- aggregate(interested_geneset_meth_probe[, -ncol(interested_geneset_meth_probe)], by = list(group = interested_geneset_annotation$gene), median)
    rownames(interested_geneset_meth_mat) <- interested_geneset_meth_mat$group
    interested_geneset_meth_mat <- interested_geneset_meth_mat[, -1]

    consist_res_T1T2 <- consist_judge_fun("T1", "T2", interested_geneset_meth_mat, SampleSheet)
    consist_res_T2T3 <- consist_judge_fun("T2", "T3", interested_geneset_meth_mat, SampleSheet)

    consist_index <- which(consist_res_T1T2 != "Inconsistent" & consist_res_T2T3 != "Inconsistent")
    filter_interested_geneset_meth_mat <- interested_geneset_meth_mat[consist_index, ]

    if ( nrow(filter_interested_geneset_meth_mat) > 1 ) {
      r <- apply( filter_interested_geneset_meth_mat, 1, function (x) {max(x) - min(x)} )
      filter_interested_geneset_meth_mat <- filter_interested_geneset_meth_mat[r > 0.05,]

      if ( nrow(filter_interested_geneset_meth_mat) > 1 ) {
        mean_meth_mat <- apply(filter_interested_geneset_meth_mat, 1, function(x) tapply(x, SampleSheet$Sample_Group, mean))
        mean_meth_mat <- t(mean_meth_mat)
        print(dim(mean_meth_mat))
        trend_p <- trend_plot(mean_meth_mat)
        plot_list_flight <- c(plot_list_flight, list(trend_p))
        term <- stringr::str_split_i(term, " - ", 1)
        names(plot_list_flight)[length(plot_list_flight)] <- term
      }
    }
  }


  saveRDS( plot_list_flight, file.path(Multi_com_res_path, paste0("plot_list_flight", exp$name, ".Rds")) )
  for (i in 1: length(plot_list_flight)) {
    filename <- paste("flight", exp$name, "_", names(plot_list_flight)[i], ".pdf", sep = "")
    pdf(file.path(Multi_com_res_path, filename), width=5, height=3)
    trend_p <- plot_list_flight[i]
    print(trend_p)
    dev.off()
  }
}


