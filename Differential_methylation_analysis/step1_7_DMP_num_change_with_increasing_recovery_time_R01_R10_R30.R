### description: using R01 vs. T1 DMP as reference, 
#  - compute consistency ratios (directional and strict) against other recovery comparisons for M13, M15, M33; 
#  - plot Δβ distributions for each comparison in a 1×3 PDF layout. ###

# Working directory need to be adjusted according to your path.
setwd("/home/lqwang/Program")

# Command line args.
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-f", "--flag"), type="character", help="flag"),
  make_option(c("-m", "--mission"), type="character", help="mission")
)
parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)
flag <- opts$flag
mission <- opts$mission

library(ggplot2)

message("step1_7_DMP_num_change_with_increasing_recovery_time_R01_R10_R30.R", flag)

DMP_result_path = paste0("./DMP_result", flag)

width = 8
height = 4


df_DMP_T1_to_T2 <- readRDS(file = file.path(DMP_result_path, mission, paste0("DMP_res_", mission, ".Rds")))[["T1 to T2"]]
df_DMP_T1_to_D10_pro <- readRDS(file = file.path(DMP_result_path, mission, paste0("DMP_res_", mission, ".Rds")))[["T1 to D10_pro"]]
df_DMP_T1_to_T3 <- readRDS(file = file.path(DMP_result_path, mission, paste0("DMP_res_", mission, ".Rds")))[["T1 to T3"]]


sig_df_DMP_T1_to_T2 <- df_DMP_T1_to_T2[df_DMP_T1_to_T2$change_deltaBeta0.05 %in% c("up", "down"), ]
sigDMPsT1_to_T2_DMP_T1_to_D10_pro <- df_DMP_T1_to_D10_pro[rownames(df_DMP_T1_to_D10_pro) %in% rownames(sig_df_DMP_T1_to_T2), ]
sigDMPsT1_to_T2_DMP_T1_to_T3 <- df_DMP_T1_to_T3[rownames(df_DMP_T1_to_T3) %in% rownames(sig_df_DMP_T1_to_T2), ]


consis_percent_df <- c()
for ( exp in c("T1_to_D10_pro", "T1_to_T3") ) {
     sig_df_DMP_T1_to_T2_x <- eval(parse(text = paste0("sigDMPsT1_to_T2_DMP_", exp)))
     df_DMP <- eval(parse(text = paste0("df_DMP_", exp)))
     totol_num <- table(sig_df_DMP_T1_to_T2_x[rownames(sig_df_DMP_T1_to_T2_x), "deltaBeta"] > 0)
     both_up_DMPs_more <- which(sig_df_DMP_T1_to_T2_x[rownames(sig_df_DMP_T1_to_T2_x),"deltaBeta"] > 0 & sig_df_DMP_T1_to_T2_x$deltaBeta > 0 & sig_df_DMP_T1_to_T2_x[rownames(sig_df_DMP_T1_to_T2_x), "deltaBeta"] > sig_df_DMP_T1_to_T2_x$deltaBeta)
     both_up_DMPs_less <- which(sig_df_DMP_T1_to_T2_x[rownames(sig_df_DMP_T1_to_T2_x),"deltaBeta"] > 0 & sig_df_DMP_T1_to_T2_x$deltaBeta > 0 & sig_df_DMP_T1_to_T2_x[rownames(sig_df_DMP_T1_to_T2_x), "deltaBeta"] < sig_df_DMP_T1_to_T2_x$deltaBeta)
     if (is.null(df_DMP$change_deltaBeta0.05)) {
          both_up_DMPs <- intersect(rownames(sig_df_DMP_T1_to_T2_x)[sig_df_DMP_T1_to_T2_x$change == "up"], rownames(df_DMP)[df_DMP$change == "up"])
     }else{
          both_up_DMPs <- intersect(rownames(sig_df_DMP_T1_to_T2_x)[sig_df_DMP_T1_to_T2_x$change_deltaBeta0.05=="up"],rownames(df_DMP)[df_DMP$change_deltaBeta0.05=="up"])
     }

     both_down_DMPs_more <- which(sig_df_DMP_T1_to_T2_x[rownames(sig_df_DMP_T1_to_T2_x), "deltaBeta"] < 0 & sig_df_DMP_T1_to_T2_x$deltaBeta < 0 & sig_df_DMP_T1_to_T2_x[rownames(sig_df_DMP_T1_to_T2_x), "deltaBeta"] < sig_df_DMP_T1_to_T2_x$deltaBeta)
     both_down_DMPs_less <- which(sig_df_DMP_T1_to_T2_x[rownames(sig_df_DMP_T1_to_T2_x),"deltaBeta"]<0 & sig_df_DMP_T1_to_T2_x$deltaBeta<0 & sig_df_DMP_T1_to_T2_x[rownames(sig_df_DMP_T1_to_T2_x), "deltaBeta"] > sig_df_DMP_T1_to_T2_x$deltaBeta)
     if (is.null(df_DMP$change_deltaBeta0.05)) {
          both_down_DMPs <- intersect(rownames(sig_df_DMP_T1_to_T2_x)[sig_df_DMP_T1_to_T2_x$change == "down"], rownames(df_DMP)[df_DMP$change == "down"])
     }else{
          both_down_DMPs <- intersect(rownames(sig_df_DMP_T1_to_T2_x)[sig_df_DMP_T1_to_T2_x$change_deltaBeta0.05 == "down"], rownames(df_DMP)[df_DMP$change_deltaBeta0.05 == "down"])
     }
     both_up_percent <- (length(both_up_DMPs_less) + length(both_up_DMPs_more)) / totol_num["TRUE"]
     both_sig_up_percent <- length(both_up_DMPs) / totol_num["TRUE"]
     both_down_percent <- (length(both_down_DMPs_less) + length(both_down_DMPs_more)) / totol_num["FALSE"]
     both_sig_down_percent <- length(both_down_DMPs) / totol_num["FALSE"]	
     consis_ratio <- (length(both_up_DMPs_less) + length(both_up_DMPs_more) + length(both_down_DMPs_less) + length(both_down_DMPs_more)) / sum(totol_num)
     consis_less_ratio <- (length(both_up_DMPs_less) + length(both_down_DMPs_less)) / sum(totol_num)
     consis_more_ratio <- (length(both_up_DMPs_more) + length(both_down_DMPs_more)) / sum(totol_num)
     strictly_consis_ratio <- (length(both_up_DMPs) + length(both_down_DMPs)) / sum(totol_num)
     consis_percent <- c(both_up_percent, both_sig_up_percent, both_down_percent, both_sig_down_percent, consis_ratio, strictly_consis_ratio, consis_less_ratio, consis_more_ratio)
     print(consis_percent)
     consis_percent_df <- rbind(consis_percent_df, consis_percent)
}
colnames(consis_percent_df) <- c("both_up_percent", "both_sig_up_percent", "both_down_percent", "both_sig_down_percent", "consis_ratio", "strictly_consis_ratio", "consis_less_ratio", "consis_more_ratio")
rownames(consis_percent_df) <- paste0( "sigDMPsT2_to_D10_pro_DMP_", c("T1_to_D10_pro", "T1_to_T3") )
print(consis_percent_df)


p <- ggplot(sig_df_DMP_T1_to_T2, aes(x = deltaBeta, fill = change_deltaBeta0.05, alpha = 0.01)) + geom_histogram(aes(y = ..density..), colour = "black", binwidth = 0.02) +
     ggtitle('DMPs R01 in R01') + theme(plot.title = element_text(size = 9)) +
     geom_density() + scale_fill_manual(values = c(up = "#FFFFCC", down = "#99CCFF"))
p <- p + geom_vline(aes(xintercept = 0.05, colour = "#BB0000", linetype = "dashed")) + geom_vline(aes(xintercept = -0.05, colour = "#BB0000", linetype = "dashed"))
plot1 <- p + scale_x_continuous(limits = c(-0.25, 0.25)) + scale_y_continuous(limits = c(0, 30)) + coord_flip()

p <- ggplot(sigDMPsT1_to_T2_DMP_T1_to_D10_pro, aes(x = deltaBeta, alpha = 0.01)) + geom_histogram(aes(y = ..density..), colour = "black", binwidth = 0.02) +
     ggtitle('DMPs R01 in R10') + theme(plot.title = element_text(size = 9)) +
     geom_density()
p <- p + geom_vline(aes(xintercept = 0.05, colour = "#BB0000", linetype = "dashed")) + geom_vline(aes(xintercept = -0.05, colour = "#BB0000", linetype="dashed"))
plot2 <- p + scale_x_continuous(limits = c(-0.25, 0.25)) + scale_y_continuous(limits = c(0, 30)) + coord_flip()

p <- ggplot(sigDMPsT1_to_T2_DMP_T1_to_T3, aes(x = deltaBeta, alpha = 0.01)) + geom_histogram(aes(y = ..density..), colour = "black", binwidth = 0.02) +
     ggtitle('DMPs R01 in R30') + theme(plot.title = element_text(size = 9)) +
     geom_density()
p <- p + geom_vline(aes(xintercept = 0.05, colour = "#BB0000", linetype = "dashed")) + geom_vline(aes(xintercept = -0.05, colour = "#BB0000", linetype="dashed"))
plot3 <- p + scale_x_continuous(limits = c(-0.25, 0.25)) + scale_y_continuous(limits = c(0, 30)) + coord_flip()


projection <- ggpubr::ggarrange(plot1, plot2, plot3, nrow=1,ncol=3, legend = "none")

ggsave(filename = paste0("projection_recovery_", mission, ".pdf"), plot = projection, path = DMP_result_path, width = width, height = height, device = "pdf")
