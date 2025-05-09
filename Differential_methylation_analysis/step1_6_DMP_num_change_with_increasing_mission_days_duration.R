### description: using M180‑1 (T2 vs. T1 DMP) as reference, compute consistency ratios (directional and strict) against other missions; plot Δβ distributions for each comparison in a 1×6 PDF layout. ###

# Working directory need to be adjusted according to your path.
setwd("/home/lqwang/Program")

# Command line args.
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-f", "--flag"), type="character", help="flag")
)
parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)
flag <- opts$flag

library(ggplot2)

message("step1_6_DMP_num_change_with_increasing_mission_days_duration.R", flag)

DMP_result_path = paste0("./DMP_result", flag)

width = 16
height = 4

df_DMP_M180_1_T1_to_T2 <- readRDS(file = file.path(DMP_result_path, "M180-1/DMP_res_M180-1.Rds"))[["T1 to T2"]]

df_DMP_M180_2_T1_to_T2 <- readRDS(file = file.path(DMP_result_path, "M180-2/DMP_res_M180-2.Rds"))[["T1 to T2"]]

df_DMP_M90_T1_to_T2 <- readRDS(file = file.path(DMP_result_path, "M90/DMP_res_M90.Rds"))[["T1 to T2"]]

df_DMP_M13_T1_to_T2 <- readRDS(file = file.path(DMP_result_path, "M13/DMP_res_M13.Rds"))[["T1 to T2"]]

df_DMP_M15_T1_to_T2 <- readRDS(file = file.path(DMP_result_path, "M15/DMP_res_M15.Rds"))[["T1 to T2"]]

df_DMP_M33_T1_to_T2 <- readRDS(file = file.path(DMP_result_path, "M33/DMP_res_M33.Rds"))[["T1 to T2"]]

sig_df_DMP_M180_1_T1_to_T2 <- df_DMP_M180_1_T1_to_T2[df_DMP_M180_1_T1_to_T2$change_deltaBeta0.05 != "stable", ]
sigDMPsM180_1_DMP_M180_2_T1_to_T2 <- df_DMP_M180_2_T1_to_T2[rownames(df_DMP_M180_2_T1_to_T2) %in% rownames(sig_df_DMP_M180_1_T1_to_T2), ]
sigDMPsM180_1_DMP_M90_T1_to_T2 <- df_DMP_M90_T1_to_T2[rownames(df_DMP_M90_T1_to_T2) %in% rownames(sig_df_DMP_M180_1_T1_to_T2), ]
sigDMPsM180_1_DMP_M33_T1_to_T2 <- df_DMP_M33_T1_to_T2[rownames(df_DMP_M33_T1_to_T2) %in% rownames(sig_df_DMP_M180_1_T1_to_T2), ]
sigDMPsM180_1_DMP_M15_T1_to_T2 <- df_DMP_M15_T1_to_T2[rownames(df_DMP_M15_T1_to_T2) %in% rownames(sig_df_DMP_M180_1_T1_to_T2), ]
sigDMPsM180_1_DMP_M13_T1_to_T2 <- df_DMP_M13_T1_to_T2[rownames(df_DMP_M13_T1_to_T2) %in% rownames(sig_df_DMP_M180_1_T1_to_T2), ]


consis_percent_df <- c()
for (exp in c("M180_2","M90","M33","M15","M13")) {
     sigDMPsM180_1_DMP_T1_to_T2 <- eval(parse(text = paste0("sigDMPsM180_1_DMP_", exp, "_T1_to_T2")))
     sigDMPsM180_1_DMP_T1_to_T2$enhancer <- df_DMP_M180_1_T1_to_T2[rownames(sigDMPsM180_1_DMP_T1_to_T2), "enhancer"]
     df_DMP_T1_to_T2 <- eval(parse(text = paste0("df_DMP_", exp, "_T1_to_T2")))
     totol_num <- table(sig_df_DMP_M180_1_T1_to_T2[rownames(sigDMPsM180_1_DMP_T1_to_T2), "deltaBeta"] > 0)
     both_up_DMPs_more <- which(sig_df_DMP_M180_1_T1_to_T2[rownames(sigDMPsM180_1_DMP_T1_to_T2),"deltaBeta"] > 0 & sigDMPsM180_1_DMP_T1_to_T2$deltaBeta > 0 & sig_df_DMP_M180_1_T1_to_T2[rownames(sigDMPsM180_1_DMP_T1_to_T2), "deltaBeta"] > sigDMPsM180_1_DMP_T1_to_T2$deltaBeta)
     both_up_DMPs_less <- which(sig_df_DMP_M180_1_T1_to_T2[rownames(sigDMPsM180_1_DMP_T1_to_T2),"deltaBeta"] > 0 & sigDMPsM180_1_DMP_T1_to_T2$deltaBeta > 0 & sig_df_DMP_M180_1_T1_to_T2[rownames(sigDMPsM180_1_DMP_T1_to_T2), "deltaBeta"] < sigDMPsM180_1_DMP_T1_to_T2$deltaBeta)
     if (is.null(df_DMP_T1_to_T2$change_deltaBeta0.05)) {
          both_up_DMPs <- intersect(rownames(sig_df_DMP_M180_1_T1_to_T2)[sig_df_DMP_M180_1_T1_to_T2$change == "up"], rownames(df_DMP_T1_to_T2)[df_DMP_T1_to_T2$change == "up"])
     }else{
          both_up_DMPs <- intersect(rownames(sig_df_DMP_M180_1_T1_to_T2)[sig_df_DMP_M180_1_T1_to_T2$change_deltaBeta0.05=="up"],rownames(df_DMP_T1_to_T2)[df_DMP_T1_to_T2$change_deltaBeta0.05=="up"])
     }

     both_down_DMPs_more <- which(sig_df_DMP_M180_1_T1_to_T2[rownames(sigDMPsM180_1_DMP_T1_to_T2), "deltaBeta"] < 0 & sigDMPsM180_1_DMP_T1_to_T2$deltaBeta < 0 & sig_df_DMP_M180_1_T1_to_T2[rownames(sigDMPsM180_1_DMP_T1_to_T2), "deltaBeta"] < sigDMPsM180_1_DMP_T1_to_T2$deltaBeta)
     both_down_DMPs_less <- which(sig_df_DMP_M180_1_T1_to_T2[rownames(sigDMPsM180_1_DMP_T1_to_T2),"deltaBeta"]<0 & sigDMPsM180_1_DMP_T1_to_T2$deltaBeta<0 & sig_df_DMP_M180_1_T1_to_T2[rownames(sigDMPsM180_1_DMP_T1_to_T2), "deltaBeta"] > sigDMPsM180_1_DMP_T1_to_T2$deltaBeta)
     if (is.null(df_DMP_T1_to_T2$change_deltaBeta0.05)) {
          both_down_DMPs <- intersect(rownames(sig_df_DMP_M180_1_T1_to_T2)[sig_df_DMP_M180_1_T1_to_T2$change == "down"], rownames(df_DMP_T1_to_T2)[df_DMP_T1_to_T2$change == "down"])
     }else{
          both_down_DMPs <- intersect(rownames(sig_df_DMP_M180_1_T1_to_T2)[sig_df_DMP_M180_1_T1_to_T2$change_deltaBeta0.05 == "down"], rownames(df_DMP_T1_to_T2)[df_DMP_T1_to_T2$change_deltaBeta0.05 == "down"])
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
rownames(consis_percent_df) <- paste0("sigDMPs_M180-1_vs_", c("M180-2", "M90", "M33", "M15", "M13"))
print(consis_percent_df)


p <- ggplot(sig_df_DMP_M180_1_T1_to_T2, aes(x = deltaBeta, fill = change_deltaBeta0.05, alpha = 0.01)) + geom_histogram(aes(y = ..density..), colour = "black", binwidth = 0.02) +
     ggtitle('DMPs_M180-1 in M180-1') + theme(plot.title = element_text(size = 9)) +
     geom_density() + scale_fill_manual(values = c(up = "#FFFFCC", down = "#99CCFF"))
p <- p + geom_vline(aes(xintercept = 0.05, colour = "#BB0000", linetype = "dashed")) + geom_vline(aes(xintercept = -0.05, colour = "#BB0000", linetype = "dashed"))
plot1 <- p + scale_x_continuous(limits = c(-0.25, 0.25)) + scale_y_continuous(limits = c(0, 30)) + coord_flip()

p <- ggplot(sigDMPsM180_1_DMP_M180_2_T1_to_T2, aes(x = deltaBeta, alpha = 0.01)) + geom_histogram(aes(y = ..density..), colour = "black", binwidth = 0.02) +
     ggtitle('DMPs_M180-1 in M180-2') + theme(plot.title = element_text(size = 9)) +
     geom_density()
p <- p + geom_vline(aes(xintercept = 0.05, colour = "#BB0000", linetype = "dashed")) + geom_vline(aes(xintercept = -0.05, colour = "#BB0000", linetype="dashed"))
plot2 <- p + scale_x_continuous(limits = c(-0.25, 0.25)) + scale_y_continuous(limits = c(0, 30)) + coord_flip()

p <- ggplot(sigDMPsM180_1_DMP_M90_T1_to_T2, aes(x = deltaBeta, alpha = 0.01)) + geom_histogram(aes(y = ..density..), colour = "black", binwidth = 0.02) +
     ggtitle('DMPs_M180-1 in M90') + theme(plot.title = element_text(size = 9)) +
     geom_density()
p <- p + geom_vline(aes(xintercept = 0.05, colour="#BB0000", linetype="dashed")) + geom_vline(aes(xintercept = -0.05, colour = "#BB0000", linetype="dashed"))
plot3 <- p + scale_x_continuous(limits = c(-0.25, 0.25)) + scale_y_continuous(limits = c(0, 30)) + coord_flip()

p <- ggplot(sigDMPsM180_1_DMP_M33_T1_to_T2,aes(x=deltaBeta, alpha = 0.01))+ geom_histogram(aes(y = ..density..),colour = "black", binwidth = 0.02) +
     ggtitle('DMPs_M180-1 in M33')+ theme(plot.title = element_text(size = 9)) +
     geom_density()
p <- p + geom_vline(aes(xintercept=0.05, colour="#BB0000", linetype="dashed")) + geom_vline(aes(xintercept=-0.05, colour="#BB0000", linetype="dashed"))
plot4 <- p + scale_x_continuous(limits = c(-0.25, 0.25)) + scale_y_continuous(limits = c(0, 30))+coord_flip()

p <- ggplot(sigDMPsM180_1_DMP_M15_T1_to_T2,aes(x=deltaBeta, alpha = 0.01))+ geom_histogram(aes(y = ..density..),colour = "black", binwidth = 0.02) +
     ggtitle('DMPs_M180-1 in M15')+ theme(plot.title = element_text(size = 9)) +
     geom_density()
p <- p + geom_vline(aes(xintercept=0.05, colour="#BB0000", linetype="dashed")) + geom_vline(aes(xintercept=-0.05, colour="#BB0000", linetype="dashed"))
plot5 <- p + scale_x_continuous(limits = c(-0.25, 0.25)) + scale_y_continuous(limits = c(0, 30))+coord_flip()

p <- ggplot(sigDMPsM180_1_DMP_M13_T1_to_T2,aes(x=deltaBeta, alpha = 0.01))+ geom_histogram(aes(y = ..density..),colour = "black", binwidth = 0.02) +
     ggtitle('DMPs_M180-1 in M13')+ theme(plot.title = element_text(size = 9))+
     geom_density()
p <- p + geom_vline(aes(xintercept=0.05, colour="#BB0000", linetype="dashed")) + geom_vline(aes(xintercept=-0.05, colour="#BB0000", linetype="dashed"))
plot6 <- p + scale_x_continuous(limits = c(-0.25, 0.25)) + scale_y_continuous(limits = c(0, 30))+coord_flip()

projection_M180_1 <- ggpubr::ggarrange(plot1,plot2,plot3,plot4,plot5,plot6, nrow=1,ncol=6, legend = "none")

ggsave(filename = "projection_M180_1.pdf", plot = projection_M180_1, path = DMP_result_path, width = width, height = height, device = "pdf")
