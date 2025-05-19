# Description: Density of gene expression logFCs in T3_vs_T1 (Fig. 4b) 
# DEG_list_trends.Rds, T1_to_T2_DEG_table.txt, T1_to_T2_DEG_table.txt, T1_to_T3_DEG_table.txt and DEG_list.Rds are available in "DEG_analysis"

# Load packages
library(ggplot2)

# 01Start: logFC density of DEGs in M90
workdir <- "90-day_spaceflight_M90/RNA-seq/DE_analysis/"
setwd(workdir)
DEG_T1_T2_01 <- read.table( file = "T1_to_T2_DEG_table.txt", header = T, sep = "\t")
DEG_T2_T3_01 <- read.table( file = "T2_to_T3_DEG_table.txt", header = T, sep = "\t")
table(DEG_T1_T2_01$change)
DEG_list_trends_01 <- readRDS(file = "DEG_list_trends.Rds")
DEG_T1_to_T3_01 <- read.table( file = "T1_to_T3_DEG_table.txt", header = T, sep = "\t", row.names = 1)
sig_DEG_T1_to_T3_01 <- DEG_T1_to_T3_01[unlist(DEG_list_trends_01),]
sig_DEG_T1_to_T3_01$change_type <- unlist(lapply(names(DEG_list_trends_01), function(x){rep(x, length(DEG_list_trends_01[[x]]))}))
sig_DEG_T1_to_T3_01$change_type <- factor(sig_DEG_T1_to_T3_01$change_type, levels = c("down_stable","up_stable","stable_up","stable_down","down_up","up_down"))
write.csv(sig_DEG_T1_to_T3_01, file = "density of logFC T3_vs_T1 M90.csv")
p_T3_vs_T1_01 <- ggplot(sig_DEG_T1_to_T3_01,aes(x=logFC,fill=change_type, alpha = 0.01))+geom_density()+ 
                 scale_x_continuous(limits = c(-2, 2))  +
                 scale_fill_manual(values=c(down_up = "#83AEDC", stable_up = "#5F559B", up_down = "#C14E8B", stable_down = "#985CB8", up_stable="#4FB895", down_stable="#C88AB9"))
ggsave("Fig. 4b M90.pdf", plot = p_T3_vs_T1_01, width = 6, height = 4)
# 01End: M90 DEG logFC density 


# 02Start: logFC density of DEGs in M180-1
workdir <- "180-day_spaceflight_M180-1/RNA-seq/DE_analysis/"
setwd(workdir)
DEG_T1_T2_02 <- read.table( file = "T1_to_T2_DEG_table.txt", header = T, sep = "\t")
DEG_T2_T3_02 <- read.table( file = "T2_to_T3_DEG_table.txt", header = T, sep = "\t")
DEG_list_trends_02 <- readRDS(file = "DEG_list_trends.Rds")
table(DEG_T1_T2_02$change)
DEG_T1_to_T3_02 <- read.table( file = "T1_to_T3_DEG_table.txt", header = T, sep = "\t", row.names = 1)
sig_DEG_T1_to_T3_02 <- DEG_T1_to_T3_02[unlist(DEG_list_trends_02),]
sig_DEG_T1_to_T3_02$change_type <- unlist(lapply(names(DEG_list_trends_02), function(x){rep(x, length(DEG_list_trends_02[[x]]))}))
sig_DEG_T1_to_T3_02$change_type <- factor(sig_DEG_T1_to_T3_02$change_type, levels = c("down_stable","up_stable","stable_up","stable_down","down_up","up_down"))
write.csv(sig_DEG_T1_to_T3_02, file = "density of logFC T3_vs_T1 M180-1.csv")
p_T3_vs_T1_02 <- ggplot(sig_DEG_T1_to_T3_02,aes(x=logFC,fill=change_type, alpha = 0.01))+geom_density()+ 
                 scale_x_continuous(limits = c(-2, 2))  +
                 scale_fill_manual(values=c(down_up = "#83AEDC", stable_up = "#5F559B", up_down = "#C14E8B", stable_down = "#985CB8", up_stable="#4FB895", down_stable="#C88AB9"))
ggsave("Fig. 4b M180-1.pdf", plot = p_T3_vs_T1_02, width = 6, height = 4)
# 02End: M180-1 DEG logFC density

# Combination of diagrams
logFC_density <- ggpubr::ggarrange(p_T3_vs_T1_01, p_T3_vs_T1_02, ncol=1, nrow=2)
ggsave("Fig. 4b.pdf", plot = logFC_density, width = 6, height = 8)
