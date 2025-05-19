# Description: Differential expression analysis (DEG) consistency checks across multiple experimental groups (Extended Data Figure 7b) 
# DEG_list_trends.Rds, T1_to_T2_DEG_table.txt, T1_to_T2_DEG_table.txt, T1_to_T3_DEG_table.txt and DEG_list.Rds are available in "DEG_analysis"

# Load packages
library(ggplot2)

# 01Start: DEG consistency checks
DEG_01_T1_to_T2 <- read.table(file = "90-day spaceflight_M90/RNA-seq/DE_analysis/T1_to_T2_DEG_table.txt")
DEG_02_T1_to_T2 <- read.table(file = "180-day spaceflight_M180-1/RNA-seq/DE_analysis/T1_to_T2_DEG_table.txt")
DEG_03_T1_to_T2 <- read.table(file = "180-day spaceflight_M180-2/RNA-seq/DE_analysis/T1_to_T2_DEG_table.txt")

head(DEG_02_T1_to_T2)
sig_DEG_02_T1_to_T2 <- DEG_02_T1_to_T2[DEG_02_T1_to_T2$change!="Stable",]
sigDEGs02_DEG_03_T1_to_T2 <- DEG_03_T1_to_T2[rownames(DEG_03_T1_to_T2)%in%rownames(sig_DEG_02_T1_to_T2),]
sigDEGs02_DEG_01_T1_to_T2 <- DEG_01_T1_to_T2[rownames(DEG_01_T1_to_T2)%in%rownames(sig_DEG_02_T1_to_T2),]

p <- ggplot(sig_DEG_02_T1_to_T2,aes(x=logFC,fill=change, alpha = 0.01))+ geom_histogram(aes(y = ..density..),colour = "black", binwidth = 0.1) +
     ggtitle('DEGs_M180-1 in M180-1')+ theme(plot.title = element_text(size = 9)) +
     geom_density()+ scale_fill_manual(values=c(Up = "#FFFFCC", Down = "#99CCFF"))
p <- p + geom_vline(aes(xintercept=log2(1.5), colour="#BB0000", linetype="dashed")) + geom_vline(aes(xintercept=-log2(1.5), colour="#BB0000", linetype="dashed"))
plot1 <- p + scale_x_continuous(limits = c(-2, 2)) + scale_y_continuous(limits = c(0, 5))+coord_flip()

p <- ggplot(sigDEGs02_DEG_03_T1_to_T2,aes(x=logFC, alpha = 0.01))+ geom_histogram(aes(y = ..density..),colour = "black", binwidth = 0.1) +
     ggtitle('DEGs_M180-1 in M180-2')+ theme(plot.title = element_text(size = 9)) +
     geom_density()
p <- p + geom_vline(aes(xintercept=log2(1.5), colour="#BB0000", linetype="dashed")) + geom_vline(aes(xintercept=-log2(1.5), colour="#BB0000", linetype="dashed"))
plot2 <- p + scale_x_continuous(limits = c(-2, 2)) + scale_y_continuous(limits = c(0, 5))+coord_flip()

p <- ggplot(sigDEGs02_DEG_01_T1_to_T2,aes(x=logFC, alpha = 0.01))+ geom_histogram(aes(y = ..density..),colour = "black", binwidth = 0.1) +
     ggtitle('DEGs_M180-1 in M90')+ theme(plot.title = element_text(size = 9)) +
     geom_density()
p <- p + geom_vline(aes(xintercept=log2(1.5), colour="#BB0000", linetype="dashed")) + geom_vline(aes(xintercept=-log2(1.5), colour="#BB0000", linetype="dashed"))
plot3 <- p + scale_x_continuous(limits = c(-2, 2)) + scale_y_continuous(limits = c(0, 5))+coord_flip()

plot_change_compare <- ggpubr::ggarrange(plot1,plot2,plot3,ncol = 3)
ggsave("Extended Data Figure 7b.pdf", plot = plot3, width = 4, height = 4)
# 01End: DEG consistency checks
