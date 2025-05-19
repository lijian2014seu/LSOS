# Description: Batch removal and Dimensionality reduction for data visualization for RNA-seq data (Extended Data Fig. 5)

# Load packages
library(Rtsne)
library(ggplot2)
library(ggpubr)

# 01Start: Dimensionality reduction for data visualization (M90/01, M180-1/02) before batch removal
# Load count data of 90-day experiment (M90),count_data_symbol.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 90-day spaceflight_M90, RNA-seq)
my_dir <- "Data/RNA-seq/"
setwd(paste0(my_dir, "01/"))
count_data_01 <- read.table("90-day spaceflight_M90/RNA-seq/count_data_symbol.txt", sep="\t", as.is=T, header=T, row.names = 1)
SampleSheet_01 <- read.csv("90-day spaceflight_M90/RNA-seq/SampleSheet.csv")
sum(colnames(count_data_01)==SampleSheet_01$Sample_Name)

# Load count data of 180-day experiment (M180-1),count_data_symbol.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 180-day spaceflight_M180-1, RNA-seq)
my_dir <- "Data/RNA-seq/"
setwd(paste0(my_dir, "02/"))
count_data_02 <- read.table("180-day spaceflight_M180-1/RNA-seq/count_data_symbol.txt", sep="\t", as.is=T, header=T, row.names = 1)
SampleSheet_02 <- read.csv("180-day spaceflight_M180-1/RNA-seq/SampleSheet.csv")
sum(colnames(count_data_02)==SampleSheet_02$Sample_Name)

# Merging of expression profiles
genes <- intersect(rownames(count_data_01), rownames(count_data_02))
count_data_all <- cbind(count_data_01[genes,], count_data_02[genes,])
count_data_all_log <- log2(count_data_all + 1)
SampleSheet_all <- rbind(SampleSheet_01, SampleSheet_02)
SampleSheet_all$Sample_Name <- paste0(SampleSheet_all$Sample_Name, c(rep("_90",9), rep("_180-1",9)))
SampleSheet_all$Subject <- c(paste0(SampleSheet_01$Sample_Well, "_90"), paste0(SampleSheet_02$Sample_Well,"_180-1"))
SampleSheet_all$Experiment <- c(rep("M90",9), rep("M180-1",9))

# Dimensionality reduction for data visualization
set.seed(123123)
tsne_out <- Rtsne(t(count_data_all_log),pca=TRUE,perplexity=2.5,theta=0.0)
tsnes <- tsne_out$Y
colnames(tsnes) <- c("tSNE1", "tSNE2")
tsnes <- as.data.frame(tsnes)
tsnes$Subject <- SampleSheet_all$Subject 
tsnes$Subject <- factor(tsnes$Subject, levels = c(
  "H01_180-1","H02_180-1","H03_180-1",
  "H01_90","H02_90","H03_90"
))
p1 <- ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col=Subject))+
  coord_cartesian(xlim = c(-150, 150)) + 
  theme_minimal()
tsnes$Time <- SampleSheet_all$Sample_Group
tsnes$Experiment=SampleSheet_all$Experiment
p2 <- ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col=Time, shape=Experiment)) +
  scale_color_discrete(name = "Time") +
  scale_shape_discrete(name = "Experiment") +
  coord_cartesian(xlim = c(-150, 150), ylim = c(-100, 120)) +  
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme_minimal()

tsne_plot_all <- ggpubr::ggarrange(p1, p2, nrow = 1)

ggsave(paste0(my_dir,"Extended Data Fig. 5a.pdf"), plot = tsne_plot_all, width = 7, height = 3)
write.csv(tsnes, file=paste0(my_dir,"Gene expression tsne data_Source data for Extended Data Fig. 5_a.csv"), row.names=F)

# Load count data of 180-day experiment (M180-2),count_data_symbol.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 180-day spaceflight_M180-2, RNA-seq)
my_dir <- "Data/RNA-seq/"
setwd(paste0(my_dir, "03/"))
count_data_03 <- read.table("180-day spaceflight_M180-2/RNA-seq/count_data_symbol.txt", sep="\t", as.is=T, header=T, row.names = 1)
SampleSheet_03 <- read.csv("180-day spaceflight_M180-2/RNA-seq/SampleSheet.csv")
sum(colnames(count_data_03)==SampleSheet_03$Sample_Name)

# Dimensionality reduction for data visualization
SampleSheet_03$Subject <- SampleSheet_03$Sample_Name
SampleSheet_03$Experiment <- ifelse(SampleSheet_03$Date %in% c("M113"),"M113","M180_2")
set.seed(123123)
tsne_out <- Rtsne(t(count_data_03),pca=TRUE,perplexity=2.5,theta=0.0)
tsnes <- tsne_out$Y
colnames(tsnes) <- c("tSNE1", "tSNE2")
tsnes <- as.data.frame(tsnes)
tsnes$Subject <- SampleSheet_03$Subject 
p1 <- ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col=Subject))
tsnes$time <- SampleSheet_03$Sample_Group
tsnes$Experiment=SampleSheet_03$Experiment
p2 <- ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col=time, shape=Experiment))
tsne_plot_all <- ggpubr::ggarrange(p1, p2, nrow = 1)
# 01End: Dimensionality reduction for data visualization before batch removal 


# 02Start: Removal of individual differences and batch effect by combat
# Removal of individual differences
combat_mat_01 <- sva::ComBat_seq(as.matrix(count_data_01), batch = SampleSheet_01$Sample_Well, group = SampleSheet_01$Sample_Group)
combat_mat_02 <- sva::ComBat_seq(as.matrix(count_data_02), batch = SampleSheet_02$Sample_Well, group = SampleSheet_02$Sample_Group)
combat_mat_03 <- sva::ComBat_seq(as.matrix(count_data_03), batch = SampleSheet_03$Sample_Well, group = SampleSheet_03$Sample_Group)
write.table(combat_mat_03, file=paste0(my_dir,"180-day spaceflight_M180-2/RNA-seq/nobatch_expression_profile.txt"),quote=F)

# Removal of batch effect 
genes <- intersect(rownames(combat_mat_01), rownames(combat_mat_02))
combat_mat_all <- cbind(combat_mat_01[genes,], combat_mat_02[genes,])
colnames(combat_mat_all) <- paste0(c(rep("90_",9), rep("180_",9)), colnames(combat_mat_all))
combat_mat_all <- sva::ComBat_seq(as.matrix(combat_mat_all), batch = c(rep("90",9), rep("180",9)))
write.table(combat_mat_all[,1:9], file=paste0(my_dir,"90-day spaceflight_M90/RNA-seq/nobatch_expression_profile.txt"),quote=F)
write.table(combat_mat_all[,10:18], file=paste0(my_dir,"180-day spaceflight_M180-1/RNA-seq/nobatch_expression_profile.txt"),quote=F)
# 02End: Removal of individual differences and batch effect by combat


# 03Start: Dimensionality reduction for data visualization after batch removal
combat_mat_all_log <- log2(combat_mat_all+1)
set.seed(123123)
tsne_out <- Rtsne(t(combat_mat_all_log),pca=TRUE,perplexity=2.5,theta=0.0)
str(tsne_out)
tsnes <- tsne_out$Y
colnames(tsnes) <- c("tSNE1", "tSNE2")
tsnes <- as.data.frame(tsnes)
tsnes$Subject <- SampleSheet_all$Subject 
p1 <- ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col=Subject))
tsnes$time <- SampleSheet_all$Sample_Group
tsnes$Experiment=SampleSheet_all$Experiment
p2 <- ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col=time, shape=Experiment))
noBatch_tsne_plot_all <- ggpubr::ggarrange(p1, p2, nrow = 1)

ggsave(paste0(my_dir,"Extended Data Fig. 5b.pdf"), plot = noBatch_tsne_plot_all, width = 7, height = 3)
write.csv(tsnes, file=paste0(my_dir,"Gene expression tsne data_Source data for Extended Data Fig. 5_b.csv"), row.names=F)

transcriptome_batch_nobatch <- ggpubr::ggarrange(tsne_plot_all, noBatch_tsne_plot_all, ncol = 1)
#03End: Dimensionality reduction for data visualization after batch removal
