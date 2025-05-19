my_dir <- "/home/momoki/HTmicro/mrna/microrna"
setwd(my_dir)
SampleSheet <- read.csv("Samplesheet.csv")
mature_miRNA_counts <- read.csv("mature_counts.csv", row.names=1)
mature_miRNA_counts_log <- log2(mature_miRNA_counts+1)

library(Rtsne)
library(ggplot2)

# Visualize raw data using t-SNE
set.seed(123123)
tsne_out <- Rtsne(mature_miRNA_counts_log,pca=FALSE,perplexity=2.5,theta=0.0)
str(tsne_out)
tsnes <- tsne_out$Y
colnames(tsnes) <- c("tSNE1", "tSNE2") 
tsnes <- as.data.frame(tsnes)
tsnes$person <- SampleSheet$Sample_Well
p1 <- ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col=person))
tsnes$time <- SampleSheet$Sample_Group
p2 <- ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col=time))
tsne_plot <- ggpubr::ggarrange(p1, p2, nrow = 1)
ggsave("tSNE_plot_pre.PDF", tsne_plot, width = 10, height = 5, dpi = 300)


# Remove control samples & correct for batch effects
SampleSheet <- SampleSheet[SampleSheet$Sample_Group!="control",]
mature_miRNA_counts <- mature_miRNA_counts[SampleSheet$Sample_Name,]
table(SampleSheet$Sample_Name==colnames(mature_miRNA_counts))
mature_miRNA_counts <- mature_miRNA_counts[,colSums(mature_miRNA_counts) > 0]
dim(mature_miRNA_counts)
mature_miRNA_counts <- t(mature_miRNA_counts)
combat_mat_01 <- sva::ComBat_seq(as.matrix(mature_miRNA_counts), batch = SampleSheet$Sample_Well, group = SampleSheet$Sample_Group)
write.table(combat_mat_01, file=paste0("nobatch_expression_profile_micro.txt"),quote=F)

#  Visualize data after batch correction
combat_mat_all_log <- log2(combat_mat_01+1)
set.seed(123123)
tsne_out <- Rtsne(t(combat_mat_all_log),pca=TRUE,perplexity=2.5,theta=0.0)
str(tsne_out)
tsnes <- tsne_out$Y
colnames(tsnes) <- c("tSNE1", "tSNE2")
tsnes <- as.data.frame(tsnes)
tsnes$Subject <- SampleSheet$Sample_Name 
p1 <- ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col=Subject))
tsnes$time <- SampleSheet$Sample_Group
tsnes$Experiment=SampleSheet$Sample_Group
p2 <- ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col=time, shape=Experiment))
noBatch_tsne_plot_all <- ggpubr::ggarrange(p1, p2, nrow = 1)
ggsave(paste0("tSNE_plot_nbatch.pdf"), plot = noBatch_tsne_plot_all, width = 7, height = 3)


#  Differential expression analysis
SampleSheet_case <- SampleSheet[!(SampleSheet$Sample_Group %in% c("control", "T3")), ]
combat_mat_01 <- combat_mat_01[,SampleSheet_case$Sample_Name]
table(SampleSheet_case$Sample_Name==colnames(combat_mat_01))

library(edgeR) 
get_DEG(exp_mat=(combat_mat_01), sample_df=SampleSheet_case,cut_off_pvalue = 0.05,
        cut_off_logFC = log2(1.5), ourdir = "/home/momoki/HTmicro/mrna/microrna/result")

get_DEG <- function(exp_mat, sample_df, cut_off_pvalue,
                    cut_off_logFC,	ourdir){
  cpm <- cpm(exp_mat)
  lcpm <- cpm(exp_mat, log=TRUE, prior.count=2)
  time_group <- sample_df$Sample_Group
  person_group <- sample_df$Sample_Well
  keep.exprs <- filterByExpr(exp_mat, group=time_group)
  exp_mat <- exp_mat[keep.exprs,]
  d0 <- DGEList(exp_mat)
  d0 <- calcNormFactors(d0, method = "TMM")
  
  par(mfrow=c(1,2))
  plotMDS(lcpm, labels=time_group)
  title(main="A. Sequencing times")
  plotMDS(lcpm, labels=person_group)
  title(main="B. Sequencing persons")
  
  time_design <- model.matrix(~0+time_group)
  colnames(time_design) <- gsub("time_group", "", colnames(time_design))
  
  contr.matrix <- makeContrasts(
    TorbitvsT1 = Torbit - T1,
    T2vsTorbit = T2 - Torbit,
    T2vsT1 = T2 - T1,
    levels = colnames(time_design))
  
  par(mfrow=c(1,3))
  v <- voom(d0, time_design, plot=TRUE)
  corfit <- duplicateCorrelation(v, block= person_group)
  v2 <- voom(d0, time_design, correlation=corfit$consensus.correlation, 
             block=person_group, plot=TRUE)
  vfit <- lmFit(v2, time_design, correlation=corfit$consensus.correlation, 
                block=person_group)
  vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- eBayes(vfit)
  plotSA(efit, main="Final model: Mean-variance trend")

  print(summary(decideTests(efit)))  
  
  DEG_T1_to_Torbit <- topTable(efit, 
                           coef = "TorbitvsT1",
                           n = Inf,
  )
  DEG_Torbit_to_T2 <- topTable(efit, 
                           coef = "T2vsTorbit",
                           n = Inf,
  )
  DEG_T1_to_T2 <- topTable(efit, 
                           coef = "T2vsT1",
                           n = Inf,
  )				
  
  
  setwd(ourdir)
  DEG_T1_to_Torbit$change = ifelse(DEG_T1_to_Torbit$adj.P.Val < cut_off_pvalue & abs(DEG_T1_to_Torbit$logFC) >= cut_off_logFC, 
                               ifelse(DEG_T1_to_Torbit$logFC> cut_off_logFC ,'Up','Down'),
                               'Stable')
  table(DEG_T1_to_Torbit$change)
  write.table(DEG_T1_to_Torbit, file = file.path(ourdir, "T1_to_Torbit_DEG_table.txt"), sep = "\t", quote = F)

  DEG_Torbit_to_T2$change = ifelse(DEG_Torbit_to_T2$adj.P.Val < cut_off_pvalue & abs(DEG_Torbit_to_T2$logFC) >= cut_off_logFC, 
                               ifelse(DEG_Torbit_to_T2$logFC> cut_off_logFC ,'Up','Down'),
                               'Stable')
  table(DEG_Torbit_to_T2$change)
  write.table(DEG_Torbit_to_T2, file = file.path(ourdir, "Torbit_to_T2_DEG_table.txt"), sep = "\t", quote = F)
  
  DEG_T1_to_T2$change = ifelse(DEG_T1_to_T2$adj.P.Val < cut_off_pvalue & abs(DEG_T1_to_T2$logFC) >= cut_off_logFC, 
                               ifelse(DEG_T1_to_T2$logFC> cut_off_logFC ,'Up','Down'),
                               'Stable')
  table(DEG_T1_to_T2$change)
  write.table(DEG_T1_to_T2, file = file.path(ourdir, "T1_to_T2_DEG_table.txt"), sep = "\t", quote = F)
  
  DEG_list <- list(up_T1_to_Torbit = rownames(DEG_T1_to_Torbit)[DEG_T1_to_Torbit$change=="Up"],
                   down_T1_to_Torbit = rownames(DEG_T1_to_Torbit)[DEG_T1_to_Torbit$change=="Down"],
                   stable_T1_to_Torbit = rownames(DEG_T1_to_Torbit)[DEG_T1_to_Torbit$change=="Stable"],
                   up_Torbit_to_T2 = rownames(DEG_Torbit_to_T2)[DEG_Torbit_to_T2$change=="Up"],
                   down_Torbit_to_T2 = rownames(DEG_Torbit_to_T2)[DEG_Torbit_to_T2$change=="Down"],
                   stable_Torbit_to_T2 = rownames(DEG_Torbit_to_T2)[DEG_Torbit_to_T2$change=="Stable"])
  saveRDS(DEG_list, file = file.path(ourdir, "DEG_list_03.Rds"))
}
