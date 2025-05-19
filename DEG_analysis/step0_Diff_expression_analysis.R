# Description: Differential expression analysis

setwd("/home/lqwang/Program")
# Load packages
library(edgeR) 

# Start: Differential expression analysis 
# Load count data of 90-day experiment (M90),nobatch_expression_profile.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 90-day spaceflight_M90, RNA-seq)
combat_mat_01 <- read.table("90-day_spaceflight_M90/RNA-seq/nobatch_expression_profile.txt", header=T, row.names=1, as.is=T, sep="\t")
SampleSheet_01 <- read.csv("90-day_spaceflight_M90/RNA-seq/SampleSheet.csv")
get_DEG(exp_mat=combat_mat_01, sample_df=SampleSheet_01, cut_off_pvalue = 0.05,
        cut_off_logFC = log2(1.5), ourdir = "90-day_spaceflight_M90/RNA-seq/DE_analysis/")

# Load count data of 180-day experiment (M180-1),nobatch_expression_profile.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 180-day spaceflight_180-1, RNA-seq)
combat_mat_02 <- read.table("180-day_spaceflight_M180-1/RNA-seq/nobatch_expression_profile.txt", header=T, row.names=1, as.is=T, sep="\t")
SampleSheet_02 <- read.csv("180-day_spaceflight_M180-1/RNA-seq/SampleSheet.csv")
get_DEG(exp_mat=combat_mat_02, sample_df=SampleSheet_02,cut_off_pvalue = 0.05,
        cut_off_logFC = log2(1.5), ourdir = "180-day_spaceflight_M180-1/RNA-seq/DE_analysis/")

# Load count data of 180-day experiment (M180-2),nobatch_expression_profile.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 180-day spaceflight_M180-2, RNA-seq)
combat_mat_03 <- read.table("180-day_spaceflight_M180-2/RNA-seq/nobatch_expression_profile.txt", header=T, row.names=1, as.is=T, sep="\t")
SampleSheet_03 <- read.csv("180-day_spaceflight_M180-2/RNA-seq/SampleSheet.csv")
SampleSheet_03 <- SampleSheet_03[SampleSheet_03$Sample_Group != "Torbit", ]
combat_mat_03 <- combat_mat_03[,SampleSheet_03$Sample_Name]
get_DEG(exp_mat=combat_mat_03, sample_df=SampleSheet_03,cut_off_pvalue = 0.05,
        cut_off_logFC = log2(1.5), ourdir = "180-day_spaceflight_M180-2/RNA-seq/DE_analysis/")

# Function for differential expression analysis 
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
    T2vsT1 = T2 - T1,
    T3vsT2 = T3 - T2,
    T3vsT1 = T3 - T1,
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
  
  DEG_T1_to_T2 <- topTable(efit, 
                           coef = "T2vsT1",
                           n = Inf,
  )
  DEG_T2_to_T3 <- topTable(efit, 
                           coef = "T3vsT2",
                           n = Inf,
  )
  DEG_T1_to_T3 <- topTable(efit, 
                           coef = "T3vsT1",
                           n = Inf,
  )				
  
  DEG_T1_to_T2$change = ifelse(DEG_T1_to_T2$adj.P.Val < cut_off_pvalue & abs(DEG_T1_to_T2$logFC) >= cut_off_logFC, 
                               ifelse(DEG_T1_to_T2$logFC> cut_off_logFC ,'Up','Down'),
                               'Stable')
  table(DEG_T1_to_T2$change)
  write.table(DEG_T1_to_T2, file = "./T1_to_T2_DEG_table.txt", sep = "\t", quote = F)
  
  DEG_T2_to_T3$change = ifelse(DEG_T2_to_T3$adj.P.Val < cut_off_pvalue & abs(DEG_T2_to_T3$logFC) >= cut_off_logFC, 
                               ifelse(DEG_T2_to_T3$logFC> cut_off_logFC ,'Up','Down'),
                               'Stable')
  table(DEG_T2_to_T3$change)
  write.table(DEG_T2_to_T3, file = "./T2_to_T3_DEG_table.txt", sep = "\t", quote = F)
  
  DEG_T1_to_T3$change = ifelse(DEG_T1_to_T3$adj.P.Val < cut_off_pvalue & abs(DEG_T1_to_T3$logFC) >= cut_off_logFC, 
                               ifelse(DEG_T1_to_T3$logFC> cut_off_logFC ,'Up','Down'),
                               'Stable')
  table(DEG_T1_to_T3$change)
  write.table(DEG_T1_to_T3, file = "./T1_to_T3_DEG_table.txt", sep = "\t", quote = F)
  
  DEG_list <- list(`up_T1_to_T2` = rownames(DEG_T1_to_T2)[DEG_T1_to_T2$change=="Up"],
                   `down_T1_to_T2` = rownames(DEG_T1_to_T2)[DEG_T1_to_T2$change=="Down"],
                   `stable_T1_to_T2` = rownames(DEG_T1_to_T2)[DEG_T1_to_T2$change=="Stable"],
                   `up_T2_to_T3` = rownames(DEG_T2_to_T3)[DEG_T2_to_T3$change=="Up"],
                   `down_T2_to_T3` = rownames(DEG_T2_to_T3)[DEG_T2_to_T3$change=="Down"],
                   `stable_T2_to_T3` = rownames(DEG_T2_to_T3)[DEG_T2_to_T3$change=="Stable"])
  saveRDS(DEG_list, file = file.path(ourdir, "DEG_list.Rds"))
}
# T1_to_T2_DEG_table.txt, T1_to_T2_DEG_table.txt, T1_to_T3_DEG_table.txt and DEG_list.Rds are available in "DE_analysis"


# DEG trend in T1_to_T2 and T2_to_T3 (M90)
DEG_list_01 <- readRDS(file = "90-day_spaceflight_M90/RNA-seq/DE_analysis/DEG_list.Rds")

state <- c("up", "stable", "down")
state_groups <- unique(t(combn(c(1:3,1:3),2)))
trends <- t(apply(state_groups, 1, function(x){return(c(state[x]))}))
trends <- as.data.frame(trends)
colnames(trends) <- c("T1_to_T2","T2_to_T3")
DEG_list_trends <- apply(trends, 1, function(x){
  group1 <- paste0(x[1], "_", colnames(trends)[1])
  group2 <- paste0(x[2], "_", colnames(trends)[2])
  trend_genes <- intersect(DEG_list_01[[group1]], DEG_list_01[[group2]])
  print(length(trend_genes))
  return(trend_genes)
})
names(DEG_list_trends) <- paste0(trends$T1_to_T2, "_", trends$T2_to_T3)
mapply(length, DEG_list_trends )
DEG_list_trends <- DEG_list_trends[-grep("stable_stable", names(DEG_list_trends))]
DEG_list_trends <- DEG_list_trends[mapply(length, DEG_list_trends )!=0]
saveRDS(DEG_list_trends, file = "90-day spaceflight_M90/RNA-seq/DE_analysis/DEG_list_trends.Rds")

# DEG trend in T1_to_T2 and T2_to_T3 (M180-1)
DEG_list_02 <- readRDS(file = "180-day_spaceflight_M180-1/RNA-seq/DE_analysis/DEG_list.Rds")
trends <- t(apply(state_groups, 1, function(x){return(c(state[x]))}))
trends <- as.data.frame(trends)
colnames(trends) <- c("T1_to_T2","T2_to_T3")
DEG_list_trends <- apply(trends, 1, function(x){
  group1 <- paste0(x[1], "_", colnames(trends)[1])
  group2 <- paste0(x[2], "_", colnames(trends)[2])
  trend_genes <- intersect(DEG_list_02[[group1]], DEG_list_02[[group2]])
  print(length(trend_genes))
  return(trend_genes)
})
names(DEG_list_trends) <- paste0(trends$T1_to_T2, "_", trends$T2_to_T3)
mapply(length, DEG_list_trends )
DEG_list_trends <- DEG_list_trends[-grep("stable_stable", names(DEG_list_trends))]
DEG_list_trends <- DEG_list_trends[mapply(length, DEG_list_trends )!=0]
saveRDS(DEG_list_trends, file = "180-day_spaceflight_M180-1/RNA-seq/DE_analysis/DEG_list_trends.Rds")

# DEG trend in T1_to_T2 and T2_to_T3 (M180-2)
DEG_list_03 <- readRDS(file = "180-day_spaceflight_M180-2/RNA-seq/DE_analysis/DEG_list.Rds")
trends <- t(apply(state_groups, 1, function(x){return(c(state[x]))}))
trends <- as.data.frame(trends)
colnames(trends) <- c("T1_to_T2","T2_to_T3")
DEG_list_trends <- apply(trends, 1, function(x){
  group1 <- paste0(x[1], "_", colnames(trends)[1])
  group2 <- paste0(x[2], "_", colnames(trends)[2])
  trend_genes <- intersect(DEG_list_03[[group1]], DEG_list_03[[group2]])
  print(length(trend_genes))
  return(trend_genes)
})
names(DEG_list_trends) <- paste0(trends$T1_to_T2, "_", trends$T2_to_T3)
mapply(length, DEG_list_trends )
DEG_list_trends <- DEG_list_trends[-grep("stable_stable", names(DEG_list_trends))]
DEG_list_trends <- DEG_list_trends[mapply(length, DEG_list_trends )!=0]
saveRDS(DEG_list_trends, file = "180-day_spaceflight_M180-2/RNA-seq/DE_analysis/DEG_list_trends.Rds")
