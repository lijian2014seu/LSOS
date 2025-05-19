##### description: Division of genes into modules by WGCNA 

#load packages
library(WGCNA)
library(reshape2)
library(stringr)

# load the expression matrix, nobatch_expression_profile.txt are available for download at https://www.spacelifescience.cn/search/ (search for 90-day spaceflight_M90 / 180-day spaceflight_180-1, RNA-seq)

combat_mat_01 <- read.table("90-day spaceflight_M90/RNA-seq/nobatch_expression_profile.txt", header=T, row.names=1, as.is=T, sep="\t")
SampleSheet_01 <- read.csv("90-day spaceflight_M90/RNA-seq/SampleSheet.csv")

combat_mat_02 <- read.table("180-day spaceflight_M180-1/RNA-seq/nobatch_expression_profile.txt", header=T, row.names=1, as.is=T, sep="\t")
SampleSheet_02 <- read.csv("180-day spaceflight_M180-1/RNA-seq/SampleSheet.csv")

dataExpr <- cbind(combat_mat_02, combat_mat_01)
dataExpr_log <- log2(dataExpr+1)
dataExpr <- as.data.frame(t(dataExpr_log))
dataExpr  <- rbind(dataExpr[10:18,], dataExpr[1:9,])
SampleSheet_all <- rbind(SampleSheet_02, SampleSheet_01)
SampleSheet_all$Sample_Name <- paste0(SampleSheet_all$Sample_Name, c(rep("_180-1",9), rep("_90",9)))
SampleSheet_all$Subject <- c(paste0(SampleSheet_02$Sample_Well, "_180-1"), paste0(SampleSheet_01$Sample_Well,"_90"))
SampleSheet_all$Experiment <- c(rep("M180-1",9), rep("M90",9))

type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

powers = c(c(1:10), seq(from = 12, to=40, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers,
                        networkType=type, verbose=5)

par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")

# Soft threshold
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")

power = sft$powerEstimate
# choose 3
power = 3
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))      
                 )
  )
}


cor <- WGCNA::cor
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.1,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType,
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0("nobatch", ".tom"),
                       verbose = 3) 


sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Modulecolors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# module eigengene
save(net, file="WGCNA_PPI_analysis/result/net_nobatch_all.Rdata")

load("WGCNA_PPI_analysis/result/net_nobatch_all.Rdata")
MEs <- net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
write.table(MEs_col, "WGCNA_PPI_analysis/result/MEs_col.txt", sep="\t", quote=F)

plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T,
                      xLabelsAngle = 90)
# Aov test between module eigengene and time				  
  aov_fun <- function(x, group){
     df <- data.frame(value=unlist(x), time=group)
     aov_res <- aov(value~time, df)
     aov_res <- summary(aov_res)
     return(aov_res[[1]][1,])
  }

MEs_col <- read.table("WGCNA_PPI_analysis/result/MEs_col.txt",sep="\t",header=T,row.names=1, as.is=T)
aov_p_all <- do.call(rbind, apply(MEs_col, 2, aov_fun, group=SampleSheet_all$Sample_Group))
aov_p_01 <- do.call(rbind, apply(MEs_col[10:18,], 2, aov_fun, group=SampleSheet_all$Sample_Group[c(10:18)]))
aov_p_02 <- do.call(rbind, apply(MEs_col[1:9,], 2, aov_fun, group=SampleSheet_all$Sample_Group[c(1:9)]))
aov_p_df <- data.frame(pvalue_all=aov_p_all[,"Pr(>F)"], pvalue_M90 = aov_p_01[,"Pr(>F)"], pvalue_M180_1 = aov_p_02[,"Pr(>F)"], row_names = rownames(aov_p_all))
aov_p_df <- reshape::melt(aov_p_df, id.vars = 'row_names', variable_name = 'pvalue_type')
colnames(aov_p_df) <- c("Module", "pvalue_type", "P_value")
ggplot(aov_p_df, aes(x = pvalue_type, y = P_value, fill = Module)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, binwidth=0.05, position = position_dodge(width=0.75), alpha=0.4) +
  labs(title = "ANOVA-based p-value to assess the correlation between gene set modules and flight",
       x = "P-value Type",
       y = "P-value",
       fill = "Gene Set Modules") +
  scale_fill_manual(values = c("#00ADEC", "#7F3A0C", "#6FAB46", "grey", "#EC0000", "turquoise", "yellow")) +
  theme_minimal()+
  theme(axis.line=element_line(color="black", size = 0.5))
aov_res_df <- cbind(aov_p_all,aov_p_01, aov_p_02)
colnames(aov_res_df) <- c(paste0(colnames(aov_p_all),"_all"), paste0(colnames(aov_p_all),"_M90"), paste0(colnames(aov_p_all),"_M180-1"))
write.table(aov_res_df, "WGCNA_PPI_analysis/result/Module_Time_aov_P.txt", sep="\t", quote=F)

# Correlation matrix of modules and genes

if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
             as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}


if (corType=="pearson") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}




load("WGCNA_PPI_analysis/result/net_nobatch_all.Rdata")
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
exprMat <- "ME7_net"
dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong 
power <- 3
plotTOM = dissTOM^power
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
probes = names(net$colors)
dimnames(TOM) <- list(probes, probes)

# Export the network into edge and node list files Cytoscape can read
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste(exprMat, ".edges.txt", sep=""),
                               nodeFile = paste(exprMat, ".nodes.txt", sep=""),
                               weighted = TRUE, threshold = 0.01,
                               nodeNames = probes, nodeAttr = moduleColors)

edges <- read.table("WGCNA_PPI_analysis/result/ME7_net.edges.txt", header = T, sep = "\t", as.is = T)
nodes <- read.table("WGCNA_PPI_analysis/result/ME7_net.nodes.txt", header = T, sep = "\t", as.is = T)
nodes_brown <- nodes[nodes[,3]=="brown",1]
write.table(nodes_brown, file = "WGCNA_PPI_analysis/result/nodes_brown.txt", quote = F, row.names = F, col.names = F)
edges_brown <- edges[(edges$fromNode %in% nodes_brown) & (edges$toNode %in% nodes_brown), ]
save(nodes_brown, edges_brown, file = "WGCNA_PPI_analysis/result/met_module_brown.RData")

nodes_blue <- nodes[nodes[,3]=="blue",1]
write.table(nodes_blue, file = "WGCNA_PPI_analysis/result/nodes_blue.txt", quote = F, row.names = F, col.names = F)
edges_blue <- edges[(edges$fromNode %in% nodes_blue) & (edges$toNode %in% nodes_blue), ]
save(nodes_blue, edges_blue, file = "WGCNA_PPI_analysis/result/met_module_blue.RData")

MEbrown_genes <- nodes_brown
MEbrown_genes_filter <- MEbrown_genes[abs(geneModuleMembership[MEbrown_genes,"MEbrown"])>0.5&MMPvalue[MEbrown_genes,"MEbrown"]<0.05]
write.table(as.data.frame(MEbrown_genes_filter), file="WGCNA_PPI_analysis/result/MEbrown_genes_filter.txt", row.names = F, quote = F)
MEblue_genes <- nodes_blue
MEblue_genes_filter <- MEblue_genes[abs(geneModuleMembership[MEblue_genes,"MEblue"])>0.5&MMPvalue[MEblue_genes,"MEblue"]<0.05]
write.table(as.data.frame(MEblue_genes_filter), file="WGCNA_PPI_analysis/result/MEblue_genes_filter.txt", row.names = F, quote = F)


#enrichment analysis for module genes
source("enrich_function.R")
library(simplifyEnrichment)
MEbrown_genes <- read.table("WGCNA_PPI_analysis/result/nodes_brown.txt")
go_enrich_res_brown <- go_enrich(gene_symbols=unlist(MEbrown_genes), gene_type="MEbrown", outdir="D:/PhD/major_program/analysis/RNA-seq/WGCNA/sorted_data/WGCNA")
write.table(go_enrich_res_brown, file = "WGCNA_PPI_analysis/result/go_enrich_res_MEbrown.txt", quote = F, sep="\t")

set.seed(888)
mat = GO_similarity(go_enrich_res_brown$ID)
df = simplifyGO(mat)

MEblue_genes <- read.table("WGCNA_PPI_analysis/result/nodes_blue.txt")
go_enrich_res_blue <- go_enrich(gene_symbols=unlist(MEblue_genes), gene_type="MEblue", outdir="D:/PhD/major_program/analysis/RNA-seq/WGCNA/sorted_data/WGCNA")
write.table(go_enrich_res_blue, file = "WGCNA_PPI_analysis/result/go_enrich_res_MEblue.txt", quote = F, sep="\t")

set.seed(888)
mat = GO_similarity(go_enrich_res_blue$ID)
df = simplifyGO(mat)



