# Description: T1->T2, T2->T3 DEG number statistic (Fig. 2d) and expression (Fig.4a and Extended Data Fig. 6a ) 

# Load packages
library(UpSetR)

# 01Start: DEG number statistic and bar plot
DEG_list_01 <- readRDS(file = "90-day spaceflight_M90/RNA-seq/DE_analysis/DEG_list.Rds") 
DEG_num_01 <- mapply(length, DEG_list_01)[-c(3,6)]
DEG_list_02 <- readRDS(file = "180-day spaceflight_M180-1/RNA-seq/DE_analysis/DEG_list.Rds") 
DEG_num_02 <- mapply(length, DEG_list_02 )[-c(3,6)]
DEG_list_03 <- readRDS(file = "180-day spaceflight_M180-2/RNA-seq/DE_analysis/DEG_list.Rds") 
DEG_num_03 <- mapply(length, DEG_list_03 )[-c(3,6)]

cor.test(DEG_num_01, DEG_num_02)
DEG_num_02/DEG_num_01
DEG_num_df <- data.frame(Expr_01=DEG_num_01, Expr_02=DEG_num_02)
DEG_num_df$ratio <- DEG_num_02/DEG_num_01
DEG_num_df$geneset <- factor(rownames(DEG_num_df), levels=rownames(DEG_num_df))
DEG_num_df <- reshape::melt(DEG_num_df, id.vars=c("geneset","ratio"))
colnames(DEG_num_df)[c(3,4)] <- c("expriment","number")
p <- ggplot(DEG_num_df, aes(x=geneset, y=number, fill=expriment))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c(Expr_01 = "#FFFFCC", Expr_02 = "#99CCFF"))
annotation <- data.frame(
  x = c(1:4),
  y = 2100,
  label = round(DEG_num_02/DEG_num_01,3)
)
annotation$expriment <- DEG_num_df$expriment[1:4]
bar_plot <- p + geom_label(data=annotation, aes( x=x, y=y, label=label),                
               color="blue", size=5 , angle= 60, fontface="bold" ) +
    geom_text(aes(label = number), position = position_dodge(0.9), vjust=0.1 )

# Number of intersections between DEGs of 90-day and 180-day
DEG_list_all <- c(DEG_list_01[c(1:2,4:5)],DEG_list_02[c(1:2,4:5)],DEG_list_03[c(1:2,4:5)])
names(DEG_list_all) <- paste(names(DEG_list_all), rep(c("M90","M180-1","M180-2"),each=4), sep="_")
DEG_matrix_all <- matrix(0,length(unique(unlist(DEG_list_all))), length(DEG_list_all))
rownames(DEG_matrix_all) <- unique(unlist(DEG_list_all))
colnames(DEG_matrix_all) <- names(DEG_list_all)
for(i in 1:length(DEG_list_all)){
   DEG_matrix_all[DEG_list_all[[i]],i] <- 1
}
DEG_df_all <- as.data.frame(DEG_matrix_all)
DEG_df_all$gene_name <- rownames(DEG_matrix_all)

intersections <- as.list(as.data.frame(combn(names(DEG_list_all),2)))
intersections_num <- lapply(intersections, function(x){length(intersect(DEG_list_all[[x[1]]],DEG_list_all[[x[2]]]))})
intersections <- intersections[intersections_num>0]

DEG_df_all_T2vsT1 <- DEG_df_all[,grep("T1_to_T2",colnames(DEG_df_all))]
colnames(DEG_df_all_T2vsT1)<- gsub("T1_to_T2","T2_vs_T1", colnames(DEG_df_all_T2vsT1))
p_T2vsT1 <- upset(DEG_df_all_T2vsT1, nset = nrow(DEG_df_all_T2vsT1), nintersects = 1000, sets= sort(colnames(DEG_df_all_T2vsT1)), decreasing = c(FALSE, TRUE), keep.order = T, main.bar.color = "#5F86A8", matrix.color = "#246D9C", sets.bar.color="#A49FC3")

DEG_df_all_T3vsT2 <- DEG_df_all[,grep("T2_to_T3",colnames(DEG_df_all))]
colnames(DEG_df_all_T3vsT2)<- gsub("T2_to_T3","T3_vs_T2", colnames(DEG_df_all_T3vsT2))
p_T3vsT2 <- upset(DEG_df_all_T3vsT2, nset = nrow(DEG_df_all_T3vsT2), nintersects = 1000, sets= sort(colnames(DEG_df_all_T3vsT2)), decreasing = c(FALSE, TRUE),keep.order = T, main.bar.color = "#5F86A8", matrix.color = "#246D9C", sets.bar.color="#A49FC3")

pdf("Fig.2d DEGs in T2_vs_T1.pdf",width=8,height=6)	  
 upset(DEG_df_all[,c(1:2,5:6)], nset = 4, nintersects = 2000, decreasing = c(TRUE, TRUE), 
    keep.order = T, main.bar.color = "#5F86A8", matrix.color = "#246D9C", sets.bar.color="#A49FC3")
dev.off()

pdf("Fig.2d DEGs in T3_vs_T2.pdf",width=8,height=6)	  
 upset(DEG_df_all[,c(3,4,7,8)], nset = 4, nintersects = 2000, decreasing = c(TRUE, TRUE), 
    keep.order = T, main.bar.color = "#5F86A8", matrix.color = "#246D9C", sets.bar.color="#A49FC3")
dev.off()
# 01End: DEG number statistic and bar plot

# 02Start: DEG heatmap plot
# For M90
DEG_list_trends <- readRDS("90-day spaceflight_M90/RNA-seq/DE_analysis/DEG_list_trends.Rds")
SampleSheet_01 <- read.csv("90-day spaceflight_M90/RNA-seq/SampleSheet.csv")
combat_mat_01 <- read.table(paste0(my_dir,"90-day spaceflight_M90/RNA-seqnobatch_expression_profile.txt"), header=T, row.names=1, as.is=T, sep=" ")

sample_order <- order(SampleSheet_01$Sample_Group, SampleSheet_01$Sample_Well)
DEG_exp_count <- combat_mat_01[unlist(DEG_list_trends), sample_order]	
rownames(DEG_exp_count) <- make.unique(rownames(DEG_exp_count))
annotation_row <- data.frame(type = unlist(lapply(1:length(DEG_list_trends), function(i) rep(names(DEG_list_trends)[i], length(DEG_list_trends[[i]])))))
rownames(annotation_row) <- rownames(DEG_exp_count)

annotation_col <- data.frame(time=SampleSheet_01$Sample_Group, person=SampleSheet_01$Sample_Well)
rownames(annotation_col) <- paste0("X90_", SampleSheet_01$Sample_Name)

gaps_row <- cumsum(mapply(length, DEG_list_trends))
my_palette <- c(colorRampPalette(c("#536BE1", "white"))(1000),colorRampPalette(c("white", "#FE4C5B"))(1000))
expression_plot_heatmap <- pheatmap::pheatmap(DEG_exp_count, cluster_rows = F, cluster_cols = F, show_rownames = F, color=my_palette,
          scale = "row", annotation_col = annotation_col, annotation_row=annotation_row, gaps_row=gaps_row)	

ggsave("Fig.4a M90.pdf", plot = expression_plot_heatmap, width = 6, height = 4)
DEG_exp_count <- as.data.frame(DEG_exp_count)
DEG_exp_count$change_type <- annotation_row$type
write.table(DEG_exp_count, file="90-day spaceflight_M90/RNA-seq/DE_analysis/DEG_expression.txt", quote=F, sep="\t")


# For M180-1
DEG_list_trends <- readRDS("180-day spaceflight_M180-1/RNA-seq/DE_analysis/DEG_list_trends.Rds")
SampleSheet_02 <- read.csv("180-day spaceflight_M180-1/RNA-seq/SampleSheet.csv")
combat_mat_02 <- read.table("180-day spaceflight_M180-1/RNA-seq/nobatch_expression_profile.txt", header=T, row.names=1, as.is=T, sep=" ")

sample_order <- order(SampleSheet_02$Sample_Group, SampleSheet_02$Sample_Well)
DEG_exp_count <- combat_mat_02[unlist(DEG_list_trends), sample_order]	
rownames(DEG_exp_count) <- make.unique(rownames(DEG_exp_count))
annotation_row <- data.frame(type = unlist(lapply(1:length(DEG_list_trends), function(i) rep(names(DEG_list_trends)[i], length(DEG_list_trends[[i]])))))
rownames(annotation_row) <- rownames(DEG_exp_count)

annotation_col <- data.frame(time=SampleSheet_02$Sample_Group, person=SampleSheet_02$Sample_Well)
rownames(annotation_col) <- paste0("X180_", SampleSheet_02$Sample_Name)

gaps_row <- cumsum(mapply(length, DEG_list_trends))
my_palette <- c(colorRampPalette(c("#536BE1", "white"))(1000),colorRampPalette(c("white", "#FE4C5B"))(1000))
expression_plot_heatmap <- pheatmap::pheatmap(DEG_exp_count, cluster_rows = F, cluster_cols = F, show_rownames = F, color=my_palette,
          scale = "row", annotation_col = annotation_col, annotation_row=annotation_row, gaps_row=gaps_row)	
ggsave("Fig.4a M180-1.pdf", plot = expression_plot_heatmap, width = 4, height = 4)
DEG_exp_count <- as.data.frame(DEG_exp_count)
DEG_exp_count$change_type <- annotation_row$type
write.table(DEG_exp_count, file="180-day spaceflight_M180-1/RNA-seq/DE_analysis/DEG_expression.txt", quote=F, sep="\t")


# For M180-2
DEG_list_trends <- readRDS("180-day spaceflight_M180-2/RNA-seq/DE_analysis/DEG_list_trends.Rds")
SampleSheet_03 <- read.csv("180-day spaceflight_M180-2/RNA-seq/SampleSheet.csv")
combat_mat_03 <- read.table(paste0(my_dir,"180-day spaceflight_M180-2/RNA-seq/nobatch_expression_profile.txt"), header=T, row.names=1, as.is=T, sep=" ")
SampleSheet_03 <- SampleSheet_03[SampleSheet_03$Sample_Group != "Torbit", ]
combat_mat_03 <- combat_mat_03[,SampleSheet_03$Sample_Name]

sample_order <- order(SampleSheet_03$Sample_Group, SampleSheet_03$Sample_Well)
DEG_exp_count <- combat_mat_03[unlist(DEG_list_trends), sample_order]	
rownames(DEG_exp_count) <- make.unique(rownames(DEG_exp_count))
annotation_row <- data.frame(type = unlist(lapply(1:length(DEG_list_trends), function(i) rep(names(DEG_list_trends)[i], length(DEG_list_trends[[i]])))))
rownames(annotation_row) <- rownames(DEG_exp_count)

annotation_col <- data.frame(time=SampleSheet_03$Sample_Group, person=SampleSheet_03$Sample_Well)
rownames(annotation_col) <-  SampleSheet_03$Sample_Name

gaps_row <- cumsum(mapply(length, DEG_list_trends))
my_palette <- c(colorRampPalette(c("#253e8f", "white"))(1000),colorRampPalette(c("white", "#c73e47"))(1000))
expression_plot_heatmap <- pheatmap::pheatmap(DEG_exp_count, cluster_rows = F, cluster_cols = F, show_rownames = F, color=my_palette,
          scale = "row", annotation_col = annotation_col, annotation_row=annotation_row, gaps_row=gaps_row)	
ggsave("Extended Data Fig. 6.pdf", plot = expression_plot_heatmap, width = 4, height = 4)
DEG_exp_count <- as.data.frame(DEG_exp_count)
DEG_exp_count$change_type <- annotation_row$type
write.table(DEG_exp_count, file="180-day spaceflight_M180-2/RNA-seq/DE_analysis/DEG_expression.txt", quote=F, sep="\t")
