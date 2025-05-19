# Description: Statistics on the number of DEGs (Fig. 2c and Extended Data Fig. 6b) and presentation of the changing trend (Fig. 4c and Extended Data Fig. 6c)
# DEG_list_trends.Rds, T1_to_T2_DEG_table.txt, T2_to_T3_DEG_table.txt, T1_to_T3_DEG_table.txt and DEG_list.Rds are available in "DEG_analysis"

setwd("/home/lqwang/Program")

# Load packages
library(foreign)
library(ggplot2)
library(ggalluvial)

# 01Start: Statistics on the number of DEGs 
# Example: For M90
workdir <- "90-day_spaceflight_M90/RNA-seq/DE_analysis/"
#Results for different experiments were obtained by changing the workdir and repeating the following code
setwd(workdir)
DEG_list_trend <- readRDS("DEG_list_trends.Rds")
DEG_num <- mapply(length, DEG_list_trend)
DEG_num_df <- data.frame(trend=names(DEG_num), number=DEG_num)
DEG_num_df$T1_to_T2 <- stringr::str_split_i(DEG_num_df$trend, "_", 1)
DEG_num_df$T2_to_T3 <- stringr::str_split_i(DEG_num_df$trend, "_", 2)
write.table(DEG_num_df, file = "DEG_number_statistic.txt", sep = "\t", quote = F, row.names = F)
DEG_num_df <- DEG_num_df[DEG_num_df$T1_to_T2!="stable",]
DEG_num_df$T1_to_T2<-factor(DEG_num_df$T1_to_T2,
                       levels = c("up","stable","down"))
DEG_num_df$T2_to_T3<-factor(DEG_num_df$T2_to_T3,
                               levels = c("up","stable","down"))						   
print(DEG_num_df)

# Creating an alluvium map					   
alluvium_map <- ggplot(data=DEG_num_df,
       aes(axis1=T1_to_T2,axis2=T2_to_T3,
           y=number))+
  geom_alluvium(aes(fill=T1_to_T2),
                #size=3,
                #color="white",
                width = 0.1,
                aes.bind = "flows")+
  geom_stratum(fill=c("#73B9CE", "#C1548E","#73B9CE", "darkgray", "#C1548E"),
               #color="white",
               #size=3,
               width=0.1)+
  scale_fill_manual(breaks = c("up","stable","down"),
                    values = c("#C1548E", "darkgray", "#73B9CE"),
                    labels = c("up","stable","down")) +
  scale_color_manual(breaks = c("up","stable","down"),
                     values = c("#C1548E", "darkgray", "#73B9CE"),
                     labels = c("up","stable","down")) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)))+
  scale_x_continuous(breaks = c(1,2),
                     labels = c("T1_to_T2","T2_to_T3"),
                     expand = expansion(mult = c(0,0)))+
  coord_cartesian(clip="off") +
  theme_classic()
ggsave("Fig. 2c.pdf", plot = alluvium_map, width = 6, height = 4)
# 01End: Statistics on the number of DEGs 


# 02Start: Changing trend of DEGs
# Load packages
library(reshape2) 
library(ggplot2)  
library(RColorBrewer)

T1_to_T2_DEG_table <- read.delim("T1_to_T2_DEG_table.txt",row.names=1)
T2_to_T3_DEG_table <- read.delim("T2_to_T3_DEG_table.txt",row.names=1)
T1_to_T3_DEG_table <- read.delim("T1_to_T3_DEG_table.txt",row.names=1)
genes <- rownames(T1_to_T2_DEG_table)
logFC_df <- data.frame(gene_name=genes, T1_value=0, T2_value=T1_to_T2_DEG_table$logFC,
                          T3_value=T1_to_T3_DEG_table[genes,"logFC"], change_type=T1_to_T2_DEG_table$change)
logFC_df <- logFC_df[logFC_df$change_type!="Stable",]
write.csv(logFC_df, file = "Changing trend of DEGs M90.csv", row.names=F)
# Changing trend of DEGs M180-1.csv
# Changing trend of DEGs M180-2.csv

# Displaying logFCs by time
logFC_df_up <- logFC_df[logFC_df$change_type=="Up",2:4]
my_col <- c("#F0D8E5", "#C1548E")
up_plot <- trend_plot(logFC_df_up, my_col, title = "Up-regulated genes in T2 vs T1")
logFC_df_down <- logFC_df[logFC_df$change_type=="Down",2:4]
my_col <- c("#C6E4F5", "#6BACD1")
down_plot <- trend_plot(logFC_df_down, my_col, title = "Down-regulated genes in T2 vs T1")
up_down_trend_plot <- ggpubr::ggarrange(up_plot, down_plot, ncol = 1, nrow = 2)
ggsave("Fig. 4c.pdf", plot = up_down_trend_plot, width = 4, height = 4)
trend_plot <- function(trend_matrix, my_col, title){
  time<-colnames(trend_matrix)
  trend_matrix<-t(apply(trend_matrix, 1, function(x)scale(x, center = F)))
  
  t<-1:length(time)
  data<-as.data.frame(t(trend_matrix))
  data$t<-t
  data_long<-reshape::melt(data,id.vars = "t")
  x<-1:length(time)
  quantile_value<-apply(trend_matrix,2,function(x){
    temp<-round(quantile(x,seq(0,1,0.01)),3)
    return(temp)
  })
  pdata.list<-list()
  for(k in 1:50){
    pdata<-data.frame(x,lower = quantile_value[k,],upper = quantile_value[102-k,])
    pdata.list[[k]]<-pdata
  }
  myPalette <- colorRampPalette(my_col)(45)
  plot.trend<-ggplot()+
    geom_rect(aes(xmin = -Inf,xmax = 1,ymin = -Inf,ymax = Inf),fill = "white",alpha = 0.4)+
    geom_rect(aes(xmin = 1,xmax = length(time),ymin = -Inf,ymax = Inf),fill = "#DEF6F3",alpha = 0.4)+
    geom_rect(aes(xmin = length(time),xmax = Inf,ymin = -Inf,ymax = Inf),fill = "white",alpha = 0.4)+
    theme(panel.grid=element_blank(),panel.border=element_rect(fill=NA,color="black", size=1))
  
  for(k in 6:50){
    pdata<-pdata.list[[k]]
    plot.trend<-plot.trend+geom_ribbon(data = pdata,aes(ymin=lower, ymax=upper, x=x), fill = myPalette[k-5], alpha = 1)
  }
  plot.trend<-plot.trend+
    theme_bw()+
    xlab("Time")+
    ylab("logFC(T2 vs T1)")+ 
    labs(title = title)+
    theme(panel.grid =element_blank()) + 
    scale_x_continuous(breaks=seq(1, length(time), 1), labels = c("T1","T2","T3"))
  return(plot.trend)
}
# 02End: Changing trend of DEGs
