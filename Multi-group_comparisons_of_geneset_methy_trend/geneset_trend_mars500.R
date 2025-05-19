##### description: Gene set methylation change trends in Mars500

#load packages
library(readxl)
library(RColorBrewer)
library(ggplot2)

#Load genesets enriched by PPI nodes 
go_enrich_res <- read.table( file = "WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_GOenrich_res.txt", sep="\t")
kegg_enrich_res <- read.table( file = "WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_KEGGenrich_res.txt", sep="\t")
geneset_all_df_PPI <- rbind(go_enrich_res, kegg_enrich_res)

#Load DNA methylation data of Mars-500,nobatch_beta_Mars.txt and SampleSheet.xlsx are available for download at https://www.spacelifescience.cn/search/ (search for Mars-500, DNA methylation)
myCombat500 <- read.table("Mars-500/DNA methylation/nobatch_beta_Mars.txt", header=T, row.names=1, as.is=T, sep="\t")
SampleSheet500 <- read_excel("Mars-500/DNA methylation/SampleSheet.xlsx")


library(ChAMP)
library(ggplot2)
times <- factor(SampleSheet500$Sample_Group, levels=unique(SampleSheet500$Sample_Group))
times <- unique(times)
compare.group <- times[1:2]
DMP_res_500 <- champ.DMP(beta = myCombat500,
                         pheno = as.character(SampleSheet500$Sample_Group),
                         compare.group = compare.group,
                         adjPVal = 1,
                         adjust.method = "BH",
                         arraytype = "450K")
DMP_res_500_df <- DMP_res_500[[1]]

  aov_fun <- function(x, group){
     df <- data.frame(value=unlist(x), time=group)
     aov_res <- aov(value~time, df)
     aov_res <- summary(aov_res)
     p <- aov_res[[1]][["Pr(>F)"]][1]
     return(p)
  }
plot_list_mars500 <- c()
for(term in geneset_all_df_PPI$Description){
  interested_geneset <- geneset_all_df_PPI[geneset_all_df_PPI$Description==term, "geneID"]
  interested_geneset <- unique(unlist(strsplit(interested_geneset,"/")))

  promoter_region<-c("TSS1500","TSS200","5'UTR","1stExon")
  interested_geneset_annotation <- DMP_res_500_df[DMP_res_500_df$gene%in%interested_geneset & DMP_res_500_df$feature%in%promoter_region, ]
  interested_geneset_meth_probe <-  as.data.frame(myCombat500[rownames(interested_geneset_annotation), ])
  interested_geneset_meth_probe$gene <- interested_geneset_annotation$gene

  interested_geneset_meth_mat <- aggregate(interested_geneset_meth_probe[,-ncol(interested_geneset_meth_probe)], by = list(group = interested_geneset_annotation$gene), median)
  rownames(interested_geneset_meth_mat) <- interested_geneset_meth_mat$group
  interested_geneset_meth_mat <- interested_geneset_meth_mat[,-1]


  aov_p <- apply(interested_geneset_meth_mat, 1, aov_fun, group=SampleSheet500$Sample_Group)
  f <- which(aov_p<0.05)
  filter_interested_geneset_meth_mat <- interested_geneset_meth_mat[f, ]
 
  if(nrow(filter_interested_geneset_meth_mat)>1){
     r <- apply(filter_interested_geneset_meth_mat, 1, function(x){max(x)-min(x)})
     filter_interested_geneset_meth_mat <- filter_interested_geneset_meth_mat[r>0.05,]

     if(nrow(filter_interested_geneset_meth_mat)>1){
        mean_exp_mat <- apply(filter_interested_geneset_meth_mat, 1, function(x) tapply(x, factor(SampleSheet500$Sample_Group, levels=unique(SampleSheet500$Sample_Group)), mean))
        mean_exp_mat <- t(mean_exp_mat)

        trend_p <- trend_plot(mean_exp_mat)

        plot_list_mars500 <- c(plot_list_mars500, list(trend_p))
		term <- stringr::str_split_i(term, " - ", 1)
        names(plot_list_mars500)[length(plot_list_mars500)] <- term
     }

    }
}


saveRDS(plot_list_mars500, "Multi-group comparisons of geneset methy trend/result/plot_list_mars500.Rds")
for(i in 1:length(plot_list_mars500)){
  	filename<-paste("Multi-group comparisons of geneset methy trend/result/mars500_",names(plot_list_mars500)[i],".pdf",sep = "")
    pdf(filename,width=5,height=3)
	trend_p <- plot_list_mars500[i]
    print(trend_p)
    dev.off()
}

trend_plot <- function(trend_matrix){
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
  my_col <- c("#C6E4F5", "#2BAAF3")
  myPalette <- colorRampPalette(my_col)(45)
  plot.trend<-ggplot()+
    geom_rect(aes(xmin = -Inf,xmax = 1,ymin = -Inf,ymax = Inf),fill = "#FFFFCC",alpha = 0.4)+
    geom_rect(aes(xmin = 1,xmax = length(time),ymin = -Inf,ymax = Inf),fill = "#DEF6F3",alpha = 0.4)+
    geom_rect(aes(xmin = length(time),xmax = Inf,ymin = -Inf,ymax = Inf),fill = "#FFFFCC",alpha = 0.4)+
    theme(panel.grid=element_blank(),panel.border=element_rect(fill=NA,color="black", size=1))

  for(k in 6:50){
    pdata<-pdata.list[[k]]
    plot.trend<-plot.trend+geom_ribbon(data = pdata,aes(ymin=lower, ymax=upper, x=x), fill = myPalette[k-5], alpha = 1)
  }
  plot.trend<-plot.trend+
    theme_bw()+
    xlab("")+
    theme(panel.grid =element_blank()) + 
    scale_x_continuous(breaks=seq(1, length(time), 1), labels = time)
  return(plot.trend)
}