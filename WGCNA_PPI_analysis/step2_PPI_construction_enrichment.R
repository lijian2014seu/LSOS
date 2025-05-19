##### description: The construction of PPI network, as well as GO and KEGG enrichment analysis for the PPI nodes

#Search for Brown module genes in Cytoscape software to obtain gene interaction information, which exists in PPI_nodes_interactions_in_string_Database_short.tsv; 
#then take the intersection with the WGCNA network for the Brown module

load("WGCNA_PPI_analysis/result/met_module_brown.RData")
head(edges_brown)

string_interactions_short <- read.csv("WGCNA_PPI_analysis/result/PPI_analysis/PPI_nodes_interactions_in_string_Database_short.tsv",sep="\t")
string_interactions_short <- string_interactions_short[string_interactions_short$combined_score>0.4,]
library(stringr)
sub_edges_brown <- edges_brown[edges_brown$weight>0.5^3,]

pairs <- apply(string_interactions_short, 1, function(x){
  pair <- x;
  index1 <- which(sub_edges_brown$fromNode==pair[1]);
  index2 <- which(sub_edges_brown$toNode==pair[2]);
  index <- intersect(index1, index2)
  if(length(index)!=0){
    return(c(unlist(sub_edges_brown[index, 1:3]), x[13]))
  }else{
    index1 <- which(sub_edges_brown$fromNode==pair[2]);
    index2 <- which(sub_edges_brown$toNode==pair[1]);
    index <- intersect(index1, index2)
    if(length(index)!=0){
      return(c(unlist(sub_edges_brown[index, 1:3]), x[13]))
    }
  }
} )
pairs <- pairs[!mapply(is.null, pairs)]
pairs <- do.call(rbind, pairs)
dim(pairs)
colnames(pairs)[4] <- "string_score"
write.table(pairs, "WGCNA_PPI_analysis/result/PPI_analysis/PPI_net_edge.txt", row.names = F, sep="\t", quote = F)

PPI_net_edge <- read.delim("WGCNA_PPI_analysis/result/PPI_analysis/PPI_net_edge.txt")
#Import the PPI network into Cytoscape, calculate the network topological properties, and then export the node attribute matrix.
PPI_net_node <- read.csv("WGCNA_PPI_analysis/result/PPI_analysis/PPI_net_node.csv")

#GO and KEGG enrichment analysis for PPI nodes
source("WGCNA_PPI_analysis/enrich_function.R")
go_enrich_res <- go_enrich(PPI_net_node$name, "PPI_node")
kegg_enrich_res <- kegg_enrich(PPI_net_node$name, "PPI_node")
save(go_enrich_res, kegg_enrich_res, file = "PPI_node_enrich_res_new.RData")
write.table(go_enrich_res, file = "WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_GOenrich_res.txt", quote = F, sep="\t")
write.table(kegg_enrich_res, file = "WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_KEGGenrich_res.txt", quote = F, sep="\t")
geneset_all_df_PPI <- rbind(go_enrich_res, kegg_enrich_res)


#Extract the subnetwork between two gene sets.
subnet_extraction <- function(geneset_name1, geneset_name2, geneset_df){
  geneset1 <- geneset_df[geneset_df$Description==geneset_name1, "geneID"]
  geneset1 <- unique(unlist(strsplit(geneset1,"/")))
  geneset2 <- geneset_df[geneset_df$Description==geneset_name2, "geneID"]
  geneset2 <- unique(unlist(strsplit(geneset2,"/")))
  interaction_index <- c(which(PPI_net_edge$fromNode%in%c(geneset1,geneset2) & PPI_net_edge$toNode%in%c(geneset1,geneset2)))
  sub_PPI_net_edge <- unique(PPI_net_edge[interaction_index, ])
  common_gene <- intersect(geneset1, geneset2)
  geneset1 <- geneset1[!geneset1%in%common_gene]
  geneset2 <- geneset2[!geneset2%in%common_gene]
  node_type <- data.frame(node=c(geneset1,geneset2,common_gene),type=c(rep(geneset_name1,length(geneset1)),rep(geneset_name2,length(geneset2)),rep("common_gene",length(common_gene))))
  write.table(node_type, paste0("WGCNA_PPI_analysis/result/gene_set_interact/",paste(geneset_name1, geneset_name2, "node_type.txt", sep="_")), row.names = F, quote = F, sep = "\t")
  sub_PPI_net_edge$cor_score <- (sub_PPI_net_edge$weight+sub_PPI_net_edge$string_score)/2
  sub_PPI_net_edge$edge_type <- "between"
  sub_PPI_net_edge$edge_type[sub_PPI_net_edge$fromNode%in%geneset1 & sub_PPI_net_edge$toNode%in%geneset1] <- "inner"
  sub_PPI_net_edge$edge_type[sub_PPI_net_edge$fromNode%in%geneset2 & sub_PPI_net_edge$toNode%in%geneset2] <- "inner"
  return(sub_PPI_net_edge)
}


term <- "MAPK signaling pathway"
term2 <- "Platelet activation"
subnet <- subnet_extraction(term, term2, geneset_all_df_PPI)
subnet[subnet$edge_type=="between",]
write.table(subnet, paste0("WGCNA_PPI_analysis/result/gene_set_interact/",paste(term, term2, "interaction_pairs.txt", sep="_")), row.names = F, quote = F, sep = "\t")
