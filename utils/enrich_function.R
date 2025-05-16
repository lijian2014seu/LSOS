### description: Perform GO Biological Process enrichment on a vector of gene symbols ###

library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(stringr)

go_enrich <- function(gene_symbols, gene_type) {
  empty_result <- data.frame(
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    Count = integer(),
    ID = character(),
    gene_type = character(),
    stringsAsFactors = FALSE
  )
  gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
  if (length(gene_symbols) == 0) {
    message("No genes provided; skipping GO enrichment.")
    return(empty_result)
  }
  cols <- c("SYMBOL", "ENSEMBL", 'ENTREZID')
  gene = AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbols, columns = cols, keytype = "SYMBOL")
  print( length(gene$ENTREZID) )
  IDs <- unique( na.omit( gene$ENTREZID ) )
  if (length(IDs) == 0) {
    message("No valid ENTREZID found; skipping GO enrichment.")
    return(empty_result)
  }
  ego_BP <- enrichGO(gene = IDs, OrgDb= org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", 
                      minGSSize = 10, pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE)
  if (is.null(ego_BP) || nrow(ego_BP) == 0) {
    message("No significant GO terms enriched.")
    return(empty_result)
  }
  ego_BP <- clusterProfiler::simplify(ego_BP)
  go.df_all <- data.frame(ego_BP)
  go.df_all$gene_type <- rep(gene_type, nrow(go.df_all))
  return( go.df_all )
}