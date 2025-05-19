#Geneset interaction analysis for genesets enriched by PPI nodes 
 
  ###Parallelism to speed up calculations (R environment in linux)
  ############
  #' Returns the total system memory in Mb.
  #' @returnType 
  #' @return 
  sysmem <- function() {
    t1 <- system("cat /proc/meminfo", intern = TRUE);
    t1 <- gsub("[[:alpha:][:punct:][:space:]]", "", t1[1])
    return(as.numeric(t1) / 1024) 
  }

  #' Returns the current memory occupied by R in linux, in Mb
  memLinux <- function() {
    temp <- as.numeric(system(paste("ps -p ", Sys.getpid(), " -o rss", sep = ""), intern = TRUE)[2])
    return(temp / 1024)
  }
  #' Calculating the right number of cores to run parallel operations
  #' @param multiCores Logical value, whether to perform multicore arithmetic
  #' @param autoCores  Logical value, whether the program automatically allocates cores based on reserved memory, if FALSE, all available thread cores will be used
  #' @param keepMem    What percentage of memory is reserved to prevent the program from running out of memory?
  workers <- function(multiCores = TRUE, autoCores = TRUE, keepMem = 0.1) {
    library(future)
    gc()
    if (multiCores) {
      if (autoCores) {
        t1 <- memLinux() ## R当前占用内存
        t2 <- sysmem() ## 系统总内存数
        works <- floor(t2 * (1 - keepMem) / t1) ## Number of cores available
        if (works < future::availableCores()) {
          works <- works
        } else {
          works <- future::availableCores()
        }
      } else {
        works <- future::availableCores()
      }
    } else {
      works <- 1
    }
    return(works)
  }

  options(future.globals.maxSize = 100.00 * 1024 ^ 3)
  worker.cores <- workers(multiCores = T, autoCores = T, keepMem = 0.3)
  future::plan(multicore, workers = worker.cores)

#Load PPI_net and genesets enriched by PPI nodes 
PPI_net_edge <- read.delim("WGCNA_PPI_analysis/result/PPI_analysis/PPI_net_edge.txt")
PPI_net_node <- read.csv("WGCNA_PPI_analysis/result/PPI_analysis/PPI_net_node.csv")
go_enrich_res <- read.table( file = "WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_GOenrich_res.txt", sep="\t")
kegg_enrich_res <- read.table( file = "WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_KEGGenrich_res.txt", sep="\t")
geneset_all_df_PPI <- rbind(go_enrich_res, kegg_enrich_res)

#geneset interaction analysis 
geneset_interact_scores_mat <- matrix(0, nrow(geneset_all_df_PPI), nrow(geneset_all_df_PPI))
rownames(geneset_interact_scores_mat) <- colnames(geneset_interact_scores_mat) <- geneset_all_df_PPI$Description
geneset_interact_scores_df <- reshape::melt(geneset_interact_scores_mat)
geneset_interact_scores_df <- geneset_interact_scores_df[geneset_interact_scores_df$X1!=geneset_interact_scores_df$X2,]
 
#res2 <- furrr::future_map(.x = 1:n2, get_cepair)
fun <- function(x) {x <- as.character(as.matrix(geneset_interact_scores_df[x,])); y <- geneset_interaction_analysis(geneset_name1=x[1], geneset_name2=x[2], geneset_df=geneset_all_df_PPI); return(y)}
#test
#system.time(interact_scores_df <- furrr::future_map(.x = 1:100, fun))
#for all
system.time(interact_scores_df <- furrr::future_map(.x = 1:nrow(geneset_interact_scores_df), fun))
interact_scores_df1 <- do.call(rbind,interact_scores_df)
geneset_interact_scores_df$value <- interact_scores_df1[,1]
geneset_interact_scores_df$p <- interact_scores_df1[,2]   
geneset_interact_scores_df <- geneset_interact_scores_df[!is.na(geneset_interact_scores_df$value), ]
geneset_interact_scores_df1 <- apply(geneset_interact_scores_df, 1, function(x){x <- as.character(as.matrix(x)); c(sort(x[1:2]),x[3:4])})
geneset_interact_scores_df <- as.data.frame(t(geneset_interact_scores_df1))
geneset_interact_scores_df <- geneset_interact_scores_df[!duplicated(geneset_interact_scores_df[,1:2]),]
geneset_interact_scores_df$V3 <- as.numeric(geneset_interact_scores_df$V3)
geneset_interact_scores_df$V4 <- as.numeric(geneset_interact_scores_df$V4)
sig_geneset_interact_scores_df <- geneset_interact_scores_df[which(geneset_interact_scores_df$V3>0 & geneset_interact_scores_df$V4<0.05),]

write.table(geneset_interact_scores_df, file = "WGCNA_PPI_analysis/result/gene_set_interact/geneset_interact_scores.txt", quote = F, sep="\t")
write.table(sig_geneset_interact_scores_df, file = "WGCNA_PPI_analysis/result/gene_set_interact/sig_geneset_interact_scores.txt", quote = F, sep="\t")


#Function for geneset interaction analysis   
geneset_interaction_analysis <- function(geneset_name1, geneset_name2, geneset_df, random_time=100){
  geneset1 <- geneset_df[geneset_df$Description==geneset_name1, "geneID"]
  geneset1 <- unique(unlist(strsplit(geneset1,"/")))
  geneset2 <- geneset_df[geneset_df$Description==geneset_name2, "geneID"]
  geneset2 <- unique(unlist(strsplit(geneset2,"/")))
  if(length(geneset1)>5 & length(geneset2)>5){
    interaction_index <- c(which(PPI_net_edge$fromNode%in%geneset1 & PPI_net_edge$toNode%in%geneset2),
                       which(PPI_net_edge$fromNode%in%geneset2 & PPI_net_edge$toNode%in%geneset1))
    sub_PPI_net_edge <- unique(PPI_net_edge[interaction_index, ])
    interaction_score_0 <- sum(sub_PPI_net_edge$weight+sub_PPI_net_edge$string_score)/(length(geneset1)*length(geneset2))
    n1 <- length(geneset1)
    n2 <- length(geneset2)
    nodes <- unique(c(PPI_net_edge$toNode,PPI_net_edge$fromNode))
    N <- length(nodes)
    ranadom_mat <- mapply(function(x){c(nodes[sample(1:N,n1)],nodes[sample(1:N,n2)])}, 1:random_time)
    random_scores <- apply(ranadom_mat, 2, function(x){
      ran_geneset1 <- x[1:n1]; ran_geneset2 <- x[1:n2]
      interaction_index <- c(which(PPI_net_edge$fromNode%in%ran_geneset1 & PPI_net_edge$toNode%in%ran_geneset2),
                             which(PPI_net_edge$fromNode%in%ran_geneset2 & PPI_net_edge$toNode%in%ran_geneset1))
      sub_PPI_net_edge <- unique(PPI_net_edge[interaction_index, ])
      interaction_score <- sum(sub_PPI_net_edge$weight+sub_PPI_net_edge$string_score)/(n1*n2)
      return(interaction_score)
    })
    p <- sum(random_scores>interaction_score_0)/random_time
    return(c(interaction_score_0,p))
  }else{return(NA)}
}
