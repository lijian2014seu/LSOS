### description:
# Lift over probe coordinates from hg19 to hg38 and save annotated mapping:
#  - Load 450K and EPIC probe annotations (hg19) and merge into unified set
#  - Convert to GRanges, apply UCSC chain file to liftOver to hg38
#  - Expand multi-mapping results, retain original probeIDs
#  - Assemble final data frame with hg19 and hg38 coords and save as RDS

library(minfi)
library(GenomicRanges)
library(rtracklayer)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(liftOver)

setwd("/home/lqwang/Program")

anno450K <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annoEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
saveRDS(rownames(anno450K), "./main/Differential_methylation_analysis/genome_conver/all_450K_probes.Rds")
saveRDS(rownames(annoEPIC), "./main/Differential_methylation_analysis/genome_conver/all_EPIC_probes.Rds")


df450K <- data.frame(
  probeID = rownames(anno450K),
  CHR     = anno450K$chr,
  MAPINFO = anno450K$pos,
  stringsAsFactors=FALSE
)

dfEPIC <- data.frame(
  probeID = rownames(annoEPIC),
  CHR     = annoEPIC$chr,
  MAPINFO = annoEPIC$pos,
  stringsAsFactors=FALSE
)

df_hg19 <- merge(
  df450K, dfEPIC,
  by="probeID", all=TRUE, suffixes=c("_450K","_EPIC")
)

df_hg19 <- df_hg19 %>% 
           mutate(
                  CHR = coalesce(CHR_EPIC, CHR_450K), 
                  MAPINFO = coalesce(MAPINFO_EPIC, MAPINFO_450K)
                 )
df_hg19 <- df_hg19[, c("probeID", "CHR", "MAPINFO")]

gr_hg19 <- GRanges(
  seqnames = df_hg19$CHR,
  ranges   = IRanges(start=df_hg19$MAPINFO, width=1),
  strand   = "*",
  probeID  = df_hg19$probeID
)
chain <- import.chain("./main/Differential_methylation_analysis/genome_conver/hg19ToHg38.over.chain")

gr_list <- liftOver(gr_hg19, chain)


gr_hg38 <- unlist(gr_list)
probe_map <- rep(gr_hg19$probeID, elementNROWS(gr_list))
mcols(gr_hg38)$probeID <- probe_map


df_out <- data.frame(
  probeID      = mcols(gr_hg38)$probeID,
  CHR          = as.character(seqnames(gr_hg19))[match(mcols(gr_hg38)$probeID, df_hg19$probeID)],
  MAPINFO      = start(gr_hg19)[match(mcols(gr_hg38)$probeID, df_hg19$probeID)],
  CHR_hg38     = as.character(seqnames(gr_hg38)),
  MAPINFO_hg38 = start(gr_hg38),
  stringsAsFactors = FALSE
)

rownames(df_out) <- df_out$probeID

saveRDS(df_out, "./main/Differential_methylation_analysis/genome_conver/probe_annotation_hg19tohg38.Rds")
