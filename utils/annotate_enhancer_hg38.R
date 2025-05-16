### description:
# Annotate CpG probes with overlapping hg38 enhancer regions:
#  - Import enhancer BED from UCSC-style file (BED6+, with “name” field)  
#  - Build GRanges for probes and find overlaps against enhancers  
#  - Flag probes in enhancers and record enhancer names per probe

# Working directory — adjust to your environment
setwd("/home/lqwang/Program")

library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(readr)
library(annotatr)

# Path to enhancer annotations in hg38
bed_file_path = "./main/Differential_methylation_analysis/genome_conver/F5.hg38.enhancers.bed"
annotation = import.bed(con = bed_file_path)


annotate_enhancer_for_probes = function(df) {
    df = df[!is.na(df$CHR_hg38), ]

    # Create GRanges of probe locations (point ranges)
    probe_gr = GRanges(
        seqnames = df$CHR_hg38,
        ranges = IRanges(start = df$MAPINFO_hg38, end = df$MAPINFO_hg38),
        strand = "*"
    )

    # Find overlaps between probes and enhancer regions
    hits = findOverlaps(probe_gr, annotation, ignore.strand = TRUE)

    df$in_enhancer = FALSE
    df$enhancer_name = NA_character_

    probe_idx = queryHits(hits)
    enhancer_idx = subjectHits(hits)

    # Summarize multiple enhancer hits per probe into semicolon-separated names
    overlap_tbl = 
        data.frame(
            probe = probe_idx,
            enh_name = mcols(annotation)$name[enhancer_idx],
            stringsAsFactors = FALSE
        ) %>%
        group_by(probe) %>%
        summarise(enhancers = paste(unique(enh_name), collapse = ";"))

    df$in_enhancer[unique(probe_idx)] = TRUE
    df$enhancer_name[overlap_tbl$probe] = overlap_tbl$enhancers

    return(df)
}