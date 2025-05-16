### description:
# Annotate CpG probes with gene and feature context across 450K and EPIC arrays:
#  - Generate placeholder beta matrices for “Case” vs “Ctrl” samples
#  - Run ChAMP’s DMP function to retrieve per‑probe gene/feature annotations
#  - Merge 450K and EPIC results, prioritizing EPIC when both exist
#  - Return a unified table of CpG → gene and genomic feature mappings

library(GenomicRanges)
library(annotatr)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ChAMP)
library(dplyr)

annotate_gene_feature_for_cpg <- function (CpG_list_450k = readRDS("./main/Differential_methylation_analysis/genome_conver/all_450K_probes.Rds"), 
                                           CpG_list_EPIC = readRDS("./main/Differential_methylation_analysis/genome_conver/all_EPIC_probes.Rds"), 
                                           n1 = 3, n2 = 3, case_val = 0.2, ctrl_val = 0.8) {
    case_samples <- paste0("Case", seq_len(n1))
    ctrl_samples <- paste0("Ctrl", seq_len(n2))
    samples <- c(case_samples, ctrl_samples)
    pheno <- data.frame(
        Sample_Group = factor(c(rep("Case", n1), rep("Ctrl", n2))),
        Sample_Name = samples,
        stringsAsFactors = FALSE
    )

    gen_anno <- function (CpG_list, arraytype, samples, pheno) {
        beta_vals <- c(rep(case_val, length(CpG_list) * n1),
                        rep(ctrl_val, length(CpG_list) * n2))
        beta <- matrix(
            beta_vals,
            nrow = length(CpG_list),
            ncol = n1 + n2,
            dimnames = list(CpG_list, samples)
        )
        dmp_list <- champ.DMP(
            beta = beta,
            pheno = pheno$Sample_Group,
            compare.group = NULL,
            adjPVal = 1,
            arraytype = arraytype
        )
        raw_tbl <- dmp_list[[1]]
        raw_tbl[, c("gene", "feature")]
    }
    tbl_450k <- gen_anno(CpG_list_450k, "450K", samples, pheno)
    tbl_EPIC <- gen_anno(CpG_list_EPIC, "EPIC", samples, pheno)

    tbl_450k$CpG <- rownames(tbl_450k)
    tbl_EPIC$CpG <- rownames(tbl_EPIC)
    tbl_merge <- merge(
        tbl_450k, tbl_EPIC,
        by = "CpG",
        suffixes = c("_450K", "_EPIC"),
        all = TRUE
    )
    tbl_final <- tbl_merge %>%
    mutate(
        gene = coalesce(gene_EPIC, gene_450K),
        feature = coalesce(feature_EPIC, feature_450K)
    )
    rownames(tbl_final) = tbl_final$CpG
    tbl_final[, c("gene", "feature")]
}

filter_DMP_in_a_particular_genome_region <- function (DMP_result_path, DMR_list, g1, g2, condition) {
    DMP_df = readRDS(DMP_result_path)[[paste0(g1, " to ", g2)]]
    DMP_df = DMP_df[!is.na(DMP_df$CHR_hg38) & DMP_df$change_deltaBeta0.05 %in% condition, ]
    if (is.null(DMR_list) || nrow(DMR_list) == 0) {
        message("DMR_list is empty. Returning an empty data.frame.")
        return(DMP_df[FALSE, , drop = FALSE])
    }

    gr_DMP <- GRanges(
            seqnames = DMP_df$CHR_hg38,
            ranges = IRanges(start = DMP_df$MAPINFO_hg38,
                                end = DMP_df$MAPINFO_hg38),
            strand = "*"
        )
    names(gr_DMP) <- rownames(DMP_df)
    gr_DMR <- GRanges(
            seqnames = DMR_list$chr_hg38,
            ranges = IRanges(start = DMR_list$start_hg38,
                                end = DMR_list$end_hg38),
            strand = "*"
        )
    mcols(gr_DMR)$DMR_index <- DMR_list$DMR_index
    
    hits <- findOverlaps(gr_DMP, gr_DMR, ignore.strand = TRUE)
    qry <- queryHits(hits)
    sub <- subjectHits(hits)

    keep_first <- !duplicated(qry)
    qry_unique <- qry[keep_first]
    sub_unique <- sub[keep_first]
    res <- DMP_df[qry_unique, , drop = FALSE]
    res$DMR_index <- mcols(gr_DMR)$DMR_index[sub_unique]
    return(res)
}


filter_probes_in_genome_regions <- function(probe_df, DMR_list) {
    if (!all(c("CHR_hg38", "MAPINFO_hg38") %in% colnames(probe_df))) {
        stop("!all(c(CHR_hg38, MAPINFO_hg38) %in% colnames(probe_df))")
    }
    if (!all(c("chr_hg38", "start_hg38", "end_hg38") %in% colnames(DMR_list))) {
        stop("!all(c(chr_hg38, start_hg38, end_hg38) %in% colnames(DMR_list))")
    }

    gr_probes <- GRanges(
        seqnames = probe_df$CHR_hg38,
        ranges = IRanges(
            start = probe_df$MAPINFO_hg38,
            end = probe_df$MAPINFO_hg38
        ),
        strand = "*"
    )
    gr_regions <- GRanges(
        seqnames = DMR_list$chr_hg38,
        ranges = IRanges(
            start = DMR_list$start_hg38,
            end = DMR_list$end_hg38
        ),
        strand = "*"
    )

    hits <- findOverlaps(gr_probes, gr_regions, ignore.strand = TRUE)
    probe_hits <- unique(queryHits(hits))

    data.frame(
        in_region = seq_len(nrow(probe_df)) %in% probe_hits,
        stringsAsFactors = FALSE,
        row.names = rownames(probe_df)
    )
}
