step1_setup.sh: Small RNA-seq Analysis Pipeline with nf-core/smrnaseq
This repository provides a ready-to-use setup for analyzing small RNA sequencing data using the [`nf-core/smrnaseq`](https://nf-co.re/smrnaseq/2.4.0/) pipeline. The workflow handles upstream preprocessing and quantification steps.

step2_analysis.R: This script performs t-SNE visualization, batch effect correction, and differential expression analysis (DEA) on mature miRNA count data using edgeR and limma.
