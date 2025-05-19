System requirements:

softwares:
* R (v4.2.1)

R packages:
* WGCNA_1.72.5
* reshape2_1.4.4
* ggplot2_3.5.0
* stringr_1.5.1
* org.Hs.eg.db_3.16.0
* clusterProfiler_4.6.2
 * future_1.34.0

Installation guide:
#Install R software
conda install -c r r-base=4.2.1
#Install R packages in R (v4.2.1)
install.packages("WGCNA")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("stringr")
install.packages("org.Hs.eg.db")
install.packages("clusterProfiler")
install.packages("future")


Demo:
The expression matrix, nobatch_expression_profile.txt are available for download at https://www.spacelifescience.cn/search/ (search for 90-day spaceflight_M90 / 180-day spaceflight_180-1, RNA-seq) (saved in "90-day spaceflight_M90/RNA-seq" and "180-day spaceflight_M180-1/RNA-seq").  
The expected results can be found  in the result/ folder. The estimated time for this task is about 2-3 hours.


Instructions for use:
1. Download the datasets from the provided link and  save them in a folder with the same name.
2. Install the required R packages as outlined in the installation guide.
3. Set the directory where you saved the datasets as the working directory in R.
4. source("enrich_function.R")
5. Execute the R scripts in order. Step 1: Division of genes into modules by WGCNA  and functional enrichment analysis for module genes. Step 2: The construction of PPI network, as well as GO and KEGG enrichment analysis for the PPI nodes. Step 3: Geneset interaction analysis for genesets enriched by PPI nodes . Step 4: Gene set annotation (categorization) and importance ranking.

Note: Make sure to adjust the file paths in the script to match the location of your downloaded datasets.