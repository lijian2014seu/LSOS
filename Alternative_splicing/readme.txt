To do alternative splicing analysis, you need to first finish DNA methylation raw data processing and differential methylation analysis.
The commands are presented as below:
nohup Rscript ./main/Alternative_splicing/DMPs_on_alt_spliced_gene.R -f _2_1 -d up_down > ./log/nohup_asDMP_DMPs_on_alt_spliced_gene.R-f_2_1-dup_down.txt
nohup Rscript ./main/Alternative_splicing/DMPs_on_alt_spliced_gene.R -f _2_1 -d up > ./log/nohup_asDMP_DMPs_on_alt_spliced_gene.R-f_2_1-dup.txt
nohup Rscript ./main/Alternative_splicing/DMPs_on_alt_spliced_gene.R -f _2_1 -d down > ./log/nohup_asDMP_DMPs_on_alt_spliced_gene.R-f_2_1-ddown.txt
