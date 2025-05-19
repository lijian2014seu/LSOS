To do expression-methylation correlation analysis, you need to first finish DNA methylation raw data processing and differential methylation analysis.
The commands are presented as below:
nohup Rscript ./main/exp_meth_corr/exp_meth_corr.R -f _2_1 -e M90 -E M180-1 -r promoter > ./log/nohup_corr_exp_meth_corr.R-f_2_1-eM90-EM180-1-rpromoter.txt
nohup Rscript ./main/exp_meth_corr/exp_meth_corr.R -f _2_1 -e M90 -E M180-1 -r enhancer > ./log/nohup_corr_exp_meth_corr.R-f_2_1-eM90-EM180-1-renhancer.txt
nohup Rscript ./main/exp_meth_corr/exp_meth_corr.R -f _2_1 -e M90 -E M180-1 -r genebody > ./log/nohup_corr_exp_meth_corr.R-f_2_1-eM90-EM180-1-rgenebody.txt
nohup Rscript ./main/exp_meth_corr/exp_meth_corr_M180_2.R -f _2_1 -e M180-2 -r promoter > ./log/nohup_corr_exp_meth_corr_M180_2.R-f_2_1-eM180-2-rpromoter.txt
nohup Rscript ./main/exp_meth_corr/exp_meth_corr_M180_2.R -f _2_1 -e M180-2 -r enhancer > ./log/nohup_corr_exp_meth_corr_M180_2.R-f_2_1-eM180-2-renhancer.txt
nohup Rscript ./main/exp_meth_corr/exp_meth_corr_M180_2.R -f _2_1 -e M180-2 -r genebody > ./log/nohup_corr_exp_meth_corr_M180_2.R-f_2_1-eM180-2-rgenebody.txt
