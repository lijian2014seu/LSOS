### Differential methylation analysis:
Conduct DMP and DMR analysis for each mission in each comparison groups.
There are 7 steps to finish the DMP analysis (step1_1 to step 1_7) and 4 steps to finish the DMR analysis (step2_1 to step2_4), finally generating raw tables and figures for the paper.
Detailed input and output information can be accessed from the annotation in the scripts.
If the code of any figure or table in the published paper is not included, please feel free to contact us.

Please run the commands in DNA_methylation_raw_data_processing first, then commands below for differential methylation analysis:
### Differential_methylation_analysis ###

# DMP analysis
nohup Rscript ./main/Differential_methylation_analysis/step1_1_DMP_analysis.R -f _2_1 > ./log/nohup_DNAm_analysis_step1_1_DMP_analysis21.txt
nohup Rscript ./main/Differential_methylation_analysis/step1_1_DMP_analysis.R -f _3_3 > ./log/nohup_DNAm_analysis_step1_1_DMP_analysis33.txt
nohup Rscript ./main/Differential_methylation_analysis/step1_2_DMP_num_statistic_and_trend_plot.R -f _2_1 > ./log/nohup_DNAm_analysis_step1_2_DMP_num_statistic_and_trend_plot21.txt
nohup Rscript ./main/Differential_methylation_analysis/step1_2_DMP_num_statistic_and_trend_plot.R -f _3_3 > ./log/nohup_DNAm_analysis_step1_2_DMP_num_statistic_and_trend_plot33.txt

nohup Rscript ./main/Differential_methylation_analysis/step1_3_DMP_growth_calculate.R -f _3_3 --group1 T1 --group2 T2 -w 1000000 -d up > ./log/nohup_DNAm_analysis_step1_3_DMP_growth_calculate.R.txt
nohup Rscript ./main/Differential_methylation_analysis/step1_3_DMP_growth_calculate.R -f _3_3 --group1 T1 --group2 T2 -w 1000000 -d down > ./log/nohup_DNAm_analysis_step1_3_DMP_growth_calculate.R.txt
nohup Rscript ./main/Differential_methylation_analysis/step1_4_DMP_growth_cal_core.R -f _3_3 --group1 T1 --group2 T2 -w 1000000 -d up > ./log/nohup_DNAm_analysis_step1_4_DMP_growth_cal_core.R.txt
nohup Rscript ./main/Differential_methylation_analysis/step1_4_DMP_growth_cal_core.R -f _3_3 --group1 T1 --group2 T2 -w 1000000 -d down > ./log/nohup_DNAm_analysis_step1_4_DMP_growth_cal_core.R.txt
nohup Rscript ./main/Differential_methylation_analysis/step1_5_DMP_growth_plot.R -f _3_3 --group1 T1 --group2 T2 -w 1000000 -l 0.15 -d up > ./log/nohup_DNAm_analysis_step1_5_DMP_growth_plot.txt
nohup Rscript ./main/Differential_methylation_analysis/step1_5_DMP_growth_plot.R -f _3_3 --group1 T1 --group2 T2 -w 1000000 -l 0.15 -d down > ./log/nohup_DNAm_analysis_step1_5_DMP_growth_plot.txt

nohup Rscript ./main/Differential_methylation_analysis/step1_6_DMP_num_change_with_increasing_inflight_time_T2vs.T1.R -f _3_3 > ./log/nohup_DNAm_analysis_step1_6_DMP_num_change_with_increasing_inflight_time_T2vs.T1.R.txt
nohup Rscript ./main/Differential_methylation_analysis/step1_7_DMP_num_change_with_increasing_recovery_time_R01_R10_R30.R -f _3_3 -m M13 > ./log/nohup_DNAm_analysis_step1_7_DMP_num_change_with_increasing_recovery_time_R01_R10_R30.R_M13.txt
nohup Rscript ./main/Differential_methylation_analysis/step1_7_DMP_num_change_with_increasing_recovery_time_R01_R10_R30.R -f _3_3 -m M15 > ./log/nohup_DNAm_analysis_step1_7_DMP_num_change_with_increasing_recovery_time_R01_R10_R30.R_M15.txt
nohup Rscript ./main/Differential_methylation_analysis/step1_7_DMP_num_change_with_increasing_recovery_time_R01_R10_R30.R -f _3_3 -m M33 > ./log/nohup_DNAm_analysis_step1_7_DMP_num_change_with_increasing_recovery_time_R01_R10_R30.R_M33.txt
nohup Rscript ./main/Differential_methylation_analysis/step1_8_DMP_num_change_with_increasing_recovery_time_R01_R60.R -f _2_1 -m M90 > ./log/nohup_DNAm_analysis_step1_8_DMP_num_change_with_increasing_recovery_time_R01_R60.R_M90.txt
nohup Rscript ./main/Differential_methylation_analysis/step1_8_DMP_num_change_with_increasing_recovery_time_R01_R60.R -f _2_1 -m M180-1 > ./log/nohup_DNAm_analysis_step1_8_DMP_num_change_with_increasing_recovery_time_R01_R60.R_M180-1.txt
nohup Rscript ./main/Differential_methylation_analysis/step1_8_DMP_num_change_with_increasing_recovery_time_R01_R60.R -f _2_1 -m M180-2 > ./log/nohup_DNAm_analysis_step1_8_DMP_num_change_with_increasing_recovery_time_R01_R60.R_M180-2.txt
nohup Rscript ./main/Differential_methylation_analysis/step1_9_promoter_region_DMP_proportion.R -f _2_1 > ./log/nohup_DNAm_analysis_step1_9_promoter_region_DMP_proportion.R.txt

# DMR analysis
nohup Rscript ./main/Differential_methylation_analysis/step2_1_DMR_analysis.R -f _2_1 -g 1000 -G 300 > ./log/nohup_DNAm_analysis_step2_1_DMR_analysis21.txt
nohup Rscript ./main/Differential_methylation_analysis/step2_1_DMR_analysis.R -f _3_3 -g 1000 -G 300 > ./log/nohup_DNAm_analysis_step2_1_DMR_analysis33.txt

nohup Rscript ./main/Differential_methylation_analysis/step2_2_check_DMR_DMP_result_match.R -m M90 --group1 T1 --group2 T2 -f _2_1 -g _1000 -G _300 > ./log/nohup_DNAm_analysis_step2_2_check_DMR_DMP_result_match.R_M90.txt
nohup Rscript ./main/Differential_methylation_analysis/step2_2_check_DMR_DMP_result_match.R -m M180-1 --group1 T1 --group2 T2 -f _2_1 -g _1000 -G _300 > ./log/nohup_DNAm_analysis_step2_2_check_DMR_DMP_result_match.R_M180-1.txt
nohup Rscript ./main/Differential_methylation_analysis/step2_2_check_DMR_DMP_result_match.R -m M180-2 --group1 T1 --group2 T2 -f _2_1 -g _1000 -G _300 > ./log/nohup_DNAm_analysis_step2_2_check_DMR_DMP_result_match.R-M180-2.txt

nohup Rscript ./main/Differential_methylation_analysis/step2_3_DMR_circlize_plot.R -f _2_1 -g _1000 -G _300 > ./log/nohup_DNAm_analysis_step2_3_DMR_circlize_2_1.txt
nohup Rscript ./main/Differential_methylation_analysis/step2_3_DMR_circlize_plot.R -f _3_3 -g _1000 -G _300 > ./log/nohup_DNAm_analysis_step2_3_DMR_circlize_3_3.txt
nohup Rscript ./main/Differential_methylation_analysis/step2_4_DMR_intersection.R -f _2_1 -g _1000 -G _300 > ./log/nohup_DNAm_analysis_step2_4_DMR_intersection.R_2_1.txt
