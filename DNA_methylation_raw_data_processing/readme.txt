### DNA methylation data preprocessing
For all the six missions named M13, M15, M33, M90, M180-1, M180-2, read the raw idat files, normalize with champ.norm, remove individual differences and across-mission batch effect, and finally generate beta.Rds, norm_beta.Rds and nobatch_beta.Rds.

Directory structure in DNA methylation analysis part should strictly follow the one below to ensure successful execution:   
--Program  
  --13day_spaceflight_M13  
    --DNA_methylation  
      --rawdata  
          *.idat  
          Samplesheet.csv  
  --15day_spaceflight_M15  
    --DNA_methylation  
      --rawdata  
          *.idat  
          Samplesheet.csv  
  --33day_spaceflight_M33  
    --DNA_methylation  
      --rawdata  
          *.idat  
          Samplesheet.csv  
  --90day_spaceflight_M90  
    --DNA_methylation  
      --rawdata  
          *.idat  
          Samplesheet.csv  
  --180day_spaceflight_M180_1  
    --DNA_methylation  
      --rawdata  
          *.idat  
          Samplesheet.csv  
  --180day_spaceflight_M180_2  
    --DNA_methylation  
      --rawdata  
          *.idat  
          Samplesheet.csv  
  --log
    *nohup*.txt
  --main  
    --Alternative_splicing  
      --result  
      *.R  
    --Differential_methylation_analysis  
      --genome_conver  
      *.R  
    --DNA_methylation_raw_data_processing  
      *.R  
    --Expression_methylation_correlation  
      *.R  
    --Multi-group_comparisons_of_geneset_methy_trend  
      *.R  
    --utils  
      *.R  

Aftering preparing the required environment and packages, run the following commands:
### DNA_methylation_raw_data_processing ###
nohup Rscript ./main/DNA_methylation_raw_data_processing/step0_samplesheet_sorting.R > ./log/nohup_DNAm_raw_step0_samplesheet_sorting.txt
nohup Rscript ./main/DNA_methylation_raw_data_processing/step1_DNA_meth_raw_data_process.R > ./log/nohup_DNAm_raw_step1_DNA_meth_raw_data_process.txt
nohup Rscript ./main/DNA_methylation_raw_data_processing/step2_DNA_meth_batch_remove.R -f _2_1 -p 2.5 -s 123123 > ./log/nohup_DNAm_raw_step2_DNA_meth_batch_remove21.txt
nohup Rscript ./main/DNA_methylation_raw_data_processing/step2_DNA_meth_batch_remove.R -f _6 -p 2.5 -s 123123 > ./log/nohup_DNAm_raw_step2_DNA_meth_batch_remove6.txt
nohup Rscript ./main/DNA_methylation_raw_data_processing/step2_DNA_meth_batch_remove.R -f _3_3 -p 2.5 -s 123123 > ./log/nohup_DNAm_raw_step2_DNA_meth_batch_remove33.txt

