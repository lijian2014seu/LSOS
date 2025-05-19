java -version
# openjdk version "1.8.0_412"

# Create and activate a new conda environment
conda create -n nf-core -y
conda init bash
source ~/.bashrc
conda activate nf-core

# Add channels and configure priority
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Install Nextflow
conda install -c bioconda nextflow

# Check version
nextflow -version

# Test the Pipeline
nextflow run nf-core/smrnaseq -profile test,docker --outdir ./results

#Run with your own data
Prepare a samplesheet.csv in the format required by nf-core. See example:
https://nf-co.re/smrnaseq/2.4.0/

conda activate nf-core

nextflow run nf-core/smrnaseq \
  -profile docker \
  --igenomes_base smrnaseq-2.4.0 \
  --input samplesheet.csv \
  --genome GRCh38 \
  --mirtrace_species hsa \
  --protocol illumina \
  --outdir <OUTDIR>
