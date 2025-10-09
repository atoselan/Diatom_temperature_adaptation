#!/bin/bash
#SBATCH --job-name=SOC_job
#SBATCH --partition=compute-64-512
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=6G
#SBATCH --time=8:00:00
#SBATCH --output=SOC_OG_output.log
#SBATCH --error=SOC_OG_error.log

module load hmmer/3.3

# List of OG prefixes
prefixes=(OG0000015_a OG0000015_b OG0000120_a OG0000238_a OG0000238_b OG0000238_c OG0000238_d
          OG0000249_a OG0000249_b OG0000249_c OG0000822 OG0001774 OG0001985 OG0002054 OG0002118
          OG0002271 OG0002593 OG0003601 OG0004069 OG0004106 OG0004405 OG0004972 OG0006595
          OG0006603 OG0015328 OG0029185)
# Define metaT and output directories
SOC_DIR="/gpfs/data/mock_lab/reference-data/SEA_OF_CHANGE/METAT/CONTIGS"
OUT_DIR="/gpfs/data/mock_lab/andrew-temp/BRCA2_hmm_model/HDR_orthogroups/SOC_HMM2_RESULTS"

# Loop through OGs and query SOC samples
for P in "${prefixes[@]}"; do
  HMM="/gpfs/data/mock_lab/andrew-temp/BRCA2_hmm_model/HDR_orthogroups/${P}_hmm2/${P}.hmm"
  mkdir -p "${OUT_DIR}/${P}"
  for f in "${SOC_DIR}"/*.faa.gz; do
    prefix=$(basename "${f%.assembled.faa.gz}")
    echo "Processing: ${prefix}..."
    hmmsearch --cpu ${SLURM_CPUS_PER_TASK} --tblout "${OUT_DIR}/${P}/${prefix}_${P}.domtblout" \
      --noali "${HMM}" <(gunzip -c "${f}") > /dev/null
  done
done
