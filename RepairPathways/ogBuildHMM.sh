#!/bin/bash
set -euo pipefail

# Usage: ./build_hmm_pipeline.sh input.fa [threads]
# Example: ./build_hmm_pipeline.sh OG0000015_a.fa 16

seq_file=${1:?Input fasta file required}
threads=${2:-16}

prefix=${seq_file%.fa}
log_file="${prefix}_pipeline.log"

echo "Starting pipeline for ${seq_file} at $(date)" | tee "$log_file"

function log_and_run {
    echo -e "\n>>> $* at $(date)" | tee -a "$log_file"
    eval "$@" 2>&1 | tee -a "$log_file"
    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        echo "ERROR: Command failed: $*" | tee -a "$log_file"
        exit 1
    fi
}

# Create working directory and link input
mkdir -p "${prefix}_hmm2"
cd "${prefix}_hmm2"
if [[ ! -e "${seq_file}" ]]; then
    ln -s "../${seq_file}" ./
fi

# Length filtering
log_and_run "seqkit fx2tab -l -n ${seq_file} > ${prefix}_lengths.txt"

MEDIAN=$(cut -f2 "${prefix}_lengths.txt" | sort -n | awk '{
    a[NR] = $1
}
END {
    if (NR % 2) {
        print a[int(NR/2)+1]
    } else {
        print int((a[NR/2] + a[NR/2 + 1]) / 2)
    }
}')

THRESHOLD=$(echo "$MEDIAN * 0.5" | bc | awk '{print int($1)}')
echo "Median seq length: $MEDIAN, threshold (50% median): $THRESHOLD" | tee -a "$log_file"

log_and_run "seqkit seq -g -m $THRESHOLD ${seq_file} > ${prefix}_filtered.fasta"
rm -f "${prefix}_lengths.txt"

# Loose clustering 20% id, 50% coverage
log_and_run "mkdir -p tmp_${prefix}"
log_and_run "mmseqs easy-cluster ${prefix}_filtered.fasta ${prefix}_mmclst tmp_${prefix} --min-seq-id 0.2 -c 0.5"
log_and_run "rm -rf tmp_${prefix}"

if [[ ! -f "${prefix}_mmclst_cluster.tsv" ]]; then
    echo "ERROR: MMseqs2 clustering output not found!" | tee -a "$log_file"
    exit 1
fi

BIGGEST_ID=$(cut -f1 "${prefix}_mmclst_cluster.tsv" | sort | uniq -c | sort -nr | head -1 | awk '{print $2}')
echo "Largest cluster representative: $BIGGEST_ID" | tee -a "$log_file"

if [[ ! -x ../get_seqs_from_cluster.sh ]]; then
    echo "ERROR: get_seqs_from_cluster.sh not found or not executable" | tee -a "$log_file"
    exit 1
fi

log_and_run "../get_seqs_from_cluster.sh $BIGGEST_ID ${prefix}_mmclst_cluster.tsv ${prefix}_mmclst_all_seqs.fasta ${prefix}_i20_c50.fa"

# More stringent clustering 95% id, 75% coverage
log_and_run "mkdir -p tmp_${prefix}"
log_and_run "mmseqs easy-cluster ${prefix}_i20_c50.fa ${prefix}_mmclst_i95_c75 tmp_${prefix} --min-seq-id 0.95 -c 0.75"
log_and_run "rm -rf tmp_${prefix}"

# Check clustered rep seq file exists
if [[ ! -f "${prefix}_mmclst_i95_c75_rep_seq.fasta" ]]; then
    echo "ERROR: MMseqs2 rep seqs not found!" | tee -a "$log_file"
    exit 1
fi

# Mafft alignment
log_and_run "mafft --thread $threads --maxiterate 1000 --localpair ${prefix}_mmclst_i95_c75_rep_seq.fasta > ${prefix}_mmclst_i95_c75.faa.mafft"

# Trim alignment
log_and_run "trimal -in ${prefix}_mmclst_i95_c75.faa.mafft -out ${prefix}_mmclst_i95_c75.faa.mafft.trimal -automated1"

# Filter gappy seqs (>75% gaps)
if [[ ! -f ~/Tools/python_scripts/filter_gappy_seqs.py ]]; then
    echo "ERROR: Python gap filter script not found at ~/Tools/python_scripts/filter_gappy_seqs.py" | tee -a "$log_file"
    exit 1
fi
log_and_run "python3 ~/Tools/python_scripts/filter_gappy_seqs.py ${prefix}_mmclst_i95_c75.faa.mafft.trimal ${prefix}_mmclst_i95_c75.faa.mafft.trimal.gapTrim 0.75"
rm -f ${prefix}_mmclst_i95_c75.faa.mafft.trimal

# Build HMM
log_and_run "hmmbuild ${prefix}.hmm ${prefix}_mmclst_i95_c75.faa.mafft.trimal.gapTrim"

# Query original sequences with hmm
log_and_run "hmmsearch \
  --cpu $threads \
  --tblout ${prefix}_hmmsearch.tbl \
  ${prefix}.hmm \
  ${seq_file} \
  > ${prefix}_hmmsearch.out"

# Get distribution of bitscores
grep -v '^#' "${prefix}_hmmsearch.tbl" | awk '{print $6}' | sort -nr > "${prefix}_bitscores.txt"


# Capture the output of the Python script
output=$(python3 ~/Tools/python_scripts/detect_bs_elbow.py ${prefix}_bitscores.txt)

# Extract the first two lines (cutoffs)
elbow_line=$(echo "$output" | sed -n '1p')
tenth_pct_line=$(echo "$output" | sed -n '2p')

# Log the cutoffs
echo "$elbow_line" | tee -a "$log_file"
echo "$tenth_pct_line" | tee -a "$log_file"

echo "Pipeline completed successfully for ${prefix} at $(date)" | tee -a "$log_file"

# Tidy up intermediate files
# Cleanup intermediate files
rm -f \
  ${prefix}_mmclst_cluster.tsv \
  ${prefix}_mmclst_all_seqs.fasta \
  ${prefix}_i20_c50.fa \
  ${prefix}_mmclst_i95_c75_cluster.tsv \
  ${prefix}_mmclst_i95_c75_all_seqs.fasta \
  ${prefix}_mmclst_i95_c75_rep_seq.fasta \
  ${prefix}_mmclst_rep_seq.fasta \
  ${prefix}_mmclst_i95_c75.faa.mafft \
  ${prefix}_hmmsearch.out \
  ${prefix}_bitscores.txt
