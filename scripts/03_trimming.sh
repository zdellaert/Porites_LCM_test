#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=80GB
#SBATCH -t 03:59:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --error=../scripts/outs_errs/%x_error.%j #if your job fails, the error report will be put in this file
#SBATCH --output=../scripts/outs_errs/%x_output.%j #once your job is completed, any final job report comments will be put in this file

# load modules needed
module load parallel/20240822

# make and define directories needed
data_dir="/project/pi_hputnam_uri_edu/zdellaert/Porites_LCM_test/data_RNA/"
out_dir="/scratch4/workspace/zdellaert_uri_edu-shared_TimeSeries/Porites_LCM_test/trimmed/"
qc_dir="/project/pi_hputnam_uri_edu/zdellaert/Porites_LCM_test/output_RNA/trimmed_qc/"

mkdir -p ${out_dir}
mkdir -p ${qc_dir}

# create an list of fastq files to process
R1_files=(${data_dir}*_R1_001.fastq.gz)

echo "There are ${#R1_files[@]} samples to process"
echo "Starting trimming at $(date)"

# define flexbar function to allow for parallel processing

run_flexbar() {
  R1_file=$1
  data_dir=$2
  out_dir=$3
  qc_dir=$4

  # extract sample name
  sample_name=$(basename "$R1_file" "_R1_001.fastq.gz")

  # define R2 file
  R2_file="${data_dir}${sample_name}_R2_001.fastq.gz"

  flexbar --reads "$R1_file" \
        --reads2 "$R2_file" \
        --stdout-reads \
        --adapters ${data_dir}oligo_fastas/tso_g_wo_hp.fasta \
        --adapter-trim-end LEFT \
        --adapter-revcomp ON \
        --adapter-revcomp-end RIGHT \
        --htrim-left GT \
        --htrim-right CA \
        --htrim-min-length 3 \
        --htrim-max-length 5 \
        --htrim-max-first \
        --htrim-adapter \
        --min-read-length 2 \
        --threads 4 | \
    flexbar \
        --reads - \
        --interleaved \
        --target "${out_dir}${sample_name}_flexbar" \
        --adapters ${data_dir}oligo_fastas/ilmn_20_2_seqs.fasta \
        --adapter-trim-end RIGHT \
        --min-read-length 2 \
        --threads 4

  echo "trimming of "${sample_name}" complete at $(date)"
}

export -f run_flexbar

# run flexbar in parallel
parallel -j 6 run_flexbar {} "$data_dir" "$out_dir" "$qc_dir" ::: "${R1_files[@]}"

# now move onto qc
echo "Starting fastqc on trimmed files at" $(date)

# load modules needed
module load fastqc/0.12.1
module load uri/main
module load MultiQC/1.12-foss-2021b

# create an list of fastq files to process
trimmed_files=(${out_dir}*trim.fastq.gz)

# Run fastqc in parallel
parallel -j 6 "fastqc {} -o "${qc_dir}" && echo 'Processed {}'" ::: "${trimmed_files[@]}"
echo "fastQC done." $(date)

#Compile MultiQC report from FastQC files
echo "Running MultiQC"
cd "${qc_dir}"
multiqc --interactive .

echo "QC of trimmed RNA-seq data complete." $(date)
