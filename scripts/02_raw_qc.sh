#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=200GB
#SBATCH -t 12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --error=../scripts/outs_errs/%x_error.%j #if your job fails, the error report will be put in this file
#SBATCH --output=../scripts/outs_errs/%x_output.%j #once your job is completed, any final job report comments will be put in this file

# load modules needed
module load parallel/20240822
module load fastqc/0.12.1
module load uri/main
module load MultiQC/1.12-foss-2021b

#go into directory with raw data (symlinks)
cd /project/pi_hputnam_uri_edu/zdellaert/Porites_LCM_test/data_RNA

#make raw_qc output folder
mkdir -p ../output_RNA/raw_qc/

# Create an list of fastq files to process
files=( *.fastq.gz )

# Run fastqc in parallel
echo "Starting fastqc..." $(date)
parallel -j 6 "fastqc {} -o ../output_RNA/raw_qc/ && echo 'Processed {}'" ::: "${files[@]}"
echo "fastQC done." $(date)

#Compile MultiQC report from FastQC files
cd ../output_RNA/raw_qc/
multiqc --interactive .

echo "Initial QC of RNA-seq data complete." $(date)
