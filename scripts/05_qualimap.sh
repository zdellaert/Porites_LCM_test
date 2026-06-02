#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --time=04:00:00
#SBATCH --error=../scripts/outs_errs/%x_%A_%a_error #if your job fails, the error report will be put in this file
#SBATCH --output=../scripts/outs_errs/%x_%A_%a_output #once your job is completed, any final job report comments will be put in this file
#SBATCH --array=0-41
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --no-requeue

species=$1
genome=$2
gtf_path=$3

# load modules needed
module load qualimap/2.2.1

# list and make required directories
scratch_dir="/scratch4/workspace/zdellaert_uri_edu-shared_TimeSeries/Porites_LCM_test"
alignments_dir="${scratch_dir}/aligned/${species}_${genome}"

qc_dir="/project/pi_hputnam_uri_edu/zdellaert/Porites_LCM_test/output_RNA/alignment_qc/${species}_${genome}"

# make the output directory if it does not exist (-p checks for this)
mkdir -p "${qc_dir}"

# make list of BAM files

bam_files=("${alignments_dir}"/*Aligned.sortedByCoord.out.bam)

# get the BAM for this array task
f="${bam_files[$SLURM_ARRAY_TASK_ID]}"
sample_name=$(basename "$f" | sed -E 's/_Aligned.*//')

echo "Running Qualimap on ${sample_name}..."

	qualimap rnaseq \
    --java-mem-size=8G \
    -gtf "${gtf_path}" \
    -pe \
    --sequencing-protocol non-strand-specific \
    -bam "${f}"  \
    -outdir "${qc_dir}"/"${sample_name}" 

