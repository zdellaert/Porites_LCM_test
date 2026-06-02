#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --no-requeue
#SBATCH --mem=16GB
#SBATCH -t 03:59:00 
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --error=../scripts/outs_errs/%x_error.%j #if your job fails, the error report will be put in this file
#SBATCH --output=../scripts/outs_errs/%x_output.%j #once your job is completed, any final job report comments will be put in this file

species=$1
genome=$2
gtf_path=$3

# load required modules
module load uri/main StringTie/2.2.1-GCC-11.2.0

# list and make required directories
scratch_dir="/scratch4/workspace/zdellaert_uri_edu-shared_TimeSeries/Porites_LCM_test"
alignments_dir="${scratch_dir}/aligned/${species}_${genome}"

out_dir="${scratch_dir}/stringtie/${species}_${genome}"

# make the output directory if it does not exist (-p checks for this)
mkdir -p "${out_dir}"

cd "${alignments_dir}"

# call the STAR bam files into an array
bams=(*Aligned.sortedByCoord.out.bam)

for f in "${bams[@]}"; do 
    sample_name=$(echo "$f" | sed -E 's/_Aligned.*//')

    # -p 16 : use 16 cores
    # --rf : library is reverse-forward stranded
    # -e : exclude novel genes
    # -B : create Ballgown input files for downstream analysis
    # -v : enable verbose mode
    # -G : gtf annotation file
    # -A : output name for gene abundance estimate files
    # -o : output name for gtf file

    stringtie -p 16 --rf -e -B -v \
        -G "${gtf_path}" \
        -A "${out_dir}"/"${sample_name}".gene_abund.tab \
        -o "${out_dir}"/"${sample_name}".gtf \
        "$f" #input bam file

    echo "StringTie assembly for seq file ${f}" $(date)
done

