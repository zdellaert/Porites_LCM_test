#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --ntasks=1 --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --time=2:00:00 --qos=short
#SBATCH --error=../scripts/outs_errs/%x_error.%j #if your job fails, the error report will be put in this file
#SBATCH --output=../scripts/outs_errs/%x_output.%j #once your job is completed, any final job report comments will be put in this file
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --no-requeue

species=$1
genome=$2
genome_path=$3
gff_path=$4
makeindex=$5

scratch_dir="/scratch4/workspace/zdellaert_uri_edu-shared_TimeSeries/Porites_LCM_test"
data_dir="${scratch_dir}/trimmed/"

genome_index_dir="${scratch_dir}/STAR_indexes/${genome}"
out_dir="${scratch_dir}/aligned/${species}_${genome}"

mkdir -p "${genome_index_dir}"
mkdir -p "${out_dir}"

cd "${scratch_dir}"

# load modules 
module load uri/main STAR/2.7.11b-GCC-12.3.0

# genome index generation
if [ "${makeindex}" = "T" ]; then
  STAR --runMode genomeGenerate \
      --runThreadN 20 \
      --genomeDir "${genome_index_dir}" \
      --genomeFastaFiles "${genome_path}" \
      --sjdbGTFfile "${gff_path}" \
      --sjdbGTFtagExonParentTranscript Parent \
      --genomeSAindexNbases 13
fi

trimmed=( "${data_dir}"*"${species}"*"_flexbar_1.fastq" )

# run star

for R1_file in "${trimmed[@]}"; do

  # extract sample name
  sample_name=$(basename "${R1_file}" "_flexbar_1.fastq")

  # define R2 file
  R2_file="${data_dir}${sample_name}_flexbar_2.fastq"

  STAR --runMode alignReads \
       --genomeDir "${genome_index_dir}" \
       --runThreadN 10 \
       --readFilesIn "${R1_file}" "${R2_file}" \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMunmapped Within \
       --outSAMattributes Standard \
       --outFileNamePrefix "${out_dir}/${sample_name}_" \
       --quantMode GeneCounts
done
