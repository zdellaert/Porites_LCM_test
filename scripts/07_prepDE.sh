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

# load required modules
module load uri/main StringTie/2.2.1-GCC-11.2.0

# list and make required directories
scratch_dir="/scratch4/workspace/zdellaert_uri_edu-shared_TimeSeries/Porites_LCM_test"
stringtie_dir="${scratch_dir}/stringtie/${species}_${genome}"
out_dir="/project/pi_hputnam_uri_edu/zdellaert/Porites_LCM_test/output_RNA/count_matrices"
script_dir="/project/pi_hputnam_uri_edu/zdellaert/Porites_LCM_test/scripts"

# make the output directory if it does not exist (-p checks for this)
mkdir -p "${out_dir}"

# move into stringtie directory
cd "${stringtie_dir}"

# make input file
for filename in *.gtf; do
    sample_name=$(basename "$filename" .gtf)

    echo $sample_name $PWD/$filename
done > listGTF.txt

#Compile the gene count matrix
python "${script_dir}"/prepDE.py3 -g "${out_dir}"/"${species}"_"${genome}"_gene_count_matrix.csv -i listGTF.txt

echo "Gene count matrix compiled." $(date)
