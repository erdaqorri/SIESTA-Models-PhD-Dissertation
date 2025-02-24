# Runs the AlphaFold2 predictions on the GPU partition
# Assumes Parafold is installed and the conda env activated

#!/bin/bash

#SBATCH -A project_title
#SBATCH --job-name=sample_name 
#SBATCH --output=sample_name_%j_.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=55GB
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

#SBATCH --array=1-n

# Define the input file containing protein names
protein_names_file="n.txt"

# Read protein file names into an array
file_names=()
while IFS= read -r line; do
  file_names+=("$line")
done < "$protein_names_file"


# Use the array task ID to select the corresponding file name
name=${file_names[SLURM_ARRAY_TASK_ID-1]}
fasta_file="${name}.fasta"


bash /hpath/to/alphafold/dir/run_alphafold.sh \
-d /path/to/alphafold/database/monomer_af2_db \
-u $CUDA_VISIBLE_DEVICES \
-o /path/to/output/dir \
-m model_1,model_2,model_3,model_4,model_5 \
-p monomer_ptm \
-i "path/to/input/fasta/files/$fasta_file" \
-t 2023-09-20 \
