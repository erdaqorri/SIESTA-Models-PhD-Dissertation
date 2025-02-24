# Sample code utilized to run the structure predictions in the CPU-only partition of Komondor
# Assumes Parafold is installed and the conda env activated

#!/bin/bash

#SBATCH -A project_title
#SBATCH --job-name=sample_name
#SBATCH --output=sample_name_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --mem=35GB
#SBATCH --partition=cpu

#SBATCH --array=1-n%n

# Define the input file containing protein names
protein_names_file="n.txt"

# Read protein file names into an array
file_names=()
while IFS= read -r line; do
  file_names+=("$line")
done < "$protein_names_file"


# Use the array task ID to select the corresponding file name
name=${file_names[SLURM_ARRAY_TASK_ID-1]}
fasta_file_name="${name}.fasta"

bash /path/to/parafold/dir/run_alphafold.sh \
-d /path/to/alphafold/databases/monomer_af2_db \
-o /path/to/output/dir \
-p monomer_ptm \
-i "/path/to/input/fasta/files/$fasta_file_name" \
-t 2023-09-20 \
-m model_1 \
-f
