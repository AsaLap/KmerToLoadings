#!/bin/sh
# Antoine Laporte 2025
#SBATCH --job-name=ponder
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-normal
#SBATCH --time=02:00:00
#SBATCH -o logs/ponder."%j".out
#SBATCH -e logs/ponder."%j".err

echo "Running on:$SLURM_NODELIST"

module purge
module load python/3.7.2

scripts_dir=$1
matrix=$2
names=$3
reads_count=$4
round=$5

echo "$matrix"

python "$scripts_dir"/cluster_ponder_by_coverage.py \
--matrix "$matrix" \
--names "$names" \
--reads_count "$reads_count" \
--round "$round"
