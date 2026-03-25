#!/bin/sh
# Antoine Laporte 2025
#SBATCH --job-name=pca
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-normal
#SBATCH -o logs/pca."%j".out
#SBATCH -e logs/pca."%j".err

echo "Running on:$SLURM_NODELIST"

module purge
module load python/3.7.2

scripts_dir=$1
matrix=$2
n_components=$3
n_loadings=$4

echo "$matrix"

python "$scripts_dir"/cluster_pca.py \
--matrix "$matrix" \
--n_components "$n_components" \
--n_loadings "$n_loadings"

echo "Done"
