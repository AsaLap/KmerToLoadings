#!/bin/sh
# Antoine Laporte 2025
#SBATCH --job-name=levenshtein
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad-normal
#SBATCH -o logs/levenshtein."%j".out
#SBATCH -e logs/levenshtein."%j".err

echo "Running on:$SLURM_NODELIST"

module purge
module load python/3.7.2

scripts_dir=$1
matrix=$2
remove_specific=$3
a_percentage=$4
t_percentage=$5
c_percentage=$6
g_percentage=$7
sep=$8
len_kmer=$9
p_value=${10}

echo "$matrix"

python "$scripts_dir"/cluster_levenshtein.py \
--matrix "$matrix" \
--remove_specific "$remove_specific" \
--a_percentage "$a_percentage" \
--t_percentage "$t_percentage" \
--c_percentage "$c_percentage" \
--g_percentage "$g_percentage" \
--sep "$sep" \
--len_kmer "$len_kmer" \
--p_value "$p_value"
