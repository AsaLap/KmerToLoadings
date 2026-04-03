#!/bin/sh
# Antoine Laporte 2025
#SBATCH --job-name=KmerToLoadings
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=cpu-dedicated
#SBATCH --account=dedicated-cpu@cirad
#SBATCH --time=00:10:00
mkdir logs
#SBATCH --output logs/Kmer."%j".out
#SBATCH --error logs/Kmer."%j".err

echo "Running on:$SLURM_NODELIST"
echo "Starting date: $(date +%d/%m/%y-%HH%M)"

# Loading parameters from config file
scripts_dir=$(pwd)
. "$scripts_dir"/run_parameters.config

printf "\nParameters:\n"
echo "do_splitting: $do_splitting"
echo "do_pondering: $do_pondering"
echo "do_levenshtein: $do_levenshtein"
echo "do_pca: $do_pca"
echo "do_concat: $do_concat"
echo "matrix: $matrix"
echo "split_factor: $split_factor"
echo "reads_count: $reads_count"
echo "names: $names"
#echo "sep: $sep"
echo "round: $round"
echo "remove_specific: $remove_specific"
echo "a_percentage: $a_percentage"
echo "t_percentage: $t_percentage"
echo "c_percentage: $c_percentage"
echo "g_percentage: $g_percentage"
echo "len_kmer: $len_kmer"
echo "p_value: $p_value"
echo "n_components: $n_components"
echo "n_loadings: $n_loadings"
echo "draw_graphs: $draw_graphs"

#Prepare filenames and output_directory
#base_dir=$(dirname "$(realpath "$matrix")")
base_name=$(basename "${matrix}")
name=${base_name%%.*}
directory_output=${scripts_dir}/${name}_tmp
logs_directory=${scripts_dir}/logs

#making output directories
mkdir -p "${directory_output}"
#mkdir -p "${logs_directory}"

if $do_splitting;then
  printf "\n------SPLIT------\n"
  #Counting number of lines to prepare the split
  nb_kmer=$(wc -l "$matrix" | awk '{ print $1 }')
  lines=$((nb_kmer / split_factor + 1)) #+1 top avoid having a file with the rest of the euclidean division

  echo "Original matrix contains $nb_kmer lines"
  echo "Splitting in $split_factor by matrices of $lines lines..."
  START_TIME=$(date +%s)
  split "$matrix" "${directory_output}/${name}_" -l $lines -d --additional-suffix .tsv

  ELAPSED=$(($(date +%s) - START_TIME))
  echo "...splitting done!"
  printf "Splitting time: %s\n\n" "$(date -d@$ELAPSED -u +%H\ hours\ %M\ min\ %S\ sec)"
else
  printf "\nSkipping splitting\n"
fi

# Estimating memory use?: 8 (bytes per float) * nb columns (=individuals) * nb lines (=kmer or line by matrix)
# Example for 1/1000 matrix PN : 8 * 610 * 561404 = 2739651520 = 2.7GB
# Two matrix will be loaded/created in the python file so previous result * 2 : 5.4GB
# But first matrix is int and not float so 4 bytes per int instead of 8 bytes per float.
# We can truncate memory need as first matrix needs half the memory of the second: (2.7/2 = 1.35) so 2.7 + 1.35 = 4.05.
# Let's put 5GB, which is first matrix * 2, truncated.
# But this calculs seems not sufficient : Ex : chardonnay with roughly 400MB files needed 10GB...
nb_individuals=$(wc -l "$names" | awk '{ print $1 }')
memory_estimation=$((8 * nb_individuals * nb_kmer * 2))
echo "Individuals: $nb_individuals"
echo "Minimum memory estimation: $memory_estimation"

if $do_pondering;then
  printf "\n------PONDER------\n"
  echo "Pondering each submatrix by coverage values..."
  START_TIME=$(date +%s)
  for matrix in $(ls "$directory_output"/* | grep -v pondered | grep -v pca);
  do
#    echo "$matrix"
    sbatch "$scripts_dir"/launch_ponder.sh \
    "$scripts_dir" \
    "$matrix" \
    "$names" \
    "$reads_count" \
    "$round";
  done
  ELAPSED=$(($(date +%s) - START_TIME))
  echo "...pondering done!"
  printf "Pondering time: %s\n\n" "$(date -d@$ELAPSED -u +%H\ hours\ %M\ min\ %S\ sec)"
else
  printf "\nSkipping pondering\n"
fi

if "$do_levenshtein";then
  printf "\n------LEVENSHTEIN------\n"
  echo "Selecting kmer with Levenshtein process for each submatrix..."
  START_TIME=$(date +%s)
  for matrix in $(ls "$directory_output"/* | grep pondered | grep -v leven | grep -v pca);
  do
#    echo "$matrix";
    sbatch "$scripts_dir"/launch_levenshtein.sh \
    "$scripts_dir" \
    "$matrix" \
    "$remove_specific" \
    "$a_percentage" \
    "$t_percentage" \
    "$c_percentage" \
    "$g_percentage" \
    "$sep" \
    "$len_kmer" \
    "$p_value";
  done
  ELAPSED=$(($(date +%s) - START_TIME))
  echo "...Levenshtein selection done!"
  printf "Levenshtein time: %s\n\n" "$(date -d@$ELAPSED -u +%H\ hours\ %M\ min\ %S\ sec)"
else
  printf "\nSkipping levenshtein\n"
fi

if "$do_pca";then
  printf "\n------PCA------\n"
  echo "Selecting pca loadings for each submatrix..."
  START_TIME=$(date +%s)
  for matrix in "$directory_output"/*pondered_levenshtein.tsv;
  do
#    echo "$matrix";
    sbatch "$scripts_dir"/launch_pca.sh \
    "$scripts_dir" \
    "$matrix" \
    "$n_components" \
    "$n_loadings" \
    "$draw_graphs"
  done
  ELAPSED=$(($(date +%s) - START_TIME))
  echo "...PCA loadings selection done!"
  printf "PCA time: %s\n\n" "$(date -d@$ELAPSED -u +%H\ hours\ %M\ min\ %S\ sec)"
else
  printf "\nSkipping PCA\n"
fi

if "$do_concat";then
  rm -f "$directory_output"/"$name"_pca_selection_"$split_factor"_chunks_"$n_loadings"_loadings_"$n_components"_components.tsv
  tail -q -n +2 "$directory_output"/*"$n_loadings"_pca_loadings_"$n_components"_components.tsv >> "$directory_output"/"$name"_pca_selection_"$split_factor"_chunks_"$n_loadings"_loadings_"$n_components"_components.tsv
#  str(n_loadings) + "_pca_loadings_" + str(n_components) + "_components.tsv"
  echo "Global PCA matrix: ${name}_pca_selection_${split_factor}_chunks_${n_loadings}_loadings_${n_components}_components.tsv"
fi
