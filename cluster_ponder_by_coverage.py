#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Antoine Laporte 2025

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Ponder kmer count by overall genome coverage')
parser.add_argument('--matrix', action="store", dest='matrix', required=True, help="Matrix of kmer frequencies, TSV file")
parser.add_argument('--names', action="store", dest='names', required=True, help="Names of lines in same order as they are in the matrix, each name on a new line type file")
parser.add_argument('--reads_count', action="store", dest='reads_count', required=True, help="Reads count value for each line (reverse and forward) with 'line_fr' and 'nb_reads' headers, TSV file")
parser.add_argument('--round', action="store", dest='round', default=2, help="Rounding number of the matrix")

args = parser.parse_args()

# --- 1. Loading data
matrix = pd.read_csv(args.matrix, sep='\t', header=None, index_col=0, engine='python')

# --- 2. Charger les noms de colonnes à partir du fichier texte
with open(args.names, 'r') as f:
    noms_colonnes = [ligne.strip() for ligne in f]

# --- 3. Affecter les noms aux colonnes
matrix.columns = noms_colonnes

# --- 4. Import coverage values dataframe
nb_reads = pd.read_csv(args.reads_count, sep='\t')
nb_reads['line'] =  nb_reads['line_fr'].apply(lambda x:  str(x).split("_R")[0]) # getting line names only
nb_reads = nb_reads.groupby('line').sum() # suming all lines to get R1 and R2 on the same line, line is now index
nb_reads["coverage"] = nb_reads["nb_reads"] * 150 / 500000000 # applying nb reads * 150 (size of a read) / 500 000 000 (size of genome)
print("nb_read", "\n", nb_reads)

# --- 5. Dividing kmer values by mean coverage
matrix_pondered = pd.DataFrame()
for i in nb_reads.index:
    try:
        matrix_pondered[i] = round(matrix[i] / nb_reads.loc[i]["coverage"], int(args.round))
    except KeyError:
        print("Key Error : ", i)
        pass

# --- 6. Transposing to get every line on a row and kmer frequencies by columns and saving
# matrix_pondered = matrix_pondered.T
matrix_pondered.to_csv(args.matrix.split(".")[0] + '_pondered.tsv', sep='\t')

# # --- 7. extracting group names to predict classes from
# # --- CHENIN STYLE
# dataframe_chenin_pondered = dataframe_chenin_pondered.T
# def extraire_groupe(nom):
#     match = re.match(r".*CHENIN(.*)$", nom)
#     return match.group(1) if match else nom
# groupes = pd.DataFrame(dataframe_chenin_pondered.index.to_series().apply(extraire_groupe), columns=["Groupe"])
# groupes.to_csv(args.matrix.split(".")[0] + '_group_list.tsv', sep='\t', index=False)


# --- 8. Remove <1 values if asked
# if args.remove_under_one:
#     dataframe_chenin_pondered_nan = dataframe_chenin_pondered.where(dataframe_chenin_pondered >= 0.5, np.NaN)
#     dataframe_chenin_pondered_nan = dataframe_chenin_pondered_nan.dropna(axis=1)
#     dataframe_chenin_pondered_nan.to_csv(args.matrix.split(".")[0] + '_no_under_one.tsv', sep='\t')
# dataframe_chenin_pondered_nan.head(26)
