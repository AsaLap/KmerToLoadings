#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Antoine Laporte - 2025

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import re
import seaborn as sns
from matplotlib import patches
from scipy.stats import chisquare
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def clustering(data, title=None, save=None, show=True, hue=None, hue_legend=None, figsize=(10,5), dpi=150, method="ward", metric="euclidean"):
    """
    Function to do and draw a clustering of lines
    :param data: pandas DataFrame, transposed dataframe with values to cluster
    :param title: str, graph title
    :param save: str, save path to save figure
    :param show: bool, show plot or not
    :param hue: list[dict, series], dictionary of color palette and series of correspondance between sample and color palette
    :param hue_legend: str, title legend of the hue variables
    :param figsize: tuple of int, size of figure
    :param dpi: int, quality of image to save
    :param method: the method of clustering, see scipy.cluster.hierarchy.linkage
    :param metric: metric of clustering, see scipy.cluster.hierarchy.linkage
    :return: None
    """
    fig, ax = plt.subplots(figsize=figsize)
    linkage_data = linkage(data, method=method, metric=metric)
    dendrogram(linkage_data, labels=data.index, ax=ax)
    plt.title(title)
    plt.xticks(rotation=90)
    if hue:
        for label in ax.get_xticklabels():
            text = label.get_text()
            label.set_color(hue[0][hue[1][text]])
        patchList = []
        #Add color to legend corresponding to line names
        for key in hue[0]:
            data_key = patches.Patch(color=hue[0][key], label=key)
            patchList.append(data_key)
        plt.legend(handles=patchList, title=hue_legend)
    if save:
        plt.savefig(save+".png", dpi=dpi, bbox_inches='tight')
    if show:
        plt.show()


def extract_domaine(name):
    first_split = name.split("_")
    second_split = first_split[-1].split("-")
    if len(second_split) == 1 or len(second_split) == 2:
        return 'Clone'
    else:
        return second_split[0]


def extract_row(name):
    return name.split("_")[-2]

def extract_type_chardonnay(name):
    return extract_type(name, "python_chardonnay_no_vintage.csv")

def extract_type_pinot(name):
    return extract_type(name, "python_pinots_no_vintage.csv")

def extract_type(name, file):
    """Function to extract type from name
    :param name: str, name of the line/accession
    :param file: str, name of the file in which type information is stored
    """
    finesse_df = pd.read_csv(file)[["ID", "Line", "Finesse", "Type", "Degree_mean", "HL_Ha_mean"]]
    id = name.split("_")[3]
    try:
        return finesse_df.loc[finesse_df['Line'] == id, 'Type'].iloc[0]
    except IndexError:
        if len(id.split("-")) == 1 or len(id.split("-")) == 2:
            return "Clone"
        else:
            print("Name not found", name)
            return "Unknown"


def extraire_groupe_chenin(name):
    """Function to extract group from which Chenin comes from. Function by Gautier Sarah
    :param name: str - name of a chenin
    :return: str - name of the group found
    """
    match = re.match(r".*CHENIN(.*)$", name)
    return match.group(1) if match else name
    # Changing pacbio chenin which doesn't go in the correct group
    # groupes = chenin_1_percent_pondered.index.to_series().apply(extraire_groupe) # .replace("BPacBio", "B")


def levenshtein_select(list_kmer, atcg_per, len_kmer=21, p_value=0.05):
    """
        Function to select kmer based on given genome ATCG composition.
        :param list_kmer: List of kmers
        :param atcg_per: dict of float - percentage of ATCG found in the genome
        :param len_kmer: int - length of kmers to be found in data
        :param p_value: float - threshold for rejection of H0 hypothesis (rejecting independence), higher value = less kmers selected
        # :param eco_mode: boolean - parameter to limit the number of dictionary and data gathered to keep it straight for the expected result, use False for debugging
        :return: list of selected kmers only if eco_mode = True, lists of selected and unselected kmers, dictionary of overall bases and kmer bases
    """
    list_selected_kmer = []
    for kmer in list_kmer:
        g_obs = kmer.count('G')
        c_obs = kmer.count('C')
        t_obs = kmer.count('T')
        a_obs = kmer.count('A')
        g_expected = atcg_per["g"] * len_kmer / 100
        c_expected = atcg_per["c"] * len_kmer / 100
        t_expected = atcg_per["t"] * len_kmer / 100
        a_expected = atcg_per["a"] * len_kmer / 100
        # performe chi² test between numbers of bases inside kmer and expected number based on genome information provided
        if chisquare([g_obs, c_obs, t_obs, a_obs], f_exp=[g_expected, c_expected, t_expected,
                                                          a_expected]).pvalue > p_value:  # take those who are dependent and do not respect H0 hypothesis
            list_selected_kmer.append(kmer)
        elif chisquare([c_obs, g_obs, a_obs, t_obs], f_exp=[g_expected, c_expected, t_expected,
                                                            a_expected]).pvalue > p_value:  # changing C in G and T in A:
            list_selected_kmer.append(kmer)
    return list_selected_kmer

# def levenshtein_kmer_selection(list_kmer, gc_percentage=None, len_kmer=21, p_value=0.05):
#     # eco_mode = True
#     """
#     Function to select kmer based on composition.
#     :param list_kmer: list of kmers
#     :param gc_percentage: float - percentage of GC found in the genome, if None, proportion is calculated on given genome. Info: Chenin: 27.1
#     :param len_kmer: int - length of kmers to be found in data
#     :param p_value: float - threshold for rejection of H0 hypothesis (rejecting independence), higher value = less kmers selected
#     # :param eco_mode: boolean - parameter to limit the number of dictionary and data gathered to keep it straight for the expected result, use False for debugging
#     :return: list of selected kmers only if eco_mode = True, lists of selected and unselected kmers, dictionary of overall bases and kmer bases
#     """
#     #Chenin: https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR11574165&display=metadata
#
#     if not gc_percentage:
#         list_all_bases = [i for x in list_kmer for i in x]
#         len_list_all_bases = len(list_all_bases)
#         dict_kmer_base = {'c': list_all_bases.count('C'),
#                           'g': list_all_bases.count('G'),
#                           't': list_all_bases.count('T'),
#                           'a': list_all_bases.count('A'),
#                           'c_proportion': list_all_bases.count('C') / len_list_all_bases,
#                           'g_proportion': list_all_bases.count('G') / len_list_all_bases,
#                           't_proportion': list_all_bases.count('T') / len_list_all_bases,
#                           'a_proportion': list_all_bases.count('A') / len_list_all_bases}
#         g_expected = dict_kmer_base['g_proportion'] * len_kmer
#         c_expected = dict_kmer_base['c_proportion'] * len_kmer
#         t_expected = dict_kmer_base['t_proportion'] * len_kmer
#         a_expected = dict_kmer_base['a_proportion'] * len_kmer
#         print("G expected by kmer (calculated from kmer) :", g_expected)
#         print("C expected by kmer (calculated from kmer) :", c_expected)
#         print("T expected by kmer (calculated from kmer) :", t_expected)
#         print("A expected by kmer (calculated from kmer) :", a_expected)
#     else:
#         g_expected = c_expected = gc_percentage * len_kmer / 100 / 2
#         t_expected = a_expected = (len_kmer - (g_expected + c_expected)) / 2
#         print("G expected by kmer (from given value) :", g_expected)
#         print("C expected by kmer (from given value) :", c_expected)
#         print("T expected by kmer (from given value) :", t_expected)
#         print("A expected by kmer (from given value) :", a_expected)
#     #--- ECO MODE ---
#     # if eco_mode: #memory economic calculation
#     list_selected_kmer = []
#     for kmer in list_kmer:
#         g_obs = kmer.count('G')
#         c_obs = kmer.count('C')
#         t_obs = kmer.count('T')
#         a_obs = kmer.count('A')
#         #performe chi² test between numbers of bases inside kmer and expected number based on genome information provided
#         chi = chisquare([g_obs, c_obs, t_obs, a_obs], f_exp=[g_expected, c_expected, t_expected, a_expected])
#         if chi.pvalue > p_value: # take those who are dependent and do not respect H0 hypothesis
#             list_selected_kmer.append(kmer)
#     return list_selected_kmer
    # #---GREEDY MODE---
    # else: #memory greedy calculation with huge dict return
    #     dict_kmer_base = {}
    #     list_selected_kmer = []
    #     list_unselected_kmer = []
    #     for kmer in list_kmer:
    #         dict_kmer_base[kmer] = {'c': kmer.count('C'),
    #                                 'g': kmer.count('G'),
    #                                 't': kmer.count('T'),
    #                                 'a': kmer.count('A'),
    #                                 'c_proportion': kmer.count('C') / len_kmer,
    #                                 'g_proportion': kmer.count('G') / len_kmer,
    #                                 't_proportion': kmer.count('T') / len_kmer,
    #                                 'a_proportion': kmer.count('A') / len_kmer}
    #         g_obs = dict_kmer_base[kmer]['g_proportion'] * len_kmer
    #         c_obs = dict_kmer_base[kmer]['c_proportion'] * len_kmer
    #         t_obs = dict_kmer_base[kmer]['t_proportion'] * len_kmer
    #         a_obs = dict_kmer_base[kmer]['a_proportion'] * len_kmer
    #         #performe chi² test between a proportion of bases inside kmer and expected proportions calculated from all kmers
    #         dict_kmer_base[kmer]["chi"] = chisquare([g_obs, c_obs, t_obs, a_obs], f_exp=[g_expected, c_expected, t_expected, a_expected])
    #         if dict_kmer_base[kmer]["chi"].pvalue > p_value: # take those who are dependent and do not respect H0 hypothesis
    #             list_selected_kmer.append(kmer)
    #         else:
    #             list_unselected_kmer.append(kmer)
    #     return list_selected_kmer, list_unselected_kmer, dict_kmer_base


#------ OLD FUNCTION ------
# def levenshtein_kmer_selection(list_kmer, len_kmer=21, p_value=0.05, eco_mode=True):
#     """
#     Function to select kmer based on composition.
#     :param list_kmer: list of kmers
#     :param len_kmer: int - length of kmers to be found in data
#     :param p_value: float - threshold for rejection of H0 hypothesis (rejecting independence), higher value = less kmers selected
#     :param eco_mode: boolean - parameter to limit the number of dictionary and data gathered to keep it straight for the expected result, use False for debugging
#     :return: list of selected kmers only if eco_mode = True, lists of selected and unselected kmers, dictionary of overall bases and kmer bases
#     """
#     list_all_bases = [i for x in list_kmer for i in x]
#     len_list_all_bases = len(list_all_bases)
#     dict_kmer_base = {'c': list_all_bases.count('C'),
#                  'g': list_all_bases.count('G'),
#                  't': list_all_bases.count('T'),
#                  'a': list_all_bases.count('A'),
#                  'c_proportion': list_all_bases.count('C') / len_list_all_bases,
#                  'g_proportion': list_all_bases.count('G') / len_list_all_bases,
#                  't_proportion': list_all_bases.count('T') / len_list_all_bases,
#                  'a_proportion': list_all_bases.count('A') / len_list_all_bases}
#     # doing calculus once to avoid repetition
#     c_expected = dict_kmer_base['c_proportion']*len_kmer
#     g_expected = dict_kmer_base['g_proportion']*len_kmer
#     t_expected = dict_kmer_base['t_proportion']*len_kmer
#     a_expected = dict_kmer_base['a_proportion']*len_kmer
#     #--- ECO MODE ---
#     if eco_mode: #memory economic calculation
#         list_selected_kmer = []
#         for kmer in list_kmer:
#             #performe chi² test between proportion of bases inside kmer and expected proportions calculated from all kmers
#             chi = stats.chisquare([kmer.count('C'), kmer.count('G'), kmer.count('T'), kmer.count('A')],
#                                   f_exp=[c_expected, g_expected, t_expected, a_expected])
#             if chi.pvalue > p_value: # take those who are dependent and do not respect H0 hypothesis
#                 list_selected_kmer.append(kmer)
#         return list_selected_kmer
#     #---GREEDY MODE---
#     else: #memory greedy calculation with huge dict return
#         dict_kmer_base = {}
#         list_selected_kmer = []
#         list_unselected_kmer = []
#         for kmer in list_kmer:
#             dict_kmer_base[kmer] = {'c': kmer.count('C'),
#                                     'g': kmer.count('G'),
#                                     't': kmer.count('T'),
#                                     'a': kmer.count('A'),
#                                     'c_proportion': kmer.count('C') / len_kmer,
#                                     'g_proportion': kmer.count('G') / len_kmer,
#                                     't_proportion': kmer.count('T') / len_kmer,
#                                     'a_proportion': kmer.count('A') / len_kmer}
#             #performe chi² test between proportion of bases inside kmer and expected proportions calculated from all kmers
#             dict_kmer_base[kmer]["chi"] = stats.chisquare([dict_kmer_base[kmer]['c'],
#                                                            dict_kmer_base[kmer]['g'],
#                                                            dict_kmer_base[kmer]['t'],
#                                                            dict_kmer_base[kmer]['a']],
#                                                           f_exp=[c_expected,
#                                                                  g_expected,
#                                                                  t_expected,
#                                                                  a_expected])
#             if dict_kmer_base[kmer]["chi"].pvalue > p_value: # take those who are dependent and do not respect H0 hypothesis
#                 list_selected_kmer.append(kmer)
#             else:
#                 list_unselected_kmer.append(kmer)
#         return list_selected_kmer, list_unselected_kmer, dict_kmer_base


# def levenshtein_kmer_reselection(dict_kmer_base, list_kmer, p_value):
#     """
#     Function to reselect kmers with another p_value if necessary
#     :param dict_kmer_base: dictionary output of "levenshtein_kmer_selection" with eco_mode=False
#     :param list_kmer: list of kmers
#     :param p_value: float - threshold for rejection of H0 hypothesis (rejecting independence), higher value = less kmers selected
#     :return: list of selected kmers
#     """
#     list_selected_kmer = []
#     for kmer in list_kmer:
#         if dict_kmer_base[kmer]["chi"].pvalue > p_value:
#             list_selected_kmer.append(kmer)
#     return list_selected_kmer


def network(values_dataframe, title=None, save=None, node_labeling=True, edge_labeling=True, show=True, node_size=800,
            font_size_nodes=10, font_size_edge=8, seed=42, figsize=(15, 15), dpi=150, layout="spring", metric="jaccard",
            **kwargs):
    """
    Function to draw a network graph of a pandas dataframe with networkx.
    :param values_dataframe: transposed pandas dataframe, values to graph with column name as kmer IDs
    :param title: str, graph title
    :param layout: str, "spring" (default), "kamada", "circular" or "spectral"
    :param metric: str, "jaccard" (default), "euclidean" or else (see sklearn.metrics.pairwise)
    :param save: str, path to save figure
    :param node_labeling: bool, toggle the labeling of nodes
    :param edge_labeling: bool, toggle the labeling of nodes
    :param show: bool, show plot or not
    :param node_size: int, size of nodes
    :param font_size_nodes: int, size of node font
    :param font_size_edge:  int, size of node edge
    :param seed: int, default = 42, for reproducibility
    :param figsize: list of tuple, default = (10,10)
    :param dpi: int - quality of image to save
    :param kwargs: parameters for 'create_coloring' function = dict_of_things_to_color and series of value \
    corresponding to dict_keys to the length of pandas dataframe (values_dataframe), ex: \
    [{'AFS': "#435F55",'B': "#EE55AA",}, ['AFS', 'AFS', 'B', 'B', 'AFS', 'B']] \
    will color nodes names AFS and B differently in a dataframe of size 6.
    :return: None
    """

    ### --- Étape 1 : matrice binaire (présence/absence de k-mers)
    # Supposons que df_binaire a les échantillons en lignes et les k-mers en colonnes
    # Si ce n'est pas le cas, transpose : df_binaire = df_binaire.T

    ### --- Étape 2 : calcul de la distance entre échantillons
    X = values_dataframe.to_numpy()
    dist_matrix = pairwise_distances(X, metric=metric)
    samples = values_dataframe.index

    ### --- Étape 3 : construire un graphe complet pondéré
    g = nx.Graph()
    for i, sample_i in enumerate(samples):
        for j in range(i + 1, len(samples)):
            g.add_edge(sample_i, samples[j], weight=dist_matrix[i][j])

    ### --- Étape 4 : calcul du Minimum Spanning Tree (MSN)
    mst = nx.minimum_spanning_tree(g)

    ### --- Étape 4,5 : customization of node color
    def create_coloring(factors, palette):
        """
        :param factors: series of factors to color and corresponding samples
        :param palette: dict of color palette for each unique element in 'factors'
        :return: a list of colors
        """
        dict_group_colors = dict(factors.map(palette))
        nodes = g.nodes()
        return [dict_group_colors.get(i) for i in nodes]

    if kwargs:
        coloring = create_coloring(**kwargs)
    else:
        coloring = "tab:blue"

    ### --- Étape 5 : dessin du graphe
    plt.figure(figsize=figsize)
    # if sys.version_info.minor >=10: # match case not implemented before 3.10
    #     match layout:
    #         case "kamada":
    #             pos = nx.kamada_kawai_layout(mst)
    #         case "spring":
    #             pos = nx.spring_layout(mst, seed=seed)
    #         case "circular":
    #             pos = nx.circular_layout(mst)
    #         case "spectral":
    #             pos = nx.spectral_layout(mst)
    #         case _:
    #             pos = nx.spring_layout(mst)
    # else:
    if layout=="kamada":
        pos = nx.kamada_kawai_layout(mst)
    if layout=="spring":
        pos = nx.spring_layout(mst, seed=seed)
    if layout=="circular":
        pos = nx.circular_layout(mst)
    if layout=="spectral":
        pos = nx.spectral_layout(mst)
    else:
        pos = nx.spring_layout(mst)
    nx.draw(mst, pos, with_labels=node_labeling, node_color=coloring, node_size=node_size, font_size=font_size_nodes)

    # Edge labeling
    if edge_labeling:
        nx.draw_networkx_edge_labels(mst, pos, edge_labels=nx.get_edge_attributes(mst, 'weight'),
                                     font_size=font_size_edge)

    ### --- Étape 6 : Titre et sauvegarde
    if title:
        plt.title(title)
    if save:
        plt.savefig(save, dpi=dpi)
    if show:
        plt.show()
    return g


def pca_computing(data, n_components=5):
    """
    :param data: dataframe with values to compute pca on. DataFrame must be untransposed (e.g. samples (k-mer 1, k-mer 2....) by lines as index, variables (Chenin220, Chenin278...) by columns)
    :param n_components: number of components to calculate
    :return: pandas dataframe with pca values for 5 principal components and index of data as index
    """
    pipeline = make_pipeline(StandardScaler(), PCA(n_components=n_components))
    data_pca = pipeline.fit_transform(data.T)

    return pd.DataFrame(data_pca, index=data.T.index,
                        columns=['PC' + str(i) for i in range(1, n_components + 1)]), pipeline


def pca_loadings_selection(pca_data, index, nb_loadings=50):
    """
    Function to get the variables with the most significant loadings of a PCA.
    :param pca_data: scikit-learn PCA dataset, untransposed
    :param index: index to put back on the new dataframe = index of the original dataframe on which PCA is computed, untransposed
    :param nb_loadings: int - number of the best loadings to select (negative and positive of every ax)
    :return: list of variables with the most significant loadings
    """
    #TODO: search a way to maximize the different values (>50 here = all same values in 1 percent kmer Chenin dataset)
    #Calculating loadings
    loadings = pd.DataFrame(pca_data.components_.T * np.sqrt(pca_data.explained_variance_.T),
                            index=index,
                            columns=[i + 1 for i in range(len(pca_data.explained_variance_.T))])
    #NOTE: "first" is used because if values are the same, it takes them all. Here, at 38, PC2 are all the same values,
    #making a non-homologous df that causes an error when transposing.
    highest_loadings = pd.DataFrame(np.array([loadings[c].nlargest(nb_loadings, keep='first').index.values for c in loadings]).T)
    highest_loadings.columns = ["PC" + str(i + 1) for i in highest_loadings.columns.values.tolist()]
    lowest_loadings = pd.DataFrame(np.array([loadings[c].nsmallest(nb_loadings, keep='first').index.values for c in loadings]).T)
    lowest_loadings.columns = ["PC" + str(i + 1) for i in lowest_loadings.columns.values.tolist()]
    list_highest_pca_loadings_variables = [i for x in [list(highest_loadings[i].values) for i in highest_loadings.columns.values] for i in x]
    list_lowest_pca_loadings_variables = [i for x in [list(lowest_loadings[i].values) for i in lowest_loadings.columns.values] for i in x]
    list_selected_pca_loading_variables = list(set(list_lowest_pca_loadings_variables + list_highest_pca_loadings_variables))
    return list_selected_pca_loading_variables


def pca_scatterplot_seaborn(df_pca, pca_pipeline, save=None, title=None, hue_legend=None, show=True, axes=("PC1", "PC2"),
                            figsize=(10, 7), dpi=150, text_label=True, **kwargs):
    """
    Function to draw a PCA scatterplot using seaborn scatterplot.
    :param df_pca: pandas dataframe containing PCA data
    :param pca_pipeline: sklearn.pipeline makepipeline object containing sklearn.preprocessing StandardScaler and sklearn.decomposition PCA results
    :param save: str, save path to save figure
    :param title: str, graph title
    :param hue_legend: str, legend on the graph about the hue variables
    :param show: bool, show plot or not
    :param axes: tuple of PCA axes = ("PC1","PC2")
    :param figsize: list of tuple = (10,7)
    :param dpi: int - quality of image to save
    :param text_label: str - show labels or not on the plot
    :param kwargs: parameters of seaborn.scatterplot
    :return: array of PCA explained variance
    """

    plt.figure(figsize=figsize)
    sns.scatterplot(data=df_pca, x=axes[0], y=axes[1], **kwargs)

    if text_label:
        for i, row in df_pca.iterrows():
            plt.text(row[axes[0]] + 0.2, row[axes[1]], i, fontsize=8)
    plt.xlabel(
        f"{axes[0]} ({pca_pipeline[1].explained_variance_ratio_[int(axes[0][2]) - 1] * 100:.1f} %)")  # pca_pipeline[1] = pca
    plt.ylabel(
        f"{axes[1]} ({pca_pipeline[1].explained_variance_ratio_[int(axes[1][2]) - 1] * 100:.1f} %)")  # pca_pipeline[1] = pca
    if title:
        plt.title(title)
    if hue_legend:
        plt.legend(title=hue_legend)
    plt.grid(True)
    plt.tight_layout()
    if save:
        plt.savefig(save + ".png", dpi=dpi)
    if show:
        plt.show()
    return pca_pipeline[1].explained_variance_ratio_


def pca_screeplot_seaborn(pca, save=None, show=True, figsize=(10,7), sum_var=0.8, dpi=150):
    """
    Function to draw the scree plot of a PCA.
    :param pca: pca matrix
    :param figsize: tuple - figure size
    :param sum_var: float - where to draw the line of cumulative variance
    :param save: str - path to save graph
    :param show: bool, show plot or not
    :param dpi: int - quality of image to save
    :return: none
    """
    fig, ax = plt.subplots(figsize=figsize)
    pc_numbers = np.arange(pca.n_components_) + 1
    sns.barplot(x=pc_numbers, y=pca.explained_variance_ratio_, linewidth=2, color='#005b96', ax=ax)
    plt.plot(range(0, len(pc_numbers)), pca.explained_variance_ratio_.cumsum(), marker='o', linestyle='--',
             color='#383434')
    plt.grid(axis="y")
    plt.axhline(y=sum_var, color='#03c03c')
    plt.title('Scree Plot')
    plt.xlabel('Principal Component')
    plt.ylabel('Variance Explained')

    if save:
        plt.savefig(save + "_screeplot.png", dpi=dpi)
    if show:
        plt.show()


def remove_specific(data, threshold=1, get_removed_values=False):
    """
    Function to remove values specific of individuals
    :param data: pandas DataFrame with individuals in columns
    :param threshold: int - threshold number of individuals to consider the value as specific. Default = 1, will remove everything present only in 0 or 1 individual.
    :param get_removed_values: True if you want the list of removed values from index to be returned as second argument
    :return: updated pandas DataFrame
    """
    binary_df = data.where(data == 0, 1) # passing dataframe binary to process to a simple removing of 0 values
    no_specific_df = data.loc[(binary_df != 0).sum(axis=1) > threshold, :] # remove if not present in more than threshold value
    if get_removed_values:
        return no_specific_df, data.index.difference(no_specific_df.index) #way of getting the list of removed values
    return no_specific_df
