# Ask if the dataset generated will be with full view, so with missing values
try_again = True
while try_again:
    full_version = input("Full version of dataset ? (y/n): ")
    # If the answer is yes, stop the loop and set variables to determine later loops
    if full_version == "y":
        full_version = True
        suffix_NoNA = ""
        suffix_full = "_full"
        try_again = False
    # If the answer is no, stop the loop and set variables in the inverse way
    elif full_version == "n":
        full_version = False
        suffix_NoNA = "_NoNA"
        suffix_full = ""
        try_again = False
    # If the answer is another sting, restart the loop
    else:
        print("You didn't enter \"y\" or \"n\"")

# Ask if the dataset generated will be in light_version
try_again = True
while try_again:
    light_version = input("Light version of dataset ? (y/n): ")
    # If the answer is yes, stop the loop and create a variable suffix to load
    # light version of views
    if light_version == "y":
        suffix_light = "_light"
        try_again = False
    # If the answer is no, stop the loop and create a empty variable suffix to
    # load normal version of views
    elif light_version == "n":
        suffix_light = ""
        try_again = False
    # If the answer is another sting, restart the loop
    else:
        print("You didn't enter \"y\" or \"n\"")

# Ask which label to use
labels_list = ["EMF", "COMPLEXES"]
try_again = True
while try_again:
    which_labels = input("Do you want to use EMF or 3'UTR complexes labels ? (EMF/complexes): ")
    which_labels = which_labels.upper()
    # If the answer is in the labels list, stop the loop
    if which_labels in labels_list:
        try_again = False
    # If it's not, restart the loop
    else:
        print("You didn't enter one of the two answers requested")

# Import of paths
RAWDATA_DIR = "../../data/rawData/"
DATA_DIR = "../../data/"
FUNCTIONS_DIR = '../Functions/'

# Import librairies
import pandas as pd
import numpy as np
import h5py
import sys

from script.Functions.functions import *

###############################################
# Load datas
###############################################

# Load the list of PPI_network_proteins
with open(RAWDATA_DIR + "PPInetwork_proteins.txt", "r") as f:
    PPInetwork_proteins = f.read().split('\n')
    PPInetwork_proteins.remove('')

# Load the labels
if which_labels == labels_list[0]:
    labels = pd.read_csv(RAWDATA_DIR + "EMF_labels.txt", sep='\t', header=0,
                         index_col=0)

elif which_labels == labels_list[1]:
    labels = pd.read_csv(RAWDATA_DIR + "3UTRcomplexes_labels.txt", sep='\t',
                         header=0, index_col=0)

# Load the all view in a dict
# We add a suffix in the files names for the light version
df_dict = {}

df_dict["PPInetwork_topology"] = pd.read_csv(DATA_DIR +
                                             'PPInetwork_topology_view/PPInetwork_topology.txt',
                                             sep='\t', header=0, index_col=0)
df_dict["Subcell_Location"] = pd.read_csv(DATA_DIR +
                                          'Subcell_Location_view/Subcell_Location_numerised.txt',
                                          sep='\t', header=0, index_col=0)
df_dict["Tissue_Expression"] = pd.read_csv(DATA_DIR +
                                           'Tissue_Expression_view/Tissue_Expression{}.txt'.format(suffix_NoNA),
                                           sep='\t', header=0, index_col=0)
df_dict["SAGE_PPInetwork"] = pd.read_csv(DATA_DIR +
                                         'PPInetwork_Embedding_view/SAGE_PPInetwork.txt',
                                         sep=' ', header=None, index_col=0)
df_dict["Gene_Ontology_BP"] = pd.read_csv(DATA_DIR +
                                          'Gene_Ontology_view/Gene_Ontology_BP{}.txt'.format(suffix_light),
                                          sep='\t', header=0, index_col=0)
df_dict["Gene_Ontology_CC"] = pd.read_csv(DATA_DIR +
                                          'Gene_Ontology_view/Gene_Ontology_CC{}.txt'.format(suffix_light),
                                          sep='\t', header=0, index_col=0)
df_dict["Gene_Ontology_MF"] = pd.read_csv(DATA_DIR +
                                          'Gene_Ontology_view/Gene_Ontology_MF{}.txt'.format(suffix_light),
                                          sep='\t', header=0, index_col=0)
df_dict["BP_PPInetwork_embed"] = pd.read_csv(DATA_DIR +
                                             'GO_PPInetwork_view/GO-BP_PPInetwork_embed.txt',
                                             sep=' ', header=None, index_col=0)
df_dict["CC_PPInetwork_embed"] = pd.read_csv(DATA_DIR +
                                             'GO_PPInetwork_view/GO-CC_PPInetwork_embed.txt',
                                             sep=' ', header=None, index_col=0)
df_dict["Phenotype_Ontology"] = pd.read_csv(DATA_DIR +
                                            'Phenotype_Ontology_view/Phenotype_Ontology{}.txt'.format(suffix_light),
                                            sep='\t', header=0, index_col=0)
df_dict["Protein_Domains"] = pd.read_csv(DATA_DIR +
                                         'Protein_Domains_view/Protein_Domains{}.txt'.format(suffix_light),
                                         sep='\t', header=0, index_col=0)
df_dict["PTM"] = pd.read_csv(DATA_DIR +
                             'Post_Traductionnal_Modifications_view/Post_Traductionnal_Modifications.txt',
                             sep='\t', header=0, index_col=0)

# If the labels are nascent/not_nascent in 3utr complexes, don't load the
# 3UTR_Complexes view
if which_labels != labels_list[1]:
    df_dict["3UTR_Complexes"] = pd.read_csv(DATA_DIR +
                                            '3UTR_Complexes_view/3UTR_Complexes{}.txt'.format(suffix_NoNA),
                                            sep='\t', header=0, index_col=0)

df_dict["Linear_Motifs"] = pd.read_csv(DATA_DIR +
                                       'Linear_Motifs_view/Linear_Motifs_SlimProb.txt',
                                       sep='\t', header=0, index_col=0)

# For full version of dataset, load the Cancer_Mutations view
# if full_version:
#    df_dict["Cancer_Mutations"] = pd.read_csv(DATA_DIR +
#                                              'Cancer_Mutations_view/Cancer_Mutations_numerised.txt',
#                                              sep='\t', header=0, index_col=0)

###############################################
# Prepare datasets
###############################################

# Change the hdf5 file name in function of the labels and the light version
h5_path = DATA_DIR + "dataset_compilation_to_hdf5/MuPPI_2025_{}{}view{}_{}.hdf5".format(
    len(df_dict), suffix_full, suffix_light, which_labels)

# Create the hdf5 file where the df will be stocked
dataset_file = h5py.File(h5_path, 'w')

# Create a copy of the PPI protein list
index_intesect = PPInetwork_proteins

# For not full version, keep only indexes that are present in all views
if not full_version:
    # Keep in the protein list only the proteins present in all df indexes
    for key in df_dict.keys():
        index_intesect = [protein for protein in index_intesect if
                          protein in df_dict[key].index]

    # For each df, select only the common indexes and sort them
    for key in df_dict.keys():
        df_dict[key] = df_dict[key][df_dict[key].index.isin(index_intesect)].sort_index()

    # For not full version, in the label df select the common indexes defined earlier
    labels = labels[labels.index.isin(index_intesect)].sort_index()

# For full version add missing proteins in each views and fill those rows by NaN
else:
    for key in df_dict.keys():
        df_dict[key] = MissingExamples(df_dict[key], zero_or_nan="nan")

        # CHANGED: Convert all the index to string
        df_dict[key].index = df_dict[key].index.astype(str)
        df_dict[key] = df_dict[key].sort_index()

    # Sort the label index
    labels.sort_index()

###############################################
# labels
###############################################

# List the labels names
if which_labels == labels_list[0]:
    labels_names = ['mono_clustered', 'multi_clustered', 'EMF']
elif which_labels == labels_list[1]:
    labels_names = ['not_nascent', 'nascent']

# In th labels column of the df, replace strings by their index in the
# labels names list
labels.labels = labels.labels.apply(lambda x: labels_names.index(x))

# Turn labels integers into a numpy array
labels = labels.labels.values

# Stock the labels array in a hdf5 file dataset "Labels"
labels_dataset = dataset_file.create_dataset("Labels", shape=labels.shape,
                                             data=labels)

# In attribute of this dataset, stock the labels names
labels_dataset.attrs["names"] = [label_name.encode() if not isinstance(label_name, bytes)
                                 else label_name for label_name in labels_names]

###############################################
# Put the datasets in hdf5 file
###############################################

# Set the view number at zero
view_nb = 0

# For each view in the dict
for key in df_dict.keys():
    # Create an array of the current view values and delete the view df to
    # free memory
    view_data = df_dict[key].values[:, :].astype(float)
    df_dict[key] = []

    # Replace np.nan by -1 because they're not supported in hdf5 format
    view_data = np.where(np.isnan(view_data), -1, view_data)

    # Create a hdf5 file dataset with the view array, its shape and
    # its view number
    view_dataset = dataset_file.create_dataset("View" + str(view_nb),
                                               view_data.shape, data=view_data)
    # In attribute of this dataset, stock the current dict key, which is the
    # view name, and a boolean to say if the view is sparse
    view_dataset.attrs["name"] = key
    view_dataset.attrs["sparse"] = False
    # Iterate the view number of 1
    view_nb += 1

# Create hdf5 file dataset "Metadata"
meta_data_grp = dataset_file.create_group("Metadata")

# In attribute of this dataset, stock the number of view, number of class and
# number of exemples
meta_data_grp.attrs["nbView"] = len(df_dict)
meta_data_grp.attrs["nbClass"] = len(labels_names)
meta_data_grp.attrs["datasetLength"] = len(labels)

# In metadata group, create a dataset which contain the view indexes
index_intesect = np.array(index_intesect).astype(np.dtype("S10"))
meta_data_grp.create_dataset("example_ids", index_intesect.shape,
                             data=index_intesect, dtype=np.dtype("S10"))

# Close the hdf5 file
dataset_file.close()
