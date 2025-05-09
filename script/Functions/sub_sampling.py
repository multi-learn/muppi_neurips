# Script used for sub-sampling
# It generates n_subsampling datasets with n_emf*ratio randomly chosen
# multi_clustered and all the EMFs, and a config file for each of these datasets.

import h5py
import numpy as np
import os
import yaml

random_state = np.random.RandomState(42)

# Define the number of sub-sampled datasets
n_subsampling = 10

# The ratio of multi_clustered to EMFs
ratio = 1

# Define the paths
base_dataset_file_path = "../../../summit/config_files/MuPPI_2025_14view_EMF.hdf5"
output_path = "../../data/dataset_compilation_to_hdf5/balanced/"
config_base_path = "../../config_files/multi_vs_emf_base.yml"
config_saving_path = "../../config_files/balanced/"

# Create directories
if not os.path.isdir( output_path[ :-1]):
    os.mkdir( output_path[ :-1])
if not os.path.isdir( config_saving_path[ :-1]):
    os.mkdir( config_saving_path[ :-1])

# Load the base configuration file
with open(config_base_path, "r") as stream_in:
    base_config = yaml.safe_load(stream_in)

# For each sub sampled dataset
for rs in range(n_subsampling):

    # Define a new random_state
    random_state = np.random.RandomState(rs+random_state.randint(1,100))

    # Load the base dataset
    database_unbalanced = h5py.File(base_dataset_file_path, 'r')

    # Get the labels
    labels = database_unbalanced["Labels"][...]

    # Get the indices of the EMF labelled proteins
    emf_indices = np.where(labels == 2 )[0]
    n_emf = emf_indices.shape[0]

    # Get the indices of the multi_clustered labelled proteins
    neg_indices = np.where(labels == 1)[0]

    # Choose randomly n_emf*ratio multi_clustered proteins
    chosen_neg_indices = random_state.choice(neg_indices, int(n_emf*ratio), replace=False)

    # Get the concatenation of the EMFs and the chosen multi_clustered
    chosen_indices = np.concatenate((emf_indices, chosen_neg_indices))

    # Reset the labels
    labels[emf_indices] = 1
    labels[chosen_neg_indices] = 0

    # Build a configuration file for the curent sub sampled dataset
    base_config["name"] = ["MuPPI_2025_14view_EMF_sub_{}".format(rs)]
    with open(os.path.join(config_saving_path,
                           "multi_vs_emf_subsampling_{}.yml".format(rs)),
              "w") as stream_out:
        yaml.dump(base_config, stream_out)

    # Build the HDF5 file for the sub-sampled dataset, compatible wit SuMMIT
    balanced_database = h5py.File(output_path + "MuPPI_14view_EMF_sub_" + str(rs) + ".hdf5", 'w')

    labels_dset = balanced_database.create_dataset("Labels", data=labels[chosen_indices])
    for att, val in dict(database_unbalanced["Labels"].attrs).items():
        labels_dset.attrs[att] = val
    labels_dset.attrs["names"] = ["multi_clustered".encode(), "EMF".encode()]

    for i in range(database_unbalanced["Metadata"].attrs["nbView"]):
        dset = balanced_database.create_dataset("View"+str(i), data=database_unbalanced["View"+str(i)][...][chosen_indices, :])
        for att, val in dict(database_unbalanced["View"+str(i)].attrs).items():
            dset.attrs[att] = val

    meta_data_grp = balanced_database.create_group("Metadata")

    for att, val in dict(database_unbalanced["Metadata"].attrs).items():
        meta_data_grp.attrs[att] = val

    meta_data_grp.attrs["datasetLength"] = len(chosen_indices)
    example_ids_dset =  meta_data_grp.create_dataset("example_ids", data=database_unbalanced["Metadata"]["example_ids"][...][chosen_indices], dtype=np.dtype('S10'))

    balanced_database.close()
    database_unbalanced.close()
