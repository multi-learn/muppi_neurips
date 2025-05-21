"""
This script displays information about a MuPPI release HDF5 file.

"""

import h5py
import numpy as np
import pandas as pd

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

# Check if your file is stored here
DATA_DIR = "../../data/dataset_compilation_to_hdf5/"


def analyze_hdf5_file(file_name):
    print(f"\nAnalyse du fichier : {file_name}")
    print("-" * 50)
    file_path = DATA_DIR + file_name
    with h5py.File(file_path, 'r') as f:
        nb_views = f['Metadata'].attrs['nbView']
        nb_classes = f['Metadata'].attrs['nbClass']
        dataset_length = f['Metadata'].attrs['datasetLength']
        print(f"Number of views : {nb_views}")
        print(f"Number of classes : {nb_classes}")
        print(f"Total number of proteins : {dataset_length}")

        # Labels distribution
        labels = f['Labels'][:]
        label_counts = np.bincount(labels)
        label_percentages = 100 * label_counts / len(labels)
        label_names = [name.decode('utf-8') if isinstance(name, bytes) else name for name in f['Labels'].attrs['names']]

        label_df = pd.DataFrame({
            "Classe": label_names,
            "Number of samples": label_counts,
            "Percent": label_percentages
        })

        print("\nLabels distribution :")
        print(label_df)

        # Views informations
        view_data = []
        for i in range(nb_views):
            view = f[f'View{i}']
            view_name = view.attrs['name']
            data = view[:]
            nb_features = data.shape[1]
            nb_samples = data.shape[0]

            proteins_all_missing = np.sum(np.all(data == -1, axis=1))
            missing_percentage = (proteins_all_missing / nb_samples) * 100
            proteins_with_values = nb_samples - proteins_all_missing

            view_data.append([i, view_name, nb_features, proteins_with_values, missing_percentage])

        view_df = pd.DataFrame(view_data, columns=[
            "View ID", "View name", "# Features",
            "# Proteins", "# Missing values (-1 on all the features)"
        ])

        print("\nView information")
        print(view_df)


if __name__ == "__main__":
    file_name = "MuPPI_2025_14_fullview_light_EMF.hdf5"
    analyze_hdf5_file(file_name)
