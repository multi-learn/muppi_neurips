
# Define paths
RAWDATA_DIR = "../../data/rawData/"

# Import modules
import pandas as pd
import numpy as np

###############################################
# Load datas
###############################################

# Load the PPI proteins list
with open( RAWDATA_DIR + "PPInetwork_proteins.txt", "r") as f:
    PPInetwork_proteins = f.read().split( '\n')
    PPInetwork_proteins.remove( '')

# Load EMFs
emf = pd.read_csv( RAWDATA_DIR + "human_EMF_list.tsv", "\t", header = 0,
                   usecols = [0])

# Load the modules list in PPI network
Clusters = pd.read_csv( RAWDATA_DIR + "PPI_modules.csv", ',', header = 0,
                        index_col = None)

###############################################
# Format datas
###############################################

# Selects PPI proteins (there are some errors due to dataset versions)
Clusters = Clusters[ Clusters.protein_id.isin( PPInetwork_proteins)]

# List multi clustered proteins
mf = Clusters[ Clusters.protein_id.duplicated()].protein_id.unique()

# Create a df with a PPI proteins column
labels = pd.DataFrame( PPInetwork_proteins)
# Create a labels column in this df and store for each protein the mono
# clustered, multi clustered or EMF label.
labels[ "labels"] = np.where( labels[ 0].isin( mf), "multi_clustered",
                              "mono_clustered")
labels[ "labels"] = np.where( labels[ 0].isin( emf.UniprotKB_AC), "EMF",
                              labels.labels)

# Rename column names
labels = labels.rename({ 0:"protein_id"}, axis = 1)

# Export
labels.to_csv( RAWDATA_DIR + "EMF_labels.txt", sep = '\t',
               header = True, encoding = "utf-8", index = False)
