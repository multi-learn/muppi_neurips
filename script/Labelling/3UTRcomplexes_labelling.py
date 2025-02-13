
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

# Load the 3UTR complexes list
UTRcomplex_df = pd.read_csv( RAWDATA_DIR + "list_3UTRcomplexes.txt", '\t',
                             header = 0)

###############################################
# Format datas
###############################################

# Create a df with a PPI proteins column
labels = pd.DataFrame( PPInetwork_proteins)

# Stock the list of unique nascents of UTRcomplex_df
nascent = UTRcomplex_df.nascent_uniprotkb_ac.unique()

# Create a labels column in this df and store for each protein the label
# nascent, or not nascent.
labels[ "labels"] = np.where( labels[ 0].isin( nascent), "nascent", "not_nascent")

# Rename column names
labels = labels.rename({ 0:"protein_id"}, axis = 1)

# Export
labels.to_csv( RAWDATA_DIR + "3UTRcomplexes_labels.txt", sep = '\t',
               header = True, encoding = "utf-8", index = False)
