# Define paths
RAWDATA_DIR = "../../data/rawData/"
VIEW_DIR = "../../data/Post_Traductionnal_Modifications_view/"
FUNCTIONS_DIR = '../Functions/'

# Import modules
import pandas as pd
import numpy as np
import sys

sys.path.insert(1, FUNCTIONS_DIR)
from script.Functions.functions import *

###############################################
# Load and format data
###############################################

# Load dbPTM data
COLNAMES_dbPTM = ["UniProtAC", "PTM_type"]
PTM_df = pd.read_csv(RAWDATA_DIR + "dbPTM_PPIproteins.txt",
                     sep='\t', header=None, names=COLNAMES_dbPTM)

# Create a count table from the dbPTM df
PTM_table = createTable(PTM_df, col_to_indexes="UniProtAC",
                        col_to_features="PTM_type")

# Export
PTM_table.to_csv(VIEW_DIR + "Post_Traductionnal_Modifications.txt", sep="\t",
                 encoding="utf-8", index=True)
