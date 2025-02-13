
# Define paths
RAWDATA_DIR = "../../data/rawData/"
VIEW_DIR = "../../data/Subcell_Location_view/"
FUNCTIONS_DIR = '../Functions/'

# Import modules
import pandas as pd
import numpy as np
import sys
sys.path.insert(1, FUNCTIONS_DIR)
from functions import *

#########################################
# Load datas
#########################################

# Load HPA data
subloc_df = pd.read_csv( RAWDATA_DIR + "subcellular_location.tsv", sep = '\t')

# Load the cross reference table
crossref = pd.read_csv( RAWDATA_DIR + "HUMAN_9606_idmapping_selected.tab", "\t",
                        header = 0, low_memory = False)

# Load the PPI proteins list
with open( RAWDATA_DIR + "PPInetwork_proteins.txt", "r") as f:
    PPInetwork_proteins = f.read().split( '\n')
    PPInetwork_proteins.remove( '')

#########################################
# Merge df
#########################################

# Create a df with only the uniprotAC and separated ENSG of the crossreference table
CR_uniprotAC_ENSG = explodeValues( crossref, "UniProtAC", "Ensembl", "; ")
# Keep only the UniProtAC presents in the PPI network proteins list
CR_uniprotAC_ENSG = CR_uniprotAC_ENSG[ CR_uniprotAC_ENSG.UniProtAC.isin(
                                       PPInetwork_proteins)]

# Merge the df with the crossref to get the UniprotACs
subloc_df = subloc_df.merge( CR_uniprotAC_ENSG, left_on = "Gene",
                             right_on = "Ensembl")


#########################################
# Format the df
#########################################

# List the confidence levels in the df colnames
confidence_levels = [ "Enhanced", "Supported", "Approved", "Uncertain",
                     "Extracellular location"]

# From the df make a coordinates list with 3 columns, protein ids,
# subcell compartment and confidence level.
# For each selected columns, which are confidence levels
subloc_coord = pd.DataFrame()
for col in confidence_levels:
     # In a tempory df with protein ids and current column,
     # separate the concatenated location into multiple rows
     tmp = explodeValues( subloc_df, "UniProtAC", col, ";")
     # Create a confidence column with the current colname at every rows
     tmp = tmp.assign( confidence = col)
     # Rename the current colname
     tmp = tmp.rename({ col:"subcell_compartment"}, axis = 1)
     # Remove rows without subcell compartment
     tmp = tmp[ tmp.subcell_compartment != ""]
     # Append this df to a final one
     subloc_coord = pd.concat([ subloc_coord, tmp])

# Define confidence score for confidence levels
confidence_scores = [ 1, 0.85, 0.7, 0.5, 1]

# In the coord list, replace confidence levels by numerical scores
subloc_coord = subloc_coord.replace( confidence_levels, confidence_scores)

# Create a table from the coord list, with subcell_compartment as columns,
# UniprotAC as indexes, and confidence score as values
subloc_df = pd.pivot_table( subloc_coord, values = "confidence",
                            index = "UniProtAC", columns = "subcell_compartment")

# Replace nan by zeroes
subloc_df = subloc_df.replace({ np.nan : 0})

# Export numerised version
subloc_df.to_csv( VIEW_DIR + "Subcell_Location_numerised.txt", sep = "\t",
                  encoding = "utf-8", index = True)

# Replace values by initial strings
confidence_scores.append( 0)
confidence_levels.append( np.nan)
subloc_df = subloc_df.replace( confidence_scores, confidence_levels)


# Export the verbose version of df
subloc_df.to_csv( VIEW_DIR + "Subcell_Location.txt", sep = "\t",
                   encoding = "utf-8", index = True)
