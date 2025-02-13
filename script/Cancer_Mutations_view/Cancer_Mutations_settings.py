
# Define paths
RAWDATA_DIR = "../../data/rawData/"
VIEW_DIR = "../../data/Cancer_Mutations_view/"

# Import modules
from scipy import stats
import pandas as pd
import numpy as np


#########################################
# Load datas
#########################################

# Define column names
Cosmic_COLNAMES = [ "ENST", "Primary_site", "Site_subtype_1", "Site_subtype_2",
                    "Site_subtype_3", "Primary_histology", "Histology_subtype_1",
                    "Histology_subtype_2", "Histology_subtype_3", "SNP",
                    "FATHMM_score", "Mutation_somatic_status"]

# Load cosmic data
Cosmic_df = pd.read_csv( RAWDATA_DIR + "CosmicMutant.tsv", sep = '\t',
                         header = None, names = Cosmic_COLNAMES)

# Load cross reference ID table
crossref_ENST = pd.read_csv( RAWDATA_DIR + "PPIproteins_ENST.txt", "\t",
                             header = 0)

# Load the PPI proteins list
with open( RAWDATA_DIR + "PPInetwork_proteins.txt", "r") as f:
    PPInetwork_proteins = f.read().split( '\n')
    PPInetwork_proteins.remove( '')

#########################################
# Merge and format the df
#########################################

# Replace "NS" strings by numpy nan
Cosmic_df = Cosmic_df.replace( "^NS$", np.nan, regex = True)
# Replace "n" and "y" in SNP column by 0 and 1
Cosmic_df[ "SNP"] = Cosmic_df.SNP.replace({ "y":1, "n":0})

# Merge the df with the crossref to get the UniprotACs
Cosmic_df = Cosmic_df.merge( crossref_ENST, left_on = "ENST",
                             right_on = "Ensembl_TRS")

# Places UniProtACs in the index, removes ENST columns and duplicated rows.
Cosmic_df = Cosmic_df.set_index( "UniProtAC").drop([ "ENST", "Ensembl_TRS"],
                                                     axis = 1).drop_duplicates()

# Export the string version of df
Cosmic_df.to_csv( VIEW_DIR + "Cancer_Mutations.txt", sep = "\t",
                  encoding = "utf-8", index = True)

#########################################
# Numerisation
#########################################

# For all columns containing strings
for str_col in Cosmic_df.drop( ["SNP", "FATHMM_score"], axis = 1).columns:

    # Create a list of unique values in the column
    str_col_values = Cosmic_df[ str_col].unique().tolist()
    # In the column replaces the values with their index in the list
    Cosmic_df[ str_col] = Cosmic_df[ str_col].apply( lambda x: str_col_values.index( x))
    # If the column contained nan, put in the -1's instead
    if np.nan in str_col_values:
        Cosmic_df[ str_col] = Cosmic_df[ str_col].replace( str_col_values.index( np.nan), -1)

# Export numerised version of df
Cosmic_df.to_csv( VIEW_DIR + "Cancer_Mutations_numerised.txt", sep = "\t", encoding = "utf-8", index = True)
