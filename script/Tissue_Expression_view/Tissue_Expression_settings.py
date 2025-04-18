# Define paths
RAWDATA_DIR = "../../data/rawData/"
VIEW_DIR = "../../data/Tissue_Expression_view/"
FUNCTIONS_DIR = '../Functions/'

# Import modules
import pandas as pd
import numpy as np
import sys

sys.path.insert(1, FUNCTIONS_DIR)
from script.Functions.functions import *

#########################################
# Load datas
#########################################

# Load HPA data
RNAtissue_df = pd.read_csv(RAWDATA_DIR + "rna_tissue_consensus.tsv", sep='\t',
                           header=0)

# Load cross reference ID table
crossref = pd.read_csv(RAWDATA_DIR + "HUMAN_9606_idmapping_selected.tab", sep="\t",
                       header=0, low_memory=False)

# Load the PPI proteins list
with open(RAWDATA_DIR + "PPInetwork_proteins.txt", "r") as f:
    PPInetwork_proteins = f.read().split('\n')
    PPInetwork_proteins.remove('')

#########################################
# Merge and format the df
#########################################

# Create a df with only the uniprotAC and separated ENSG of the crossreference table
CR_uniprotAC_ENSG = explodeValues(crossref, "UniProtAC", "Ensembl", "; ")
# Keep only the UniProtAC presents in the PPI network proteins list
CR_uniprotAC_ENSG = CR_uniprotAC_ENSG[CR_uniprotAC_ENSG.UniProtAC.isin(
    PPInetwork_proteins)]

# CHANGED: Ignore the annotation version
CR_uniprotAC_ENSG["Ensembl"] = CR_uniprotAC_ENSG["Ensembl"].apply(lambda s: s.split(".")[0])

# Merge the df with the crossref to get the UniprotACs
RNAtissue_df = RNAtissue_df.merge(CR_uniprotAC_ENSG, left_on="Gene",
                                  right_on="Ensembl")

# Create a table from this df, with tissues as columns, UniprotAC as indexes,
# and nTPM (NX) as values
# HACK: Think nTPM is NX (previously used as values)
RNAtissue_df = pd.pivot_table(RNAtissue_df, values="nTPM", index="UniProtAC",
                              columns="Tissue")

# Export it
RNAtissue_df.to_csv(VIEW_DIR + "Tissue_Expression.txt", sep="\t",
                    encoding="utf-8", index=True)

#########################################
# Remove missing values
#########################################

# Create a new df without rows with NaN
RNAtissue_NoNA_df = RNAtissue_df[~pd.isnull(RNAtissue_df).any(axis=1)]

# Export it
RNAtissue_NoNA_df.to_csv(VIEW_DIR + "Tissue_Expression_NoNA.txt", sep="\t",
                         encoding="utf-8", index=True)
