# Define paths
RAWDATA_DIR = "../../data/rawData/"
FUNCTIONS_DIR = '../Functions/'

# Import modules
import pandas as pd
import numpy as np
import sys

# sys.path.insert(1, FUNCTIONS_DIR)
from script.Functions.functions import *

#########################################
# Load and Format datas
#########################################

# Load PPI network edge list
PPInetwork = pd.read_csv(RAWDATA_DIR + "PPInetwork_edgelist.csv",
                         sep=",", header=0)

# List the sorted unique proteins in the network
PPIproteins = PPInetwork.interactorA_id._append(
    PPInetwork.interactorB_id).sort_values().unique()

# Export it
with open(RAWDATA_DIR + "PPInetwork_proteins.txt", "w") as f:
    for prot in PPIproteins:
        f.write(prot + "\n")

# Define colnames for uniprot IDmap
IDMAP_COLNAMES = ["UniProtAC", "UniProtID", "GeneID", "RefSeq", "GI", "PDB",
                  "GO", "UniRef100", "UniRef90", "UniRef50", "UniParc", "PIR",
                  "NCBI-taxon", "MIM", "UniGene", "PubMed", "EMBL", "EMBL-CDS",
                  "Ensembl", "Ensembl_TRS", "Ensembl_PRO", "Additional_PubMed"]

# Load uniprot IDmap with new colnames
crossref = pd.read_csv(RAWDATA_DIR + "HUMAN_9606_idmapping_selected.tab",
                       sep="\t", header=None, low_memory=False, names=IDMAP_COLNAMES)

# Export it
crossref.to_csv(RAWDATA_DIR + "HUMAN_9606_idmapping_selected.tab",
                sep="\t", index=False, encoding="utf-8")

# Extract UniProtAC and ENST, and explode ENST into multiple rows
crossref_ENST = explodeValues(crossref,
                              index_col="UniProtAC",
                              col_to_explode="Ensembl_TRS",
                              separator="; ")

# HACK: PPInetwork_proteins is not defined in this script so I'm going to use PPIproteins instead
# Select only PPI proteins, and sort the df
crossref_ENST = crossref_ENST[crossref_ENST.UniProtAC.isin(PPIproteins)  # PPI proteins not PPInetwork proteins I think
                              & (crossref_ENST.Ensembl_TRS != '')]

# Export
crossref_ENST.to_csv(RAWDATA_DIR + "PPIproteins_ENST.txt",
                     sep='\t', encoding="utf-8", index=False)
