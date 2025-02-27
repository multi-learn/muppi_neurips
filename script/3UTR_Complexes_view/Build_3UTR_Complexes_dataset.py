# Import du path
RAWDATA_DIR = "../../data/rawData/"
VIEW_DIR = "../../data/3UTR_Complexes_view/"
FUNCTIONS_DIR = '../Functions/'

# Import librairies
import pandas as pd
import numpy as np
import re
import sys

sys.path.insert(1, FUNCTIONS_DIR)
from script.Functions.functions import *

###############################################
# Load datas
###############################################

# Load the list of PPI_network_proteins
with open(RAWDATA_DIR + "PPInetwork_proteins.txt", "r") as f:
    PPInetwork_proteins = f.read().split('\n')
    PPInetwork_proteins.remove('')

# Load the list of 3UTR complexes and change the colnames
UTRcomplex_df = pd.read_csv(RAWDATA_DIR + "list_3UTRcomplexes.txt",
                            sep='\t', header=0)

# Rename df columns
COL_COMPLEX = ["Nascent", "Intermediate", "RBP"]
UTRcomplex_df.columns = COL_COMPLEX

###############################################
# Format datas
###############################################

# Creat a count table with proteins as entries, and role in complex as columns
zeros_table = np.zeros((len(PPInetwork_proteins), 3))
UTRcomplex_table = pd.DataFrame(zeros_table, index=PPInetwork_proteins,
                                columns=COL_COMPLEX)

# Count the occurence of each proteins in each columns and put the number in the
# table created earlier
for col in COL_COMPLEX:
    for prot in UTRcomplex_df[col].unique():
        UTRcomplex_table.loc[prot, col] = sum(UTRcomplex_df[col] == prot)

# Create to table columns of zeroes
UTRcomplex_table["Inter_linked"] = 0
UTRcomplex_table["RBP_linked"] = 0

# For each protein predicted as nascent in complex
for nas in UTRcomplex_df.Nascent.unique():
    # List the different intermediates and RBP linked to them
    inter_linked = UTRcomplex_df[UTRcomplex_df.Nascent == nas].Intermediate.unique()
    rbp_linked = UTRcomplex_df[UTRcomplex_df.Nascent == nas].RBP.unique()

    # Then stock their number in the table
    UTRcomplex_table.loc[nas, "Inter_linked"] = len(inter_linked)
    UTRcomplex_table.loc[nas, "RBP_linked"] = len(rbp_linked)

###############################################
# 3UTR altrenatives
###############################################

# Load the alternative 3UTR data for all PPI proteins
ALTER_3UTR_COL = ["RBP", "Nascent", "multi_or_single_UTR", "binding_class"]
alter_3UTR = pd.read_csv(RAWDATA_DIR + "RBP-ALL-nascent-complexes_3UTR-alternative.txt",
                         sep='\t', names=ALTER_3UTR_COL, header=None)
# Load the alternative 3UTR data for only EMF proteins, that had been manually
# verified
emf_alter_3UTR = pd.read_csv(RAWDATA_DIR + "RBP-EMF-nascent-complexes_3UTR-alternative.txt",
                             sep='\t', names=ALTER_3UTR_COL, header=None)

# Replace the rows of EMF in the df for all proteins by the df of verified data
alter_3UTR = alter_3UTR[~alter_3UTR.Nascent.isin(emf_alter_3UTR.Nascent)]
alter_3UTR = pd.concat([alter_3UTR, emf_alter_3UTR], axis=0)

# Load the CrossReference table GeneSymbol / UniprotAC
# HACK: Replace the file GeneSymbole_rainet.csv by GeneSymbolCrossReference.csv
# TODO: Find the correct file
crossref = pd.read_csv(RAWDATA_DIR + "GeneSymbolCrossReference.csv",
                       sep=',', header=0)

# Identify the nascents which are identified as multi_utr and single_utr
wrong = alter_3UTR[["Nascent", "multi_or_single_UTR"]].drop_duplicates()
wrong = wrong[wrong.Nascent.duplicated()].Nascent

# Delete rows where those nascents are single_utr
alter_3UTR = alter_3UTR[~((alter_3UTR.Nascent.isin(wrong)) &
                          (alter_3UTR.multi_or_single_UTR == "single_utr"))]

# Translate the nascent column from GeneSymbol to UniProtAC
alter_3UTR = alter_3UTR.merge(crossref, left_on="Nascent",
                              right_on="uniprotGeneSymbol")
alter_3UTR = alter_3UTR[alter_3UTR.protein_id.isin(UTRcomplex_df.Nascent)]

# Create a counting table with proteins in index and binding_classes in columns
alter_3UTR_table = createTable(alter_3UTR,
                               col_to_indexes="protein_id",
                               col_to_features="binding_class")

# Add rows full of zeroes for the PPI proteins that are not in the table
alter_3UTR_table = MissingExamples(alter_3UTR_table, zero_or_nan="zero")

# Join these table to the UTRcomplex_table
UTRcomplex_table = UTRcomplex_table.join(alter_3UTR_table)

# Load the fasta file of alternative 3UTR
with open(RAWDATA_DIR + "Ensemblv90_3UTR-sequences.txt", "r") as f:
    UTRfasta_list = f.read().split('\n>')

# For each fasta
ENSG_list = []
for i in range(len(UTRfasta_list)):
    # Extract the protein id and stock it in a list
    ENSG_list.append("ENSG" + re.search("ENSG(.*)\n", UTRfasta_list[i]).group(1))

ENSG_list = pd.Series(ENSG_list)

# Create an empty df to put the ENSG ids and their alternative 3UTR number
UTRfasta_df = pd.DataFrame(ENSG_list.value_counts(), columns=["count"])
UTRfasta_df = UTRfasta_df.reset_index()

# Load the cross reference id table
# FIXME: The file rainetCrossReference.csv is missing
crossref_ENSG = pd.read_csv(RAWDATA_DIR + "rainetCrossReference.csv",
                            sep="\t", header=0)

# In the id mapping table, select only the GeneSymbol to UniProtAC translations
crossref_ENSG = crossref_ENSG[(crossref_ENSG.sourceDB == "Ensembl") &
                              (crossref_ENSG.protein_id.isin(PPInetwork_proteins))][["protein_id", "crossReferenceID"]]

# Merge the id mapping table with the UTRfasta df to translate the ids
UTRfasta_df = UTRfasta_df.merge(crossref_ENSG, left_on="index",
                                right_on="crossReferenceID")

# List the UniprotAC that are translated into several ENSG and delete their
# rows in the df
dupliENSG = UTRfasta_df[UTRfasta_df.protein_id.duplicated()].protein_id.unique()
UTRfasta_df = UTRfasta_df[["protein_id", "count"]][~UTRfasta_df.protein_id.isin(dupliENSG)]

# Join this df to the UTRcomplex_table
UTRcomplex_table = UTRcomplex_table.join(UTRfasta_df.set_index("protein_id"))

# Export the table
UTRcomplex_table.to_csv(VIEW_DIR + "3UTR_Complexes.txt",
                        sep="\t", encoding="utf-8", index=True)

# Create a new table without rows with NaN
UTRcomplex_table_NoNA = UTRcomplex_table[~pd.isnull(UTRcomplex_table).any(axis=1)]

# Export it
UTRcomplex_table_NoNA.to_csv(VIEW_DIR + "3UTR_Complexes_NoNA.txt",
                             sep="\t", encoding="utf-8", index=True)
