# Define paths
RAWDATA_DIR = "../../data/rawData/"
VIEW_DIR = "../../data/Linear_Motifs_view/"
FUNCTIONS_DIR = '../Functions/'

# Import modules
import pandas as pd
import numpy as np
import re
import sys

sys.path.insert(1, FUNCTIONS_DIR)
from script.Functions.functions import *

###############################################
# Load datas
###############################################

# Load the PPI proteins list
with open(RAWDATA_DIR + "PPInetwork_proteins.txt", "r") as f:
    PPInetwork_proteins = f.read().split('\n')
    PPInetwork_proteins.remove('')

# Gets the list of predicted elm's
ELM_SLIMProb = pd.read_csv(RAWDATA_DIR + "ELM_interactome-proteins_predicted-motifs.txt",
                           sep='\t', header=0)

# Get elm classes data with their probabilities
ELM_classes = pd.read_csv(RAWDATA_DIR + "elm_classes.txt", sep='\t', header=0)

# Load the ELM classes list that have been experimentaly observed in human
with open(RAWDATA_DIR + "elm_instances_TP_HomoSapiens.txt", "r") as f:
    verified_ELM_classes = f.read()
    verified_ELM_classes = re.sub("\"", "", verified_ELM_classes).split("\n")
    verified_ELM_classes.remove('')

# Load the protein disorder predictions
PPI_IUPred = pd.read_csv(RAWDATA_DIR + "PPIproteins_iupred.txt",
                         sep='\t', header=None,
                         names=["protein_id", "windows_pred", "IUP"])

# Load the protein lengths data
protLen = pd.read_csv(RAWDATA_DIR + "protein_length.csv",
                      sep=',', header=0)

###############################################
# Format datas
###############################################

# Replace 2 class names
ELM_SLIMProb.Motif = ELM_SLIMProb.Motif.replace("DOC_CYCLIN_RxL_1", "DOC_CyclinA_RxL_1")
ELM_SLIMProb.Motif = ELM_SLIMProb.Motif.replace("LIG_SH2_GRB2", "LIG_SH2_GRB2like")

# Select the ELM classes that have been experimentaly observed in human
# and have an occurence probability under 0.001
ELM_classes = ELM_classes[(ELM_classes.ELMIdentifier.isin(verified_ELM_classes)) &
                          (ELM_classes.Probability < 0.001)]
ELM_SLIMProb = ELM_SLIMProb[ELM_SLIMProb.Motif.isin(ELM_classes.ELMIdentifier)]

# Removes lines in the df that do not correspond to the interactome and elm with
# a disorder (IUP) less than 0.4
ELM_SLIMProb = ELM_SLIMProb[(ELM_SLIMProb.AccNum.isin(PPInetwork_proteins)) &
                            (ELM_SLIMProb.IUP > 0.4)]

# Create an occurrence table with prot in index and ipr_id in column
ELM_table = createTable(ELM_SLIMProb, col_to_features="Motif",
                        col_to_indexes="AccNum")

# Add for each column in the table a density column
count_cols = ELM_table.columns
density_cols = list(ELM_table.columns + "_density")
ELM_table = ELM_table.reindex(columns=ELM_table.columns.tolist() + density_cols)

# Add missing proteins in indexes and fill their features with zeroes
ELM_table = MissingExamples(ELM_table, zero_or_nan="zero")

# Add a protein lenght column in the table
ELM_table = ELM_table.join(protLen.set_index("uniprotAC"))

# Fill the density columns by the count columns divided by the protein length
for density_col in density_cols:
    count_col = density_col[:-8]
    ELM_table[density_col] = ELM_table[count_col].values / ELM_table.proteinLength.values

# Add the disorder values of each protein in a new table column
PPI_IUPred = pd.pivot_table(PPI_IUPred, index="protein_id", values="IUP",
                            columns="windows_pred")
PPI_IUPred = PPI_IUPred.rename(columns={"long": "Disorder_long",
                                        "short": "Disorder_short"})
ELM_table = ELM_table.join(PPI_IUPred, how="inner").sort_index()

# Convert th count columns of the df in integer type
ELM_table[count_cols] = ELM_table[count_cols].astype(int)

# Export
ELM_table.to_csv(VIEW_DIR + "Linear_Motifs_SlimProb.txt", sep='\t',
                 header=True, index=True, encoding="utf-8")
