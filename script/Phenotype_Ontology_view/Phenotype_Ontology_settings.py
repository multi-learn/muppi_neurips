# Define paths
RAWDATA_DIR = "../../data/rawData/"
VIEW_DIR = "../../data/Phenotype_Ontology_view/"
FUNCTIONS_DIR = '../Functions/'

# Import modules
import pandas as pd
import numpy as np
import os
from goatools.obo_parser import GODag
from goatools.godag.go_tasks import get_go2parents
import sys

sys.path.insert(1, FUNCTIONS_DIR)
from script.Functions.functions import *

#########################################
# Load datas
#########################################

# Load HPO data in a df
HPOA_df = pd.read_csv(RAWDATA_DIR + "HPO_genes_to_phenotype.txt",
                      sep='\t', header=0)

# Load the human phenotype DAG
HPO_tree = GODag(RAWDATA_DIR + "hp.obo")

# Load the PPI proteins list
with open(RAWDATA_DIR + "PPInetwork_proteins.txt", "r") as f:
    PPInetwork_proteins = f.read().split('\n')
    PPInetwork_proteins.remove('')

# Load the cross reference table GeneSymbol / protein_id
crossref = pd.read_csv(RAWDATA_DIR + "GeneSymbolCrossReference.csv",
                       sep=',', header=0)

#########################################
# Merge and format the df
#########################################

crossref = crossref[crossref.protein_id.isin(PPInetwork_proteins)]

# Selects the proteins present in the PPInetwork and removes duplicates
HPOA_df = HPOA_df.merge(crossref, left_on="gene_symbol",
                        right_on="uniprotGeneSymbol")

# CHANGED: Using the right column name
HPOA_df = HPOA_df[["protein_id", "hpo_id"]].drop_duplicates()

# Test hypergeometric sur les emf annotées HPO
# with open( RAWDATA_DIR + "238_EMFs.txt", "r") as f:
#     emf = f.read().split( '\n')
# import scipy.stats as stats
# emf_annotated = len( HPOA_df[ HPOA_df.protein_id.isin( emf)].protein_id.unique())
# PPIprot = len( PPInetwork_proteins)
# emf_all = len( emf)
# prot_annotated = len(HPOA_df["protein_id"].unique())
# non_emf_non_annotated = ( PPIprot - emf_all) - ( prot_annotated - emf_annotated)
# 1-stats.hypergeom.cdf( emf_annotated, PPIprot-emf_all, emf_all, prot_annotated)
# prot_annotated/PPIprot*100
# emf_annotated/emf_all*100


# Create a dict containing the hpo_id parents of each hpo_id
HPO_parents = get_go2parents(HPO_tree, relationships={'is_a'})

# Create an occurrence table with prot in index and hpo_id in column.
HPOA_table = createTable(HPOA_df, col_to_features="hpo_id",
                         col_to_indexes="protein_id")

# Define root hpo term
ROOT_HPO_ID = ['HP:0000001']

# For all hpo_id in the graph except the root terms
for hpo_id in HPOA_df.hpo_id.unique():

    # Create an empty hpo_id ancestor dict and store in value the parent terms
    # of the hpo_id with hpo_id as key
    hpo_ancestors = {}
    hpo_ancestors[hpo_id] = list(HPO_parents[hpo_id])

    # As long as the values in the ancestor's dictum do not contain the term root
    while ROOT_HPO_ID not in list(hpo_ancestors.values()):

        # Retrieves in a list the values of the dict that are not present in the
        # keys. Either the terms that we don't yet have the parents of
        new_parents = [x for x in list(hpo_ancestors.values()) if
                       x not in list(hpo_ancestors.keys())]
        # Each of these terms are stored in the dict as keys, and their values
        # are their parent terms
        for parents in new_parents:
            for i in range(len(parents)):
                hpo_ancestors[parents[i]] = list(HPO_parents[parents[i]])

    # Adds in the table the ancestors columns with the annotations of the
    # heir column
    HPOA_table = fillAncestorsTerms(HPOA_table, hpo_id, hpo_ancestors.keys())

# Replace values other than 0 with 1's to remove additions.
HPOA_table[HPOA_table != 0] = 1

# Supersedes redundant columns by keeping their names concatenated in the name
# of the remaining column
HPOA_table = delRedundantCols(HPOA_table)

###############################################
# Export
###############################################

# Add missing proteins in indexes and fill these rows with zeroes
HPOA_table = MissingExamples(HPOA_table, zero_or_nan="zero")

# Convert the df in integer type
HPOA_table = HPOA_table.astype(int)

# Export
HPOA_table.to_csv(VIEW_DIR + "Phenotype_Ontology.txt", sep="\t",
                  encoding="utf-8", index=True)

# Removes uninformative features and exports this light version of df
HPOA_table = CutSmallLeafs(HPOA_table)
HPOA_table.to_csv(VIEW_DIR + "Phenotype_Ontology_light.txt", sep="\t",
                  encoding="utf-8", index=True)
