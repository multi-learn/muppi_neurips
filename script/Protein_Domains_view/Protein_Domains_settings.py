
# Define paths
RAWDATA_DIR = "../../data/rawData/"
VIEW_DIR = "../../data/Protein_Domains_view/"
FUNCTIONS_DIR = '../Functions/'

# Import modules
import pandas as pd
import numpy as np
import re
import sys
sys.path.insert(1, FUNCTIONS_DIR)
from functions import *

###############################################
# Load datas
###############################################

# Retrieves the Interpro term tree
with open( RAWDATA_DIR + "ParentChildTreeFile.txt", "r") as f:
           domains = f.read().split( '\n')

# Load Interpro data
COLNAMES_IPR = [ "UniProtAC", "IPR_ID"]
IPR_df = pd.read_csv( RAWDATA_DIR + "PPIprotein2ipr.txt", '\t',
                      header = None, names = COLNAMES_IPR)

###############################################
# Construire l'arbre des termes
###############################################

# Remove the domain names in strings to keep only the IDs
for i in range( len( domains)):
    domains[ i] = re.sub( "::.*::", "", domains[ i])

# For each terms
domain_parents = {}
for i in reversed( range( len( domains))):
    y = i - 1
    parent = ''
    # Look at the previous terms in the list up to the parent term
    while y >= 0 and domains[ y].count( '--') - domains[ i].count( '--') <= 0:
        # Once the parent term is found, store it in a variable and stop the
        # while loop
        if domains[ y].count( '--') - domains[ i].count( '--') == -1:
            parent = re.sub( "[-]*", '', domains[ y])
            y = 0
        y -= 1
    # If the variable contains a term, we store it in the dict with the term son
    # as the key.
    if parent:
        domain_parents[ re.sub( "[-]*", '', domains[ i])] = parent

###############################################
# Formatter le dataset
###############################################

# Create an occurrence table with proteins in index and ipr_id in column
IPR_table = createTable( IPR_df, col_to_features = "IPR_ID",
                         col_to_indexes = "UniProtAC")

# For all ipr_id annotated with proteins and having a parent term (so that the
# ipr_id is present in the keys of the domain_parents dict)
for ipr_id in IPR_df.IPR_ID.unique():
    if ipr_id in domain_parents.keys():
        # Store in a list the parent term of the current term
        ipr_ancestors = []
        ipr_ancestors.append( domain_parents[ ipr_id])
        y = 0
        # Then search and store the term grandparent, then great-grandparent to
        # the last one that is not present in the domain parents keys
        while ipr_ancestors[ y] in domain_parents.keys():
            ipr_ancestors.append( domain_parents[ ipr_ancestors[ y]])
            y += 1

        # Adds in the table the ancestors columns with the annotations of the
        # heir column
        IPR_table = fillAncestorsTerms( IPR_table, ipr_id, ipr_ancestors)

###############################################
# Export
###############################################

# Remove redundant columns keeping their names concatenated in the remaining
# column name
IPR_table = delRedundantCols( IPR_table)

# Export
IPR_table.to_csv( VIEW_DIR + "Protein_Domains.txt", sep = "\t",
                  encoding = "utf-8", index = True)

# Removes uninformative features and exports this light version of df
IPR_table = CutSmallLeafs( IPR_table)
IPR_table.to_csv( VIEW_DIR + "Protein_Domains_light.txt", sep = "\t",
                  encoding = "utf-8", index = True)
