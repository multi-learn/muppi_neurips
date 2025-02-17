# Define paths
RAWDATA_DIR = "../../data/rawData/"
VIEW_DIR = "../../data/Gene_Ontology_view/"
FUNCTIONS_DIR = '../Functions/'

# Import modules
import pandas as pd
import numpy as np
import os
from goatools.obo_parser import GODag
from goatools.godag.go_tasks import get_go2parents
from goatools.gosubdag.gosubdag import GoSubDag
import sys

sys.path.insert(1, FUNCTIONS_DIR)
from script.Functions.functions import *

#########################################
# Load datas
#########################################

# Define column names
GOA_COLNAMES = ["DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID",
                "DB_Reference", "Evidence_Code", "From", "GO_class",
                "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type",
                "Taxon", "Date", "Assigned_By", "Annotation_Extension",
                "Annotation_Properties"]

# Load the gene Ontology annotations data in a df
GOA_df = pd.read_csv(RAWDATA_DIR + "goa_uniprot_human.gaf", sep='\t',
                     header=None, low_memory=False, names=GOA_COLNAMES)

# Load the Gene Ontology graph
GO_tree = GODag(RAWDATA_DIR + "go-basic.obo", load_obsolete=True, prt=None,
                optional_attrs={'relationship', 'consider', 'replaced_by'})

# Load the PPI proteins list
with open(RAWDATA_DIR + "PPInetwork_proteins.txt", "r") as f:
    PPInetwork_proteins = f.read().split('\n')
    PPInetwork_proteins.remove('')

#########################################
# Format the df and manage obsolete terms
#########################################

# Keep only usefull columns
GOA_df = GOA_df[["DB_Object_ID", "GO_ID", "GO_class"]]
# Keep only annotations about PPI proteins, and drop dupplicated ones
GOA_df = GOA_df[GOA_df.DB_Object_ID.isin(PPInetwork_proteins)].drop_duplicates()

# Rename GO classes
GOA_df.GO_class = GOA_df.GO_class.replace({'P': 'BP', 'C': 'CC', 'F': 'MF'})

# List obsolete go terms
obsolete_terms = [term for term in GO_tree.keys() if GO_tree[term].is_obsolete]

# Create a dict with obsolete terms in keys
obsolete_terms_dict = {}
for term in obsolete_terms:
    # If the obsolete term has considered term, concatenate them and put them
    # in dict value
    # if GO_tree[ term].consider:
    #     obsolete_terms_dict[ term] = '|'.join( GO_tree[ term].consider)
    #
    # If it has replaced_by term, put it in the dict value
    if GO_tree[term].replaced_by[:3] == "GO:":
        obsolete_terms_dict[term] = GO_tree[term].replaced_by

    # If it has none of those, put an empty string in the dict value
    else:
        obsolete_terms_dict[term] = ''

# Replace obsolete terms in the GOA df by their replaced term
GOA_df[GOA_df.GO_ID.isin(obsolete_terms)] = GOA_df[GOA_df.GO_ID.isin(
    obsolete_terms)].replace({
    "GO_ID": obsolete_terms_dict})

# Drop df rows with no GO_ID
GOA_df = GOA_df.drop(GOA_df[GOA_df.GO_ID == ''].index, axis=0).drop_duplicates()

# Create a variable that will indicate which optional relationships to consider
# when searching for the ancestors of a GO_ID.
optional_relationships = {'regulates', 'negatively_regulates',
                          'positively_regulates'}
# Create a dict containing the parent GO_ID of each GO_ID
go2parents = get_go2parents(GO_tree, optional_relationships)

# Define the terms roots cellular component, biological process, and molecular
# function that have no parent terms.
ROOT_GO_ID = ["GO:0005575", "GO:0008150", "GO:0003674"]


###############################################
# Function to build the ontology table
###############################################

### Create a protein / sub-Ontology terms table ###
def get_subGOA_table(GO_subClass):
    # Create a sub df with the annotations of a single GO_class
    subGOA_df = GOA_df[GOA_df.GO_class == GO_subClass].drop("GO_class", axis=1).reset_index()

    # Create an occurrence table with proteins in indexes and GO_ID in columns
    subGOA_table = createTable(subGOA_df, col_to_features="GO_ID",
                               col_to_indexes="DB_Object_ID")

    # For all GO_IDs of the sub graph except the root terms
    for go in subGOA_df.GO_ID.unique():
        if go not in ROOT_GO_ID:
            # Search and store the previous terms of the current GO_ID
            gosubdag = GoSubDag([go], GO_tree, prt=None,
                                relationship=optional_relationships)
            GO_ancestors = gosubdag.rcntobj.go2parents[go]

            # Adds in the table the ancestors columns with the annotations of
            # the heir column
            fillAncestorsTerms(subGOA_table, go, GO_ancestors)

    # Replace values other than 0 with 1's to remove additions.
    subGOA_table[subGOA_table != 0] = 1

    # Removes redundant columns by keeping their names concatenated in the
    # name of the remaining column
    subGOA_table = delRedundantCols(subGOA_table)
    # Convert the df in integer type
    subGOA_table = subGOA_table.astype(int)

    # Return the table
    return (subGOA_table)


###############################################
# Creating and exporting subGOA tables
###############################################

# Store the sub classes of Gene_Ontology
subClasses = GOA_df.GO_class.unique()

# For each sub class
for subClass in subClasses:
    # Generate the sub GOA table and export it
    GOA_table = get_subGOA_table(subClass)
    GOA_table.to_csv(VIEW_DIR + "Gene_Ontology_{}.txt".format(subClass),
                     sep="\t", encoding="utf-8", index=True)
    # Removes uninformative features and exports this light version of df
    GOA_table = CutSmallLeafs(GOA_table)
    GOA_table.to_csv(VIEW_DIR + "Gene_Ontology_{}_light.txt".format(subClass),
                     sep="\t", encoding="utf-8", index=True)
    del (GOA_table)
