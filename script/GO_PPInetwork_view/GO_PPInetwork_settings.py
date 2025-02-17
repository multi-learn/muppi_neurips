# Define paths
RAWDATA_DIR = "../../data/rawData/"
GO_VIEW_DIR = "../../data/Gene_Ontology_view/"
VIEW_DIR = "../../data/GO_PPInetwork_view/"

# Import modules
import pandas as pd
import numpy as np

#########################################
# Load datas
#########################################

# Load Gene_Ontology datasets
GOA_BP_df = pd.read_csv(GO_VIEW_DIR + 'Gene_Ontology_BP.txt',
                        sep='\t', header=0, index_col=0)
GOA_CC_df = pd.read_csv(GO_VIEW_DIR + 'Gene_Ontology_CC.txt',
                        sep='\t', header=0, index_col=0)

# Load the edgelist of PPI network
PPInet_edges_df = pd.read_csv(RAWDATA_DIR + "human_binary_network.txt",
                              sep='\t',
                              header=None,
                              names=["Interactor_A", "Interactor_B"])

# Load the PPI proteins list
with open(RAWDATA_DIR + "PPInetwork_proteins.txt", "r") as f:
    PPInetwork_proteins = f.read().split('\n')
    PPInetwork_proteins.remove('')


#########################################
# Merge and format the df
#########################################

def jaccard_index(edge, GOA_df):
    """Calculate the similarity of the GO annotations of the two proteins on an
    edge by taking into account the number of total annotations of the two
    proteins."""

    GOA_prot1 = GOA_df.loc[edge[0]]
    GOA_prot2 = GOA_df.loc[edge[1]]
    common_GOterms = sum(GOA_prot1 & GOA_prot2)
    jaccard_simmilarity = common_GOterms / (sum(GOA_prot1) + sum(GOA_prot2)
                                            - common_GOterms)
    return (jaccard_simmilarity)


def get_edge_weights(GOA_df):
    """This function takes the list of edges of the PPI network and adds to each
    one a weight according to the common GO annotations. Takes as input the
    table of GO annotations by proteins"""

    # Create a sub df from the list of edges of the PPI, by selecting the edges
    # between the proteins that are present in the df GOA used to calculate the
    # weights.
    mask = PPInet_edges_df.Interactor_A.isin(GOA_df.index) & \
           PPInet_edges_df.Interactor_B.isin(GOA_df.index)
    weighted_edges_df = PPInet_edges_df[mask]

    # Compute the number of common GOterms between the proteins of each edge and
    # the stocks in a new column of the edge df.
    weighted_edges_df["weight"] = weighted_edges_df.apply(lambda edge:
                                                          jaccard_index(edge, GOA_df), axis=1)

    # Returns the list of edges with their weights
    return (weighted_edges_df)


# Create 2 df with the list of edges of the PPI network with a column of weights
# related to the similarity of GO annotations of proteins
PPInet_CCweighted_df = get_edge_weights(GOA_CC_df)
PPInet_BPweighted_df = get_edge_weights(GOA_BP_df)

# Export
PPInet_CCweighted_df.to_csv(VIEW_DIR + "GO-CC_PPInetwork.txt", index=False,
                            header=False, sep="\t", encoding="utf-8")
PPInet_BPweighted_df.to_csv(VIEW_DIR + "GO-BP_PPInetwork.txt", index=False,
                            header=None, sep="\t", encoding="utf-8")
