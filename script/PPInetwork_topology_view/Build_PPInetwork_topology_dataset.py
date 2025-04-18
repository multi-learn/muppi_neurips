# Define paths
RAWDATA_DIR = "../../data/rawData/"
VIEW_DIR = "../../data/PPInetwork_topology_view/"

# Import module
import pandas as pd

# Load node table in a df
PPInetwork_df = pd.read_csv(RAWDATA_DIR + 'rawPPI_network_nodeTable.csv',
                            sep=',')

# Define wanted columns
SELCTED_COLUMNS = ['name', 'AverageShortestPathLength', 'BetweennessCentrality',
                   'ClusteringCoefficient', 'ClosenessCentrality', 'Eccentricity',
                   'Stress', 'Degree', 'NeighborhoodConnectivity',
                   'Radiality', 'TopologicalCoefficient']

# Select them and set the index on proteins ids
PPInetwork_df = PPInetwork_df[PPInetwork_df.columns.intersection(SELCTED_COLUMNS)].set_index('name')

# CHANGED: Adding indexes
# Export
PPInetwork_df.to_csv(VIEW_DIR + "PPInetwork_topology.txt", sep="\t",
                     encoding="utf-8", index=True)
