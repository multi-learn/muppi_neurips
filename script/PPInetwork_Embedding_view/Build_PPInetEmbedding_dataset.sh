#! /bin/bash

# Generate the graph embedding
cd ../../data/
python -m openne --method sdne --input rawData/PPI_network.txt --graph-format edgelist --output PPInetwork_Embedding_view/SDNE_PPInetwork.txt

# Remove the first row of embedding file
sed -i '1d' PPInetwork_Embedding_view/SDNE_PPInetwork.txt
