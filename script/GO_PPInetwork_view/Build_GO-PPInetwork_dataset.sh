#! /bin/bash

# Créer des edgelists avec des poids liés aux GOA communes
python GO_PPInetwork_settings.py

# Génère les embeddings de graphs pondérés
cd ../../data/GO_PPInetwork_view/
python -m openne --method sdne --input GO-CC_PPInetwork.txt --graph-format edgelist --output GO-CC_PPInetwork_embed.txt --weighted
python -m openne --method sdne --input GO-BP_PPInetwork.txt --graph-format edgelist --output GO-BP_PPInetwork_embed.txt --weighted

# Supprime la première ligne
sed -i '1d' GO-CC_PPInetwork_embed.txt
sed -i '1d' GO-BP_PPInetwork_embed.txt
