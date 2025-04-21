#! /bin/bash

# Créer des edgelists avec des poids liés aux GOA communes
python GO_PPInetwork_settings.py

# Génère les embeddings de graphs pondérés
python Build_Weighted_Embeddings.py
