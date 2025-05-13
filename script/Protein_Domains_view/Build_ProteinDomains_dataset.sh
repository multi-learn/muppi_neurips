#!/bin/bash

set -e  # Stop on error

# Créer le dossier de destination
mkdir -p ../../data/rawData/Interpro_data
cd ../../data/rawData/

# Télécharger les fichiers si non présents
[ ! -f ParentChildTreeFile.txt ] && wget https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/ParentChildTreeFile.txt
cd Interpro_data/
[ ! -f protein2ipr.dat.gz ] && wget https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz

# Vérifier que le fichier de protéines PPI existe
if [ ! -f ../PPInetwork_proteins.txt ]; then
    echo "Fichier ../PPInetwork_proteins.txt introuvable."
    exit 1
fi

# Extraire les protéines du réseau PPI directement depuis l’archive compressée
gzcat protein2ipr.dat.gz \
    | cut -f 1,2 \
    | awk -F'\t' 'NR==FNR{c[$1]; next} $1 in c' ../PPInetwork_proteins.txt - \
    > ../PPIprotein2ipr.txt

# Nettoyage
cd ..
rm -r Interpro_data

# Exécuter le script Python de formatage
cd ../../script/Protein_Domains_view/
python Protein_Domains_settings.py