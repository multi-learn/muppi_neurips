#! /bin/bash

# Connect to https://cancer.sanger.ac.uk/cosmic then
# https://cancer.sanger.ac.uk/cosmic/download and download the file
# CosmicMutantExport.tsv.gz

# Create a working directory for cosmic datas
mkdir ../../data/rawData/Cosmic_data
cd ../../data/rawData/Cosmic_data

# Transfer and unzip the downloaded file in this directory
cp ~/Téléchargements/CosmicMutantExport.tsv.gz .
gunzip CosmicMutantExport.tsv.gz

# Extract the ENST column
cut -f 2 CosmicMutantExport.tsv > cosmicENST.txt

# Remove version number of ENST ids (numbers following the point)
sed -i 's/\.[0-9]*$//g' cosmicENST.txt

# Concatenate modified ENST column with the initial file
paste -d"\t" cosmicENST.txt CosmicMutantExport.tsv > 2CosmicMutantExport.tsv

# Extract columns of interest in a new file
cut -f 1,9-16,29,32-33 2CosmicMutantExport.tsv > 3CosmicMutant.tsv

# echo Créer un nouveau tableau ou seules les lignes contenant les ENST presentes dans un autre fichier sont sélectionnées
awk -F'\t' 'NR==FNR{c[$2]++;next};c[$1]' ../PPIproteins_ENST.txt 3CosmicMutant.tsv > ../CosmicMutant.tsv

# Remove useless files
cd ..
rm -r Cosmic_data

# Run the python script to format the dataset
cd ../../script/Cancer_Mutations_view/
python Cancer_Mutations_settings.py
