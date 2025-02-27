#! /bin/bash

cd ../../data/rawData

# Download HPA data
wget https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip
unzip rna_tissue_consensus.tsv.zip
rm rna_tissue_consensus.tsv.zip

# Run the python script to format the dataset
cd ../../script/Tissue_Expression_view/
python Tissue_Expression_settings.py
