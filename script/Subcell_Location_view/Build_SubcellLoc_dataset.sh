#! /bin/bash

cd ../../data/rawData

# Download HPA data
wget https://www.proteinatlas.org/download/tsv/subcellular_location.tsv.zip
unzip subcellular_location.tsv.zip
rm subcellular_location.tsv.zip

# Run the python script to format the dataset
cd ../../script/Subcell_Location_view/
python Subcell_Location_settings.py
