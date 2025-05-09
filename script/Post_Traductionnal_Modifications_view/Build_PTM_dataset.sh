#! /bin/bash

# Create a folder where the dbPTM data will be put during processing
mkdir ../../data/rawData/PTM_data
cd ../../data/rawData/PTM_data/

# CHANGED: Download links of all dbPTM files
# Download all dbPTM files
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Phosphorylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Acetylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Ubiquitination.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Succinylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Methylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Malonylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/N-linked\ Glycosylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/O-linked\ Glycosylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Sumoylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/S-nitrosylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Glutathionylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Amidation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Hydroxylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Pyrrolidone\ carboxylic\ acid.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Glutarylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/S-palmitoylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/O-palmitoylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/N-palmitoylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Gamma-carboxyglutamic\ acid.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Crotonylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Oxidation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Myristoylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/C-linked\ Glycosylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Sulfation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Formylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Citrullination.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/GPI-anchor.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Nitration.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/S-diacylglycerol.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Carboxylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Lipoylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Carbamidation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Neddylation.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/Pyruvate.gz
wget https://biomics.lab.nycu.edu.tw/dbPTM/download/experiment/S-linked\ Glycosylation.gz

# Unzip all files
gunzip *

# CHANGED: How we extract the content of the files
for file in ./*; do
  tar -xOf "$file" > "${file}.txt"
done

# Concatenate them
cat *.txt > dbPTM_allSpecies.txt

# Selects the two columns of interest and rows relating to human proteins
grep HUMAN dbPTM_allSpecies.txt | cut -f 2,4 > dbPTM_human.txt

# Keep only PPI network proteins
awk -F'\t' 'NR==FNR{c[$1]++;next};c[$1]' ../PPInetwork_proteins.txt dbPTM_human.txt > ../dbPTM_PPIproteins.txt

# Remove useless files
cd ..
rm -r PTM_data

# Run the python script to format the dataset
cd ../../script/Post_Traductionnal_Modifications_view/
python PTM_settings.py
