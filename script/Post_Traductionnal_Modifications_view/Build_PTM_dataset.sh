#! /bin/bash

# Create a folder where the dbPTM data will be put during processing
mkdir ../../data/rawData/PTM_data
cd ../../data/rawData/PTM_data/

# Download all dbPTM files
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Phosphorylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Acetylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Ubiquitination.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Succinylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Methylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Malonylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/N-linkedGlycosylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/O-linkedGlycosylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Sumoylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/S-nitrosylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Glutathionylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Amidation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Hydroxylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Pyrrolidonecarboxylicacid.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Glutarylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Palmitoylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Gamma-carboxyglutamicacid.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Crotonylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Oxidation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Myristoylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/C-linkedGlycosylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Sulfation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Formylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Citrullination.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/GPI-anchor.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Nitration.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/S-diacylglycerol.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Carboxylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Lipoylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Carbamidation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Neddylation.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/Pyruvate.txt.gz
wget http://dbptm.mbc.nctu.edu.tw/download/experiment/S-linkedGlycosylation.txt.gz

# Unzip all files
gunzip *

# Concatenate them
cat *.txt > dbPTM_allSpecies.txt

# Selects the two columns of interest and rows relating to human proteins
grep HUMAN dbPTM_allSpecies.txt | cut -f 2,4 > dbPTM_human.txt

# Keep only PPI network proteins
awk -F'\t' 'NR==FNR{c[$1]++;next};c[$2]' ../PPInetwork_proteins.txt dbPTM_human.txt > ../dbPTM_PPIproteins.txt

# Remove useless files
cd ..
rm -r PTM_data

# Run the python script to format the dataset
cd ../../script/Post_Traductionnal_Modifications_view/
python PTM_settings.py
