#! /bin/bash

# Create a file where the Interpro data will be put in process
mkdir ../../data/rawData/Interpro_data
cd ../../data/rawData/

# Download interpro data
wget ftp://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt
cd Interpro_data/
wget ftp://ftp.ebi.ac.uk/pub/databases/interpro/protein2ipr.dat.gz
gunzip protein2ipr.dat.gz

# Extract all the columns except the third into a new file
cut -f 1,2 protein2ipr.dat > protein2ipr.txt

# Divide the file in segments of 100M
split -l 100000000 protein2ipr.txt segment_

# In all segments, select the proteins of PPI network
for segment in segment_*
do
	awk -F'\t' 'NR==FNR{c[$1]++;next};c[$1]' ../PPInetwork_proteins.txt $segment > PPIprotein2ipr_$segment.txt
	echo $segment \done
done

# Concatenate extracted rows
cat PPIprotein2ipr_* > ../PPIprotein2ipr.txt

# Remove usless files
cd ..
rm -r Interpro_data

# Run the python script to format the dataset
cd ../../script/Protein_Domains_view/
python Protein_Domains_settings.py
