#! /bin/bash

cd ../../data/rawData

# Get elm classes data
wget http://elm.eu.org/elms/elms_index.tsv

# Delete comments in the firsts rows
sed '1,5d' elms_index.tsv | cut -f 2,6 > elm_classes.txt
rm elms_index.tsv

# Get elm instances data
wget http://elm.eu.org/instances.tsv?q=*

# Select elm instances that have been experimentaly observed in human
# And extract only the list of elm
sed '1,5d' instances.tsv?q=* | grep sapiens | grep 'true positive' | cut -f 3 | sort | uniq > elm_instances_TP_HomoSapiens.txt
rm instances.tsv?q=*

# Execute the python script
cd ../../script/Linear_Motifs_view/
python Linear_Motifs_SLIMProb.py
