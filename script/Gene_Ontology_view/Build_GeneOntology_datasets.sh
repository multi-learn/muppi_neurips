#! /bin/bash

# Go to the rawData directory
cd ../../data/rawData

# Download and unzip GO annotations file
wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz
gunzip goa_uniprot_all.gaf.gz

# Download the hierarchical tree of GO terms
wget http://purl.obolibrary.org/obo/go-basic.obo

# Selects from the annotation file the lines concerning human proteins.
egrep -w 'protein\staxon:9606\|?' goa_uniprot_all.gaf > goa_uniprot_human.gaf

# Delete useless files
rm goa_uniprot_all*

# Run the python script to format the dataset
cd ../../script/Gene_Ontology_view/
python Gene_Ontology_settings.py
