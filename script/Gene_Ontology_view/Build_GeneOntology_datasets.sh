#! /bin/bash

# Go to the rawData directory
cd ../../data/rawData

# Download and unzip GO annotations file
wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz

# Download the hierarchical tree of GO terms
wget https://current.geneontology.org/ontology/go-basic.obo

# Selects from the annotation file the lines concerning human proteins.
# HACK: Instead of unzipping the entire goa_uniprot_all.gaf.gz file, we read it line by line
gzcat goa_uniprot_all.gaf.gz | egrep -w 'protein\staxon:9606\|?' | tee goa_uniprot_human.gaf | pv > /dev/null

# Delete useless files
rm goa_uniprot_all*

# Run the python script to format the dataset
cd ../../script/Gene_Ontology_view/
python Gene_Ontology_settings.py
