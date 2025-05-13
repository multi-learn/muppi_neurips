#! /bin/bash

cd ../../data/rawData

# Download the HPO ontology and HPO annotation file
wget - 0 HPO_genes_to_phenotype.txt "https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2025-05-06/genes_to_phenotype.txt" HPO_genes_to_phenotype.txt
wget http://purl.obolibrary.org/obo/hp.obo

# Run the python script to format the dataset
cd ../../script/Phenotype_Ontology_view/
python Phenotype_Ontology_settings.py
