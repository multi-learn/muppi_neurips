#! /bin/bash

cd ../../data/rawData

# Download the HPO ontology and HPO annotation file
wget http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt
wget http://purl.obolibrary.org/obo/hp.obo

# Change the name of the annotation file
mv ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt HPO_genes_to_phenotype.txt

# Run the python script to format the dataset
cd ../../script/Phenotype_Ontology_view/
python Phenotype_Ontology_settings.py
