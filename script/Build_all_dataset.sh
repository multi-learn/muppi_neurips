#! /bin/bash

# CHANGED: Replacing sh by python for python scripts

cd 3UTR_Complexes_view
python Build_3UTR_Complexes_dataset.py

cd ../Cancer_Mutations_view
sh Build_CancerMutations_dataset.sh

cd ../Gene_Ontology_view
sh Build_GeneOntology_datasets.sh

cd ../GO_PPInetwork_view
sh Build_GO-PPInetwork_dataset.sh

cd ../Linear_Motifs_view
sh Build_LinearMotifs_dataset.sh

cd ../Phenotype_Ontology_view
sh Build_PhenotypeOntology_dataset.sh

cd ../Post_Traductionnal_Modifications_view
sh Build_PTM_dataset.sh

cd ../PPInetwork_Embedding_view
sh Build_PPInetEmbedding_dataset.sh

cd ../PPInetwork_topology_view
python Build_PPInetwork_topology_dataset.py

cd ../Protein_Domains_view
sh Build_ProteinDomains_dataset.sh

cd ../Subcell_Location_view
sh Build_SubcellLoc_dataset.sh

cd ../Tissue_Expression_view
sh Build_TissueExpression_dataset.sh
