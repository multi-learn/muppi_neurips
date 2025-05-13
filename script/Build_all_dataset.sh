#! /bin/bash
set -e

cd PPInetwork_topology_view
python Build_PPInetwork_topology_dataset.py
echo "===== 1. Finished PPInetwork_topology_view ======"

cd ../PPInetwork_Embedding_view
python Build_PPInetEmbedding_dataset.py
echo "===== 2. Finished PPInetwork_Embedding_view ====="

cd ../3UTR_Complexes_view
python Build_3UTR_Complexes_dataset.py
echo "===== 3. Finished 3UTR_Complexes_view ====="

cd ../Gene_Ontology_view
sh Build_GeneOntology_datasets.sh
echo "===== 4. Finished Gene_Ontology_view ====="

cd ../GO_PPInetwork_view
python GO_PPInetwork_settings.py
python Build_Weighted_Embeddings.py
echo "===== 5. Finished GO_PPInetwork_view ====="

cd ../Linear_Motifs_view
sh Build_LinearMotifs_dataset.sh
echo "===== 6. Finished Linear_Motifs_view ====="

cd ../Phenotype_Ontology_view
sh Build_PhenotypeOntology_dataset.sh
echo "===== 7. Finished Phenotype_Ontology_view ====="

cd ../Post_Traductionnal_Modifications_view
sh Build_PTM_dataset.sh
echo "===== 8. Finished Post_Traductionnal_Modifications_view ====="

cd ../Protein_Domains_view
sh Build_ProteinDomains_dataset.sh
echo "===== 9. Finished Protein_Domains_view ====="

cd ../Subcell_Location_view
sh Build_SubcellLoc_dataset.sh
echo "===== 10. Finished Subcell_Location_view ====="

cd ../Tissue_Expression_view
sh Build_TissueExpression_dataset.sh
echo "===== 11. Finished Tissue_Expression_view ====="
