########################################################################################################################
###########Raimondi et al., Insights into cancer severity from biomolecular interaction mechanisms######################
##############Analysis of COSMIC GENOMES (v74) somatic missense mutations with Mechismo - Raw data############################
########################################################################################################################
###
###Each cancer type-specific directory contains a "sample_id.txt" file with Mechismo predictions for each cancer sample
###
###Fields description:

#1 Cosmic Sample ID
#2 gene name of query
#3 primary identifier of query (usually UniProt accession)
#4 primary id/variant (wild-type residue,position in the query sequence, mutated residue / PTM)
#5 name of the interactor:
   -protein-protein interactions: gene name of interactor
   -protein-chemical interactions: '[CHEM:type:id]'
   -DNA/RNA interactions: '[DNA/RNA]'
   -single-structure matches (i.e. no interaction): '[PROT]'
#6 primary identifier of interactor (usually UniProt accession for protein-protein interactions)
#7 mechismo score (if the mechismo predicted effect (col.#8) is predicted to be disabling(Weak), the score is multiplied by -1)
#8 mechismo predicted effect
#9 confidence in the structural template
#10 intEvLTP - interaction known from low-throughput evidence: true (1) or false (0)
#11 ntEvHTP - interaction known from high-throughput evidence: true (1) or false (0)
#12 intEvStructure - interaction known from structure: true (1) or false (0)
#13 position in the query sequence
#14 Query Gene name(Pfam domain at the corresponding query position)
#15 Pfam domain position:representative amino acid (upper case is for highly conserved)
#16 position in interactor sequence
#17 Interactor Gene name(Pfam domain at the corresponding query position) in case of protein-protein interactions
#18 Pfam domain position:representative amino acid (upper case is for highly conserved)
#19 pdb identifier of template structure
#20 position in the template structure sequence
#21 PDB residue sequence identifier
#22 PDB chain identifier
#23 Pfam domain at the corresponding PDB  query position
#24 interacting residue position in the template structure sequence
#25 interacting PDB residue sequence identifier
