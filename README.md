https://github.com/JBLED/methylobacterium-phylogenomics.git

##########################################
IF USING ANY INFORMATION RELATED TO THIS PROJECT, PLEASE READ AND CITE:

Leducq et al. (2022) Comprehensive phylogenomics of Methylobacterium reveals four evolutionary distinct groups and underappreciated phyllosphere diversity.

PRE-PRINT VERSION: xxx and updated versions

##########################################
DESCRIPTION

###1-Filtering-24-Methylobacterium-genome-assemblies

- 1-assemblies.R (R script)

For each genome (24), this folder contains: 
1) Fasta files with megahit assemblies of 24 Methylobacterium genomes (xxx_final.contigs.fa)
2) Filtering summary file (xxx_summary.txt)
3) Fasta files without filtered scaffolds (xxx_filtered.fa)
4) A summary of assembly statistics before filtering (ummary-megahit.txt ; summary-megahit.pdf)
5) Fasta and summary files for 3 contaminated assemblies, before and after removing contaminating scaffolds (E-016, E-025, J-059; clean assemblies in xxx_filtered-clean.fa; contaminant genomes in xxx_filtered-contaminant.fa)

###2-Determination-Methylobacterium-core-genome

- 2-Determination-Core-Genome-Original-Dataset.R

###3-Genome-Annotation-Phylogeny-Synteny

- 3-ADD-Alessa-genomes.R

###4-Synteny-visualisation

- 4-Formating-Input-Cytoscape.R
