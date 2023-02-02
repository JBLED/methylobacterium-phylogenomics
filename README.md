https://github.com/JBLED/methylobacterium-phylogenomics.git

##########################################
IF USING ANY INFORMATION RELATED TO THIS PROJECT, PLEASE READ AND CITE:

Leducq et al. (2022) Comprehensive phylogenomics of Methylobacterium reveals four evolutionary distinct groups and underappreciated phyllosphere diversity.
Genome Biology and Evolution, Volume 14, Issue 8, August 2022

https://doi.org/10.1093/gbe/evac123

##########################################
DESCRIPTION

###0-Original-code-and-ressources
This section contains the original R code that was used in Leducq et al. 2022 (doi.org/10.1093/gbe/evac123) and major outputs.
As this code is a bit messy and I keep adding new Methylobacteriaceae genomes to the workflow, I decided to provide it as raw. 
Updated and cleaner versions of the code and outputs will be progressively added to the following sections.

###1-Filtering-24-Methylobacterium-genome-assemblies
Filtering and cleaning of 24 Methylobacterium genomes assemblies

- 1-assemblies.R (R script)

- summary of assembly statistics before assembly filtering (ummary-megahit.txt ; summary-megahit.pdf)

- For each assembly (24), this folder also contains: 
1) Fasta file with megahit assembly (xxx_final.contigs.fa)
2) Filtering summary file (xxx_summary.txt)
3) Fasta file without filtered scaffolds (xxx_filtered.fa)

- Fasta and summary files for 3 contaminated assemblies (E-016, E-025, J-059), before and after removing contaminating scaffolds (clean assemblies in xxx_filtered-clean.fa; contaminant genomes in xxx_filtered-contaminant.fa)

###2-Determination-Methylobacterium-core-genome

- 2-Determination-Core-Genome-Original-Dataset.R

###3-Genome-Annotation-Phylogeny-Synteny

- 3-ADD-Alessa-genomes.R

###4-Synteny-visualisation

- 4-Formating-Input-Cytoscape.R
