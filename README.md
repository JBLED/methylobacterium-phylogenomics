https://github.com/JBLED/methylobacterium-phylogenomics.git

##########################################
## IF USING ANY INFORMATION RELATED TO THIS PROJECT, PLEASE READ AND CITE:

## Leducq et al. (2022b) Comprehensive phylogenomics of Methylobacterium reveals four evolutionary distinct groups and underappreciated phyllosphere diversity.
Genome Biology and Evolution, Volume 14, Issue 8, August 2022

https://doi.org/10.1093/gbe/evac123

##########################################
# DESCRIPTION

## 0-Original-codes-and-ressources (AVAILABLE)

This section contains the original R code that was used in Leducq et al. 2022b (doi.org/10.1093/gbe/evac123), annotation files from the original 213 Methylobacteriaceae genomes (myRAST outputs), major outputs from the R codes (alignments and statistics) and tree files (outputs from RAxML, ASTRAL and SVDquartet).
As the R code is a bit messy, and as I keep adding new Methylobacteriaceae genomes to the workflow, I decided to provide it as raw. 
Updated and cleaner versions of the code and outputs will be progressively added to the following sections.

## 1-Filtering-24-Methylobacterium-genome-assemblies (AVAILABLE)

Filtering and cleaning of 24 Methylobacterium genomes assemblies (sequenced in Leducq et al. 2022b)

- 1-assemblies.R (R script)

- summary of assembly statistics before assembly filtering (summary-megahit.txt ; summary-megahit.pdf)

- For each assembly (24), this folder also contains: 
1) Fasta file with megahit assembly (xxx_final.contigs.fa)
2) Filtering summary file (xxx_summary.txt)
3) Fasta file without filtered scaffolds (xxx_filtered.fa)

- Fasta and summary files for 3 contaminated assemblies (E-016, E-025, J-059), before and after removing contaminating scaffolds (clean assemblies in xxx_filtered-clean.fa; contaminant genomes in xxx_filtered-contaminant.fa)

## 2-Determination-Methylobacterium-core-genome (AVAILABLE)

Methylobacterium core genome determination and production of alignments - Starting for 213 Methylobacteriaceae genomes used in Leducq et al. 2022b (doi.org/10.1093/gbe/evac123) and enriched with 21 newly released Methylobacterium genomes. Unlike in Leducq et al. 2022b, the core genome is here defined only for Methylobacterium (Microvirga and Enterovirga excluded, but still nucleotide sequences were retrieved when available), resulting in 616 core genes (384 in the study).

- 2-Determination-Core-Genome-Feb-2023.R (R-script)

1) Annotation statistics for each genome (Annotation-Summary.txt)
2) Raw gene abundance per genome and per group (Genome-summary.pdf; Gene-copy-number-raw.txt)
3) Proportion of genes with 1,2,3,4 or more copies in function of the number of scaffolds (Gene-copy.txt; GeneCopy-function-scaffold.pdf)
4) Estimate of Methylobacterium core genome size by K random ressampling of genomes within groups (Annot-Raw-Rarefaction_K=K.pdf; Annot-Raw-Rarefaction_K=K.txt)
5) Look at the number of core genes per group in function of different threshold (Annot-Numb-Core-Gene-per-group.txt)
6) Refine the copy number estimation for these genes in all genomes using the sequence size normalized by size in complete genomes - look at distributions and define threshold (NCx) to keep only 1:1 genes within groups AND within Methylobacterium (Gene-modified-names.txt; RawFcoreCount_TH>=NCx.txt; RawFcoreSize_TH>=NCx.txt; NormalizedSize-vs-CopyNumber_TH>=NCx.pdf; Filtered-NormalizedSize-vs-CopyNumber_TH>=NCx.pdf; Gene-occurence-after-filtering.pdf; CoreGene-per-genome-N50.pdf). Also produce a fasta file for each genome, containing nucleotide sequences of core genes (231 files: Core-Genes-stXXX.fas zipped in two files in 2-Determination-Methylobacterium-core-genome/OUT-FASTA)
7) For each core gene, produce a fasta file (geneXXXX.fas; not provided here) from genomes for which sequences are available (Summary-per-core-gene.txt; summary-core-genes-final.pdf; FinalFcoreCount_TH>=NCx.txt)
8) Align sequences for each core gene (Summary-per-core-gene-after-alignment.txt; 616 fasta files: GeneXXXX_alignedNF.fas zipped in two files in 2-Determination-Methylobacterium-core-genome/OUT-FASTA; intermediate alignement files are not provided)

## 3-Genome-Annotation-Phylogeny-Synteny (TO DO)

## 4-Synteny-visualisation (TO DO)
