# Assessing mutation rates in putatively constrained sites in the human genome
This repository contains pre-processed data (mutation rates and constraint scores) and scripts required to generate the results presented in Barroso & Ragsdale 2026b.<br>
The 'example' directory holds intermediate representations between raw data files (chromosome 22, mutation maps and GPN-Star scores for the Mammalian phylogeny) and the compressed tables used to generate plots. Due to their large size, the original files cannot be stored in this repository.

# Hierarchy of scripts to handle input data
Here is an overview of the R pipeline used to generate our results. We use the example of GPN-Star scores because they require one extra (upstream) step (level 0):

## Level 0
process_gpn-star.R<br>
This script is not needed for the other annotation sets. It reads the original files (parquet format) for a given phylogenetic depth and creates simplified TSV files.
This setup was chosen for compatibility reasons.

## Level 1
parse_muts.R, parse_gpn-x.R<br>
These scripts read the Roulette files (VCF format) and TSV score files for a given chromosome and further compress them, calculating the average rate per site.

## Level 2
quant_bins_gpn-x.R<br>
This script takes single-nucleotide constraint scores and assign the top 30% to 12 equal-size constraint classes.

## Level 3
summarize_mut_gpn-x.R<br>
This script merges genomic maps of constraint, mutation and B-values and generates compressed summaries used for plotting.

## Level 4
plot_gpn-x.R<br>
This script produces intermediate plots.

## Level 5
comb_plots.R<br>
This script combines data across annotation sets and produces the plots shown in the manuscript.
