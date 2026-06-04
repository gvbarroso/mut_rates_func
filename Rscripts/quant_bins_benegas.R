#################
#
# This script reads processed benegas scores 
# and defines the quantile bins per chr mimicking the phastCons classification
#
#################

library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
chr <- args[1]

cat(paste0("Gathering data from chr ", chr, "..."))

score_files <- list.files("split_scores/", pattern=paste0("^scores_chr", chr, "_"), full.names=T)
scores_chr <- data.table::rbindlist(lapply(score_files, fread))

cat("done.\nComputing quantiles...")

# splitting top 30% of sites into bins, treat the other 70% as neutral
threshold <- quantile(scores_chr$benegas_score, 0.30)
scores_chr[, is_benegas := benegas_score <= threshold]

cat("done.\nAssigning bins...")

nbins <- 12 # discretizing distribution of constraint scores
N_ben <- scores_chr[is_benegas == TRUE, .N]
scores_chr[is_benegas==T, rank_ben := frank(benegas_score, ties.method = "average")]
scores_chr[is_benegas==T, benegas_bin := pmin(nbins, ceiling(rank_ben / (N_ben / nbins)))]
scores_chr[is_benegas==F, benegas_bin := nbins + 3] # neutral bin
scores_chr[, c("is_benegas", "rank_ben", "benegas_score") := NULL] 

cat("done.\nWriting to file...")

fwrite(scores_chr, paste0("score_bins/benegas_bins_chr", chr, ".csv.gz"))

cat("Finished!\n")
