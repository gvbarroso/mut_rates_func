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

score_files <- list.files("transfer/", pattern=paste0("scores_chr", chr, "*"), full.names=T)
scores_chr <- data.table::rbindlist(lapply(score_files, fread))

# splitting top 30% of sites into bins, treat the other 70% as neutral
threshold <- quantile(scores_chr$benegas_score, 0.30)
scores_chr[, is_benegas := ifelse(benegas_score <= threshold, T, F)] # more negative -> more constrained

nbins <- 12 # discretizing distribution of constraint scores
scores_chr[is_benegas==T, benegas_bin := as.numeric(cut(frank(benegas_score, ties.method="average") / .N,
                                                 breaks=seq(0, 1, by = 1 / nbins), labels = 1:nbins, include.lowest=T))]
scores_chr[is_benegas==F, benegas_bin := nbins + 3] # neutral bin
scores_chr[, c("is_benegas", "benegas_score") := NULL] 

fwrite(scores_chr, paste0("score_bins/benegas_bins_chr", chr, ".csv.gz"))
