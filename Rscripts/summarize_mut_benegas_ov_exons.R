#########################
#
# Setting up
#
########################

library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly=T)
chr <- args[1]

#########################
#
# Summarizing tables
#
########################

cat(paste0("Reading mutation maps for chr ", chr, "..."))
mut_files <- list.files("/../../media/gvbarroso/extradrive1/mut_rates_func/split_muts/", pattern=paste0("^mut_map_chr", chr, "_"), full.names=T)
#mut_files <- list.files("../../../Data/transfer/", pattern=paste0("^mut_map_chr", chr, "_"), full.names=T)
mut_map <- data.table::rbindlist(lapply(mut_files, fread))

cat("done.\nReading score map...") # reads all scores (not just those overlapping exons)
scores_chr <- fread(paste0("/../../media/gvbarroso/extradrive1/mut_rates_func/benegas/score_bins/benegas_bins_chr", chr, ".csv.gz"))
#scores_chr <- fread(paste0("../../../Data/transfer/benegas_bins_chr", chr, ".csv.gz"))
setnames(scores_chr, old="benegas_bin", new="benegas_class")

cat("done.\nReading overlap benegas-exons...") # Benegas && Exons from benegas_exons_overlaps.R
ov_chr <- fread(paste0("summary_tbls/overlaps_benegas_exons_chr", chr, ".csv.gz"))

cat("done.\nJoining maps...")
scores_chr <- mut_map[scores_chr, on=.(chrom, pos), nomatch=0] 
scores_chr <- merge(scores_chr, ov_chr, by=c("chrom", "pos"), all=T) 
# sanitizing
scores_chr[, benegas_class := fcoalesce(benegas_class.x, benegas_class.y)]
scores_chr[, benegas_class.x := NULL]
scores_chr[, benegas_class.y := NULL]
# setting status of sites the do not overlap exons to FALSE
set(scores_chr, which(is.na(scores_chr[["overlap_benegas_exon"]])), "overlap_benegas_exon", F)

cat("done.\nSummarizing tables...")

# average mutation rates per class of constraint
scores_chr[, mean_class_roulette := mean(roulette, na.rm=T), by=.(benegas_class, overlap_benegas_exon)]
scores_chr[, mean_class_gnomad := mean(gnomad, na.rm=T), by=.(benegas_class, overlap_benegas_exon)]
scores_chr[, mean_class_carlson := mean(carlson, na.rm=T), by=.(benegas_class, overlap_benegas_exon)]
# NOTE skipping these summaries ti keep table clean
#scores_chr[, se_class_roulette := sd(roulette, na.rm=T) / sqrt(.N), by=.(benegas_class, overlap_benegas_exon)]
#scores_chr[, se_class_gnomad := sd(gnomad, na.rm=T) / sqrt(.N), by=.(benegas_class, overlap_benegas_exon)]
#scores_chr[, se_class_carlson := sd(carlson, na.rm=T) / sqrt(.N), by=.(benegas_class, overlap_benegas_exon)]

# average mutation rates per class of constraint per triplet context
#scores_chr[, mean_class_triplet_roulette := mean(roulette, na.rm=T), by=.(benegas_class, triplet, overlap_benegas_exon)]
#scores_chr[, mean_class_triplet_carlson := mean(carlson, na.rm=T), by=.(benegas_class, triplet, overlap_benegas_exon)]
#scores_chr[, mean_class_triplet_gnomad := mean(gnomad, na.rm=T), by=.(benegas_class, triplet, overlap_benegas_exon)]
#scores_chr[, se_class_triplet_roulette := sd(roulette, na.rm=T) / sqrt(.N), by=.(benegas_class, triplet, overlap_benegas_exon)]
#scores_chr[, se_class_triplet_carlson := sd(carlson, na.rm=T) / sqrt(.N), by=.(benegas_class, triplet, overlap_benegas_exon)]
#scores_chr[, se_class_triplet_gnomad := sd(gnomad, na.rm=T) / sqrt(.N), by=.(benegas_class, triplet, overlap_benegas_exon)]

# storing number of sites per class per triplet per overlap status (used to compute the means for focal chr)
scores_chr[, `:=`(num_sites_roulette=sum(!is.na(roulette)),
                  num_sites_carlson=sum(!is.na(carlson)),
                  num_sites_gnomad=sum(!is.na(gnomad))),
           by=.(benegas_class, overlap_benegas_exon)]

unique_vals <- unique(scores_chr, by=c("benegas_class", "decile", "overlap_benegas_exon"))
setorder(unique_vals, overlap_benegas_exon, benegas_class, triplet)
cat("done.\nWriting to file...")

fwrite(unique_vals[,-c("pos", "triplet")], paste0("summary_tbls/overlap_summaries_benegas_exons_chr", chr, ".csv.gz"))
cat("Finished.\n")
