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
mut_files <- list.files("../split_muts/", pattern=paste0("^mut_map_chr", chr, "_"), full.names=T)
mut_map <- data.table::rbindlist(lapply(mut_files, fread))

cat("done.\nReading score map...")
scores_chr <- fread(paste0("/../../media/gvbarroso/extradrive1/mut_rates_func/benegas/score_bins/benegas_bins_chr", chr, ".csv.gz"))
setnames(scores_chr, old="benegas_bin", new="benegas_class") # reserving "bin" to '1 kb bins'

cat("done.\nJoining maps...")
scores_chr <- mut_map[scores_chr, on=.(chrom, pos), nomatch=0] 

plot_df <- scores_chr[, .N, by=.(chrom, triplet, benegas_class)] 
fwrite(plot_df, paste0("summary_tbls/stacks_benegas_chr", chr, ".csv.gz"))

cat("done.\nSummarizing tables...")

# average mutation rates per class of constraint
scores_chr[, mean_class_roulette := mean(roulette, na.rm=T), by=.(benegas_class)]
scores_chr[, mean_class_gnomad := mean(gnomad, na.rm=T), by=.(benegas_class)]
scores_chr[, mean_class_carlson := mean(carlson, na.rm=T), by=.(benegas_class)]
scores_chr[, se_class_roulette := sd(roulette, na.rm=T) / sqrt(.N), by=.(benegas_class)]
scores_chr[, se_class_gnomad := sd(gnomad, na.rm=T) / sqrt(.N), by=.(benegas_class)]
scores_chr[, se_class_carlson := sd(carlson, na.rm=T) / sqrt(.N), by=.(benegas_class)]

# average mutation rates per class of constraint per triplet context
scores_chr[, mean_class_triplet_roulette := mean(roulette, na.rm=T), by=.(benegas_class, triplet)]
scores_chr[, mean_class_triplet_carlson := mean(carlson, na.rm=T), by=.(benegas_class, triplet)]
scores_chr[, mean_class_triplet_gnomad := mean(gnomad, na.rm=T), by=.(benegas_class, triplet)]
scores_chr[, se_class_triplet_roulette := sd(roulette, na.rm=T) / sqrt(.N), by=.(benegas_class, triplet)]
scores_chr[, se_class_triplet_carlson := sd(carlson, na.rm=T) / sqrt(.N), by=.(benegas_class, triplet)]
scores_chr[, se_class_triplet_gnomad := sd(gnomad, na.rm=T) / sqrt(.N), by=.(benegas_class, triplet)]
  
# storing number of sites per class per triplet (used to compute the means for focal chr)
scores_chr[, `:=`(num_sites_roulette=sum(!is.na(roulette)),
                  num_sites_carlson=sum(!is.na(carlson)),
                  num_sites_gnomad=sum(!is.na(gnomad))),
              by=.(benegas_class, triplet)]

unique_vals <- unique(scores_chr, by=c("benegas_class", "triplet"))
setorder(unique_vals, benegas_class, triplet)
cat("done.\nWriting to file...")

fwrite(unique_vals, paste0("summary_tbls/summaries_chr", chr, ".csv.gz"))
cat("Finished! Moving on to 1 kb windows...\n")

#########################
#
# 1 kb windows
#
########################

# keep only relevant columns
keep <- c("chrom", "pos", "roulette",  "carlson", "gnomad", "benegas_class")
tbl_chr <- scores_chr[, ..keep]

tbl_chr[, bin_1kb := pos %/% 1e3] # defining 1 kb bins

# computing summaries across 1 kb bins, stratified by benegas class
tbl_means <- tbl_chr[, .(
  mean_roulette = mean(roulette, na.rm=T),
  mean_carlson = mean(carlson, na.rm=T),
  mean_gnomad = mean(gnomad, na.rm=T),
  n_sites_class = .N), by = .(bin_1kb, benegas_class)]

# pivoting to wide format
wide <- dcast(tbl_means, bin_1kb ~ benegas_class, value.var = c("mean_roulette", "mean_carlson", "mean_gnomad", "n_sites_class"))
wide[, chrom := as.integer(chr)]

cat("Reading B-map...")
b_chr <- fread(paste0("../B_1kb_roulette/B_map_YRI_chr", chr, "_1kb.csv.gz"))
b_chr[, bin_1kb := pos %/% 1e3] # for joining

cat("done.\nJoining and sanitizing...")
tbl_chr <- wide[b_chr, on=.(chrom, bin_1kb), nomatch=0] 
tbl_chr[, bin_1kb := NULL]

# setting NaN's to NA
for(i in seq_along(tbl_chr)) {
  set(tbl_chr, which(is.nan(tbl_chr[[i]])), i, NA)
}

# NOTE since we are joining B-values (see ggplot p3 in plot_benegas.R) 
# it makes sense to compute rations within each 1 kb bin, then summarize them later
cat("done.\nComputing ratios...")

# class 15 is putatively neutral
for(i in 1:12) {
  tbl_chr[, paste0("ratio_roulette_class_", i) := get(paste0("mean_roulette_", i)) / mean_roulette_15]
  tbl_chr[, paste0("mean_roulette_", i) := NULL]
}
tbl_chr[, mean_roulette_15 := NULL]

for(i in 1:12) {
  tbl_chr[, paste0("ratio_carlson_class_", i) := get(paste0("mean_carlson_", i)) / mean_carlson_15]
  tbl_chr[, paste0("mean_carlson_", i) := NULL]
}
tbl_chr[, mean_carlson_15 := NULL]

for(i in 1:12) {
  tbl_chr[, paste0("ratio_gnomad_class_", i) := get(paste0("mean_gnomad_", i)) / mean_gnomad_15]
  tbl_chr[, paste0("mean_gnomad_", i) := NULL]
}
tbl_chr[, mean_gnomad_15 := NULL]

cat("done.\nRe-organizing table...")

# "transposing" table (onloy relevant columns)
tbl_inv <- pivot_longer(tbl_chr, cols=starts_with("ratio_")) %>% setDT()
tbl_inv[, benegas_group := as.integer(sub(".*_", "", name))]
tbl_inv[, variable := sub("_class.*", "", name)]

# getting 1 kb windows where only one class of benegas elements appear
## replace NA's with 0's
counts <- paste0("n_sites_class_", 1:12)
tbl_inv[, (counts) := lapply(.SD, function(x) fifelse(is.na(x), 0L, x)), .SDcols = counts]
tbl_inv[, sum_constrained := rowSums(.SD), .SDcols = counts]
tbl_inv[, map := sub(".*_", "", variable)]

single_benegas <- tbl_inv[tbl_inv[, do.call(pmax, c(.SD, na.rm=T)), .SDcols = counts] == sum_constrained]
single_benegas <- single_benegas[sum_constrained > 0 & !is.na(value),]

cat("done.\nSummarizing 1kb maps...")

single_benegas[, chrom := as.integer(chr)]
fwrite(single_benegas, paste0("summary_tbls/exclusive_1kb_chr", chr, ".csv.gz"))

cat("done.\nFinished unfiltered tables!\n")

#########################
#
# 1 kb windows (filtered by triplet to exclude CpG)
#
########################

# keep only relevant columns
keep <- c("chrom", "pos", "roulette",  "carlson", "gnomad", "benegas_class", "triplet")
tbl_chr <- scores_chr[, ..keep]

cat("Filtering out CpG sites...")

CpGs <- c("ACG", "CCG", "GCG", "TCG", "CGA", "CGC", "CGG", "CGT")
tbl_chr <- tbl_chr[!triplet %in% CpGs,]

cat("done.\nComputing summaries across 1 kb bins, stratified by functional class...")

tbl_chr[, triplet := NULL]
tbl_chr[, bin_1kb := pos %/% 1e3] # defining 1 kb bins

tbl_means <- tbl_chr[, .(
  mean_roulette = mean(roulette, na.rm=T),
  mean_carlson = mean(carlson, na.rm=T),
  mean_gnomad = mean(gnomad, na.rm=T),
  n_sites_class = .N), by = .(bin_1kb, benegas_class)]

cat("done.\nPivoting to wide format...")

wide <- dcast(tbl_means, bin_1kb ~ benegas_class, value.var = c("mean_roulette", "mean_carlson", "mean_gnomad", "n_sites_class"))
wide[, chrom := as.integer(chr)]

cat("done.\nReading B-map...")
b_chr <- fread(paste0("../B_1kb_roulette/B_map_YRI_chr", chr, "_1kb.csv.gz"))
b_chr[, bin_1kb := pos %/% 1e3] # for joining

cat("done.\nJoining and sanitizing...")
tbl_chr <- wide[b_chr, on=.(chrom, bin_1kb), nomatch=0] 
tbl_chr[, bin_1kb := NULL]

# setting NaN's to NA
for(i in seq_along(tbl_chr)) {
  set(tbl_chr, which(is.nan(tbl_chr[[i]])), i, NA)
}

# NOTE since we are joining B-values (see ggplot p3 in plot_phastcons.R) 
# it makes sense to compute rations within each 1 kb bin, then summarize them later
cat("done.\nComputing ratios...")

# class 15 is putatively neutral
for(i in 1:12) {
  tbl_chr[, paste0("ratio_roulette_class_", i) := get(paste0("mean_roulette_", i)) / mean_roulette_15]
  tbl_chr[, paste0("mean_roulette_", i) := NULL]
}
tbl_chr[, mean_roulette_15 := NULL]

for(i in 1:12) {
  tbl_chr[, paste0("ratio_carlson_class_", i) := get(paste0("mean_carlson_", i)) / mean_carlson_15]
  tbl_chr[, paste0("mean_carlson_", i) := NULL]
}
tbl_chr[, mean_carlson_15 := NULL]

for(i in 1:12) {
  tbl_chr[, paste0("ratio_gnomad_class_", i) := get(paste0("mean_gnomad_", i)) / mean_gnomad_15]
  tbl_chr[, paste0("mean_gnomad_", i) := NULL]
}
tbl_chr[, mean_gnomad_15 := NULL]

cat("done.\nRe-organizing table...")

# "transposing" table (only relevant columns)
tbl_inv <- pivot_longer(tbl_chr, cols=starts_with("ratio_")) %>% setDT()
tbl_inv[, benegas_group := as.integer(sub(".*_", "", name))]
tbl_inv[, variable := sub("_class.*", "", name)]

# getting 1 kb windows where only one class of benegas elements appear
counts <- paste0("n_sites_class_", 1:12)
tbl_inv[, (counts) := lapply(.SD, function(x) fifelse(is.na(x), 0L, x)), .SDcols = counts]
tbl_inv[, sum_constrained := rowSums(.SD), .SDcols = counts]
tbl_inv[, map := sub(".*_", "", variable)]

single_benegas <- tbl_inv[tbl_inv[, do.call(pmax, c(.SD, na.rm=T)), .SDcols = counts] == sum_constrained]
single_benegas <- single_benegas[sum_constrained > 0 & !is.na(value),]

cat("done.\nSummarizing 1kb maps...")

single_benegas[, chrom := as.integer(chr)]
fwrite(single_benegas, paste0("summary_tbls/exclusive_1kb_chr", chr, "_nonCpG.csv.gz"))

cat("done.\nFinished filtered tables!\n")
