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
mut_files <- list.files("chr22_muts/", pattern=paste0("^mut_map_chr22_"), full.names=T)
mut_map <- data.table::rbindlist(lapply(mut_files, fread))

# NOTE: unlike Benegas scores (which are per site), phastCons are given in BED intervals
cat("done.\nReading score map...")
phast <- vector("list", 12) # 12 phastcons classes
for(l in seq(0, 55, 5)) {
  
  u <- l + 5
  # path -> git clone https://github.com/nwcol/bgs_lmr into $HOME/Devel
  tmp <- fread(paste("~/Devel/bgs_lmr/data/annotations/phastcons/top", l, "-", u,
                     "/phastcons_top", l, "-", u, "_chr", chr, ".bed.gz", sep="")) 
  
  tmp$phast_class <- u / 5
  phast[[l %/% 5 + 1]] <- tmp
}
phast <- data.table::rbindlist(phast) %>% setDT()
phast[, chrom := as.integer(chr)]
phast[, chrom := as.integer(chrom)]

# putatively neutral sites
chr_range <- phast[,.(start=1, end=max(chromEnd)), by=chrom][, .(pos=seq(start, end)), by=chrom]
# rolling BED intervals to single-nucleotide positions to join mutation rates
phast <- phast[, .(pos=seq(chromStart, chromEnd)), by=.(chrom, phast_class, chromStart, chromEnd)][, c("chromStart", "chromEnd") := NULL]
phast <- phast[chr_range, on=.(chrom, pos)]
phast[is.na(phast_class), phast_class := 15L] # class 15 -> putatively neutral sites

cat("done.\nJoining maps...")
phast <- mut_map[phast, on=.(chrom, pos), nomatch=0] 

plot_df <- phast[, .N, by=.(chrom, triplet, phast_class)]
fwrite(plot_df, paste0("summary_tbls/stacks_phastcons_chr", chr, ".csv.gz"))

cat("done.\nSummarizing tables...")

# average mutation rates per class of constraint
phast[, mean_class_roulette := mean(roulette, na.rm=T), by=.(phast_class)]
phast[, mean_class_gnomad := mean(gnomad, na.rm=T), by=.(phast_class)]
phast[, mean_class_carlson := mean(carlson, na.rm=T), by=.(phast_class)]
phast[, se_class_roulette := sd(roulette, na.rm=T) / sqrt(.N), by=.(phast_class)]
phast[, se_class_gnomad := sd(gnomad, na.rm=T) / sqrt(.N), by=.(phast_class)]
phast[, se_class_carlson := sd(carlson, na.rm=T) / sqrt(.N), by=.(phast_class)]

# average mutation rates per class of constraint per triplet context
phast[, mean_class_triplet_roulette := mean(roulette, na.rm=T), by=.(phast_class, triplet)]
phast[, mean_class_triplet_carlson := mean(carlson, na.rm=T), by=.(phast_class, triplet)]
phast[, mean_class_triplet_gnomad := mean(gnomad, na.rm=T), by=.(phast_class, triplet)]
phast[, se_class_triplet_roulette := sd(roulette, na.rm=T) / sqrt(.N), by=.(phast_class, triplet)]
phast[, se_class_triplet_carlson := sd(carlson, na.rm=T) / sqrt(.N), by=.(phast_class, triplet)]
phast[, se_class_triplet_gnomad := sd(gnomad, na.rm=T) / sqrt(.N), by=.(phast_class, triplet)]

# storing number of sites per class per triplet (used to compute the means for each chr)
phast[, `:=`(num_sites_roulette=sum(!is.na(roulette)),
             num_sites_carlson=sum(!is.na(carlson)),
             num_sites_gnomad=sum(!is.na(gnomad))),
           by=.(phast_class, triplet)]

unique_vals <- unique(phast, by=c("phast_class", "triplet"))
setorder(unique_vals, phast_class, triplet)
cat("done.\nWriting to file...")

fwrite(unique_vals, paste0("summary_tbls/summaries_chr", chr, ".csv.gz"))
cat("Finished! Moving on to 1 kb windows...\n")

#########################
#
# 1 kb windows
#
########################

# keep only relevant columns
keep <- c("chrom", "pos", "roulette",  "carlson", "gnomad", "phast_class")
tbl_chr <- phast[, ..keep]

tbl_chr[, bin_1kb := pos %/% 1e3] # defining 1 kb bins

# computing summaries across 1 kb bins, stratified by phastcons group
tbl_means <- tbl_chr[, .(mean_roulette=mean(roulette, na.rm=T),
                         mean_carlson=mean(carlson, na.rm=T),
                         mean_gnomad=mean(gnomad, na.rm=T),
                         n_sites_class=.N), by=.(bin_1kb, phast_class)]

# pivoting to wide format
wide <- dcast(tbl_means, bin_1kb ~ phast_class, value.var=c("mean_roulette", "mean_carlson", "mean_gnomad", "n_sites_class"))
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

# NOTE since we are joining B-values (see ggplot p3 in plot_phastcons.R) 
# it makes sense to compute ratios within each 1 kb bin, then summarize them later
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
tbl_inv[, phast_group := as.integer(sub(".*_", "", name))]
tbl_inv[, variable := sub("_class.*", "", name)]

# getting 1 kb windows where only one class of phastcons elements appear
## replace NA's with 0's
counts <- paste0("n_sites_class_", 1:12)
tbl_inv[, (counts) := lapply(.SD, function(x) fifelse(is.na(x), 0L, x)), .SDcols=counts]
tbl_inv[, sum_constrained := rowSums(.SD), .SDcols=counts]
tbl_inv[, map := sub(".*_", "", variable)]

single_phast <- tbl_inv[tbl_inv[, do.call(pmax, c(.SD, na.rm=T)), .SDcols=counts] == sum_constrained]
single_phast <- single_phast[sum_constrained > 0 & !is.na(value),]

cat("done.\nSummarizing 1kb maps...")

single_phast[, chrom := as.integer(chr)]
fwrite(single_phast, paste0("summary_tbls/exclusive_1kb_chr", chr, ".csv.gz"))

cat("done.\nFinished unfiltered tables!\n")

#########################
#
# 1 kb windows (filtered by triplet to exclude CpG)
#
########################

# keep only relevant columns
keep <- c("chrom", "pos", "roulette",  "carlson", "gnomad", "phast_class", "triplet")
tbl_chr <- phast[, ..keep]

cat("Filtering out CpG sites...")

CpGs <- c("ACG", "CCG", "GCG", "TCG", "CGA", "CGC", "CGG", "CGT")
tbl_chr <- tbl_chr[!triplet %in% CpGs,]

cat("done.\nComputing summaries across 1 kb bins, stratified by functional class...")

tbl_chr[, triplet := NULL]
tbl_chr[, bin_1kb := pos %/% 1e3] # defining 1 kb bins

tbl_means <- tbl_chr[, .(
  mean_roulette=mean(roulette, na.rm=T),
  mean_carlson=mean(carlson, na.rm=T),
  mean_gnomad=mean(gnomad, na.rm=T),
  n_sites_class=.N), by=.(bin_1kb, phast_class)]

cat("done.\nPivoting to wide format...")

wide <- dcast(tbl_means, bin_1kb ~ phast_class, value.var=c("mean_roulette", "mean_carlson", "mean_gnomad", "n_sites_class"))
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
# it makes sense to compute ratios within each 1 kb bin, then summarize them later
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
tbl_inv[, phast_group := as.integer(sub(".*_", "", name))]
tbl_inv[, variable := sub("_class.*", "", name)]

# getting 1 kb windows where only one class of phastcons elements appear
counts <- paste0("n_sites_class_", 1:12)
tbl_inv[, (counts) := lapply(.SD, function(x) fifelse(is.na(x), 0L, x)), .SDcols=counts]
tbl_inv[, sum_constrained := rowSums(.SD), .SDcols=counts]
tbl_inv[, map := sub(".*_", "", variable)]

single_phast <- tbl_inv[tbl_inv[, do.call(pmax, c(.SD, na.rm=T)), .SDcols=counts] == sum_constrained]
single_phast <- single_phast[sum_constrained > 0 & !is.na(value),]

cat("done.\nSummarizing 1kb maps...")

single_phast[, chrom := as.integer(chr)]
fwrite(single_phast, paste0("summary_tbls/exclusive_1kb_chr", chr, "_nonCpG.csv.gz"))

cat("done.\nFinished filtered tables!\n")
