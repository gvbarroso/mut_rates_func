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

cat("done.\nReading functional maps...")
# exons split by decile of constraint (Zeng et al 2024 NatGen)
deciles <- vector("list", 11)
for(d in 1:11) {
  # path -> git clone https://github.com/nwcol/bgs_lmr into $HOME/Devel
  tmp <- fread(paste0("~/Devel/bgs_lmr/data/annotations/cds/lof_d", d, "/", "cds_lof_d", d, "_chr", chr, ".bed.gz"))
  
  tmp[, functional_class := paste0("decile_", d)]
  deciles[[d]] <- tmp
}
deciles <- data.table::rbindlist(deciles)
names(deciles)[1:3] <- c("chrom", "chromStart", "chromEnd")

# promoters, path -> git clone https://github.com/nwcol/bgs_lmr into $HOME/Devel
promoters <- fread(paste0("~/Devel/bgs_lmr/data/annotations/regulatory/promoters/promoters_chr", chr, ".bed.gz"))
promoters[, functional_class := "promoter"]
promoters[, chrom := as.integer(chr)]
promoters[, chrom := as.integer(chrom)]

# enhancers, # path -> git clone https://github.com/nwcol/bgs_lmr into $HOME/Devel
enhancers <- fread(paste0("~/Devel/bgs_lmr/data/annotations/regulatory/enhancers/enhancers_chr", chr, ".bed.gz"))
enhancers[, functional_class := "enhancer"]
enhancers[, chrom := as.integer(chr)]
enhancers[, chrom := as.integer(chrom)]

# rolling intervals
deciles <- deciles[, .(pos=seq(chromStart, chromEnd)), by=.(chrom, functional_class, chromStart, chromEnd)][, c("chromStart", "chromEnd") := NULL]
enhancers <- enhancers[, .(pos=seq(chromStart, chromEnd)), by=.(chrom, functional_class, chromStart, chromEnd)][, c("chromStart", "chromEnd") := NULL]
promoters <- promoters[, .(pos=seq(chromStart, chromEnd)), by=.(chrom, functional_class, chromStart, chromEnd)][, c("chromStart", "chromEnd") := NULL]

functional <- rbind.data.frame(deciles, enhancers, promoters)
setorder(functional, chrom, pos, functional_class)
setcolorder(functional, c("chrom", "pos", "functional_class"))

# after all functional functional_classents are loaded, identify putatively neutral sites
chr_range <- functional[, .(start=1, end=max(pos)), by=chrom][, .(pos=seq(start, end)), by=chrom]
functional <- functional[chr_range, on=.(chrom, pos)]
functional[is.na(functional_class), functional_class := "neutral"]

functional[, functional_class := factor(functional_class, levels=c(paste0("decile_", 1:11), "enhancer", "promoter", "neutral"))]

cat("done.\nJoining maps...")
functional <- mut_map[functional, on=.(chrom, pos), nomatch=0] 

plot_df <- functional[, .N, by=.(chrom, triplet, functional_class)]
fwrite(plot_df, paste0("summary_tbls/stacks_functional_chr", chr, ".csv.gz"))

cat("done.\nSummarizing tables...")

# average mutation rates per class of constraint
functional[, mean_class_roulette := mean(roulette, na.rm=T), by=.(functional_class)]
functional[, mean_class_gnomad := mean(gnomad, na.rm=T), by=.(functional_class)]
functional[, mean_class_carlson := mean(carlson, na.rm=T), by=.(functional_class)]
functional[, se_class_roulette := sd(roulette, na.rm=T) / sqrt(.N), by=.(functional_class)]
functional[, se_class_gnomad := sd(gnomad, na.rm=T) / sqrt(.N), by=.(functional_class)]
functional[, se_class_carlson := sd(carlson, na.rm=T) / sqrt(.N), by=.(functional_class)]

# average mutation rates per class of constraint per triplet context
functional[, mean_class_triplet_roulette := mean(roulette, na.rm=T), by=.(functional_class, triplet)]
functional[, mean_class_triplet_carlson := mean(carlson, na.rm=T), by=.(functional_class, triplet)]
functional[, mean_class_triplet_gnomad := mean(gnomad, na.rm=T), by=.(functional_class, triplet)]
functional[, se_class_triplet_roulette := sd(roulette, na.rm=T) / sqrt(.N), by=.(functional_class, triplet)]
functional[, se_class_triplet_carlson := sd(carlson, na.rm=T) / sqrt(.N), by=.(functional_class, triplet)]
functional[, se_class_triplet_gnomad := sd(gnomad, na.rm=T) / sqrt(.N), by=.(functional_class, triplet)]

# storing number of sites per class per triplet (used to compute the means for focal chr)
functional[, `:=`(num_sites_roulette=sum(!is.na(roulette)),
                  num_sites_carlson=sum(!is.na(carlson)),
                  num_sites_gnomad=sum(!is.na(gnomad))),
              by=.(functional_class, triplet)]

unique_vals <- unique(functional, by=c("functional_class", "triplet"))
setorder(unique_vals, functional_class, triplet)
cat("done.\nWriting to file...")

fwrite(unique_vals, paste0("summary_tbls/summaries_chr", chr, ".csv.gz"))
cat("Finished! Moving on to 1 kb windows...\n")

#########################
#
# 1 kb windows
#
########################

# keep only relevant columns
keep <- c("chrom", "pos", "roulette",  "carlson", "gnomad", "functional_class")
tbl_chr <- functional[, ..keep]

tbl_chr[, bin_1kb := pos %/% 1e3] # defining 1 kb bins

# computing summaries across 1 kb bins, stratified by functional class
tbl_means <- tbl_chr[, .(mean_roulette=mean(roulette, na.rm=T),
                         mean_carlson=mean(carlson, na.rm=T),
                         mean_gnomad=mean(gnomad, na.rm=T),
                         n_sites_class=.N), by=.(bin_1kb, functional_class)]

# pivoting to wide format
wide <- dcast(tbl_means, bin_1kb ~ functional_class, value.var=c("mean_roulette", "mean_carlson", "mean_gnomad", "n_sites_class"))
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

# NOTE since we are joining B-values (see ggplot p3 in plot_functional.R) 
# it makes sense to compute rations within each 1 kb bin, then summarize them later
cat("done.\nComputing ratios...")
functional_classs <- c(paste0("decile_", 1:11), "enhancer", "promoter")

for(i in functional_classs) {
  tbl_chr[, paste0("ratio_roulette_class_", i) := get(paste0("mean_roulette_", i)) / mean_roulette_neutral]
  tbl_chr[, paste0("mean_roulette_", i) := NULL]
}
tbl_chr[, mean_roulette_neutral := NULL]

for(i in functional_classs) {
  tbl_chr[, paste0("ratio_carlson_class_", i) := get(paste0("mean_carlson_", i)) / mean_carlson_neutral]
  tbl_chr[, paste0("mean_carlson_", i) := NULL]
}
tbl_chr[, mean_carlson_neutral := NULL]

for(i in functional_classs) {
  tbl_chr[, paste0("ratio_gnomad_class_", i) := get(paste0("mean_gnomad_", i)) / mean_gnomad_neutral]
  tbl_chr[, paste0("mean_gnomad_", i) := NULL]
}
tbl_chr[, mean_gnomad_neutral := NULL]

cat("done.\nRe-organizing table...")

# "transposing" table (only relevant columns)
tbl_inv <- pivot_longer(tbl_chr, cols=starts_with("ratio_")) %>% setDT()
tbl_inv[, func_class := sub("ratio_*.*_class_", "", name)]
tbl_inv[, variable := sub("_class.*", "", name)]

# getting 1 kb windows where only one class of functional functional_classents appear
## replace NA's with 0's
counts <- paste0("n_sites_class_", functional_classs)
tbl_inv[, (counts) := lapply(.SD, function(x) fifelse(is.na(x), 0L, x)), .SDcols=counts]
tbl_inv[, sum_constrained := rowSums(.SD), .SDcols=counts]
tbl_inv[, map := sub(".*_", "", variable)]

single_func <- tbl_inv[tbl_inv[, do.call(pmax, c(.SD, na.rm=T)), .SDcols=counts] == sum_constrained]
single_func <- single_func[sum_constrained > 0 & !is.na(value),]

cat("done.\nSummarizing 1kb maps...")

single_func <- single_func[func_class != "neutral",]
single_func[, chrom := as.integer(chr)]
fwrite(single_func, paste0("summary_tbls/exclusive_1kb_chr", chr, ".csv.gz"))

cat("done.\nFinished unfiltered tables!\n")

#########################
#
# 1 kb windows (filtered by triplet)
#
########################

# keep only relevant columns
keep <- c("chrom", "pos", "roulette",  "carlson", "gnomad", "functional_class", "triplet")
tbl_chr <- functional[, ..keep]

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
  n_sites_class=.N), by=.(bin_1kb, functional_class)]

cat("done.\nPivoting to wide format...")

wide <- dcast(tbl_means, bin_1kb ~ functional_class, value.var=c("mean_roulette", "mean_carlson", "mean_gnomad", "n_sites_class"))
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

# NOTE since we are joining B-values (see ggplot p3 in plot_functional.R) 
# it makes sense to compute rations within each 1 kb bin, then summarize them later
cat("done.\nComputing ratios...")
functional_classs <- c(paste0("decile_", 1:11), "enhancer", "promoter")

for(i in functional_classs) {
  tbl_chr[, paste0("ratio_roulette_class_", i) := get(paste0("mean_roulette_", i)) / mean_roulette_neutral]
  tbl_chr[, paste0("mean_roulette_", i) := NULL]
}
tbl_chr[, mean_roulette_neutral := NULL]

for(i in functional_classs) {
  tbl_chr[, paste0("ratio_carlson_class_", i) := get(paste0("mean_carlson_", i)) / mean_carlson_neutral]
  tbl_chr[, paste0("mean_carlson_", i) := NULL]
}
tbl_chr[, mean_carlson_neutral := NULL]

for(i in functional_classs) {
  tbl_chr[, paste0("ratio_gnomad_class_", i) := get(paste0("mean_gnomad_", i)) / mean_gnomad_neutral]
  tbl_chr[, paste0("mean_gnomad_", i) := NULL]
}
tbl_chr[, mean_gnomad_neutral := NULL]

cat("done.\nRe-organizing table...")

# "transposing" table (only relevant columns)
tbl_inv <- pivot_longer(tbl_chr, cols=starts_with("ratio_")) %>% setDT()
tbl_inv[, func_class := sub("ratio_*.*_class_", "", name)]
tbl_inv[, variable := sub("_class.*", "", name)]

# getting 1 kb windows where only one class of functional functional_classents appear
counts <- paste0("n_sites_class_", functional_classs)
tbl_inv[, (counts) := lapply(.SD, function(x) fifelse(is.na(x), 0L, x)), .SDcols=counts]
tbl_inv[, sum_constrained := rowSums(.SD), .SDcols=counts]
tbl_inv[, map := sub(".*_", "", variable)]

single_func <- tbl_inv[tbl_inv[, do.call(pmax, c(.SD, na.rm=T)), .SDcols=counts] == sum_constrained]
single_func <- single_func[sum_constrained > 0 & !is.na(value),]

cat("done.\nSummarizing 1kb maps...")

single_func <- single_func[func_class != "neutral",]
single_func[, chrom := as.integer(chr)]
fwrite(single_func, paste0("summary_tbls/exclusive_1kb_chr", chr, "_nonCpG.csv.gz"))

cat("done.\nFinished filtered tables!\n")
