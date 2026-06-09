#########################
#
# Setting up
#
########################

library(data.table)
library(tidyverse)
library(cowplot)

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

# NOTE: unlike Benegas scores (which are per site), phastCons are given in BED intervals
cat("done.\nReading score map...")
phast <- vector("list", 12) # 12 phastcons bins
for(l in seq(0, 55, 5)) {
  
  u <- l + 5
  
  tmp <- fread(paste("~/Devel/bgs_lmr/data/annotations/phastcons/top", l, "-", u,
                     "/phastcons_top", l, "-", u, "_chr", chr, ".bed.gz", sep="")) 
  
  tmp$phast_bin <- u / 5
  phast[[l %/% 5 + 1]] <- tmp
}
phast <- data.table::rbindlist(phast) %>% setDT()
phast[, chrom := as.integer(chr)]
phast[, chrom := as.integer(chrom)]

# putatively neutral sites
chr_range <- phast[,.(start=1, end=max(chromEnd)), by=chrom][, .(pos=seq(start, end)), by=chrom]
# rolling BED intervals to single-nucleotide positions to join mutation rates
phast <- phast[, .(pos=seq(chromStart, chromEnd)), by=.(chrom, phast_bin, chromStart, chromEnd)][, c("chromStart", "chromEnd") := NULL]
phast <- phast[chr_range, on=.(chrom, pos)]
phast[is.na(phast_bin), phast_bin := 15L]

cat("done.\nJoining maps...")
phast <- mut_map[phast, on=.(chrom, pos), nomatch=0] 

cat("done.\nPlotting stacks...")
plot_df <- phast[, .N, by=.(triplet, phast_bin)]
plot_df[, prop := N / sum(N), by=phast_bin]
plot_df[, phast_bin := factor(phast_bin, levels=sort(unique(phast_bin)))]

p <- ggplot(plot_df, aes(x=phast_bin, y=prop, fill=triplet)) +
  theme_classic() + geom_col() +
  scale_fill_viridis_d(option="C", direction=1, guide=guide_legend(nrow=4), name=NULL) +
  scale_x_discrete(breaks=sort(unique(phast$phast_bin)), 
                   labels=c(as.character(sort(unique(phast$phast_bin))[-length(sort(unique(phast$phast_bin)))]), 
                            "Outside")) +
  labs(x="PhastCons bin", y="Proportion", fill="Trinucleotide") +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=18),
        legend.position="bottom",
        legend.box="horizontal")
save_plot(paste0("plots/phastcons_triplets_chr", chr, ".pdf"), p, base_height=7, base_width=10)

cat("done.\nSummarizing tables...")

# average mutation rates per bin of constraint
phast[, mean_bin_roulette := mean(roulette, na.rm=T), by=.(phast_bin)]
phast[, mean_bin_gnomad := mean(gnomad, na.rm=T), by=.(phast_bin)]
phast[, mean_bin_carlson := mean(carlson, na.rm=T), by=.(phast_bin)]
phast[, se_bin_roulette := sd(roulette, na.rm=T) / sqrt(.N), by=.(phast_bin)]
phast[, se_bin_gnomad := sd(gnomad, na.rm=T) / sqrt(.N), by=.(phast_bin)]
phast[, se_bin_carlson := sd(carlson, na.rm=T) / sqrt(.N), by=.(phast_bin)]

# average mutation rates per bin of constraint per triplet context
phast[, mean_bin_triplet_roulette := mean(roulette, na.rm=T), by=.(phast_bin, triplet)]
phast[, mean_bin_triplet_carlson := mean(carlson, na.rm=T), by=.(phast_bin, triplet)]
phast[, mean_bin_triplet_gnomad := mean(gnomad, na.rm=T), by=.(phast_bin, triplet)]
phast[, se_bin_triplet_roulette := sd(roulette, na.rm=T) / sqrt(.N), by=.(phast_bin, triplet)]
phast[, se_bin_triplet_carlson := sd(carlson, na.rm=T) / sqrt(.N), by=.(phast_bin, triplet)]
phast[, se_bin_triplet_gnomad := sd(gnomad, na.rm=T) / sqrt(.N), by=.(phast_bin, triplet)]

phast[, `:=`(num_sites_roulette = sum(!is.na(roulette)),
                  num_sites_carlson = sum(!is.na(carlson)),
                  num_sites_gnomad = sum(!is.na(gnomad))),
           by = .(phast_bin, triplet)
]

unique_vals <- unique(phast, by = c("phast_bin", "triplet"))
setorder(unique_vals, phast_bin, triplet)
cat("done.\nWriting to file...")

fwrite(unique_vals, paste0("summary_tbls/summaries_chr", chr, ".csv.gz"))
cat("Finished! Moving on to 1 kb windows...\n")

#########################
#
# 1 kb windows
#
########################

# keep only relevant columns
keep <- c("chrom", "pos", "roulette",  "carlson", "gnomad", "phast_bin")
tbl_chr <- phast[, ..keep]

tbl_chr[, bin_1kb := pos %/% 1e3] # defining 1 kb bins

# computing summaries across 1 kb bins, stratified by phastcons group
tbl_means <- tbl_chr[, .(
  mean_roulette = mean(roulette, na.rm=T),
  mean_carlson = mean(carlson, na.rm=T),
  mean_gnomad = mean(gnomad, na.rm=T),
  n_sites_bin = .N), by = .(bin_1kb, phast_bin)]

# pivoting to wide format
wide <- dcast(tbl_means, bin_1kb ~ phast_bin, value.var = c("mean_roulette", "mean_carlson", "mean_gnomad", "n_sites_bin"))
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

cat("done.\nComputing ratios...")

for(i in 1:12) {
  tbl_chr[, paste0("ratio_roulette_bin_", i) := get(paste0("mean_roulette_", i)) / mean_roulette_15]
  tbl_chr[, paste0("mean_roulette_", i) := NULL]
}
tbl_chr[, mean_roulette_15 := NULL]

for(i in 1:12) {
  tbl_chr[, paste0("ratio_carlson_bin_", i) := get(paste0("mean_carlson_", i)) / mean_carlson_15]
  tbl_chr[, paste0("mean_carlson_", i) := NULL]
}
tbl_chr[, mean_carlson_15 := NULL]

for(i in 1:12) {
  tbl_chr[, paste0("ratio_gnomad_bin_", i) := get(paste0("mean_gnomad_", i)) / mean_gnomad_15]
  tbl_chr[, paste0("mean_gnomad_", i) := NULL]
}
tbl_chr[, mean_gnomad_15 := NULL]

cat("done.\nRe-organizing table...")

tb_inv <- pivot_longer(tbl_chr, cols=starts_with("ratio_")) %>% setDT()
tb_inv[, phast_group := as.integer(sub(".*_", "", name))]
tb_inv[, variable := sub("_bin.*", "", name)]

# getting 1 kb windows where only one class of phastcons elements appear
counts <- paste0("n_sites_bin_", 1:12)
tb_inv[, (counts) := lapply(.SD, function(x) fifelse(is.na(x), 0L, x)), .SDcols = counts]
tb_inv[, sum_constrained := rowSums(.SD), .SDcols = counts]
tb_inv[, map := sub(".*_", "", variable)]

single_phast <- tb_inv[tb_inv[, do.call(pmax, c(.SD, na.rm=T)), .SDcols = counts] == sum_constrained]
single_phast <- single_phast[sum_constrained > 0 & !is.na(value),]

cat("done.\nSummarizing 1kb maps...")

# weighted average by num_sites of each phastcons group
counts <- paste0("n_sites_bin_", 1:12)
single_phast[, weight := as.matrix(.SD)[cbind(seq_len(.N), phast_group)], .SDcols = counts]
single_phast[, chrom := as.integer(chr)]

fwrite(single_phast, paste0("summary_tbls/exclusive_1kb_chr", chr, ".csv.gz"))

cat("done.\nFinished unfiltered tables!\n")

#########################
#
# 1 kb windows (filtered by triplet)
#
########################

# keep only relevant columns
keep <- c("chrom", "pos", "roulette",  "carlson", "gnomad", "phast_bin", "triplet")
tbl_chr <- phast[, ..keep]

cat("Filtering out CpG sites...")

CpGs <- c("ACG", "CCG", "GCG", "TCG", "CGA", "CGC", "CGG", "CGT")
tbl_chr <- tbl_chr[!triplet %in% CpGs,]

tbl_chr[, triplet := NULL]
tbl_chr[, bin_1kb := pos %/% 1e3] # defining 1 kb bins

# computing summaries across 1 kb bins, stratified by phastcons bin
tbl_means <- tbl_chr[, .(
  mean_roulette = mean(roulette, na.rm=T),
  mean_carlson = mean(carlson, na.rm=T),
  mean_gnomad = mean(gnomad, na.rm=T),
  n_sites_bin = .N), by = .(bin_1kb, phast_bin)]

# pivoting to wide format
wide <- dcast(tbl_means, bin_1kb ~ phast_bin, value.var = c("mean_roulette", "mean_carlson", "mean_gnomad", "n_sites_bin"))
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

cat("done.\nComputing ratios...")

for(i in 1:12) {
  tbl_chr[, paste0("ratio_roulette_bin_", i) := get(paste0("mean_roulette_", i)) / mean_roulette_15]
  tbl_chr[, paste0("mean_roulette_", i) := NULL]
}
tbl_chr[, mean_roulette_15 := NULL]

for(i in 1:12) {
  tbl_chr[, paste0("ratio_carlson_bin_", i) := get(paste0("mean_carlson_", i)) / mean_carlson_15]
  tbl_chr[, paste0("mean_carlson_", i) := NULL]
}
tbl_chr[, mean_carlson_15 := NULL]

for(i in 1:12) {
  tbl_chr[, paste0("ratio_gnomad_bin_", i) := get(paste0("mean_gnomad_", i)) / mean_gnomad_15]
  tbl_chr[, paste0("mean_gnomad_", i) := NULL]
}
tbl_chr[, mean_gnomad_15 := NULL]

cat("done.\nRe-organizing table...")

tb_inv <- pivot_longer(tbl_chr, cols=starts_with("ratio_")) %>% setDT()
tb_inv[, phast_group := as.integer(sub(".*_", "", name))]
tb_inv[, variable := sub("_bin.*", "", name)]

# getting 1 kb windows where only one class of phastcons elements appear
counts <- paste0("n_sites_bin_", 1:12)
tb_inv[, (counts) := lapply(.SD, function(x) fifelse(is.na(x), 0L, x)), .SDcols = counts]
tb_inv[, sum_constrained := rowSums(.SD), .SDcols = counts]
tb_inv[, map := sub(".*_", "", variable)]

single_phast <- tb_inv[tb_inv[, do.call(pmax, c(.SD, na.rm=T)), .SDcols = counts] == sum_constrained]
single_phast <- single_phast[sum_constrained > 0 & !is.na(value),]

cat("done.\nSummarizing 1kb maps...")

# weighted average by num_sites of each benegas group
counts <- paste0("n_sites_bin_", 1:12)
single_phast[, weight := as.matrix(.SD)[cbind(seq_len(.N), phast_group)], .SDcols = counts]
single_phast[, chrom := as.integer(chr)]

fwrite(single_phast, paste0("summary_tbls/exclusive_1kb_chr", chr, "_nonCpG.csv.gz"))

cat("done. Finished filtered tables!\n")
