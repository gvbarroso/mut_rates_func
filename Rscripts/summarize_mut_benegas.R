#########################
#
# Setting up
#
########################

library(data.table)
library(tidyverse)
library(cowplot)

args <- commandArgs(trailingOnly=TRUE)
chr <- args[1]

#########################
#
# Summarizing tables
#
########################

cat(paste0("Reading mutation maps for chr ", chr, "..."))
#mut_files <- list.files("split_muts/", pattern=paste0("^mut_map_chr", chr, "_"), full.names=T)
mut_files <- list.files("transfer/", pattern=paste0("^mut_map_chr", chr, "_"), full.names=T) # macbook
mut_map <- data.table::rbindlist(lapply(mut_files, fread))

cat("done.\nReading score map...")
#scores_chr <- fread(paste0("score_bins/benegas_bins_chr", chr, ".csv.gz"))
scores_chr <- fread(paste0("transfer/benegas_bins_chr", chr, ".csv.gz")) # macbook
  
cat("done.\nJoining maps...")
scores_chr <- mut_map[scores_chr, on=.(chrom, pos), nomatch=0] 

cat("done.\nPlotting stacks...")
plot_df <- scores_chr[, .N, by=.(triplet, benegas_bin)]
plot_df[, prop := N / sum(N), by=benegas_bin]
plot_df[, benegas_bin := factor(benegas_bin, levels=sort(unique(benegas_bin)))]

p <- ggplot(plot_df, aes(x=benegas_bin, y=prop, fill=triplet)) +
  theme_classic() + geom_col() +
  scale_fill_viridis_d(option="C", direction=1, guide=guide_legend(nrow=4), name=NULL) +
  scale_x_discrete(breaks=sort(unique(scores_chr$benegas_bin)), 
                   labels=c(as.character(sort(unique(scores_chr$benegas_bin))[-length(sort(unique(scores_chr$benegas_bin)))]), 
                            "Outside")) +
  labs(x="Benegas score bin", y="Proportion", fill="Trinucleotide") +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=18),
        legend.position="bottom",
        legend.box="horizontal")
save_plot(paste0("plots/benegas_triplets_chr", chr, ".pdf"), p, base_height=7, base_width=10)
  
cat("done.\nSummarizing tables...")

# average mutation rates per bin of constraint
scores_chr[, mean_bin_roulette := mean(roulette, na.rm=T), by=.(benegas_bin)]
scores_chr[, mean_bin_gnomad := mean(gnomad, na.rm=T), by=.(benegas_bin)]
scores_chr[, mean_bin_carlson := mean(carlson, na.rm=T), by=.(benegas_bin)]
scores_chr[, se_bin_roulette := sd(roulette, na.rm=T) / sqrt(.N), by=.(benegas_bin)]
scores_chr[, se_bin_gnomad := sd(gnomad, na.rm=T) / sqrt(.N), by=.(benegas_bin)]
scores_chr[, se_bin_carlson := sd(carlson, na.rm=T) / sqrt(.N), by=.(benegas_bin)]

# average mutation rates per bin of constraint per triplet context
scores_chr[, mean_bin_triplet_roulette := mean(roulette, na.rm=T), by=.(benegas_bin, triplet)]
scores_chr[, mean_bin_triplet_carlson := mean(carlson, na.rm=T), by=.(benegas_bin, triplet)]
scores_chr[, mean_bin_triplet_gnomad := mean(gnomad, na.rm=T), by=.(benegas_bin, triplet)]
scores_chr[, se_bin_triplet_roulette := sd(roulette, na.rm=T) / sqrt(.N), by=.(benegas_bin, triplet)]
scores_chr[, se_bin_triplet_carlson := sd(carlson, na.rm=T) / sqrt(.N), by=.(benegas_bin, triplet)]
scores_chr[, se_bin_triplet_gnomad := sd(gnomad, na.rm=T) / sqrt(.N), by=.(benegas_bin, triplet)]
  
scores_chr[, `:=`(num_sites_roulette = sum(!is.na(roulette)),
                  num_sites_carlson = sum(!is.na(carlson)),
                  num_sites_gnomad = sum(!is.na(gnomad))),
  by = .(benegas_bin, triplet)
]

unique_vals <- unique(scores_chr, by = c("benegas_bin", "triplet"))
setorder(unique_vals, benegas_bin, triplet)
cat("done.\nWriting to file...")

fwrite(unique_vals, paste0("summary_tbls/summaries_chr", chr, ".csv.gz"))
cat("Finished!\n")

#########################
#
# Moving on to more detailed analyses
#
########################



# TODO repeat after filtering CpG sites



# keep only relevant columns
keep <- c("chrom", "pos", "roulette",  "carlson", "gnomad", "benegas_bin")
scores_chr <- scores_chr[, ..keep]

scores_chr[, bin_1kb := pos %/% 1e3] # defining 1 kb bins

# computing summaries across 1 kb bins, stratified by benegas bin
tbl_means <- scores_chr[, .(
  mean_roulette = mean(roulette, na.rm = TRUE),
  mean_carlson = mean(carlson, na.rm = TRUE),
  mean_gnomad = mean(gnomad, na.rm = TRUE),
  n_sites_bin = .N
), by = .(bin_1kb, benegas_bin)]

# pivoting to wide format
wide <- dcast(tbl_means, bin_1kb ~ benegas_bin, value.var = c("mean_roulette", "mean_carlson", "mean_gnomad", "n_sites_bin"))
wide[, chrom := as.integer(chr)]

cat("Reading B-map...")
b_chr <- fread(paste0("B_1kb_roulette/B_map_YRI_chr", chr, "_1kb.csv.gz"))
b_chr[, bin_1kb := pos %/% 1e3] # for joining

cat("done.\nJoining and sanitizing...")
scores_chr <- wide[b_chr, on=.(chrom, bin_1kb), nomatch=0] 
scores_chr[, bin_1kb := NULL]

# setting NaN's to NA
for(i in seq_along(scores_chr)) {
  set(scores_chr, which(is.nan(scores_chr[[i]])), i, NA)
}

#tbl <- dplyr::select(scores_chr, c("chrom", "pos", "B", starts_with("n_sites"))) %>%
#    pivot_longer(., cols=starts_with("n_sites"), names_to="benegas_bin", values_to="count") %>% setDT()
#tbl[, benegas_bin := as.integer(sub(".*_", "", benegas_bin))]

# hist <- tbl %>% ggplot(aes(x=count, fill=as.factor(benegas_bin))) +
#   geom_histogram(alpha=0.3, binwidth=0.05, color="black", position="identity") +
#   scale_x_log10() + theme_bw() +
#   scale_fill_viridis_d(option="C", direction=1) +
#   labs(x="Total Length (max=1kb)", y="Density", fill=NULL) +
#   theme(panel.grid.minor=element_blank(),
#         axis.text=element_text(size=12),
#         axis.title=element_text(size=16),
#         axis.text.y=element_text(hjust=1),
#         legend.position="bottom")
# TODO stack plot

cat("done.\nComputing ratios...")

for(i in 1:12) {
  scores_chr[, paste0("ratio_roulette_bin_", i) := get(paste0("mean_roulette_", i)) / mean_roulette_15]
  scores_chr[, paste0("mean_roulette_", i) := NULL]
}
scores_chr[, mean_roulette_15 := NULL]

for(i in 1:12) {
  scores_chr[, paste0("ratio_carlson_bin_", i) := get(paste0("mean_carlson_", i)) / mean_carlson_15]
  scores_chr[, paste0("mean_carlson_", i) := NULL]
}
scores_chr[, mean_carlson_15 := NULL]

for(i in 1:12) {
  scores_chr[, paste0("ratio_gnomad_bin_", i) := get(paste0("mean_gnomad_", i)) / mean_gnomad_15]
  scores_chr[, paste0("mean_gnomad_", i) := NULL]
}
scores_chr[, mean_gnomad_15 := NULL]

cat("done.\nRe-organizing table...")

tb_inv <- pivot_longer(scores_chr, cols=starts_with("ratio_")) %>% setDT()
tb_inv[, benegas_group := as.integer(sub(".*_", "", name))]
tb_inv[, variable := sub("_bin.*", "", name)]

# getting 1 kb windows where only one class of benegas elements appear
counts <- paste0("n_sites_bin_", 1:12)
tb_inv[, (counts) := lapply(.SD, function(x) fifelse(is.na(x), 0L, x)), .SDcols = counts]
tb_inv[, sum_constrained := rowSums(.SD), .SDcols = counts]
tb_inv[, map := sub(".*_", "", variable)]

single_benegas <- tb_inv[tb_inv[, do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = counts] == sum_constrained]
single_benegas <- single_benegas[sum_constrained > 0 & !is.na(value),]

cat("done.\nSummarizing 1kb maps...")

# weighted average by num_sites of each benegas group
counts <- paste0("n_sites_bin_", 1:12)
single_benegas[, weight := as.matrix(.SD)[cbind(seq_len(.N), benegas_group)], .SDcols = counts]
single_benegas[, chrom := as.integer(chr)]
  
fwrite(single_benegas, paste0("summary_tbls/exclusive_1kb_chr", chr, ".csv.gz"))

cat("Finished!")
