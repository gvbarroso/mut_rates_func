library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]

mut_map <- fread(paste0(chr, "_rate_v5.2_TFBS_correction_all.vcf.bgz"))
mut_map[, c("ID", "REF", "ALT", "QUAL", "FILTER") := NULL] 
names(mut_map) <- c("chrom", "pos", "rates")

cat(paste0("Done reading mutation maps for chr", chr, "."))

mut_map[, roulette := as.numeric(sub(".*MR=([0-9.]+).*", "\\1", rates))]
mut_map[, carlson := as.numeric(sub(".*MC=([0-9.]+).*", "\\1", rates))]
mut_map[, gnomad := as.numeric(sub(".*MG=([0-9.]+).*", "\\1", rates))]
mut_map[, triplet := sub(".*PN=([^;]+).*", "\\1", rates)] # 5-mer
mut_map[, triplet := substr(triplet, 2, nchar(triplet) - 1)] # keeps inner 3-mer
mut_map[, rates := NULL]

# https://github.com/vseplyarskiy/Roulette/tree/main/adding_mutation_rate
roulette_scale <- 1.015e-7 / 2 
gnomad_scale <- roulette_scale
carlson_scale <- 2.086e-9 / 2

mut_map[, roulette := roulette * roulette_scale]
mut_map[, carlson := carlson * carlson_scale]
mut_map[, gnomad := gnomad * gnomad_scale]

mut_map[, roulette := mean(roulette, na.rm=T), by=.(chrom, pos)]
mut_map[, carlson := mean(carlson, na.rm=T), by=.(chrom, pos)]
mut_map[, gnomad := mean(gnomad, na.rm=T), by=.(chrom, pos)]

mut_map <- mut_map[!duplicated(pos)] 

cat("Done averaging mutation rates for sites. Now reading benegas files...")

# original scores downloaded from:
# https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores/resolve/main/scores.tsv.bgz
files <- list.files("split_scores/", pattern=paste0("scores_chr", chr, "*"), full.names=T)
pb <- txtProgressBar(min=0, max=length(files), style=3)
for(i in 1:length(files)) {
  
  setTxtProgressBar(pb, i)
  f <- files[i]
  
  benegas <- fread(f)
  interval <- sub(".*_(\\d+-\\d+).*", "\\1", f)
  # joins conservation scores, mut_map coverage may reduce site count
  dat <- mut_map[benegas, on=.(chrom, pos), nomatch=0] 
  
  if(nrow(dat)==0L) {
    next
  }
  
  # splitting top 30% of sites into bins, treat the other 70% as neutral
  threshold <- quantile(dat$benegas_score, 0.30)
  dat[, is_benegas := ifelse(benegas_score <= threshold, T, F)] # more negative -> more constrained
  
  nbins <- 12 # discretizing distribution of constraint scores
  dat[is_benegas==T, benegas_bin := as.numeric(cut(frank(benegas_score, ties.method="average") / .N,
                                                   breaks=seq(0, 1, by = 1 / nbins), labels = 1:nbins, include.lowest=T))]
  dat[is_benegas==F, benegas_bin := nbins + 3] 
  dat[, c("is_benegas", "benegas_score") := NULL] 
  
  #last_pos <- max(dat$pos) # marking last position in data to process phastcons more efficiently
  # stacked plot x=benegas_bin, colors are triplets
  # plot_df <- dat[, .N, by=.(triplet, benegas_bin)]
  # plot_df[, prop := N / sum(N), by=benegas_bin]
  # plot_df[, benegas_bin := factor(benegas_bin, levels=sort(unique(benegas_bin)))]
  # 
  # s1 <- ggplot(plot_df, aes(x=benegas_bin, y=prop, fill=triplet)) +
  #   theme_classic() + geom_col() +
  #   scale_fill_viridis_d(option="C", direction=1, guide=guide_legend(nrow=4), name=NULL) +
  #   scale_x_discrete(breaks=c(1:nbins, nbins + 3), labels=c(as.character(1:nbins), "Outside")) +
  #   labs(x="Benegas score bin", y="Proportion", fill="Trinucleotide") +
  #   theme(panel.grid=element_blank(),
  #         axis.text=element_text(size=14),
  #         axis.title=element_text(size=18),
  #         legend.position="bottom",
  #         legend.box="horizontal")
  # save_plot("~/Desktop/mut_rates/benegas_triplets.pdf", s1, base_height=7, base_width=10)
  
  dat[, mean_roulette := mean(roulette, na.rm=T), by=.(benegas_bin)]
  dat[, mean_gnomad := mean(gnomad, na.rm=T), by=.(benegas_bin)]
  dat[, mean_carlson := mean(carlson, na.rm=T), by=.(benegas_bin)]
  dat[, se_roulette := sd(roulette, na.rm=T) / sqrt(.N), by=.(benegas_bin)]
  dat[, se_gnomad := sd(gnomad, na.rm=T) / sqrt(.N), by=.(benegas_bin)]
  dat[, se_carlson := sd(carlson, na.rm=T) / sqrt(.N), by=.(benegas_bin)]
  
  dat[, num_sites := .N, by=.(benegas_bin, triplet, mean_roulette, mean_carlson, mean_gnomad, se_gnomad, se_carlson, se_roulette)]
  
  unique_vals <- dat[, .(benegas_bin, num_sites, triplet, mean_roulette, mean_carlson, mean_gnomad, se_gnomad, se_carlson, se_roulette)]
  unique_vals <- unique(unique_vals)
  setorder(unique_vals, benegas_bin, triplet)
  
  fwrite(unique_vals, paste0("summary_tbls/benegas_mut_summaries_", interval, ".csv.gz"))
}
close(pb)