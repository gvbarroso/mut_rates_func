library(data.table)
library(tidyverse)
library(cowplot)

args <- commandArgs(trailingOnly=TRUE)
chr <- args[1]

mut_files <- list.files("transfer/", pattern=paste0("mut_map_chr", chr, "*"), full.names=T)
mut_map <- data.table::rbindlist(lapply(mut_files, fread))

scores_chr <- fread(paste0("score_bins/benegas_bins_chr", chr, ".csv.gz"))
# joins conservation scores, coverage may differ between muts and scores 
scores_chr <- mut_map[scores_chr, on=.(chrom, pos), nomatch=0] 
  
# stacked plot x=benegas_bin, colors are triplets
plot_df <- scores_chr[, .N, by=.(triplet, benegas_bin)]
plot_df[, prop := N / sum(N), by=benegas_bin]
plot_df[, benegas_bin := factor(benegas_bin, levels=sort(unique(benegas_bin)))]

p <- ggplot(plot_df, aes(x=benegas_bin, y=prop, fill=triplet)) +
  theme_classic() + geom_col() +
  scale_fill_viridis_d(option="C", direction=1, guide=guide_legend(nrow=4), name=NULL) +
  scale_x_discrete(breaks=c(1:nbins, nbins + 3), labels=c(as.character(1:nbins), "Outside")) +
  labs(x="Benegas score bin", y="Proportion", fill="Trinucleotide") +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=18),
        legend.position="bottom",
        legend.box="horizontal")
save_plot(paste0("plots/benegas_triplets_chr", chr, ".pdf"), p, base_height=7, base_width=10)
  
scores_chr[, mean_roulette := mean(roulette, na.rm=T), by=.(benegas_bin)]
scores_chr[, mean_gnomad := mean(gnomad, na.rm=T), by=.(benegas_bin)]
scores_chr[, mean_carlson := mean(carlson, na.rm=T), by=.(benegas_bin)]
scores_chr[, se_roulette := sd(roulette, na.rm=T) / sqrt(.N), by=.(benegas_bin)]
scores_chr[, se_gnomad := sd(gnomad, na.rm=T) / sqrt(.N), by=.(benegas_bin)]
scores_chr[, se_carlson := sd(carlson, na.rm=T) / sqrt(.N), by=.(benegas_bin)]
  
scores_chr[, num_sites := .N, by=.(benegas_bin, triplet, mean_roulette, mean_carlson, mean_gnomad, se_gnomad, se_carlson, se_roulette)]
  
unique_vals <- scores_chr[, .(benegas_bin, num_sites, triplet, mean_roulette, mean_carlson, mean_gnomad, se_gnomad, se_carlson, se_roulette)]
unique_vals <- unique(unique_vals)
setorder(unique_vals, benegas_bin, triplet)
  
fwrite(unique_vals, paste0("summary_tbls/summaries_chr", chr, ".csv.gz"))

