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
  
fwrite(unique_vals, paste0("summary_tbls/summaries_chr", chr, ".csv.gz"))

