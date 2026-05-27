library(data.table)
library(tidyverse)

# original scores downloaded from:
# https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores/resolve/main/scores.tsv.bgz

args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]

benegas <- fread(file_name)

roulette_file <- 
# https://github.com/vseplyarskiy/Roulette/tree/main/adding_mutation_rate
roulette_scale <- 1.015e-7 / 2 
gnomad_scale <- roulette_scale
carlson_scale <- 2.086e-9 / 2

mut_map <- fread(roulette_file)
mut_map[, c("V3", "V4", "V5", "V6", "V7") := NULL] 
names(mut_map) <- c("chrom", "pos", "rates")

mut_map[, roulette := as.numeric(sub(".*MR=([0-9.]+).*", "\\1", rates))]
mut_map[, carlson := as.numeric(sub(".*MC=([0-9.]+).*", "\\1", rates))]
mut_map[, gnomad := as.numeric(sub(".*MG=([0-9.]+).*", "\\1", rates))]
mut_map[, triplet := sub(".*PN=([^;]+).*", "\\1", rates)] # 5-mer
mut_map[, triplet := substr(triplet, 2, nchar(triplet) - 1)] # keeps inner 3-mer
mut_map[, rates := NULL]

mut_map[, roulette := roulette * roulette_scale]
mut_map[, carlson := carlson * carlson_scale]
mut_map[, gnomad := gnomad * gnomad_scale]

mut_map[, roulette := mean(roulette, na.rm=T), by=.(chrom, pos)]
mut_map[, carlson := mean(carlson, na.rm=T), by=.(chrom, pos)]
mut_map[, gnomad := mean(gnomad, na.rm=T), by=.(chrom, pos)]

mut_map <- mut_map[!duplicated(pos)] 

dat <- mut_map[benegas, on=.(chrom, pos), nomatch=0] # joins conservation scores

# splitting top 30% of sites into bins, treat the other 70% as neutral
threshold <- quantile(dat$benegas_score, 0.30)
dat[, is_benegas := ifelse(benegas_score <= threshold, T, F)] # more negative -> more constrained

nbins <- 12 # discretizing distribution of constraint scores
dat[is_benegas==T, benegas_bin := as.numeric(cut(frank(benegas_score, ties.method="average") / .N,
                                                 breaks=seq(0, 1, by = 1 / nbins), labels = 1:nbins, include.lowest=T))]
dat[is_benegas==F, benegas_bin := nbins + 3] 
dat[, c("is_benegas", "benegas_score") := NULL] 

last_pos <- max(dat$pos) # marking last position in data to process phastcons more efficiently

rm(mut_map, benegas)
gc()

# stacked plot x=benegas_bin, colors are triplets
plot_df <- dat[, .N, by=.(triplet, benegas_bin)]
plot_df[, prop := N / sum(N), by=benegas_bin]
plot_df[, benegas_bin := factor(benegas_bin, levels=sort(unique(benegas_bin)))]

s1 <- ggplot(plot_df, aes(x=benegas_bin, y=prop, fill=triplet)) +
  theme_classic() + geom_col() +
  scale_fill_viridis_d(option="C", direction=1, guide=guide_legend(nrow=4), name=NULL) +
  scale_x_discrete(breaks=c(1:nbins, nbins + 3), labels=c(as.character(1:nbins), "Outside")) +
  labs(x="Benegas score bin", y="Proportion", fill="Trinucleotide") +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=18),
        legend.position="bottom",
        legend.box="horizontal")
save_plot("~/Desktop/mut_rates/benegas_triplets.pdf", s1, base_height=7, base_width=10)

# plot of mean rates per benegas_bin
dat[, mean_roulette := mean(roulette, na.rm=T), by=.(benegas_bin)]
dat[, mean_gnomad := mean(gnomad, na.rm=T), by=.(benegas_bin)]
dat[, mean_carlson := mean(carlson, na.rm=T), by=.(benegas_bin)]
dat[, se_roulette := sd(roulette, na.rm=T) / sqrt(.N), by=.(benegas_bin)]
dat[, se_gnomad := sd(gnomad, na.rm=T) / sqrt(.N), by=.(benegas_bin)]
dat[, se_carlson := sd(carlson, na.rm=T) / sqrt(.N), by=.(benegas_bin)]

unique_vals <- dat[, .(benegas_bin, mean_roulette, mean_carlson, mean_gnomad, se_gnomad, se_carlson, se_roulette)]
unique_vals <- unique(unique_vals)

dt_m <- unique_vals %>%
  pivot_longer(cols=c(mean_roulette, mean_carlson, mean_gnomad, se_roulette, se_carlson, se_gnomad),
               names_to=c(".value", "source"), names_pattern="(mean|se)_(.*)")

p1 <- ggplot(dt_m, aes(x=benegas_bin, y=mean, color=source, group=paste0(source, benegas_bin < 15))) +
  geom_line(linewidth=1) + geom_point(size=3) + 
  geom_errorbar(aes(ymin=mean - se, ymax=mean + se), width=0.2, linewidth=0.6) +
  scale_x_continuous(breaks=c(1:nbins, nbins + 3), labels=c(as.character(1:nbins), "Outside")) + theme_classic() + 
  scale_color_manual(values=c("cyan3", "brown1", "green4"), name=NULL) +
  labs(x="Conservation score (percentile bin)", y="Mean Rate", color="Map",
       title="Mutation rates per conservation score") +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=14),
        strip.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.position="bottom")
save_plot("~/Desktop/mut_rates/benegas_mut_rates.pdf", p1, base_height=6, base_width=10)