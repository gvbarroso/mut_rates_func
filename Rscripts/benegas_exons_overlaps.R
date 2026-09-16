#########################
#
# Setting up
#
########################

pdf(NULL)

library(data.table)
library(tidyverse)
library(scales)
library(cowplot)

args <- commandArgs(trailingOnly=T)
chr <- args[1]

#########################
#
# Finding overlaps between Benegas scores (classes) and exonic deciles (LOF -> sh)
#
########################

cat("\nReading Benegas score map...")
scores_chr <- fread(paste0("/../../media/gvbarroso/extradrive1/mut_rates_func/benegas/score_bins/benegas_bins_chr", chr, ".csv.gz"))
setnames(scores_chr, old="benegas_bin", new="benegas_class") # reserving "bin" to '1 kb bins'
scores_chr <- scores_chr[benegas_class < 15,] # assuming a classification w/ 1-12 -> constrained, 15 -> neutral

cat("done.\nReading exonic maps...")
# exons split by decile of constraint (Zeng et al 2024 NatGen)
deciles <- vector("list", 10)
for(d in 1:10) {
  # path -> git clone https://github.com/nwcol/bgs_lmr into $HOME/Devel
  tmp <- fread(paste0("~/Devel/bgs_lmr/data/annotations/cds/lof_d", d, "/", "cds_lof_d", d, "_chr", chr, ".bed.gz"))
  
  tmp[, decile := d]
  deciles[[d]] <- tmp
}
deciles <- data.table::rbindlist(deciles)
names(deciles)[1:3] <- c("chrom", "chromStart", "chromEnd")

# rolling intervals
deciles <- deciles[, .(pos=seq(chromStart + 1L, chromEnd)), by=.(chrom, decile, chromStart, chromEnd)][, c("chromStart", "chromEnd") := NULL]
# de-duplicates keeping the lowest decile (minor exonic overlaps)
min_dec <- deciles[, .(decile = min(decile)), by = .(chrom, pos)]
deciles <- deciles[min_dec, on = .(chrom, pos, decile)]

cat("done.\nJoining maps...")
overlaps <- merge(scores_chr, deciles, by=c("chrom", "pos"), all=F) # keeping just the overlapping sites
cat("done.\nWriting and plotting...\n")

plot_df_overlaps <- overlaps[, .N, by=.(benegas_class, decile)]
plot_df_overlaps[, decile := factor(as.character(decile), levels = as.character(1:10))]
plot_df_overlaps[, benegas_class := factor(as.character(benegas_class), levels = as.character(1:12))]
plot_df_overlaps[, chr := chr]

# writing to compute stratified ratios in a dedicated script
overlaps[, overlap_benegas_exon := T] # set all to TRUE (by construction)
fwrite(overlaps, paste0("summary_tbls/overlaps_benegas_exons_chr", chr, ".csv.gz"))

# plotting distributions from both perspectives
p1 <- ggplot(plot_df_overlaps, aes(x=decile, y=N, fill=benegas_class)) +
  geom_col() + theme_classic() +
  scale_x_discrete(breaks=1:10) +
  scale_y_continuous(breaks=pretty_breaks()) +
  scale_fill_viridis_d(option="C", direction=1, guide=guide_legend(nrow=4), name=NULL) +
  labs(title=NULL, x=NULL, y=NULL, fill="Benegas") +
  guides(fill=guide_legend(nrow=1)) +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=16),
        axis.text.y=element_text(hjust=1),
        legend.position="bottom",
        legend.box="horizontal")

p2 <- ggplot(plot_df_overlaps, aes(x=benegas_class, y=N, fill=decile)) +
  geom_col() + theme_classic() +
  scale_x_discrete(breaks=1:12) +
  scale_y_continuous(breaks=pretty_breaks()) +
  scale_fill_viridis_d(option="C", direction=1, guide=guide_legend(nrow=4), name=NULL) +
  labs(title=NULL, x=NULL, y=NULL, fill="Benegas") +
  guides(fill=guide_legend(nrow=1)) +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=16),
        axis.text.y=element_text(hjust=1),
        legend.position="bottom",
        legend.box="horizontal")

p <- plot_grid(p1, p2, nrow=1, labels="AUTO")
save_plot(paste0("plots/ov_exons_benegas_chr", chr, ".pdf"), p, base_height=6, base_width=14)

cat("Finished.\n")
