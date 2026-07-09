#########################
#
# Setting up
#
########################

library(data.table)
library(tidyverse)
library(scales)

args <- commandArgs(trailingOnly=T)
chr <- args[1]

#########################
#
# 
#
########################

cat("\nReading Benegas score map...")
scores_chr <- fread(paste0("/../../media/gvbarroso/extradrive1/mut_rates_func/benegas/score_bins/benegas_bins_chr", chr, ".csv.gz"))
setnames(scores_chr, old="benegas_bin", new="benegas_class") # reserving "bin" to '1 kb bins'
scores_chr <- scores_chr[benegas_class < 15,] # assuming a classification w/ 1-12 -> constrained, 15 -> neutral

# NOTE: unlike Benegas scores (which are per site), phastCons are given in BED intervals
cat("done.\nReading phastCons scores map...")
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

# rolling BED intervals to single-nucleotide positions to join mutation rates
phast <- phast[, .(pos=seq(chromStart, chromEnd)), by=.(chrom, phast_class, chromStart, chromEnd)][, c("chromStart", "chromEnd") := NULL]

cat("done.\nJoining maps...")
overlaps <- scores_chr[phast, on=.(chrom, pos)] 

plot_df_overlaps <- overlaps[, .N, by=.(benegas_class, phast_class)]
plot_df_overlaps[, phast_class := factor(as.character(phast_class), levels = as.character(1:12))]
plot_df_overlaps[is.na(benegas_class), benegas_class := 15]
plot_df_overlaps[, benegas_class := factor(as.character(benegas_class), levels = c(as.character(1:12), "15"))]
plot_df_overlaps[, chr := chr]

p <- ggplot(plot_df_overlaps, aes(x=phast_class, y=N, fill=benegas_class)) +
  geom_col() + theme_classic() +
  scale_x_discrete(breaks=1:12) +
  scale_y_continuous(breaks=pretty_breaks()) +
  scale_fill_viridis_d(option="C", direction=1, guide=guide_legend(nrow=4), name=NULL) +
  labs(title="Overlaps between phastCons and Benegas scores", x=NULL, y=NULL, fill="Benegas") +
  guides(fill=guide_legend(nrow=1)) +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=16),
        axis.text.y=element_text(hjust=1),
        legend.position="bottom",
        legend.box="horizontal")
save_plot(paste0("benegas_phast_overlaps_chr", chr, ".pdf"), p, base_height=6, base_width=7)

