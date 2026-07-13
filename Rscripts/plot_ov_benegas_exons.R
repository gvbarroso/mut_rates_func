####################
#
# Setting up
#
####################

pdf(NULL)

library(data.table)
library(tidyverse)
library(scales)
library(cowplot)

setwd("~/Devel/mut_rates_func/benegas/")

CpGs <- c("ACG", "CCG", "GCG", "TCG", "CGA", "CGC", "CGG", "CGT")

ov_summary_files <- list.files("~/Devel/mut_rates_func/benegas/summary_tbls/", pattern=paste0("^overlap_summaries_benegas_exons_chr"), full.names=T)
dat <- data.table::rbindlist(lapply(ov_summary_files, fread))

nclasses <- length(unique(na.omit(dat$benegas_class))) - 1 # exclude last class (putatively neutral sites)

####################
#
# Comparing rates within and without benegas elements
# Here looking at benegas sites vs genome-wide background
#
####################

dat_withCpG <- dplyr::select(dat, c(mean_class_roulette, mean_class_carlson, mean_class_gnomad,
                                    num_sites_roulette, num_sites_carlson, num_sites_gnomad,
                                    overlap_benegas_exon, benegas_class)) %>% setDT()

# computing means across chromosomes, weighted by num sites in each category
dat_withCpG[, mean_roulette_gw := sum(mean_class_roulette * num_sites_roulette, na.rm=T) / sum(num_sites_roulette, na.rm = TRUE), by=.(benegas_class, overlap_benegas_exon)]
dat_withCpG[, mean_carlson_gw := sum(mean_class_carlson * num_sites_carlson, na.rm=T) / sum(num_sites_carlson, na.rm = TRUE), by=.(benegas_class, overlap_benegas_exon)]
dat_withCpG[, mean_gnomad_gw := sum(mean_class_gnomad * num_sites_gnomad, na.rm=T) / sum(num_sites_gnomad, na.rm = TRUE), by=.(benegas_class, overlap_benegas_exon)]

dat_withCpG <- unique(dat_withCpG, by=c("benegas_class", "overlap_benegas_exon")) %>% 
  dplyr::select(., c(benegas_class, overlap_benegas_exon, mean_roulette_gw, mean_carlson_gw, mean_gnomad_gw)) %>% setDT()

denom_roulette <- dat_withCpG[benegas_class==15, mean_roulette_gw]
denom_carlson <- dat_withCpG[benegas_class==15, mean_carlson_gw]
denom_gnomad <- dat_withCpG[benegas_class==15, mean_gnomad_gw]

dat_withCpG[, ratio_roulette := mean_roulette_gw / denom_roulette, by=.(benegas_class, overlap_benegas_exon)]
dat_withCpG[, ratio_carlson := mean_carlson_gw / denom_carlson, by=.(benegas_class, overlap_benegas_exon)]
dat_withCpG[, ratio_gnomad := mean_gnomad_gw / denom_carlson, by=.(benegas_class, overlap_benegas_exon)]
dat_withCpG[, CpG := T] 
dat_withCpG <- dat_withCpG[benegas_class < 15, .(benegas_class, ratio_roulette, ratio_carlson, ratio_gnomad, overlap_benegas_exon, CpG)]

m_ratios <- pivot_longer(dat_withCpG, cols=starts_with("ratio"), values_to="ratio", names_to="map") %>% setDT()
m_ratios[, annotation := "Benegas"]
fwrite(m_ratios, "gw_ratios_benegas.csv")

map_labels <- c(ratio_roulette="Roulette", ratio_gnomad="gnomAD", ratio_carlson="Carlson")

# "Ratios of mutation rates within Benegas elements w.r.t. genome-wide background, stratified by overlap with exons"
p1 <- ggplot(m_ratios, aes(x=overlap_benegas_exon, y=ratio, color=as.factor(benegas_class))) +
  facet_wrap(~map, labeller=labeller(map=map_labels)) + geom_jitter(size=3) + theme_classic() + 
  geom_hline(yintercept=1, linetype="dashed", color="grey") +
  scale_color_viridis_d(option="C", direction=1, name="Benegas class", breaks=c(1, 6, 12)) +
  labs(x="Overlap with exons", y=expression(paste(mu, " ratio")), title=NULL) +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=14),
        strip.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
save_plot("plots/benegas_ratios_ov_exons.pdf", p1, base_height=5, base_width=9)
