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

# overlap with exons
ov_summary_files <- list.files("~/Devel/mut_rates_func/benegas/summary_tbls/", 
                               pattern=paste0("^overlap_summaries_benegas_exons_chr"), full.names=T)
dat <- data.table::rbindlist(lapply(ov_summary_files, fread))
nclasses <- length(unique(na.omit(dat$benegas_class))) - 1 # exclude last class (putatively neutral sites)

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
map_labels <- c(ratio_roulette="Roulette", ratio_gnomad="gnomAD", ratio_carlson="Carlson")

# for ordering points within TRUE/FALSE groups
m_ratios$xfac <- interaction(m_ratios$overlap_benegas_exon, m_ratios$benegas_class, sep="_", drop=T)
m_ratios$xfac <- factor(m_ratios$xfac, levels=unique(m_ratios$xfac)) 

# "Ratios of mutation rates within Benegas elements w.r.t. genome-wide background, stratified by overlap with exons"
p1 <- ggplot(m_ratios, aes(x=xfac, y=ratio, color=as.factor(benegas_class), shape=overlap_benegas_exon)) +
  facet_wrap(~map, labeller=labeller(map=map_labels)) +
  geom_point(size=3) + theme_classic() + 
  scale_shape_manual(values=c(16, 17), guide="none") +
  geom_hline(yintercept=1, linetype="dashed", color="grey") +
  scale_color_viridis_d(option="C", direction=1, name="Benegas class", breaks=1:12) +
  labs(x=NULL, y=expression(paste(mu, " ratio")), title=NULL) +
  scale_x_discrete(breaks=levels(m_ratios$xfac)[c(6.5, 18.5)], labels = c("FALSE", "TRUE")) +
  guides(color=guide_legend(nrow=1)) +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=14),
        axis.text.x=element_blank(),
        strip.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="none")

# overlap with phastCons
ov_summary_files <- list.files("~/Devel/mut_rates_func/benegas/summary_tbls/", 
                               pattern=paste0("^overlap_summaries_benegas_phast_chr"), full.names=T)
dat <- data.table::rbindlist(lapply(ov_summary_files, fread))
nclasses <- length(unique(na.omit(dat$benegas_class))) - 1 # exclude last class (putatively neutral sites)

dat_withCpG <- dplyr::select(dat, c(mean_class_roulette, mean_class_carlson, mean_class_gnomad,
                                    num_sites_roulette, num_sites_carlson, num_sites_gnomad,
                                    overlap_benegas_phast, benegas_class)) %>% setDT()

# computing means across chromosomes, weighted by num sites in each category
dat_withCpG[, mean_roulette_gw := sum(mean_class_roulette * num_sites_roulette, na.rm=T) / sum(num_sites_roulette, na.rm = TRUE), by=.(benegas_class, overlap_benegas_phast)]
dat_withCpG[, mean_carlson_gw := sum(mean_class_carlson * num_sites_carlson, na.rm=T) / sum(num_sites_carlson, na.rm = TRUE), by=.(benegas_class, overlap_benegas_phast)]
dat_withCpG[, mean_gnomad_gw := sum(mean_class_gnomad * num_sites_gnomad, na.rm=T) / sum(num_sites_gnomad, na.rm = TRUE), by=.(benegas_class, overlap_benegas_phast)]

dat_withCpG <- unique(dat_withCpG, by=c("benegas_class", "overlap_benegas_phast")) %>% 
  dplyr::select(., c(benegas_class, overlap_benegas_phast, mean_roulette_gw, mean_carlson_gw, mean_gnomad_gw)) %>% setDT()

denom_roulette <- dat_withCpG[benegas_class==15, mean_roulette_gw]
denom_carlson <- dat_withCpG[benegas_class==15, mean_carlson_gw]
denom_gnomad <- dat_withCpG[benegas_class==15, mean_gnomad_gw]

dat_withCpG[, ratio_roulette := mean_roulette_gw / denom_roulette, by=.(benegas_class, overlap_benegas_phast)]
dat_withCpG[, ratio_carlson := mean_carlson_gw / denom_carlson, by=.(benegas_class, overlap_benegas_phast)]
dat_withCpG[, ratio_gnomad := mean_gnomad_gw / denom_carlson, by=.(benegas_class, overlap_benegas_phast)]
dat_withCpG[, CpG := T] 
dat_withCpG <- dat_withCpG[benegas_class < 15, .(benegas_class, ratio_roulette, ratio_carlson, ratio_gnomad, overlap_benegas_phast, CpG)]

m_ratios <- pivot_longer(dat_withCpG, cols=starts_with("ratio"), values_to="ratio", names_to="map") %>% setDT()
map_labels <- c(ratio_roulette="Roulette", ratio_gnomad="gnomAD", ratio_carlson="Carlson")

# for ordering points within TRUE/FALSE groups
m_ratios$xfac <- interaction(m_ratios$overlap_benegas_phast, m_ratios$benegas_class, sep="_", drop=T)
m_ratios$xfac <- factor(m_ratios$xfac, levels=unique(m_ratios$xfac)) 

# "Ratios of mutation rates within Benegas elements w.r.t. genome-wide background, stratified by overlap with phastCons"
p2 <- ggplot(m_ratios, aes(x=xfac, y=ratio, color=as.factor(benegas_class), shape=overlap_benegas_phast)) +
  facet_wrap(~map, labeller=labeller(map=map_labels)) +
  geom_point(size=3) + theme_classic() + 
  scale_shape_manual(values=c(16, 17), guide="none") +
  geom_hline(yintercept=1, linetype="dashed", color="grey") +
  scale_color_viridis_d(option="C", direction=1, name="Benegas class", breaks=1:12) +
  labs(x="Overlap", y=expression(paste(mu, " ratio")), title=NULL) +
  scale_x_discrete(breaks=levels(m_ratios$xfac)[c(6.5, 18.5)], labels = c("FALSE", "TRUE")) +
  guides(color=guide_legend(nrow=1)) +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=14),
        strip.text=element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom",
        legend.direction="horizontal")

p <- plot_grid(p1, p2, ncol=1, align="v", labels="AUTO")
save_plot("plots/benegas_ov.pdf", p, base_height=8, base_width=12)
