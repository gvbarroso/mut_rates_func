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

dat_withCpG <- dplyr::select(dat, c(mean_class_triplet_roulette, mean_class_triplet_carlson, mean_class_triplet_gnomad,
                                    num_sites_roulette, num_sites_carlson, num_sites_gnomad,
                                    overlap_benegas_exon, benegas_class, triplet)) %>% setDT()

# computing means across chromosomes, weighted by num sites in each category
dat_withCpG[, mean_class_triplet_roulette_gw := sum(mean_class_triplet_roulette * num_sites_roulette, na.rm=T) / sum(num_sites_roulette, na.rm = TRUE), by=.(benegas_class, overlap_benegas_exon)]
dat_withCpG[, mean_class_triplet_carlson_gw := sum(mean_class_triplet_carlson * num_sites_carlson, na.rm=T) / sum(num_sites_carlson, na.rm = TRUE), by=.(benegas_class, overlap_benegas_exon)]
dat_withCpG[, mean_class_triplet_gnomad_gw := sum(mean_class_triplet_gnomad * num_sites_gnomad, na.rm=T) / sum(num_sites_gnomad, na.rm = TRUE), by=.(benegas_class, overlap_benegas_exon)]

dat_withCpG <- unique(dat_withCpG, by=c("benegas_class", "overlap_benegas_exon")) %>% 
  dplyr::select(., c(benegas_class, overlap_benegas_exon, mean_class_triplet_roulette_gw, mean_class_triplet_carlson_gw, mean_class_triplet_gnomad_gw)) %>% setDT()

denom_roulette <- dat_withCpG[benegas_class==15, mean_class_triplet_roulette_gw]
denom_carlson <- dat_withCpG[benegas_class==15, mean_class_triplet_carlson_gw]
denom_gnomad <- dat_withCpG[benegas_class==15, mean_class_triplet_gnomad_gw]

dat_withCpG[, ratio_roulette := mean_class_triplet_roulette_gw / denom_roulette, by=benegas_class]
dat_withCpG[, ratio_carlson := mean_class_triplet_carlson_gw / denom_carlson, by=benegas_class]
dat_withCpG[, ratio_gnomad := mean_class_triplet_gnomad_gw / denom_carlson, by=benegas_class]
dat_withCpG[, CpG := T] 
dat_withCpG <- dat_withCpG[benegas_class < 15, .(benegas_class, ratio_roulette, ratio_carlson, ratio_gnomad, CpG)]

# filtering out CpG sites
dat_nonCpG <- dplyr::select(dat, c(mean_class_triplet_roulette, mean_class_triplet_carlson, mean_class_triplet_gnomad,
                                   num_sites_roulette, num_sites_carlson, num_sites_gnomad, benegas_class, triplet)) %>% 
  filter(., !triplet %in% CpGs) %>% setDT()

# computing means across chromosomes, weighted by num sites in each window
dat_nonCpG[, mean_class_triplet_roulette_gw := sum(mean_class_triplet_roulette * num_sites_roulette, na.rm=T) / sum(num_sites_roulette, na.rm = TRUE), by=.(benegas_class)]
dat_nonCpG[, mean_class_triplet_carlson_gw := sum(mean_class_triplet_carlson * num_sites_carlson, na.rm=T) / sum(num_sites_carlson, na.rm = TRUE), by=.(benegas_class)]
dat_nonCpG[, mean_class_triplet_gnomad_gw := sum(mean_class_triplet_gnomad * num_sites_gnomad, na.rm=T) / sum(num_sites_gnomad, na.rm = TRUE), by=.(benegas_class)]

dat_nonCpG <- unique(dat_nonCpG, by=c("benegas_class")) %>% 
  dplyr::select(., c(benegas_class, mean_class_triplet_roulette_gw, mean_class_triplet_carlson_gw, mean_class_triplet_gnomad_gw)) %>% setDT()

denom_roulette <- dat_nonCpG[benegas_class==15, mean_class_triplet_roulette_gw]
denom_carlson <- dat_nonCpG[benegas_class==15, mean_class_triplet_carlson_gw]
denom_gnomad <- dat_nonCpG[benegas_class==15, mean_class_triplet_gnomad_gw]

dat_nonCpG[, ratio_roulette := mean_class_triplet_roulette_gw / denom_roulette, by=benegas_class]
dat_nonCpG[, ratio_carlson := mean_class_triplet_carlson_gw / denom_carlson, by=benegas_class]
dat_nonCpG[, ratio_gnomad := mean_class_triplet_gnomad_gw / denom_carlson, by=benegas_class]
dat_nonCpG[, CpG := F] 
dat_nonCpG <- dat_nonCpG[benegas_class < 15, .(benegas_class, ratio_roulette, ratio_carlson, ratio_gnomad, CpG)]

m_ratios <- pivot_longer(rbind.data.frame(dat_withCpG, dat_nonCpG), cols=starts_with("ratio"), values_to="ratio", names_to="map") %>% setDT()
m_ratios[, annotation := "Benegas"]
fwrite(m_ratios, "gw_ratios_benegas.csv")

# "Ratios of mutation rates within Benegas elements w.r.t. genome-wide background"
p1 <- ggplot(m_ratios, aes(x=benegas_class, y=ratio, color=map, group=paste0(map, CpG))) +
  geom_line(aes(linetype=CpG), linewidth=1) + geom_point(size=3) + 
  geom_hline(yintercept=1, linetype="dashed", color="grey") +
  scale_x_continuous(breaks=1:nclasses) + theme_classic() + 
  scale_color_manual(values=c("brown1", "cyan3", "seagreen"), name=NULL,
                     labels=c("ratio_roulette"="Roulette", "ratio_gnomad"="gnomAD", "ratio_carlson"="Carlson")) +
  labs(x="Constraint class", y=expression(paste(mu, " ratio")), title=NULL) +
  scale_linetype_manual(name=NULL, values=c("FALSE"="solid", "TRUE"="dashed"),
                        labels=c("TRUE"="With CpG", "FALSE"="Without CpG")) +
  guides(linetype=guide_legend(keywidth=unit(1.5, "cm"), keyheight=unit(0.2, "cm"),
                               override.aes=list(color="black", linewidth=1.2, x=0, xend=1, y=0.5, yend=0.5))) +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=14),
        strip.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.position="bottom")
save_plot("plots/benegas_ratios.pdf", p1, base_height=5, base_width=9)

# stratifying by trinucleotide context
dat[, mean_class_triplet_gw_roulette := sum(mean_class_triplet_roulette * num_sites_roulette, na.rm=T) / sum(num_sites_roulette, na.rm = TRUE), by=.(benegas_class, triplet)]
dat[, mean_class_triplet_gw_carlson := sum(mean_class_triplet_carlson * num_sites_carlson, na.rm=T) / sum(num_sites_carlson, na.rm = TRUE),, by=.(benegas_class, triplet)]
dat[, mean_class_triplet_gw_gnomad := sum(mean_class_triplet_gnomad * num_sites_gnomad, na.rm=T) / sum(num_sites_gnomad, na.rm = TRUE),, by=.(benegas_class, triplet)]

denoms <- dat[benegas_class == 15, 
              .(den_roulette=mean_class_triplet_gw_roulette,
                den_carlson=mean_class_triplet_gw_carlson,
                den_gnomad=mean_class_triplet_gw_gnomad,
                triplet=triplet, chrom=chrom)]

nums <- dat[benegas_class %in% 1:12, .(num_roulette=mean_class_triplet_gw_roulette,
                                       num_carlson=mean_class_triplet_gw_carlson,
                                       num_gnomad=mean_class_triplet_gw_gnomad,
                                       benegas_class=benegas_class,
                                       triplet=triplet, chrom=chrom)]

ratios_benegas <- merge(nums, denoms, by=c("triplet", "chrom"), all.x = TRUE)
setorder(ratios_benegas, chrom, triplet, benegas_class)

ratios_benegas[, `:=`(ratio_roulette=num_roulette / den_roulette,
                      ratio_carlson=num_carlson / den_carlson,
                      ratio_gnomad=num_gnomad / den_gnomad)]

ratios_benegas[, c("num_roulette", "num_carlson", "num_gnomad", "den_roulette", "den_carlson", "den_gnomad") := NULL]
ratios_benegas_m <- pivot_longer(ratios_benegas, cols=starts_with("ratio_"), names_to="map", values_to="ratio") %>% setDT()

ratios_benegas_m[, map := factor(sub("^ratio_", "", map), levels = c("roulette", "carlson", "gnomad"))]
ratios_benegas_m[, annotation := "Benegas"]
fwrite(ratios_benegas_m, "ratios_triplets_benegas.csv")

# "Ratios of mutation rates within Benegas elements stratified by triplet context"
p2a <- ggplot(filter(ratios_benegas_m, map=="carlson"),
              aes(x=triplet, y=ratio, color=benegas_class)) + 
  annotate(xmin = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) - 0.5,
           xmax = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) + 0.5,
           geom="rect", ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_point(size=2.5) + theme_classic() + 
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_viridis_c(option="C", direction=1, name="Constraint Class", breaks=c(1, 12)) +
  labs(x=NULL, y=expression(paste(mu, " ratio (Carlson)")), title=NULL) +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_blank(),
        legend.text=element_text(size=16),
        legend.position="none",
        legend.title=element_text(size=18),
        legend.box="horizontal")

p2b <- ggplot(filter(ratios_benegas_m, map=="roulette"),
              aes(x=triplet, y=ratio, color=benegas_class)) + 
  annotate(xmin = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) - 0.5,
           xmax = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) + 0.5,
           geom="rect", ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_point(size=2.5) + theme_classic() + 
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_viridis_c(option="C", direction=1, name="Constraint Class", breaks=c(1, 12)) +
  labs(x=NULL, y=expression(paste(mu, " ratio (Roulette)")), title=NULL) +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12, angle=90, hjust=1, vjust=0.5),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title=element_text(size=18),
        legend.box="horizontal")
p2 <- plot_grid(p2a, p2b, ncol=1, rel_heights=c(1, 1.35))
save_plot("plots/ratios_benegas_triplet_classes.pdf", p2, base_height=8, base_width=14)
