
library(data.table)
library(tidyverse)
library(scales)

args <- commandArgs(trailingOnly=TRUE)
chr <- args[1]

file_names <- list.files("~/Devel/mut_rates_func/summary_tbls/", pattern=paste0("^summaries_chr"), full.names=T)
dat <- data.table::rbindlist(lapply(file_names, fread))

# computing genome-wide means (averaging across autosomes), keeping bin x triplet
dat[, mean_gw_roulette := mean(mean_bin_roulette), by=.(benegas_bin, triplet)]
dat[, mean_gw_carlson := mean(mean_bin_carlson), by=.(benegas_bin, triplet)]
dat[, mean_gw_gnomad := mean(mean_bin_gnomad), by=.(benegas_bin, triplet)]

# se across autosomes, keeping bin x triplet
dat[, se_gw_roulette := sd(mean_gw_roulette, na.rm=T) / sqrt(.N), by=.(benegas_bin, triplet)]
dat[, se_gw_carlson := sd(mean_gw_carlson, na.rm=T) / sqrt(.N), by=.(benegas_bin, triplet)]
dat[, se_gw_gnomad := sd(mean_gw_gnomad, na.rm=T) / sqrt(.N), by=.(benegas_bin, triplet)]

dt_m <- dat %>%
  pivot_longer(cols=c(mean_gw_roulette, mean_gw_carlson, mean_gw_gnomad, se_gw_roulette, se_gw_carlson, se_gw_gnomad),
               names_to=c(".value", "source"), names_pattern="(mean_gw|se_gw)_(.*)")

p1 <- ggplot(dt_m, aes(x=benegas_bin, y=mean_gw, color=source, group=paste0(source, benegas_bin < 15))) +
  geom_line(linewidth=1) + geom_point(size=3) + 
  geom_errorbar(aes(ymin=mean_gw - se_gw, ymax=mean_gw + se_gw), width=0.2, linewidth=0.6) +
  scale_x_continuous(breaks=c(1:nbins, nbins + 3), labels=c(as.character(1:nbins), "Outside")) + theme_classic() + 
  scale_color_manual(values=c("cyan3", "brown1", "green4"), name=NULL) +
  labs(x="5-percentile bin (Benegas score)", y="Mean Rate", color="Map",
       title="Mutation rates per group (Mean +- SE across autosomes)") +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=14),
        strip.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.position="bottom")
save_plot("plots/benegas_ratios.pdf", p1, base_height=5, base_width=10)
save_plot("plots/benegas_ratios.png", p1, base_height=5, base_width=10)

dat[, mean_bin_triplet_gw_roulette := mean(mean_bin_triplet_roulette), by=.(benegas_bin, triplet)]
dat[, mean_bin_triplet_gw_carlson := mean(mean_bin_triplet_carlson), by=.(benegas_bin, triplet)]
dat[, mean_bin_triplet_gw_gnomad := mean(mean_bin_triplet_gnomad), by=.(benegas_bin, triplet)]

denoms <- dat[benegas_bin == 15, 
             .(den_roulette=mean_bin_triplet_gw_roulette,
               den_carlson=mean_bin_triplet_gw_carlson,
               den_gnomad=mean_bin_triplet_gw_gnomad,
               triplet=triplet, chrom=chrom)]

nums <- dat[benegas_bin %in% 1:12, .(num_roulette=mean_bin_triplet_gw_roulette,
                                     num_carlson=mean_bin_triplet_gw_carlson,
                                     num_gnomad=mean_bin_triplet_gw_gnomad,
                                     benegas_bin=benegas_bin,
                                     triplet=triplet, chrom=chrom)]

ratios_benegas <- merge(nums, denoms, by=c("triplet", "chrom"), all.x = TRUE)
setorder(ratios_benegas, chrom, triplet, benegas_bin)

ratios_benegas[, `:=`(ratio_roulette=num_roulette / den_roulette,
                      ratio_carlson=num_carlson / den_carlson,
                      ratio_gnomad=num_gnomad / den_gnomad)]

ratios_benegas[, c("num_roulette", "num_carlson", "num_gnomad", "den_roulette", "den_carlson", "den_gnomad") := NULL]
ratios_benegas_m <- pivot_longer(ratios_benegas, cols=starts_with("ratio_"), names_to="map", values_to="ratio") %>% setDT()
CpGs <- c("ACG", "CCG", "GCG", "TCG", "CGA", "CGC", "CGG", "CGT")

ratios_benegas_m[, map := factor(sub("^ratio_", "", map), levels = c("roulette", "carlson", "gnomad"))]

p4 <- ggplot(filter(ratios_benegas_m, map!="gnomad"), aes(x=triplet, y=ratio, color=map)) + 
  annotate(xmin = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) - 0.5,
           xmax = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) + 0.5,
           geom="rect", ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_point(size=2.5) + theme_classic() + 
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_manual(values=c("brown1", "cyan3"), name=NULL) +
  labs(x=NULL, y="Ratio of mutation rates", color="Map", 
       title="Ratios of mutation rates across Benegas bins and triplets") +
  theme(axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12, angle=90, hjust=1),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.box="horizontal")
save_plot("plots/ratios_benegas_triplet.pdf", p4, base_height=5, base_width=14)

p5a <- ggplot(filter(ratios_benegas_m, map=="carlson"),
              aes(x=triplet, y=ratio, color=benegas_bin)) + 
  annotate(xmin = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) - 0.5,
           xmax = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) + 0.5,
           geom="rect", ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_point(size=2.5) + theme_classic() + 
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_viridis_c(option="C", direction=1, name="Bin", breaks=c(1, 12)) +
  labs(x=NULL, y=expression(paste(mu, " ratio")), color="Map") +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_blank(),
        legend.text=element_text(size=16),
        legend.position="none",
        legend.title=element_text(size=18),
        legend.box="horizontal")

p5b <- ggplot(filter(ratios_benegas_m, map=="roulette"),
              aes(x=triplet, y=ratio, color=benegas_bin)) + 
  annotate(xmin = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) - 0.5,
           xmax = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) + 0.5,
           geom="rect", ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_point(size=2.5) + theme_classic() + 
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_viridis_c(option="C", direction=1, name="Bin", breaks=c(1, 12)) +
  labs(x=NULL, y=expression(paste(mu, " ratio")), color="Map", title=NULL) +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12, angle=90, hjust=1),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title=element_text(size=18),
        legend.box="horizontal")
p5 <- plot_grid(p5a, p5b, ncol=1, rel_heights=c(1, 1.35), labels="AUTO")
save_plot("plots/ratios_benegas_triplet_bins.pdf", p5, base_height=8, base_width=14)

# TODO join B-values and Rec map
# define 1 kb windows
# split coding vs noncoding sites
# etc

# single_benegas has 1 kb windows with elements exclusively from one benegas group
#> names(single_benegas)
#[1] "n_sites_bin_1"   "n_sites_bin_2"   "n_sites_bin_3"   "n_sites_bin_4"   "n_sites_bin_5"   "n_sites_bin_6"  
#[7] "n_sites_bin_7"   "n_sites_bin_8"   "n_sites_bin_9"   "n_sites_bin_10"  "n_sites_bin_11"  "n_sites_bin_12" 
#[13] "n_sites_bin_15"  "chrom"           "pos"             "B"               "name"            "value"          
#[19] "benegas_group"   "variable"        "sum_constrained" "map"             "weight"          "mean_B" 

tbl_mean_ratios <- single_benegas[, .(avg_ratio = sum(value * weight, na.rm = TRUE) / sum(weight, na.rm = TRUE),
                                      se_ratio = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value)))), by = .(benegas_group, map)]

# joining mean B-value
b_group <- single_benegas[, .(mean_B = mean(B)), by=.(benegas_group)]
tbl <- tbl_mean_ratios[b_group, on=.(benegas_group)]

p2 <- ggplot(tbl, aes(x=benegas_group, y=avg_ratio, color=mean_B)) +
  facet_wrap(~map) + theme_classic() + 
  geom_point(size=3) + geom_line() +
  geom_errorbar(aes(ymin = avg_ratio - se_ratio, ymax = se_ratio + se_ratio), width = 0.2, linewidth = 0.7) +
  labs(x="phastCons bin", y=expression(paste(mu, " ratio"))) +
  scale_x_continuous(breaks=1:12) +
  scale_y_continuous(breaks=pretty_breaks(), limits=c(min(tbl$avg_ratio) - tbl$se_ratio[which.min(tbl$avg_ratio)], 1)) +
  scale_color_viridis_c(option="C", direction=1, name="B",
                        breaks=c(round(min(tbl$mean_B)+0.01, 2), round(max(tbl$mean_B)-0.01, 2))) +
  theme(strip.text=element_text(size = 16),
        axis.title=element_text(size = 20),
        axis.text=element_text(size = 16),
        legend.text=element_text(size = 16), 
        legend.title=element_text(size = 16, margin = margin(b = 20)), 
        legend.position="bottom")


####################
#
# mutation rates per trinucleotide
#
####################

dat[, `:=`(mean_roulette_triplets=mean(roulette, na.rm=T),
           mean_carlson_triplets=mean(carlson, na.rm=T),
           mean_gnomad_triplets=mean(gnomad, na.rm=T)), by=triplet]

tmp <- dat[!duplicated(triplet),] %>% 
  pivot_longer(., cols=ends_with("_triplets"), names_to="map", values_to="rate")

p5 <- ggplot(tmp, aes(x=triplet, y=rate, color=map)) + 
  geom_point(size=2.5) + theme_classic() + 
  scale_color_manual(values=c("cyan3", "green4", "brown1"), 
                     labels=c(mean_roulette_triplets="roulette", 
                              mean_carlson_triplets="carlson",
                              mean_gnomad_triplets="gnomad"), name=NULL) +
  labs(x=NULL, y=expression(mu), color="Map", 
       title="Mutation rates across triplets") +
  theme(axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12, angle=90, vjust=1, hjust=0.5),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.box="horizontal")
save_plot("~/Desktop/mut_rates/rates_triplets.pdf", p5, base_height=6, base_width=10)
