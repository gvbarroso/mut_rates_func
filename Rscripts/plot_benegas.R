
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

gw_summary_files <- list.files("~/Devel/mut_rates_func/benegas/summary_tbls/", pattern=paste0("^summaries_chr"), full.names=T)
dat <- data.table::rbindlist(lapply(gw_summary_files, fread))

nclasses <- length(unique(na.omit(dat$benegas_class))) - 1 # exclude last class (putatively neutral sites)

####################
#
# Mutation rates per trinucleotide (descriptive)
#
####################

dat[, `:=`(mean_roulette_triplets=mean(roulette, na.rm=T),
           mean_carlson_triplets=mean(carlson, na.rm=T),
           mean_gnomad_triplets=mean(gnomad, na.rm=T)), by=triplet]

p0 <- dat[!duplicated(triplet),] %>% pivot_longer(., cols=ends_with("_triplets"), names_to="map", values_to="rate") %>%
  ggplot(aes(x=triplet, y=rate, color=map)) + 
  geom_point(size=2.5) + theme_classic() + 
  scale_color_manual(values=c("cyan3", "green4", "brown1"), 
                     labels=c(mean_roulette_triplets="roulette", 
                              mean_carlson_triplets="carlson",
                              mean_gnomad_triplets="gnomad"), name=NULL) +
  labs(x=NULL, y=expression(mu), color="Map", title=NULL) +
  theme(axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12, angle=90, vjust=1, hjust=1),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.box="horizontal")
save_plot("plots/rates_triplets.pdf", p0, base_height=4, base_width=12)

stacked_files <- list.files("~/Devel/mut_rates_func/benegas/summary_tbls/", pattern=paste0("^stacks_benegas_chr"), full.names=T)
plot_df <- data.table::rbindlist(lapply(stacked_files, fread))
plot_df[, prop := N / sum(N), by=benegas_class]
plot_df[, benegas_class := factor(benegas_class, levels=sort(unique(benegas_class)))]

p0 <- ggplot(plot_df, aes(x=benegas_class, y=prop, fill=triplet)) +
  theme_classic() + geom_col() +
  scale_fill_viridis_d(option="C", direction=1, guide=guide_legend(nrow=4), name=NULL) +
  scale_x_discrete(breaks=sort(unique(plot_df$benegas_class)), 
                   labels=c(as.character(sort(unique(plot_df$benegas_class))[-length(sort(unique(plot_df$benegas_class)))]), 
                            "Outside")) +
  labs(x="Benegas class", y="Proportion", fill="Trinucleotide") +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=18),
        legend.position="bottom",
        legend.box="horizontal")
save_plot(paste0("plots/benegas_triplets.pdf"), p0, base_height=7, base_width=10)

####################
#
# Comparing rates within and without benegas elements
# Here looking at benegas sites vs genome-wide background
#
####################

dat_withCpG <- dplyr::select(dat, c(mean_class_triplet_roulette, mean_class_triplet_carlson, mean_class_triplet_gnomad,
                                    num_sites_roulette, num_sites_carlson, num_sites_gnomad, benegas_class, triplet)) %>% setDT()

# computing means across chromosomes, weighted by num sites in each window
dat_withCpG[, mean_class_triplet_roulette_gw := sum(mean_class_triplet_roulette * num_sites_roulette, na.rm=T) / sum(num_sites_roulette, na.rm = TRUE), by=.(benegas_class)]
dat_withCpG[, mean_class_triplet_carlson_gw := sum(mean_class_triplet_carlson * num_sites_carlson, na.rm=T) / sum(num_sites_carlson, na.rm = TRUE), by=.(benegas_class)]
dat_withCpG[, mean_class_triplet_gnomad_gw := sum(mean_class_triplet_gnomad * num_sites_gnomad, na.rm=T) / sum(num_sites_gnomad, na.rm = TRUE), by=.(benegas_class)]

dat_withCpG <- unique(dat_withCpG, by=c("benegas_class")) %>% 
  dplyr::select(., c(benegas_class, mean_class_triplet_roulette_gw, mean_class_triplet_carlson_gw, mean_class_triplet_gnomad_gw)) %>% setDT()

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

####################
#
# Comparing rates within and without benegas elements
# Here looking at benegas sites vs other sites in 1 kb windows
#
####################

summary_files_1kb <- list.files("~/Devel/mut_rates_func/benegas/summary_tbls/", pattern=paste0("^exclusive_1kb_chr"), full.names=T)
withCpG_files <- summary_files_1kb[!grepl("nonCpG", summary_files_1kb)]
nonCpG_files <- summary_files_1kb[grepl("nonCpG", summary_files_1kb)]

withCpG <- data.table::rbindlist(lapply(withCpG_files, fread), fill=T, use.names=T)
nonCpG <- data.table::rbindlist(lapply(nonCpG_files, fread), fill=T, use.names=T)

withCpG[, hasCpG := T]
nonCpG[, hasCpG := F]
dat <- rbind.data.frame(withCpG, nonCpG)

# NOTE use log(ratio)?
# using median as a summary to mitigate outliers
tbl_med_ratios <- dat[, .(med_ratio=median(value, na.rm=T), 
                          se_ratio=sd(value, na.rm=T) / sqrt(sum(!is.na(value)))),
                      by=.(benegas_group, map, hasCpG)]

b_group <- dat[, .(mean_B = mean(B)), by=.(benegas_group)] # joining mean B-value
tbl <- tbl_med_ratios[b_group, on=.(benegas_group)]
tbl[, annotation := "Benegas"]
fwrite(tbl, "ratios_1kb_benegas.csv")

map_labels <- c("carlson"="Carlson", "gnomad"="gnomAD", "roulette"="Roulette")

# separate table to plot segments because lines cannot be plotted with both color and linetype
seg_df <- tbl %>%
  arrange(map, hasCpG, benegas_group) %>%
  group_by(map, hasCpG) %>%
  mutate(x=benegas_group,
         y=med_ratio,
         xend=lead(benegas_group),
         yend=lead(med_ratio),
         mean_B_mid=(mean_B + lead(mean_B)) / 2) %>%
  filter(!is.na(xend)) %>%
  ungroup()

# "Ratios of mutation rates within (isolated) Benegas elements w.r.t. 1 kb background"
p3 <- ggplot(tbl, aes(x=benegas_group, y=med_ratio)) +
  facet_wrap(~map, labeller=labeller(map=map_labels)) +
  geom_hline(yintercept=1, linetype="dashed", color="grey") + theme_classic() +
  geom_segment(data=seg_df,
               aes(x=x, xend=xend, y=y, yend=yend, color=mean_B_mid, linetype=hasCpG, group=interaction(map, hasCpG)),
               linewidth=0.9, lineend="round", inherit.aes=FALSE) +
  geom_point(aes(color=mean_B, group=hasCpG), size=3) +
  geom_errorbar(aes(ymin=med_ratio - se_ratio, ymax=med_ratio + se_ratio, color=mean_B), width=0.2, linewidth=0.7) +
  labs(x="Constraint class", y=expression(paste(mu, " ratio")), title=NULL) +
  scale_x_continuous(breaks=1:12) +
  scale_linetype_manual(name=NULL, values=c("FALSE"="solid", "TRUE"="dashed"), labels=c("TRUE"="With CpG", "FALSE"="Without CpG")) +
  scale_y_continuous(breaks=pretty_breaks(), limits=c(min(tbl$med_ratio) - tbl$se_ratio[which.min(tbl$med_ratio)], 1)) +
  scale_color_viridis_c(option="C", direction=1, name="B-value", breaks=c(round(min(tbl$mean_B) + 0.01, 2), round(max(tbl$mean_B) - 0.01, 2))) +
  guides(linetype = guide_legend(keywidth = unit(1.5, "cm"), keyheight = unit(0.2, "cm"),
         override.aes = list(color = "black", linewidth = 1.2, x = 0, xend = 1, y = 0.5, yend = 0.5))) +
  theme(strip.text=element_text(size=16),
        axis.title=element_text(size=20),
        axis.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16, margin=margin(b=20)),
        legend.position="bottom")
save_plot("plots/benegas_ratios_isolated_1kb.pdf", p3, base_height=5, base_width=13)
