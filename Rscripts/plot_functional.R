
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

setwd("~/Devel/mut_rates_func/functional/")

CpGs <- c("ACG", "CCG", "GCG", "TCG", "CGA", "CGC", "CGG", "CGT")

gw_summary_files <- list.files("~/Devel/mut_rates_func/functional/summary_tbls/", pattern=paste0("^summaries_chr"), full.names=T)
dat <- data.table::rbindlist(lapply(gw_summary_files, fread))

nclasses <- length(unique(na.omit(dat$functional_class))) - 1 # last class -> putatively neutral sites

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

stacked_files <- list.files("~/Devel/mut_rates_func/functional/summary_tbls/", pattern=paste0("^stacks_functional_chr"), full.names=T)
plot_df <- data.table::rbindlist(lapply(stacked_files, fread))
plot_df[, prop := N / sum(N), by=functional_class]
plot_df[, functional_class := factor(functional_class, levels=sort(unique(functional_class)))]

p0 <- ggplot(plot_df, aes(x=functional_class, y=prop, fill=triplet)) +
  theme_classic() + geom_col() +
  scale_fill_viridis_d(option="C", direction=1, guide=guide_legend(nrow=4), name=NULL) +
  scale_x_discrete(breaks=sort(unique(plot_df$functional_class)), 
                   labels=c(as.character(sort(unique(plot_df$functional_class))[-length(sort(unique(plot_df$functional_class)))]), 
                            "Outside")) +
  labs(x="functional_classent", y="Proportion", fill="Trinucleotide") +
  theme(panel.grid=functional_classent_blank(),
        axis.text=element_text(size=14),
        axis.title=element_text(size=18),
        legend.position="bottom",
        legend.box="horizontal")
save_plot(paste0("plots/functional_triplets.pdf"), p0, base_height=7, base_width=10)

####################
#
# Comparing rates within and without benegas functional_classents
# Here looking at benegas sites vs genome-wide background
#
####################





# start test
# TODO compute genome-wide means weighted by num sites in each chr, e.g
dat[, mean_class_triplet_roulette_gw := sum(mean_class_triplet_roulette * num_sites_roulette, na.rm=T) / sum(num_sites_roulette, na.rm = TRUE), by=.(elem)]
# TODO analogously for gnomad and carlson

# functional_class is the new term for elem
x <- filter(dat, elem=="decile_1") # dat in RAM still has old name
nrow(x) # 64 * 22
length(unique(x$mean_class_triplet_roulette)) # differs every row
# end test




dat_withCpG <- dplyr::select(dat, c(mean_class_triplet_roulette_gw, mean_class_triplet_carlson, mean_class_triplet_gnomad, functional_class, triplet)) %>%
  unique(., by=c("functional_class", "triplet")) %>% setDT()

dat_withCpG[, `:=`(mean_roluette_group=mean(mean_class_triplet_roulette),
                   mean_carlson_group=mean(mean_class_triplet_carlson),
                   mean_gnomad_group=mean(mean_class_triplet_gnomad)), by=functional_class] 
dat_withCpG <- unique(dat_withCpG, by="functional_class") %>%
  dplyr::select(., c(functional_class, mean_roluette_group, mean_carlson_group, mean_gnomad_group)) %>% setDT()

denom_roulette <- dat_withCpG[functional_class=="neutral", mean_roluette_group]
denom_carlson <- dat_withCpG[functional_class=="neutral", mean_carlson_group]

dat_withCpG[, ratio_roulette := mean_roluette_group / denom_roulette, by=functional_class]
dat_withCpG[, ratio_carlson := mean_carlson_group / denom_carlson, by=functional_class]
dat_withCpG[, CpG := T] 
dat_withCpG <- dat_withCpG[functional_class != "neutral", .(functional_class, ratio_roulette, ratio_carlson, CpG)]

# filtering out CpG sites
dat_nonCpG <- filter(dat, !triplet %in% CpGs) %>% 
  dplyr::select(., c(mean_bin_triplet_roulette, mean_bin_triplet_carlson, mean_bin_triplet_gnomad, functional_class, triplet)) %>%
  unique(., by=c("functional_class", "triplet")) %>% setDT()

dat_nonCpG[, `:=`(mean_roluette_group=mean(mean_bin_triplet_roulette),
                  mean_carlson_group=mean(mean_bin_triplet_carlson),
                  mean_gnomad_group=mean(mean_bin_triplet_gnomad)), by=functional_class] 
dat_nonCpG <- unique(dat_nonCpG, by="functional_class") %>%
  dplyr::select(., c(functional_class, mean_roluette_group, mean_carlson_group, mean_gnomad_group)) %>% setDT()

denom_roulette <- dat_nonCpG[functional_class==15, mean_roluette_group]
denom_carlson <- dat_nonCpG[functional_class==15, mean_carlson_group]

dat_nonCpG[, ratio_roulette := mean_roluette_group / denom_roulette, by=functional_class]
dat_nonCpG[, ratio_carlson := mean_carlson_group / denom_carlson, by=functional_class]
dat_nonCpG[, CpG := F] 
dat_nonCpG <- dat_nonCpG[functional_class < 15, .(functional_class, ratio_roulette, ratio_carlson, CpG)]

m_ratios <- pivot_longer(rbind.data.frame(dat_withCpG, dat_nonCpG), cols=starts_with("ratio"), values_to="ratio", names_to="map")

p1 <- ggplot(m_ratios, aes(x=functional_class, y=ratio, color=map, group=paste0(map, CpG))) +
  geom_line(aes(linetype=CpG), linewidth=1) + geom_point(size=3) + 
  geom_hline(yintercept=1, linetype="dashed", color="grey") +
  scale_x_continuous(breaks=1:nclasses) + theme_classic() + 
  scale_color_manual(values=c("cyan3", "brown1"), name=NULL,
                     labels=c("ratio_roulette"="Roulette", "ratio_carlson"="Carlson")) +
  labs(x="Constraint class", y=expression(paste(mu, " ratio")),
       title="Ratios of mutation rates within Benegas functional_classents w.r.t. genome-wide background") +
  scale_linetype_manual(name=NULL, values=c("FALSE"="solid", "TRUE"="dashed"),
                        labels=c("TRUE"="With CpG", "FALSE"="Without CpG")) +
  guides(linetype=guide_legend(keywidth=unit(1.5, "cm"), keyheight=unit(0.2, "cm"),
                               override.aes=list(color="black", linewidth=1.2, x=0, xend=1, y=0.5, yend=0.5))) +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=14),
        strip.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.position="bottom")
save_plot("plots/benegas_ratios.pdf", p1, base_height=4, base_width=8)

dat[, mean_bin_triplet_gw_roulette := mean(mean_bin_triplet_roulette), by=.(functional_class, triplet)]
dat[, mean_bin_triplet_gw_carlson := mean(mean_bin_triplet_carlson), by=.(functional_class, triplet)]
dat[, mean_bin_triplet_gw_gnomad := mean(mean_bin_triplet_gnomad), by=.(functional_class, triplet)]

denoms <- dat[functional_class == 15, 
              .(den_roulette=mean_bin_triplet_gw_roulette,
                den_carlson=mean_bin_triplet_gw_carlson,
                den_gnomad=mean_bin_triplet_gw_gnomad,
                triplet=triplet, chrom=chrom)]

nums <- dat[functional_class %in% 1:12, .(num_roulette=mean_bin_triplet_gw_roulette,
                                     num_carlson=mean_bin_triplet_gw_carlson,
                                     num_gnomad=mean_bin_triplet_gw_gnomad,
                                     functional_class=functional_class,
                                     triplet=triplet, chrom=chrom)]

ratios_benegas <- merge(nums, denoms, by=c("triplet", "chrom"), all.x = TRUE)
setorder(ratios_benegas, chrom, triplet, functional_class)

ratios_benegas[, `:=`(ratio_roulette=num_roulette / den_roulette,
                      ratio_carlson=num_carlson / den_carlson,
                      ratio_gnomad=num_gnomad / den_gnomad)]

ratios_benegas[, c("num_roulette", "num_carlson", "num_gnomad", "den_roulette", "den_carlson", "den_gnomad") := NULL]
ratios_benegas_m <- pivot_longer(ratios_benegas, cols=starts_with("ratio_"), names_to="map", values_to="ratio") %>% setDT()

ratios_benegas_m[, map := factor(sub("^ratio_", "", map), levels = c("roulette", "carlson", "gnomad"))]

p2a <- ggplot(filter(ratios_benegas_m, map=="carlson"),
              aes(x=triplet, y=ratio, color=functional_class)) + 
  annotate(xmin = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) - 0.5,
           xmax = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) + 0.5,
           geom="rect", ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_point(size=2.5) + theme_classic() + 
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_viridis_c(option="C", direction=1, name="Bin", breaks=c(1, 12)) +
  labs(x=NULL, y=expression(paste(mu, " ratio")),
       title="Ratios of mutation rates within Benegas functional_classents stratified by triplet context") +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=functional_classent_blank(),
        legend.text=element_text(size=16),
        legend.position="none",
        legend.title=element_text(size=18),
        legend.box="horizontal")

p2b <- ggplot(filter(ratios_benegas_m, map=="roulette"),
              aes(x=triplet, y=ratio, color=functional_class)) + 
  annotate(xmin = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) - 0.5,
           xmax = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) + 0.5,
           geom="rect", ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_point(size=2.5) + theme_classic() + 
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_viridis_c(option="C", direction=1, name="Bin", breaks=c(1, 12)) +
  labs(x=NULL, y=expression(paste(mu, " ratio")), title=NULL) +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12, angle=90, hjust=1, vjust=0.5),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title=element_text(size=18),
        legend.box="horizontal")
p2 <- plot_grid(p2a, p2b, ncol=1, rel_heights=c(1, 1.35), labels="AUTO")
save_plot("plots/ratios_benegas_triplet_bins.pdf", p2, base_height=8, base_width=14)

####################
#
# Comparing rates within and without benegas functional_classents
# Here looking at benegas sites vs other sites in 1 kb windows
#
####################

summary_files_1kb <- list.files("~/Devel/mut_rates_func/summary_tbls_benegas/", pattern=paste0("^exclusive_1kb_chr"), full.names=T)
withCpG_files <- summary_files_1kb[!grepl("nonCpG", summary_files_1kb)]
nonCpG_files <- summary_files_1kb[grepl("nonCpG", summary_files_1kb)]

# NOTE: chr 7 ends up with four extra (trivial) columns; TODO: fix within summarize_mut_benegas.R 
withCpG <- data.table::rbindlist(lapply(withCpG_files, fread), fill=T, use.names=T)
nonCpG <- data.table::rbindlist(lapply(nonCpG_files, fread), fill=T, use.names=T)

# manually removing them
withCpG[, c("mean_roulette_NA", "mean_carlson_NA", "mean_gnomad_NA", "n_sites_bin_NA") := NULL]
nonCpG[, c("mean_roulette_NA", "mean_carlson_NA", "mean_gnomad_NA", "n_sites_bin_NA") := NULL]

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

map_labels <- c("carlson"="Carlson", "gnomad"="gnomAD", "roulette"="Roulette")

# TODO adapt
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

p3 <- ggplot(tbl, aes(x=benegas_group, y=med_ratio)) +
  facet_wrap(~map, labeller=labeller(map=map_labels)) +
  theme_classic() +
  geom_segment(data=seg_df,
               aes(x=x, xend=xend, y=y, yend=yend, color=mean_B_mid, linetype=hasCpG, group=interaction(map, hasCpG)),
               linewidth=0.9, lineend="round", inherit.aes=FALSE) +
  geom_point(aes(color=mean_B, group=hasCpG), size=3) +
  geom_errorbar(aes(ymin=med_ratio - se_ratio, ymax=med_ratio + se_ratio, color=mean_B), width=0.2, linewidth=0.7) +
  labs(x="Constraint class", y=expression(paste(mu, " ratio")),
       title="Ratios of mutation rates within (isolated) functional_classents w.r.t. 1 kb background") +
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
