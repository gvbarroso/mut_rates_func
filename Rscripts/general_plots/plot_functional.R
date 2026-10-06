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
library(ggnewscale)

setwd("~/Devel/mut_rates_func/functional/")

CpGs <- c("ACG", "CCG", "GCG", "TCG", "CGA", "CGC", "CGG", "CGT")

gw_summary_files <- list.files("~/Devel/mut_rates_func/functional/summary_tbls/", pattern=paste0("^summaries_chr"), full.names=T)
dat <- data.table::rbindlist(lapply(gw_summary_files, fread))

nclasses <- length(unique(na.omit(dat$functional_class))) - 1 # exclude last class (putatively neutral sites)

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
  labs(x="Functional class", y="Proportion", fill="Trinucleotide") +
  theme(panel.grid=element_blank(),
        axis.text.x=element_text(size=14, angle=45, hjust=1),
        axis.text=element_text(size=14),
        axis.title=element_text(size=18),
        legend.position="bottom",
        legend.box="horizontal")
save_plot(paste0("plots/functional_triplets.pdf"), p0, base_height=7, base_width=10)

####################
#
# Comparing rates within and without functional elements
# Here looking at functional sites vs genome-wide background
#
####################

dat_withCpG <- dplyr::select(dat, c(mean_class_triplet_roulette, mean_class_triplet_carlson, mean_class_triplet_gnomad,
                                    se_class_roulette, se_class_gnomad, se_class_carlson,
                                    num_sites_roulette, num_sites_carlson, num_sites_gnomad, functional_class, triplet)) %>% setDT()

# computing means across chromosomes, weighted by num sites in each window
dat_withCpG[, mean_class_triplet_roulette_gw := sum(mean_class_triplet_roulette * num_sites_roulette, na.rm=T) / sum(num_sites_roulette, na.rm = TRUE), by=.(functional_class)]
dat_withCpG[, mean_class_triplet_carlson_gw := sum(mean_class_triplet_carlson * num_sites_carlson, na.rm=T) / sum(num_sites_carlson, na.rm = TRUE), by=.(functional_class)]
dat_withCpG[, mean_class_triplet_gnomad_gw := sum(mean_class_triplet_gnomad * num_sites_gnomad, na.rm=T) / sum(num_sites_gnomad, na.rm = TRUE), by=.(functional_class)]

# SEs across chromosomes
gw_se_roulette <- dat_withCpG[, { m <- dat_withCpG[functional_class == .BY$functional_class, mean_class_triplet_roulette_gw]
v <- sum(num_sites_roulette * (se_class_roulette^2 + (mean_class_triplet_roulette_gw - m)^2)) / sum(num_sites_roulette) 
.(se_gw_roulette = sqrt(v)) }, by = functional_class]

gw_se_carlson <- dat_withCpG[, { m <- dat_withCpG[functional_class == .BY$functional_class, mean_class_triplet_carlson_gw]
v <- sum(num_sites_carlson * (se_class_carlson^2 + (mean_class_triplet_carlson_gw - m)^2)) / sum(num_sites_carlson) 
.(se_gw_carlson = sqrt(v)) }, by = functional_class]

gw_se_gnomad <- dat_withCpG[, { m <- dat_withCpG[functional_class == .BY$functional_class, mean_class_triplet_gnomad_gw]
v <- sum(num_sites_gnomad * (se_class_gnomad^2 + (mean_class_triplet_gnomad_gw - m)^2)) / sum(num_sites_gnomad) 
.(se_gw_gnomad = sqrt(v)) }, by = functional_class]

dat_withCpG <- unique(dat_withCpG, by=c("functional_class")) %>% 
  dplyr::select(., c(functional_class, mean_class_triplet_roulette_gw, mean_class_triplet_carlson_gw, mean_class_triplet_gnomad_gw)) %>% setDT()

dat_withCpG <- merge(dat_withCpG, gw_se_roulette, by="functional_class")
dat_withCpG <- merge(dat_withCpG, gw_se_carlson, by="functional_class")
dat_withCpG <- merge(dat_withCpG, gw_se_gnomad, by="functional_class")

denom_roulette <- dat_withCpG[functional_class=="neutral", mean_class_triplet_roulette_gw]
denom_carlson <- dat_withCpG[functional_class=="neutral", mean_class_triplet_carlson_gw]
denom_gnomad <- dat_withCpG[functional_class=="neutral", mean_class_triplet_gnomad_gw]

dat_withCpG[, ratio_roulette := mean_class_triplet_roulette_gw / denom_roulette, by=functional_class]
dat_withCpG[, ratio_carlson := mean_class_triplet_carlson_gw / denom_carlson, by=functional_class]
dat_withCpG[, ratio_gnomad := mean_class_triplet_gnomad_gw / denom_carlson, by=functional_class]
dat_withCpG[, ratio_roulette_se := se_gw_roulette / denom_roulette]
dat_withCpG[, ratio_carlson_se := se_gw_carlson / denom_carlson]
dat_withCpG[, ratio_gnomad_se := se_gw_gnomad / denom_gnomad]
dat_withCpG[, CpG := T] 
dat_withCpG <- dat_withCpG[functional_class != "neutral", .(functional_class, ratio_roulette, ratio_carlson, ratio_gnomad,
                                                            ratio_roulette_se, ratio_carlson_se, ratio_gnomad_se, CpG)]

# filtering out CpG sites
dat_nonCpG <- dplyr::select(dat, c(mean_class_triplet_roulette, mean_class_triplet_carlson, mean_class_triplet_gnomad,
                                   se_class_roulette, se_class_gnomad, se_class_carlson,
                                   num_sites_roulette, num_sites_carlson, num_sites_gnomad, functional_class, triplet)) %>% 
  filter(., !triplet %in% CpGs) %>% setDT()

# computing means across chromosomes, weighted by num sites in each window
dat_nonCpG[, mean_class_triplet_roulette_gw := sum(mean_class_triplet_roulette * num_sites_roulette, na.rm=T) / sum(num_sites_roulette, na.rm = TRUE), by=.(functional_class)]
dat_nonCpG[, mean_class_triplet_carlson_gw := sum(mean_class_triplet_carlson * num_sites_carlson, na.rm=T) / sum(num_sites_carlson, na.rm = TRUE), by=.(functional_class)]
dat_nonCpG[, mean_class_triplet_gnomad_gw := sum(mean_class_triplet_gnomad * num_sites_gnomad, na.rm=T) / sum(num_sites_gnomad, na.rm = TRUE), by=.(functional_class)]

# SEs across chromosomes
gw_se_roulette <- dat_nonCpG[, { m <- dat_nonCpG[functional_class == .BY$functional_class, mean_class_triplet_roulette_gw]
v <- sum(num_sites_roulette * (se_class_roulette^2 + (mean_class_triplet_roulette_gw - m)^2)) / sum(num_sites_roulette) 
.(se_gw_roulette = sqrt(v)) }, by = functional_class]

gw_se_carlson <- dat_nonCpG[, { m <- dat_nonCpG[functional_class == .BY$functional_class, mean_class_triplet_carlson_gw]
v <- sum(num_sites_carlson * (se_class_carlson^2 + (mean_class_triplet_carlson_gw - m)^2)) / sum(num_sites_carlson) 
.(se_gw_carlson = sqrt(v)) }, by = functional_class]

gw_se_gnomad <- dat_nonCpG[, { m <- dat_nonCpG[functional_class == .BY$functional_class, mean_class_triplet_gnomad_gw]
v <- sum(num_sites_gnomad * (se_class_gnomad^2 + (mean_class_triplet_gnomad_gw - m)^2)) / sum(num_sites_gnomad) 
.(se_gw_gnomad = sqrt(v)) }, by = functional_class]

dat_nonCpG <- unique(dat_nonCpG, by=c("functional_class")) %>% 
  dplyr::select(., c(functional_class, mean_class_triplet_roulette_gw, mean_class_triplet_carlson_gw, mean_class_triplet_gnomad_gw)) %>% setDT()

dat_nonCpG <- merge(dat_nonCpG, gw_se_roulette, by="functional_class")
dat_nonCpG <- merge(dat_nonCpG, gw_se_carlson, by="functional_class")
dat_nonCpG <- merge(dat_nonCpG, gw_se_gnomad, by="functional_class")

denom_roulette <- dat_nonCpG[functional_class=="neutral", mean_class_triplet_roulette_gw]
denom_carlson <- dat_nonCpG[functional_class=="neutral", mean_class_triplet_carlson_gw]
denom_gnomad <- dat_nonCpG[functional_class=="neutral", mean_class_triplet_gnomad_gw]

dat_nonCpG[, ratio_roulette := mean_class_triplet_roulette_gw / denom_roulette, by=functional_class]
dat_nonCpG[, ratio_carlson := mean_class_triplet_carlson_gw / denom_carlson, by=functional_class]
dat_nonCpG[, ratio_gnomad := mean_class_triplet_gnomad_gw / denom_gnomad, by=functional_class]
dat_nonCpG[, ratio_roulette_se := se_gw_roulette / denom_roulette]
dat_nonCpG[, ratio_carlson_se := se_gw_carlson / denom_carlson]
dat_nonCpG[, ratio_gnomad_se := se_gw_gnomad / denom_gnomad]

dat_nonCpG[, CpG := F] 
dat_nonCpG <- dat_nonCpG[functional_class != "neutral", .(functional_class, ratio_roulette, ratio_carlson, ratio_gnomad,
                                                          ratio_roulette_se, ratio_carlson_se, ratio_gnomad_se, CpG)]

m_means <- pivot_longer(
  rbind.data.frame(dat_withCpG, dat_nonCpG),
  cols = c(ratio_roulette, ratio_carlson, ratio_gnomad),
  names_to = "map",
  values_to = "ratio"
) %>% setDT()

m_ses <- pivot_longer(
  rbind.data.frame(dat_withCpG, dat_nonCpG),
  cols = c(ratio_roulette_se, ratio_carlson_se, ratio_gnomad_se),
  names_to = "map",
  values_to = "se"
) %>% setDT()

m_ses[, map := sub("_se$", "", map)] 

m_ratios <- m_means[m_ses, on = c("functional_class", "CpG", "map")]
dat_nonCpG <- dat_nonCpG[functional_class != "neutral", .(functional_class, ratio_roulette, ratio_carlson, ratio_gnomad, CpG)]

m_ratios <- filter(m_ratios, functional_class != "decile_11") %>% setDT() # filtering out "extra" exons
m_ratios[, functional_class := factor(gsub("^decile_", "d_", functional_class), levels=c(paste0("d_", 1:10), "enhancer", "promoter"))]

m_ratios[, annotation := "Functional"]
m_ratios <- m_ratios[, .(functional_class, CpG, map, ratio, se, annotation)]
fwrite(m_ratios, "gw_ratios_functional.csv")

# "Ratios of mutation rates within functional elements w.r.t. genome-wide background"
p1 <- ggplot(m_ratios, aes(x=functional_class, y=ratio, color=map, group=paste0(map, CpG))) +
  geom_line(aes(linetype=CpG), linewidth=1) + geom_point(size=3) + 
  geom_hline(yintercept=1, linetype="dashed", color="grey") + theme_classic() + 
  scale_color_manual(values=c("brown1", "cyan3", "seagreen"), name=NULL,
                     labels=c("ratio_roulette"="Roulette", "ratio_gnomad"="gnomAD", "ratio_carlson"="Carlson")) +
  labs(x="Constraint class", y=expression(paste(mu, " ratio")), title=NULL) +
  scale_linetype_manual(name=NULL, values=c("FALSE"="solid", "TRUE"="dashed"),
                        labels=c("TRUE"="With CpG", "FALSE"="Without CpG")) +
  guides(linetype=guide_legend(keywidth=unit(1.5, "cm"), keyheight=unit(0.2, "cm"),
                               override.aes=list(color="black", linewidth=1.2, x=0, xend=1, y=0.5, yend=0.5))) +
  theme(axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14, angle=45, vjust=0.6),
        strip.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.position="bottom")
save_plot("plots/functional_ratios.pdf", p1, base_height=5, base_width=9)

# stratifying by trinucleotide context
dat[, mean_class_triplet_gw_roulette := sum(mean_class_triplet_roulette * num_sites_roulette, na.rm=T) / sum(num_sites_roulette, na.rm=TRUE), by=.(functional_class, triplet)]
dat[, mean_class_triplet_gw_carlson := sum(mean_class_triplet_carlson * num_sites_carlson, na.rm=T) / sum(num_sites_carlson, na.rm=TRUE),, by=.(functional_class, triplet)]
dat[, mean_class_triplet_gw_gnomad := sum(mean_class_triplet_gnomad * num_sites_gnomad, na.rm=T) / sum(num_sites_gnomad, na.rm=TRUE),, by=.(functional_class, triplet)]

denoms <- dat[functional_class == "neutral", 
              .(den_roulette=mean_class_triplet_gw_roulette,
                den_carlson=mean_class_triplet_gw_carlson,
                den_gnomad=mean_class_triplet_gw_gnomad,
                triplet=triplet, chrom=chrom)]

nums <- dat[functional_class != "neutral", 
            .(num_roulette=mean_class_triplet_gw_roulette,
              num_carlson=mean_class_triplet_gw_carlson,
              num_gnomad=mean_class_triplet_gw_gnomad,
              functional_class=functional_class,
              triplet=triplet, chrom=chrom)]

ratios_functional <- merge(nums, denoms, by=c("triplet", "chrom"), all.x=TRUE)
setorder(ratios_functional, chrom, triplet, functional_class)

ratios_functional[, `:=`(ratio_roulette=num_roulette / den_roulette,
                      ratio_carlson=num_carlson / den_carlson,
                      ratio_gnomad=num_gnomad / den_gnomad)]

ratios_functional[, c("num_roulette", "num_carlson", "num_gnomad", "den_roulette", "den_carlson", "den_gnomad") := NULL]
ratios_functional_m <- pivot_longer(ratios_functional, cols=starts_with("ratio_"), names_to="map", values_to="ratio") %>% setDT()

ratios_functional_m[, map := factor(sub("^ratio_", "", map), levels=c("roulette", "carlson", "gnomad"))]
ratios_functional_m[, class := sub(".*_", "", functional_class)] # for splitting plot colors between exons and regulatory
ratios_functional_m[, annotation := "Functional"]
fwrite(ratios_functional_m, "ratios_triplets_functional.csv")

# "Ratios of mutation rates within functional classes stratified by triplet context"
p2a <- ggplot(filter(ratios_functional_m, map=="carlson")) + 
  annotate(xmin=which(levels(factor(ratios_functional_m$triplet)) %in% CpGs) - 0.5,
           xmax=which(levels(factor(ratios_functional_m$triplet)) %in% CpGs) + 0.5,
           geom="rect", ymin=-Inf, ymax=Inf, fill="grey85", alpha=0.6) +
geom_point(data=subset(ratios_functional_m, map=="carlson" & grepl("^decile_", functional_class)),
           aes(x=triplet, y=ratio, color=as.numeric(class)), size=2.5) +
  scale_color_viridis_c(option="C", direction=1, name="Decile", breaks=c(1, 10)) +
  ggnewscale::new_scale_color() +
geom_point(data=subset(ratios_functional_m, map=="carlson" & functional_class %in% c("enhancer", "promoter")),
           aes(x=triplet, y=ratio, color=functional_class), size=3) +
  scale_color_manual(name=NULL, values=c(enhancer="seagreen", promoter="cyan3")) +
geom_hline(yintercept=1, linetype="dashed", color="black") +
  theme_classic() +
  labs(x=NULL, y=expression(paste(mu, " ratio (Carlson)")), title=NULL) +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_blank(),
        legend.position="none")

p2b <- ggplot(filter(ratios_functional_m, map=="roulette")) + 
  annotate(xmin=which(levels(factor(ratios_functional_m$triplet)) %in% CpGs) - 0.5,
           xmax=which(levels(factor(ratios_functional_m$triplet)) %in% CpGs) + 0.5,
           geom="rect", ymin=-Inf, ymax=Inf, fill="grey85", alpha=0.6) +
  geom_point(data=subset(ratios_functional_m, map=="roulette" & grepl("^decile_", functional_class)),
             aes(x=triplet, y=ratio, color=as.numeric(class)), size=2.5) +
  scale_color_viridis_c(option="C", direction=1, name="Decile", breaks=c(1, 10)) +
  ggnewscale::new_scale_color() +
  geom_point(data=subset(ratios_functional_m, map=="roulette" & functional_class %in% c("enhancer", "promoter")),
             aes(x=triplet, y=ratio, color=functional_class), size=3) +
  scale_color_manual(name=NULL, values=c(enhancer="seagreen", promoter="cyan3")) +
  geom_hline(yintercept=1, linetype="dashed", color="black") +
  theme_classic() +
  labs(x=NULL, y=expression(paste(mu, " ratio (Roulette)")), title=NULL) +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12, angle=90, hjust=1, vjust=0.5),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title=element_text(size=18),
        legend.box="horizontal")
p2 <- plot_grid(p2a, p2b, ncol=1, rel_heights=c(1, 1.35), labels="AUTO")
save_plot("plots/ratios_functional_triplet_bins.pdf", p2, base_height=8, base_width=14)

####################
#
# Comparing rates within and without functional functional_classents
# Here looking at functional sites vs other sites in 1 kb windows
#
####################

summary_files_1kb <- list.files("~/Devel/mut_rates_func/functional/summary_tbls/", pattern=paste0("^exclusive_1kb_chr"), full.names=T)
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
                      by=.(func_class, map, hasCpG)]

b_group <- dat[, .(mean_B=mean(B)), by=.(func_class)] # joining mean B-value
tbl <- tbl_med_ratios[b_group, on=.(func_class)]
tbl <- tbl[func_class != "decile_11", ]
tbl[, class := factor(sub(".*_", "d_", func_class), levels=c(paste0("d_", 1:10), "enhancer", "promoter"))]
tbl[, annotation := "Functional"]
fwrite(tbl, "ratios_1kb_functional.csv")

map_labels <- c("carlson"="Carlson", "gnomad"="gnomAD", "roulette"="Roulette")

# separate table to plot segments because lines cannot be plotted with both color and linetype
seg_df <- tbl %>%
  arrange(map, hasCpG, class) %>%
  group_by(map, hasCpG) %>%
  mutate(x=class,
         y=med_ratio,
         xend=lead(class),
         yend=lead(med_ratio),
         mean_B_mid=(mean_B + lead(mean_B)) / 2) %>%
  filter(!is.na(xend)) %>%
  ungroup()

# "Ratios of mutation rates within (isolated) functional_elements w.r.t. 1 kb background"
p3 <- ggplot(tbl, aes(x=class, y=med_ratio)) +
  facet_wrap(~map, labeller=labeller(map=map_labels)) +
  geom_hline(yintercept=1, linetype="dashed", color="grey") + theme_classic() +
  geom_segment(data=seg_df,
               aes(x=x, xend=xend, y=y, yend=yend, color=mean_B_mid, linetype=hasCpG, group=interaction(map, hasCpG)),
               linewidth=0.9, lineend="round", inherit.aes=FALSE) +
  geom_point(aes(color=mean_B, group=hasCpG), size=3) +
  geom_errorbar(aes(ymin=med_ratio - se_ratio, ymax=med_ratio + se_ratio, color=mean_B), width=0.2, linewidth=0.7) +
  labs(x="Constraint class", y=expression(paste(mu, " ratio")), title=NULL) +
  scale_linetype_manual(name=NULL, values=c("FALSE"="solid", "TRUE"="dashed"), labels=c("TRUE"="With CpG", "FALSE"="Without CpG")) +
  scale_y_continuous(breaks=pretty_breaks(), limits=c(min(tbl$med_ratio) - tbl$se_ratio[which.min(tbl$med_ratio)],
                                                      max(tbl$med_ratio) + tbl$se_ratio[which.max(tbl$med_ratio)])) +
  scale_color_viridis_c(option="C", direction=1, name="B-value", breaks=c(round(min(tbl$mean_B) + 0.01, 2), round(max(tbl$mean_B) - 0.01, 2))) +
  guides(linetype=guide_legend(keywidth=unit(1.5, "cm"), keyheight=unit(0.2, "cm"),
                                 override.aes=list(color="black", linewidth=1.2, x=0, xend=1, y=0.5, yend=0.5))) +
  theme(strip.text=element_text(size=16),
        axis.title=element_text(size=20),
        axis.text=element_text(size=14),
        axis.text.x=element_text(size=14, angle=45, vjust=1, hjust=1),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16, margin=margin(b=20)),
        legend.position="bottom")
save_plot("plots/functional_ratios_isolated_1kb.pdf", p3, base_height=5, base_width=13)
