
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

CpGs <- c("ACG", "CCG", "GCG", "TCG", "CGA", "CGC", "CGG", "CGT")

gw_summary_files_M <- list.files("~/Devel/mut_rates_func/gpn-star/summaries_gpn-star/summary_tbls-M/", pattern=paste0("^summaries_chr"), full.names=T)
dat_M <- data.table::rbindlist(lapply(gw_summary_files_M, fread))
dat_M[, phylo:="M"]

gw_summary_files_V <- list.files("~/Devel/mut_rates_func/gpn-star/summaries_gpn-star/summary_tbls-V/", pattern=paste0("^summaries_chr"), full.names=T)
dat_V <- data.table::rbindlist(lapply(gw_summary_files_V, fread))
dat_V[, phylo:="V"]

gw_summary_files_P <- list.files("~/Devel/mut_rates_func/gpn-star/summaries_gpn-star/summary_tbls-P/", pattern=paste0("^summaries_chr"), full.names=T)
dat_P <- data.table::rbindlist(lapply(gw_summary_files_P, fread))
dat_P[, phylo:="P"]

dat <- rbindlist(list(dat_M, dat_V, dat_P))
nclasses <- length(unique(na.omit(dat$benegas_class))) - 1 # -1 excludes last class (putatively neutral sites)

####################
#
# Mutation rates per trinucleotide (descriptive)
#
####################

dat[, `:=`(mean_roulette_triplets=mean(roulette, na.rm=T),
           mean_carlson_triplets=mean(carlson, na.rm=T),
           mean_gnomad_triplets=mean(gnomad, na.rm=T)), by=c("triplet", "phylo")]

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
# Comparing rates within and without GPN-star elements
# Here looking at GPN-star sites vs genome-wide background
#
####################

dat_withCpG <- dplyr::select(dat, c(mean_class_triplet_roulette, mean_class_triplet_carlson, mean_class_triplet_gnomad,
                                    se_class_roulette, se_class_gnomad, se_class_carlson,
                                    num_sites_roulette, num_sites_carlson, num_sites_gnomad, benegas_class, triplet, phylo)) %>% setDT()

# computing means across chromosomes, weighted by num sites in each window
dat_withCpG[, mean_class_triplet_roulette_gw := sum(mean_class_triplet_roulette * num_sites_roulette, na.rm=T) / sum(num_sites_roulette, na.rm = TRUE), by=.(benegas_class, phylo)]
dat_withCpG[, mean_class_triplet_carlson_gw := sum(mean_class_triplet_carlson * num_sites_carlson, na.rm=T) / sum(num_sites_carlson, na.rm = TRUE), by=.(benegas_class, phylo)]
dat_withCpG[, mean_class_triplet_gnomad_gw := sum(mean_class_triplet_gnomad * num_sites_gnomad, na.rm=T) / sum(num_sites_gnomad, na.rm = TRUE), by=.(benegas_class, phylo)]

# SEs across chromosomes
gw_se_roulette <- dat_withCpG[, { m <- dat_withCpG[benegas_class == .BY$benegas_class, mean_class_triplet_roulette_gw]
v <- sum(num_sites_roulette * (se_class_roulette^2 + (mean_class_triplet_roulette_gw - m)^2)) / sum(num_sites_roulette) 
.(se_gw_roulette = sqrt(v)) }, by=.(benegas_class, phylo)]

gw_se_carlson <- dat_withCpG[, { m <- dat_withCpG[benegas_class == .BY$benegas_class, mean_class_triplet_carlson_gw]
v <- sum(num_sites_carlson * (se_class_carlson^2 + (mean_class_triplet_carlson_gw - m)^2)) / sum(num_sites_carlson) 
.(se_gw_carlson = sqrt(v)) }, by=.(benegas_class, phylo)]

gw_se_gnomad <- dat_withCpG[, { m <- dat_withCpG[benegas_class == .BY$benegas_class, mean_class_triplet_gnomad_gw]
v <- sum(num_sites_gnomad * (se_class_gnomad^2 + (mean_class_triplet_gnomad_gw - m)^2)) / sum(num_sites_gnomad) 
.(se_gw_gnomad = sqrt(v)) }, by=.(benegas_class, phylo)]

dat_withCpG <- unique(dat_withCpG, by=c("benegas_class", "phylo")) %>% 
  dplyr::select(., c(benegas_class, mean_class_triplet_roulette_gw, mean_class_triplet_carlson_gw, mean_class_triplet_gnomad_gw, phylo)) %>% setDT()

dat_withCpG <- merge(dat_withCpG, gw_se_roulette, by=c("benegas_class", "phylo"))
dat_withCpG <- merge(dat_withCpG, gw_se_carlson, by=c("benegas_class", "phylo"))
dat_withCpG <- merge(dat_withCpG, gw_se_gnomad, by=c("benegas_class", "phylo"))

denom_roulette <- dat_withCpG[benegas_class==15, mean_class_triplet_roulette_gw]
denom_carlson <- dat_withCpG[benegas_class==15, mean_class_triplet_carlson_gw]
denom_gnomad <- dat_withCpG[benegas_class==15, mean_class_triplet_gnomad_gw]

names(denom_roulette) <- c("M", "P", "V")
names(denom_carlson) <- c("M", "P", "V")
names(denom_gnomad) <- c("M", "P", "V")

dat_withCpG[, ratio_roulette := mean_class_triplet_roulette_gw / denom_roulette[phylo]]
dat_withCpG[, ratio_carlson := mean_class_triplet_carlson_gw / denom_carlson[phylo]]
dat_withCpG[, ratio_gnomad := mean_class_triplet_gnomad_gw / denom_carlson[phylo]]
dat_withCpG[, ratio_roulette_se := se_gw_roulette / denom_roulette[phylo]]
dat_withCpG[, ratio_carlson_se := se_gw_carlson / denom_carlson[phylo]]
dat_withCpG[, ratio_gnomad_se := se_gw_gnomad / denom_gnomad[phylo]]
dat_withCpG[, CpG := T] 
dat_withCpG <- dat_withCpG[benegas_class < 15, .(benegas_class, ratio_roulette, ratio_carlson, ratio_gnomad,
                                                 ratio_roulette_se, ratio_carlson_se, ratio_gnomad_se, CpG, phylo)]

# filtering out CpG sites
dat_nonCpG <- dplyr::select(dat, c(mean_class_triplet_roulette, mean_class_triplet_carlson, mean_class_triplet_gnomad,
                                   se_class_roulette, se_class_gnomad, se_class_carlson,
                                   num_sites_roulette, num_sites_carlson, num_sites_gnomad, benegas_class, triplet, phylo)) %>% 
  filter(., !triplet %in% CpGs) %>% setDT()

# computing means across chromosomes, weighted by num sites in each window
dat_nonCpG[, mean_class_triplet_roulette_gw := sum(mean_class_triplet_roulette * num_sites_roulette, na.rm=T) / sum(num_sites_roulette, na.rm = TRUE), by=.(benegas_class, phylo)]
dat_nonCpG[, mean_class_triplet_carlson_gw := sum(mean_class_triplet_carlson * num_sites_carlson, na.rm=T) / sum(num_sites_carlson, na.rm = TRUE), by=.(benegas_class, phylo)]
dat_nonCpG[, mean_class_triplet_gnomad_gw := sum(mean_class_triplet_gnomad * num_sites_gnomad, na.rm=T) / sum(num_sites_gnomad, na.rm = TRUE), by=.(benegas_class, phylo)]

# SEs across chromosomes
gw_se_roulette <- dat_nonCpG[, { m <- dat_nonCpG[benegas_class == .BY$benegas_class, mean_class_triplet_roulette_gw]
v <- sum(num_sites_roulette * (se_class_roulette^2 + (mean_class_triplet_roulette_gw - m)^2)) / sum(num_sites_roulette) 
.(se_gw_roulette = sqrt(v)) }, by=.(benegas_class, phylo)]

gw_se_carlson <- dat_nonCpG[, { m <- dat_nonCpG[benegas_class == .BY$benegas_class, mean_class_triplet_carlson_gw]
v <- sum(num_sites_carlson * (se_class_carlson^2 + (mean_class_triplet_carlson_gw - m)^2)) / sum(num_sites_carlson) 
.(se_gw_carlson = sqrt(v)) }, by=.(benegas_class, phylo)]

gw_se_gnomad <- dat_nonCpG[, { m <- dat_nonCpG[benegas_class == .BY$benegas_class, mean_class_triplet_gnomad_gw]
v <- sum(num_sites_gnomad * (se_class_gnomad^2 + (mean_class_triplet_gnomad_gw - m)^2)) / sum(num_sites_gnomad) 
.(se_gw_gnomad = sqrt(v)) }, by=.(benegas_class, phylo)]

dat_nonCpG <- unique(dat_nonCpG, by=c("benegas_class", "phylo")) %>% 
  dplyr::select(., c(benegas_class, mean_class_triplet_roulette_gw, mean_class_triplet_carlson_gw, mean_class_triplet_gnomad_gw, phylo)) %>% setDT()

dat_nonCpG <- merge(dat_nonCpG, gw_se_roulette, by=c("benegas_class", "phylo"))
dat_nonCpG <- merge(dat_nonCpG, gw_se_carlson, by=c("benegas_class", "phylo"))
dat_nonCpG <- merge(dat_nonCpG, gw_se_gnomad, by=c("benegas_class", "phylo"))

denom_roulette <- dat_nonCpG[benegas_class==15, mean_class_triplet_roulette_gw]
denom_carlson <- dat_nonCpG[benegas_class==15, mean_class_triplet_carlson_gw]
denom_gnomad <- dat_nonCpG[benegas_class==15, mean_class_triplet_gnomad_gw]

names(denom_roulette) <- c("M", "P", "V")
names(denom_carlson) <- c("M", "P", "V")
names(denom_gnomad) <- c("M", "P", "V")

dat_nonCpG[, ratio_roulette := mean_class_triplet_roulette_gw / denom_roulette[phylo]]
dat_nonCpG[, ratio_carlson := mean_class_triplet_carlson_gw / denom_carlson[phylo]]
dat_nonCpG[, ratio_gnomad := mean_class_triplet_gnomad_gw / denom_gnomad[phylo]]
dat_nonCpG[, ratio_roulette_se := se_gw_roulette / denom_roulette[phylo]]
dat_nonCpG[, ratio_carlson_se := se_gw_carlson / denom_carlson[phylo]]
dat_nonCpG[, ratio_gnomad_se := se_gw_gnomad / denom_gnomad[phylo]]

dat_nonCpG[, CpG := F] 
dat_nonCpG <- dat_nonCpG[benegas_class < 15, .(benegas_class, ratio_roulette, ratio_carlson, ratio_gnomad,
                                               ratio_roulette_se, ratio_carlson_se, ratio_gnomad_se, CpG, phylo)]

m_means <- pivot_longer(
  rbind.data.frame(dat_withCpG, dat_nonCpG),
  cols = c(ratio_roulette, ratio_carlson, ratio_gnomad),
  names_to="map",
  values_to="ratio"
) %>% setDT()

m_ses <- pivot_longer(
  rbind.data.frame(dat_withCpG, dat_nonCpG),
  cols = c(ratio_roulette_se, ratio_carlson_se, ratio_gnomad_se),
  names_to = "map",
  values_to = "se"
) %>% setDT()

m_ses[, map := sub("_se$", "", map)] 

m_ratios <- m_means[m_ses, on=c("benegas_class", "CpG", "map", "phylo")]
m_ratios[, annotation := "GNP-Star"]
m_ratios <- m_ratios[, .(benegas_class, CpG, map, ratio, se, annotation, phylo)]
fwrite(m_ratios, "gw_ratios_gpn-star.csv")

# "Ratios of mutation rates within Benegas elements w.r.t. genome-wide background"
p1 <- ggplot(m_ratios, aes(x=benegas_class, y=ratio, color=map, group=paste0(map, CpG))) +
  facet_wrap(~phylo, nrow=1) +
  geom_line(aes(linetype=CpG), linewidth=1) + geom_point(size=3) + 
  geom_errorbar(aes(ymin = ratio - se, ymax = ratio + se), width = 0.2) +
  geom_hline(yintercept=1, linetype="dashed", color="grey") +
  scale_x_continuous(breaks=1:nclasses) + theme_classic() + 
  scale_y_continuous(breaks=seq(0.5, 1.5, 0.2)) +
  scale_color_manual(values=c("brown1", "cyan3", "seagreen"), name=NULL,
                     labels=c("ratio_roulette"="Roulette", "ratio_gnomad"="gnomAD", "ratio_carlson"="Carlson")) +
  labs(x="Constraint class", y=expression(paste(mu["÷"])), title=NULL) +
  scale_linetype_manual(name=NULL, values=c("FALSE"="solid", "TRUE"="dashed"),
                        labels=c("TRUE"="With CpG", "FALSE"="Without CpG")) +
  guides(linetype=guide_legend(keywidth=unit(1.5, "cm"), keyheight=unit(0.2, "cm"),
                               override.aes=list(color="black", linewidth=1.2, x=0, xend=1, y=0.5, yend=0.5))) +
  geom_errorbar(data=m_ratios[grepl("_se$", map)], aes(ymin=ratio - ratio, ymax=ratio + ratio)) +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=14),
        strip.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.position="bottom")
save_plot("plots/gpn-star_ratios.pdf", p1, base_height=7, base_width=14)

# stratifying by trinucleotide context
dat[, mean_class_triplet_gw_roulette := sum(mean_class_triplet_roulette * num_sites_roulette, na.rm=T) / sum(num_sites_roulette, na.rm = TRUE), by=.(benegas_class, triplet, phylo)]
dat[, mean_class_triplet_gw_carlson := sum(mean_class_triplet_carlson * num_sites_carlson, na.rm=T) / sum(num_sites_carlson, na.rm = TRUE),, by=.(benegas_class, triplet, phylo)]
dat[, mean_class_triplet_gw_gnomad := sum(mean_class_triplet_gnomad * num_sites_gnomad, na.rm=T) / sum(num_sites_gnomad, na.rm = TRUE),, by=.(benegas_class, triplet, phylo)]

denoms <- dat[benegas_class == 15, 
              .(den_roulette=mean_class_triplet_gw_roulette,
                den_carlson=mean_class_triplet_gw_carlson,
                den_gnomad=mean_class_triplet_gw_gnomad,
                triplet=triplet, chrom=chrom, phylo=phylo)]

nums <- dat[benegas_class %in% 1:12, .(num_roulette=mean_class_triplet_gw_roulette,
                                       num_carlson=mean_class_triplet_gw_carlson,
                                       num_gnomad=mean_class_triplet_gw_gnomad,
                                       benegas_class=benegas_class,
                                       triplet=triplet, chrom=chrom, phylo=phylo)]

ratios_benegas <- merge(nums, denoms, by=c("triplet", "chrom", "phylo"), all.x = TRUE)
setorder(ratios_benegas, chrom, triplet, benegas_class, phylo)

ratios_benegas[, `:=`(ratio_roulette=num_roulette / den_roulette,
                      ratio_carlson=num_carlson / den_carlson,
                      ratio_gnomad=num_gnomad / den_gnomad)]

ratios_benegas[, c("num_roulette", "num_carlson", "num_gnomad", "den_roulette", "den_carlson", "den_gnomad") := NULL]
ratios_benegas_m <- pivot_longer(ratios_benegas, cols=starts_with("ratio_"), names_to="map", values_to="ratio") %>% setDT()

ratios_benegas_m[, map := factor(sub("^ratio_", "", map), levels = c("roulette", "carlson", "gnomad"))]
ratios_benegas_m[, annotation := "GPN-star"]
fwrite(ratios_benegas_m, "ratios_triplets_pgn-star.csv")

# "Ratios of mutation rates within GPN-star elements stratified by triplet context"
p2a <- ggplot(filter(ratios_benegas_m, map=="carlson"),
              aes(x=triplet, y=ratio, color=benegas_class)) + 
  annotate(xmin = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) - 0.5,
           xmax = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) + 0.5,
           geom="rect", ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_point(size=2.5) + theme_classic() + 
  facet_wrap(~phylo, ncol=1) +
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_viridis_c(option="C", direction=1, name="Constraint Class", breaks=c(1, 12)) +
  labs(x=NULL, y=expression(paste(mu["÷"])), title=NULL) +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12, angle=90, hjust=1, vjust=0.5),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title=element_text(size=18),
        legend.box="horizontal")

p2b <- ggplot(filter(ratios_benegas_m, map=="roulette"),
              aes(x=triplet, y=ratio, color=benegas_class)) + 
  annotate(xmin = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) - 0.5,
           xmax = which(levels(factor(ratios_benegas_m$triplet)) %in% CpGs) + 0.5,
           geom="rect", ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_point(size=2.5) + theme_classic() + 
  facet_wrap(~phylo, ncol=1) +
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_viridis_c(option="C", direction=1, name="Constraint Class", breaks=c(1, 12)) +
  labs(x=NULL, y=NULL, title=NULL) +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12, angle=90, hjust=1, vjust=0.5),
        legend.text=element_text(size=16),
        legend.position="none",
        legend.title=element_text(size=18),
        legend.box="horizontal")

leg <- get_legend(p2a)
p2a <- p2a + theme(legend.position="none")
p2 <- plot_grid(plot_grid(p2a, p2b, ncol=2, rel_heights=c(1, 1.35), labels="AUTO"), leg, ncol=1, rel_heights=c(1, 0.5))
save_plot("plots/ratios_gpn-star_triplet_classes.pdf", p2, base_height=16, base_width=20)

####################
#
# Comparing rates within and without GPN-Star elements
# Here looking at GPN-Star sites vs other sites in 1 kb windows
#
####################

summary_files_1kb_M <- list.files("~/Devel/mut_rates_func/gpn-star/summaries_gpn-star/summary_tbls-M/", pattern=paste0("^exclusive_1kb_chr"), full.names=T)
withCpG_files <- summary_files_1kb_M[!grepl("nonCpG", summary_files_1kb_M)]
nonCpG_files <- summary_files_1kb_M[grepl("nonCpG", summary_files_1kb_M)]

withCpG <- data.table::rbindlist(lapply(withCpG_files, fread), fill=T, use.names=T)
nonCpG <- data.table::rbindlist(lapply(nonCpG_files, fread), fill=T, use.names=T)

withCpG[, hasCpG := T]
nonCpG[, hasCpG := F]
dat_M <- rbind.data.frame(withCpG, nonCpG)
dat_M$phylo <- "M"

summary_files_1kb_P <- list.files("~/Devel/mut_rates_func/gpn-star/summaries_gpn-star/summary_tbls-P/", pattern=paste0("^exclusive_1kb_chr"), full.names=T)
withCpG_files <- summary_files_1kb_P[!grepl("nonCpG", summary_files_1kb_P)]
nonCpG_files <- summary_files_1kb_P[grepl("nonCpG", summary_files_1kb_P)]

withCpG <- data.table::rbindlist(lapply(withCpG_files, fread), fill=T, use.names=T)
nonCpG <- data.table::rbindlist(lapply(nonCpG_files, fread), fill=T, use.names=T)

withCpG[, hasCpG := T]
nonCpG[, hasCpG := F]
dat_P <- rbind.data.frame(withCpG, nonCpG)
dat_P$phylo <- "P"

summary_files_1kb_V <- list.files("~/Devel/mut_rates_func/gpn-star/summaries_gpn-star/summary_tbls-V/", pattern=paste0("^exclusive_1kb_chr"), full.names=T)
withCpG_files <- summary_files_1kb_V[!grepl("nonCpG", summary_files_1kb_V)]
nonCpG_files <- summary_files_1kb_V[grepl("nonCpG", summary_files_1kb_V)]

withCpG <- data.table::rbindlist(lapply(withCpG_files, fread), fill=T, use.names=T)
nonCpG <- data.table::rbindlist(lapply(nonCpG_files, fread), fill=T, use.names=T)

withCpG[, hasCpG := T]
nonCpG[, hasCpG := F]
dat_V <- rbind.data.frame(withCpG, nonCpG)
dat_V$phylo <- "V"

dat <- rbindlist(list(dat_M, dat_P, dat_V))
setnames(dat, old="benegas_group", new="gpn_star_group")

# brief pause to look at length distributions of these isolated elements
tbl <- dat[hasCpG==T & !duplicated(pos),] %>%
  dplyr::select(., c(starts_with("n_sites_class"), "phylo")) %>% 
  pivot_longer(., cols=starts_with("n_sites_class"), values_to="lengths", names_to="class") %>%
  dplyr::filter(., lengths > 0, class != "n_sites_class_15") %>% setDT()

tbl[, class := as.integer(sub(".*_", "", class))]
counts <- as.data.frame(table(tbl$class))
names(counts) <- c("class", "log10(# windows)")
counts$class <- as.integer(counts$class)
counts$`log10(# windows)` <- log10(counts$`log10(# windows)`)
tbl <- merge(tbl, counts, by="class")

p <- ggplot(tbl, aes(x=class, y=lengths, group=class, color=`log10(# windows)`)) +
  geom_boxplot() + facet_wrap(~phylo) + theme_classic() + 
  scale_x_continuous(breaks=1:12) +
  labs(x="Constraint class", y="Length within 1kb windows", title=NULL) +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=14),
        strip.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position="bottom")
save_plot("plots/gpn-star_lengths_1kb.pdf", p, base_height=5, base_width=12)
# done with lengths, back to the main task


# NOTE use log(ratio)?
# using median as a summary to mitigate outliers
tbl_med_ratios <- dat[, .(med_ratio=median(value, na.rm=T), 
                          se_ratio=sd(value, na.rm=T) / sqrt(sum(!is.na(value)))),
                      by=.(gpn_star_group, map, hasCpG, phylo)]

b_group <- dat[, .(mean_B = mean(B)), by=.(gpn_star_group, phylo)] # joining mean B-value
tbl <- tbl_med_ratios[b_group, on=.(gpn_star_group, phylo)]
tbl[, annotation := "GPN-star"]
fwrite(tbl, "ratios_1kb_benegas.csv")

map_labels <- c("carlson"="Carlson", "gnomad"="gnomAD", "roulette"="Roulette")

# separate table to plot segments because lines cannot be plotted with both color and linetype
seg_df <- tbl %>%
  arrange(map, hasCpG, phylo, gpn_star_group) %>%
  group_by(map, hasCpG, phylo) %>%
  mutate(x=gpn_star_group,
         y=med_ratio,
         xend=lead(gpn_star_group),
         yend=lead(med_ratio),
         mean_B_mid=(mean_B + lead(mean_B)) / 2) %>%
  filter(!is.na(xend)) %>%
  ungroup()

# "Ratios of mutation rates within (isolated) GPN-star elements w.r.t. 1 kb background"
p3 <- ggplot(tbl, aes(x=gpn_star_group, y=med_ratio)) +
  facet_grid(phylo~map, labeller=labeller(map=map_labels)) +
  geom_hline(yintercept=1, linetype="dashed", color="grey") + theme_classic() +
  geom_segment(data=seg_df,
               aes(x=x, xend=xend, y=y, yend=yend, color=mean_B_mid, linetype=hasCpG, group=interaction(map, hasCpG, phylo)),
               linewidth=0.9, lineend="round", inherit.aes=FALSE) +
  geom_point(aes(color=mean_B, group=hasCpG), size=3) +
  geom_errorbar(aes(ymin=med_ratio - se_ratio, ymax=med_ratio + se_ratio, color=mean_B), width=0.2, linewidth=0.7) +
  labs(x="Constraint class", y=expression(paste(mu["÷"])), title=NULL) +
  scale_x_continuous(breaks=1:12) +
  scale_linetype_manual(name=NULL, values=c("FALSE"="solid", "TRUE"="dashed"), labels=c("TRUE"="With CpG", "FALSE"="Without CpG")) +
  scale_y_continuous(breaks=pretty_breaks()) + 
  coord_cartesian(ylim = c(0.25, 1.75)) +
  scale_color_viridis_c(option="C", direction=1, name="B-value", breaks=c(round(min(tbl$mean_B) + 0.01, 2), round(max(tbl$mean_B) - 0.01, 2))) +
  guides(linetype = guide_legend(keywidth = unit(1.5, "cm"), keyheight = unit(0.2, "cm"),
                                 override.aes = list(color = "black", linewidth = 1.2, x = 0, xend = 1, y = 0.5, yend = 0.5))) +
  theme(strip.text=element_text(size=16),
        axis.title=element_text(size=20),
        axis.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16, margin=margin(b=20)),
        legend.position="bottom")
save_plot("plots/gpn-star_ratios_isolated_1kb.pdf", p3, base_height=8, base_width=14)
