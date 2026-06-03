
library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
chr <- args[1]

dat <- fread(paste0("summary_tbls/summaries_chr", chr, ".csv.gz"))

dt_m <- dat %>%
  pivot_longer(cols=c(mean_bin_roulette, mean_bin_carlson, mean_bin_gnomad, se_bin_roulette, se_bin_carlson, se_bin_gnomad),
               names_to=c(".value", "source"), names_pattern="(mean|se)_(.*)")

p1 <- ggplot(dt_m, aes(x=benegas_bin, y=mean, color=source, group=paste0(source, benegas_bin < 15))) +
  geom_line(linewidth=1) + geom_point(size=3) + 
  geom_errorbar(aes(ymin=mean - se, ymax=mean + se), width=0.2, linewidth=0.6) +
  scale_x_continuous(breaks=c(1:nbins, nbins + 3), labels=c(as.character(1:nbins), "Outside")) + theme_classic() + 
  scale_color_manual(values=c("cyan3", "brown1", "green4"), name=NULL) +
  labs(x="Conservation score (percentile bin)", y="Mean Rate", color="Map",
       title="Mutation rates per conservation score") +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=14),
        strip.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.position="bottom")
save_plot(paste0("plots/benegas_ratios_chr", chr, ".pdf"), p1, base_height=6, base_width=10)


denoms <- dat[benegas_bin == 15, .(den_roulette=mean_bin_triplet_roulette,
                                   den_carlson=mean_bin_triplet_carlson,
                                   den_gnomad=mean_bin_triplet_gnomad,
                                   triplet=triplet)]

nums <- dat[benegas_bin %in% 1:12, .(num_roulette=mean_bin_triplet_roulette,
                                     num_carlson=mean_bin_triplet_carlson,
                                     num_gnomad=mean_bin_triplet_gnomad,
                                     benegas_bin=benegas_bin,
                                     triplet=triplet)]

ratios_benegas <- merge(nums, denoms, by="triplet", all.x = TRUE)
setorder(ratios_benegas, triplet, benegas_bin)

ratios_benegas[, `:=`(ratio_roulette=num_roulette / den_roulette,
                      ratio_carlson=num_carlson / den_carlson,
                      ratio_gnomad=num_gnomad / den_gnomad)]

ratios_benegas[, c("num_roulette", "num_carlson", "num_gnomad", "den_roulette", "den_carlson", "den_gnomad") := NULL]
ratios_benegas_m <- pivot_longer(ratios_benegas, cols=starts_with("ratio_"), names_to="map", values_to="ratio") %>% setDT()

ratios_benegas_m[, map := factor(sub("^ratio_", "", map), levels = c("roulette", "carlson", "gnomad"))]

p4 <- ggplot(filter(ratios_benegas_m, map!="gnomad"), aes(x=triplet, y=ratio, color=map)) + 
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
save_plot("~/Desktop/mut_rates/ratios_benegas_triplet.pdf", p4, base_height=7, base_width=12)

p4b <- ggplot(filter(ratios_benegas_m, map!="gnomad"),
              aes(x=triplet, y=ratio, color=benegas_bin)) + 
  facet_wrap(~map, ncol=1) + 
  geom_point(size=2.5) + theme_classic() + 
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_viridis_c(option="C", direction=1, name="Bin", breaks=c(1, 12)) +
  labs(x=NULL, y="Ratio of mutation rates", color="Map", 
       title="Ratios of mutation rates across Benegas bins and triplets") +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12, angle=90, hjust=1),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title=element_text(size=18),
        legend.box="horizontal")


# TODO join B-values and Rec map
# define 1 kb windows
# etc

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
