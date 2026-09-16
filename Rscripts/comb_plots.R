library(data.table)
library(tidyverse)
library(cowplot)

####################
#
# genome-wide ratios
#
####################

gw_benegas <- fread("benegas/gw_ratios_benegas.csv")
gw_phast <- fread("phastCons/gw_ratios_phastCons.csv")
gw_func <- fread("functional/gw_ratios_functional.csv")

setnames(gw_benegas, old="benegas_class", new="class")
setnames(gw_phast, old="phast_class", new="class")
setnames(gw_func, old="functional_class", new="class")

gw_phast[, class := factor(class, levels=1:12)]
gw_benegas[, class := factor(class, levels=1:12)]
gw_func[, class := factor(fifelse(class == "enhancer", "e",
                          fifelse(class == "promoter", "p",
                          gsub("^d_", "", class))),
                          levels = c(as.character(1:10), "e", "p"))]

gw_ratios <- rbind.data.frame(gw_benegas, gw_phast, gw_func) %>% setDT()
gw_ratios[, plot_group :=
            fifelse(class %chin% as.character(1:10), "quant",
                    fifelse(class == "e", "enhancer",
                            fifelse(class == "p", "promoter", NA_character_)))]

p1a <- ggplot(filter(gw_ratios, annotation == "Functional"),
             aes(x=class, y=ratio, color=map, shape=CpG, group=paste0(map, CpG, plot_group))) + 
  geom_line(linewidth=1) +
  geom_point(size=4) + geom_errorbar(aes(ymin = ratio - se, ymax = ratio + se), width = 0.2, color="black") +
  geom_hline(yintercept=1, linetype="dashed", color="grey") + theme_classic() + 
  scale_color_manual(values=c("brown1", "cyan3", "seagreen"), name=NULL,
                     labels=c("ratio_roulette"="Roulette", "ratio_gnomad"="gnomAD", "ratio_carlson"="Carlson")) +
  labs(x="Constraint class", y=expression(paste(mu["÷"])), title=NULL) +
  scale_shape_manual(name=NULL, values=c("FALSE"=16, "TRUE"=17),
                     labels=c("TRUE"="With CpG", "FALSE"="Without CpG")) +
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=22),
        axis.text=element_text(size=14),
        legend.text=element_text(size=16),
        legend.position="none")

p1b <- ggplot(filter(gw_ratios, annotation == "phastCons"),
              aes(x=class, y=ratio, color=map, shape=CpG, group=paste0(map, CpG))) + theme_classic() + 
  geom_line(linewidth=1) +
  geom_point(size=4) + geom_errorbar(aes(ymin = ratio - se, ymax = ratio + se), width = 0.2, color="black") +
  geom_hline(yintercept=1, linetype="dashed", color="grey") +
  scale_color_manual(values=c("brown1", "cyan3", "seagreen"), name=NULL,
                     labels=c("ratio_roulette"="Roulette", "ratio_gnomad"="gnomAD", "ratio_carlson"="Carlson")) +
  labs(x="Constraint class", y=NULL, title=NULL) +
  scale_shape_manual(name=NULL, values=c("FALSE"=16, "TRUE"=17),
                     labels=c("TRUE"="With CpG", "FALSE"="Without CpG")) +
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=22),
        axis.text=element_text(size=14),
        legend.text=element_text(size=16),
        legend.position="none")

p1c <- ggplot(filter(gw_ratios, annotation == "Benegas"),
              aes(x=class, y=ratio, color=map, shape=CpG, group=paste0(map, CpG))) + theme_classic() + 
  geom_line(linewidth=1) +
  geom_point(size=4) + geom_errorbar(aes(ymin = ratio - se, ymax = ratio + se), width = 0.2, color="black") +
  geom_hline(yintercept=1, linetype="dashed", color="grey") +
  scale_color_manual(values=c("brown1", "cyan3", "seagreen"), name=NULL,
                     labels=c("ratio_roulette"="Roulette", "ratio_gnomad"="gnomAD", "ratio_carlson"="Carlson")) +
  labs(x="Constraint class", y=NULL, title=NULL) +
  scale_shape_manual(name=NULL, values=c("FALSE"=16, "TRUE"=17),
                        labels=c("TRUE"="With CpG", "FALSE"="Without CpG")) +
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=22),
        axis.text=element_text(size=14),
        legend.text=element_text(size=16),
        legend.position="bottom")

leg <- get_legend(p1c)

p1c <- p1c + theme(legend.position="none")

p1 <- plot_grid(plot_grid(p1a, p1b, p1c, nrow=1, labels="AUTO", label_size=14, rel_widths=c(1.2, 1, 1)), 
                leg, nrow=2, rel_heights=c(1, 0.05))
save_plot("gw_ratios.pdf", p1, base_height=6, base_width=15)

####################
#
# ratios stratified by triplets
#
####################

benegas_triplets <- fread("benegas/ratios_triplets_benegas.csv")
phast_triplets <- fread("phastcons/ratios_triplets_phastCons.csv")
func_triplets <- fread("functional/ratios_triplets_functional.csv")

p2a <- ggplot(filter(func_triplets, map=="carlson")) +
  annotate(xmin=which(levels(factor(func_triplets$triplet)) %in% CpGs) - 0.5,
           xmax=which(levels(factor(func_triplets$triplet)) %in% CpGs) + 0.5,
           geom="rect", ymin=-Inf, ymax=Inf, fill="grey85", alpha=0.6) +
  geom_point(data=subset(func_triplets, map=="carlson" & grepl("^decile_", functional_class)),
             aes(x=triplet, y=ratio, color=as.numeric(class)), size=2.5) +
  scale_color_viridis_c(option="C", direction=1, name="Decile", breaks=c(1, 10)) +
  ggnewscale::new_scale_color() +
  geom_point(data=subset(func_triplets, map=="carlson" & functional_class %in% c("enhancer", "promoter")),
             aes(x=triplet, y=ratio, color=functional_class), size=3) +
  scale_color_manual(name=NULL, values=c(enhancer="seagreen", promoter="cyan3")) +
  geom_hline(yintercept=1, linetype="dashed", color="black") +
  theme_classic() +
  labs(x=NULL, y=expression(paste(mu, " ratio")), title=NULL) +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_blank(),
        legend.position="none")

p2b <- ggplot(filter(phast_triplets, map=="carlson"),
              aes(x=triplet, y=ratio, color=phast_class)) + 
  annotate(xmin=which(levels(factor(phast_triplets$triplet)) %in% CpGs) - 0.5,
           xmax=which(levels(factor(phast_triplets$triplet)) %in% CpGs) + 0.5,
           geom="rect", ymin=-Inf, ymax=Inf, fill="grey85", alpha=0.6) +
  geom_point(size=2.5) + theme_classic() + 
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_viridis_c(option="C", direction=1, name="Class", breaks=c(1, 12)) +
  labs(x=NULL, y=expression(paste(mu, " ratio")), title=NULL) +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_blank(),
        legend.text=element_text(size=16),
        legend.position="none",
        legend.title=element_text(size=18),
        legend.box="horizontal")

p2c <- ggplot(filter(benegas_triplets, map=="carlson"),
              aes(x=triplet, y=ratio, color=benegas_class)) + 
  annotate(xmin=which(levels(factor(benegas_triplets$triplet)) %in% CpGs) - 0.5,
           xmax=which(levels(factor(benegas_triplets$triplet)) %in% CpGs) + 0.5,
           geom="rect", ymin=-Inf, ymax=Inf, fill="grey85", alpha=0.6) +
  geom_point(size=2.5) + theme_classic() + 
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_viridis_c(option="C", direction=1, name="Class", breaks=c(1, 12)) +
  labs(x=NULL, y=expression(paste(mu, " ratio")), title=NULL) +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12, angle=90, hjust=1, vjust=0.5),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title=element_text(size=18),
        legend.box="horizontal")

p2d <- ggplot(filter(func_triplets, map=="roulette")) +
  annotate(xmin=which(levels(factor(func_triplets$triplet)) %in% CpGs) - 0.5,
           xmax=which(levels(factor(func_triplets$triplet)) %in% CpGs) + 0.5,
           geom="rect", ymin=-Inf, ymax=Inf, fill="grey85", alpha=0.6) +
  geom_point(data=subset(func_triplets, map=="carlson" & grepl("^decile_", functional_class)),
             aes(x=triplet, y=ratio, color=as.numeric(class)), size=2.5) +
  scale_color_viridis_c(option="C", direction=1, name="Decile", breaks=c(1, 10)) +
  ggnewscale::new_scale_color() +
  geom_point(data=subset(func_triplets, map=="carlson" & functional_class %in% c("enhancer", "promoter")),
             aes(x=triplet, y=ratio, color=functional_class), size=3) +
  scale_color_manual(name=NULL, values=c(enhancer="seagreen", promoter="cyan3")) +
  geom_hline(yintercept=1, linetype="dashed", color="black") +
  theme_classic() +
  labs(x=NULL, y=NULL, title=NULL) +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_blank(),
        legend.position="none")

p2e <- ggplot(filter(phast_triplets, map=="roulette"),
              aes(x=triplet, y=ratio, color=phast_class)) + 
  annotate(xmin=which(levels(factor(phast_triplets$triplet)) %in% CpGs) - 0.5,
           xmax=which(levels(factor(phast_triplets$triplet)) %in% CpGs) + 0.5,
           geom="rect", ymin=-Inf, ymax=Inf, fill="grey85", alpha=0.6) +
  geom_point(size=2.5) + theme_classic() + 
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_viridis_c(option="C", direction=1, name="Class", breaks=c(1, 12)) +
  labs(x=NULL, y=NULL, title=NULL) +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_blank(),
        legend.text=element_text(size=16),
        legend.position="none",
        legend.title=element_text(size=18),
        legend.box="horizontal")

p2f <- ggplot(filter(benegas_triplets, map=="roulette"),
              aes(x=triplet, y=ratio, color=benegas_class)) + 
  annotate(xmin=which(levels(factor(benegas_triplets$triplet)) %in% CpGs) - 0.5,
           xmax=which(levels(factor(benegas_triplets$triplet)) %in% CpGs) + 0.5,
           geom="rect", ymin=-Inf, ymax=Inf, fill="grey85", alpha=0.6) +
  geom_point(size=2.5) + theme_classic() + 
  geom_hline(yintercept=1, linetype="dashed", color="green4") +
  scale_color_viridis_c(option="C", direction=1, name="Class", breaks=c(1, 12)) +
  labs(x=NULL, y=NULL, title=NULL) +
  theme(strip.text=element_text(size=18),
        axis.title=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=12, angle=90, hjust=1, vjust=0.5),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title=element_text(size=18),
        legend.box="horizontal")

p2 <- plot_grid(p2a, p2d, p2b, p2e, p2c, p2f, labels="AUTO", label_size=16, ncol=2, align="vh")
save_plot("ratios_triplets.pdf", p2, base_height=12, base_width=18)

####################
#
# ratios in 1 kb windows
#
####################

benegas_1kb <- fread("benegas/ratios_1kb_bengas.csv")
phast_1kb <- fread("phastCons/ratios_1kb_phastCons.csv")
func_1kb <- fread("functional/ratios_1kb_functional.csv") 
func_1kb[, func_class := NULL]

setnames(benegas_1kb, old="benegas_group", new="class")
setnames(phast_1kb, old="phast_group", new="class")

phast_1kb[, class := factor(class, levels=1:12)]
benegas_1kb[, class := factor(class, levels=5:12)]
func_1kb[, class := factor(fifelse(class == "enhancer", "e",
                           fifelse(class == "promoter", "p", gsub("^d_", "", class))),
                           levels = c(as.character(1:10), "e", "p"))]

func_1kb[, plot_group :=
            fifelse(class %chin% as.character(1:10), "quant",
                    fifelse(class == "e", "enhancer",
                            fifelse(class == "p", "promoter", NA_character_)))]

map_labels <- c("carlson"="Carlson", "gnomad"="gnomAD", "roulette"="Roulette")

# separate table to plot segments because lines cannot be plotted with both color and linetype
seg_benegas <- benegas_1kb %>%
  arrange(map, hasCpG, class) %>%
  group_by(map, hasCpG) %>%
  mutate(x=class,
         y=med_ratio,
         xend=lead(class),
         yend=lead(med_ratio),
         mean_B=(mean_B + lead(mean_B)) / 2) %>%
  filter(!is.na(xend)) %>%
  ungroup()

seg_phast <- phast_1kb %>%
  arrange(map, hasCpG, class) %>%
  group_by(map, hasCpG) %>%
  mutate(x=class,
         y=med_ratio,
         xend=lead(class),
         yend=lead(med_ratio),
         mean_B=(mean_B + lead(mean_B)) / 2) %>%
  filter(!is.na(xend)) %>%
  ungroup()

seg_func <- func_1kb %>%
  arrange(map, hasCpG, class) %>%
  group_by(map, hasCpG, plot_group) %>%
  mutate(x=class,
         y=med_ratio,
         xend=lead(class),
         yend=lead(med_ratio),
         mean_B=(mean_B + lead(mean_B)) / 2) %>%
  filter(!is.na(xend)) %>%
  ungroup()

p3a <- ggplot(func_1kb, aes(x=class, y=med_ratio, shape=hasCpG)) +
  facet_wrap(~map, labeller=labeller(map=map_labels)) +
  geom_hline(yintercept=1, linetype="dashed", color="grey") + theme_classic() +
  geom_segment(data=seg_func,
               aes(x=x, xend=xend, y=y, yend=yend, color=mean_B, group=interaction(map, hasCpG, plot_group)),
               linewidth=0.9, lineend="round", inherit.aes=FALSE) +
  geom_point(aes(color=mean_B), size=3) +
  geom_errorbar(aes(ymin=med_ratio - se_ratio, ymax=med_ratio + se_ratio, color=mean_B), width=0.2, linewidth=0.7) +
  labs(x=NULL, y=expression(paste(mu["÷"])), title=NULL) +
  scale_shape_manual(name=NULL, values=c("FALSE"=16, "TRUE"=17),
                     labels=c("TRUE"="With CpG", "FALSE"="Without CpG"), guide="none") +
  scale_y_continuous(breaks=pretty_breaks(), limits=c(min(func_1kb$med_ratio) - func_1kb$se_ratio[which.min(func_1kb$med_ratio)],
                                                      max(func_1kb$med_ratio) + func_1kb$se_ratio[which.max(func_1kb$med_ratio)])) +
  scale_color_viridis_c(option="C", direction=1, name="B-value", breaks=c(round(min(func_1kb$mean_B) + 0.01, 2), round(max(func_1kb$mean_B) - 0.01, 2))) +
  theme(strip.text=element_text(size=16),
        axis.title=element_text(size=20),
        axis.text=element_text(size=15),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14, margin=margin(b=20)),
        legend.position=c(0.15, 0.15),
        legend.direction="horizontal",
        legend.box="horizontal")

p3b <- ggplot(benegas_1kb, aes(x=class, y=med_ratio, shape=hasCpG)) +
  facet_wrap(~map, labeller=labeller(map=map_labels)) +
  theme_classic() + geom_hline(yintercept=1, linetype="dashed", color="grey") +
  geom_segment(data=seg_benegas,
               aes(x=x, xend=xend, y=y, yend=yend, color=mean_B, group=interaction(map, hasCpG)),
               linewidth=0.9, lineend="round", inherit.aes=FALSE) +
  geom_point(aes(color=mean_B), size=3) +
  geom_errorbar(aes(ymin=med_ratio - se_ratio, ymax=med_ratio + se_ratio, color=mean_B), width=0.2, linewidth=0.7) +
  labs(x=NULL, y=expression(paste(mu["÷"])), title=NULL) +
  scale_shape_manual(name=NULL, values=c("FALSE"=16, "TRUE"=17),
                     labels=c("TRUE"="With CpG", "FALSE"="Without CpG"), guide="none") +
  scale_y_continuous(breaks=pretty_breaks(), limits=c(min(benegas_1kb$med_ratio) - benegas_1kb$se_ratio[which.min(benegas_1kb$med_ratio)],
                                                      max(benegas_1kb$med_ratio) + benegas_1kb$se_ratio[which.min(benegas_1kb$med_ratio)])) +
  scale_color_viridis_c(option="C", direction=1, name="B-value", breaks=c(round(min(benegas_1kb$mean_B) + 0.01, 2), round(max(benegas_1kb$mean_B) - 0.01, 2))) +
  theme(strip.text=element_text(size=16),
        axis.title=element_text(size=20),
        axis.text=element_text(size=15),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14, margin=margin(b=20)),
        legend.position=c(0.15, 0.15),
        legend.direction="horizontal",
        legend.box="horizontal")

p3c <- ggplot(phast_1kb, aes(x=class, y=med_ratio, shape=hasCpG)) +
  facet_wrap(~map, labeller=labeller(map=map_labels)) +
  theme_classic() + geom_hline(yintercept=1, linetype="dashed", color="grey") +
  geom_segment(data=seg_phast,
               aes(x=x, xend=xend, y=y, yend=yend, color=mean_B, group=interaction(map, hasCpG)),
               linewidth=0.9, lineend="round", inherit.aes=FALSE) +
  geom_point(aes(color=mean_B), size=3) +
  geom_errorbar(aes(ymin=med_ratio - se_ratio, ymax=med_ratio + se_ratio, color=mean_B), width=0.2, linewidth=0.7) +
  labs(x="Constraint class", y=expression(paste(mu["÷"])), title=NULL) +
  scale_shape_manual(name=NULL, values=c("FALSE"=16, "TRUE"=17),
                     labels=c("TRUE"="With CpG", "FALSE"="Without CpG")) +
  scale_y_continuous(breaks=pretty_breaks(), limits=c(min(seg_phast$med_ratio) - seg_phast$se_ratio[which.min(seg_phast$med_ratio)],
                                                      max(seg_phast$med_ratio) + seg_phast$se_ratio[which.min(seg_phast$med_ratio)])) +
  scale_color_viridis_c(option="C", direction=1, name="B-value", 
                        breaks=c(round(min(seg_phast$mean_B) + 0.01, 2), round(max(seg_phast$mean_B) - 0.01, 2))) +
  theme(strip.text=element_text(size=16),
        axis.title=element_text(size=20),
        axis.text=element_text(size=15),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14, margin=margin(b=20)),
        legend.position=c(0.15, 0.1),
        legend.direction="horizontal",
        legend.box="horizontal")

# keeping color legends inside while bringing shape legend to bottom
p_base <- p3c + theme(legend.position="none")

p_color_legend <- p3c +
  guides(shape="none") +
  theme(legend.position="bottom",
        legend.direction="horizontal")

p_shape_legend <- p3c +
  scale_color_viridis_c(guide="none") +  
  theme(legend.position="bottom",
        legend.direction="horizontal")

g_color <- get_legend(p_color_legend)
g_shape <- get_legend(p_shape_legend)

p3c2 <- ggdraw() +
  draw_plot(p_base) +
  draw_grob(g_color, x=0.065, y=0.28, width=0.30, height=0.05) 

p3 <- plot_grid(p3a, p3b, p3c2, g_shape, labels=c("A", "B", "C"), label_size=16, ncol=1, rel_heights=c(1, 1, 1, 0.075))
save_plot("ratios_isolated_1kb_NEW.pdf", p3, base_height=11, base_width=12)
