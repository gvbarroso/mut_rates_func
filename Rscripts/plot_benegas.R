

dt_m <- unique_vals %>%
  pivot_longer(cols=c(mean_roulette, mean_carlson, mean_gnomad, se_roulette, se_carlson, se_gnomad),
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
save_plot("~/Desktop/mut_rates/benegas_mut_rates.pdf", p1, base_height=6, base_width=10)