library(ggplot2)
library(cowplot)
library(grid)
library(ggpubr)

variants <- data.frame(
  CHR = 1,
  POS = c(10221:10230, 11215:11220, 980772:980775),
  TRIPLET = c("AAA","AAC","ACG","CGT","GTT",
              "TTA","TAT","ATC","TCG","CGG",
              "GGT","GTA","TAG","AGC","GCC",
              "CCT","CTA","TAA","AAG","AGC"),
  CLASS = c(rep("2", 5), rep("Neu", 5), rep("7", 3), rep("8", 3), rep("Neu", 4)),
  ROULETTE = round(runif(20, 1, 15), 1) * 1e-9,
  GNOMAD = round(runif(20, 1, 15), 1) * 1e-9,
  CARLSON = round(runif(20, 1, 15), 1) * 1e-9)

theme_tight <- gridExtra::ttheme_minimal(
  core = list(
    fg_params = list(cex = 1.4),          # larger text
    padding.h = unit(0.8, "mm"),          # reduce horizontal padding
    padding.v = unit(0.5, "mm")           # reduce vertical padding
  ),
  colhead = list(
    fg_params = list(cex = 1.5, fontface = 2),
    padding.h = unit(0.8, "mm"),
    padding.v = unit(0.5, "mm")
  )
)

tab <- gridExtra::tableGrob(
  variants[, c("CHR","POS","TRIPLET","CLASS","ROULETTE","GNOMAD","CARLSON")],
  rows = NULL,
  theme = theme_tight
)

left_panel <- ggpubr::as_ggplot(tab)
ggsave("scheme/left_panel.pdf", left_panel, width=8.5, height=7.25) # loaded in Inkscape

