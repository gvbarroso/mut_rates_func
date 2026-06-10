library(ggplot2)
library(cowplot)
library(grid)
library(ggpubr)

variants <- data.frame(
  CHROM   = "chr1",
  POS     = c(10120, 10220, 10320, 10420, 10520,
              10620, 10720, 10820, 10920, 11020,
              11120, 11220, 11320, 11420, 11520),
  TRIPLET = c("AAA","AAC","AAG","ACA","ACC",
              "ACG","AGA","AGC","AGG","ATA",
              "ATC","ATG","CAA","CAC","CAG"),
  GROUP   = c("G1","G2","G3","G4","G5",
              "G6","G7","G8","G9","G10",
              "G11","G12","G1","G2","G3"),
  RATE    = round(runif(15, 0.5, 1.5), 3)
)

# Flag rows in the 1 kb window [10320, 11320]
variants$in_window <- with(variants, POS >= 10320 & POS <= 11320)

# Left panel
tab <- gridExtra::tableGrob(
  variants[, c("CHROM","POS","TRIPLET","GROUP","RATE")],
  rows = NULL,
  theme = gridExtra::ttheme_minimal(
    core = list(
      fg_params = list(cex = 0.7),
      bg_params = list(fill = NA, col = NA)
    ),
    colhead = list(
      fg_params = list(cex = 0.7, fontface = 2),
      bg_params = list(fill = NA, col = NA)
    )
  )
)

# Compute which rows are in the window (for bracket placement)
window_rows <- which(variants$in_window)
row_top    <- min(window_rows)
row_bottom <- max(window_rows)

# Function to draw a thin minimal bracket on the right side of the table
add_bracket <- function(table_grob, row_top, row_bottom, label = "1 kb") {
  g <- table_grob
  g <- gTree(children = gList(g), cl = "bracket_table")
  
  bracket_fun <- function(x) {
    grid.draw(x$children[[1]])
    
    # Get layout info
    lay <- x$children[[1]]$layout
    n_rows <- max(lay$t)
    
    # Right edge x position (npc)
    x_right <- unit(1, "npc") + unit(0.01, "npc")
    
    # Convert row indices to npc y positions
    y_top    <- unit(1 - (row_top - 0.5) / n_rows, "npc")
    y_bottom <- unit(1 - (row_bottom - 0.5) / n_rows, "npc")
    
    # Thin bracket lines
    grid.lines(x = unit.c(x_right, x_right),
               y = unit.c(y_top, y_bottom),
               gp = gpar(lwd = 0.5))
    
    grid.lines(x = unit.c(x_right - unit(0.01, "npc"), x_right),
               y = unit.c(y_top, y_top),
               gp = gpar(lwd = 0.5))
    
    grid.lines(x = unit.c(x_right - unit(0.01, "npc"), x_right),
               y = unit.c(y_bottom, y_bottom),
               gp = gpar(lwd = 0.5))
    
    # Label
    grid.text(label,
              x = x_right + unit(0.02, "npc"),
              y = (y_top + y_bottom) / 2,
              just = "left",
              gp = gpar(cex = 0.7))
  }
  
  class(bracket_fun) <- "grob"
  g$draw <- bracket_fun
  g
}

tab_with_bracket <- add_bracket(tab, row_top, row_bottom, label = "1 kb")

## -----------------------------
## 3. Right panel: 12 group ratios
## -----------------------------
set.seed(1)
ratios <- data.frame(
  group = factor(paste0("Group ", 1:12),
                 levels = paste0("Group ", 12:1)), # reverse for top-down
  ratio = runif(12, 0.85, 0.95)
)

p_right <- ggplot(ratios, aes(x = ratio, y = group)) +
  geom_segment(aes(x = 0.85, xend = 0.95, y = group, yend = group),
               color = "grey80", linewidth = 0.3) +
  geom_point(size = 1.2) +
  scale_x_continuous(limits = c(0.84, 0.96), breaks = c(0.85, 0.90, 0.95)) +
  labs(x = "Group / neutral ratio", y = NULL) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 5)
  )

## -----------------------------
## 4. Combine panels + arrow
## -----------------------------
# Left: table grob as a cowplot object
left_panel  <- ggpubr::as_ggplot(tab_with_bracket)
right_panel <- p_right

# Base two-panel layout
two_panel <- plot_grid(
  left_panel,
  right_panel,
  nrow = 1,
  rel_widths = c(1.2, 1),
  align = "h"
)

# Add a short arrow between panels using draw_plot + annotation_custom
arrow_layer <- ggplot() +
  xlim(0, 1) + ylim(0, 1) +
  theme_void()

arrow_layer <- arrow_layer +
  annotation_custom(
    grob = linesGrob(
      x = unit(c(0.48, 0.52), "npc"),
      y = unit(c(0.5, 0.5), "npc"),
      gp = gpar(lwd = 0.5),
      arrow = arrow(length = unit(0.02, "npc"), type = "closed")
    )
  )

final_plot <- ggdraw() +
  draw_plot(two_panel, 0, 0, 1, 1) +
  draw_plot(arrow_layer, 0, 0, 1, 1)

## -----------------------------
## 5. Save as crisp vector PDF
## -----------------------------
#ggsave("two_panel_schematic.pdf", final_plot,
#       width = 7, height = 3, units = "in", device = cairo_pdf)
