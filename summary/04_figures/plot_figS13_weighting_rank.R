#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README.
## Weighting-sensitivity, WEIGHT-BY-RANK version (manager revision, Aug 2026):
##  1. restrict to 50-90% performance weight (0% ignores performance, 100% ignores usability -> degenerate).
##  2. rank methods on performance and on usability SEPARATELY (each 1..12), then weight the RANKS:
##        weighted_rank(w) = w * perf_rank + (1 - w) * usab_rank        (lower = better)
##     This puts both axes on the same 1..12 scale, unlike weighting the raw (differently-scaled) scores.
##  3. weighted_rank is LINEAR in w, so it renders as a smooth line; a 10% grid lands exactly on it.
## Component ranks come from the sweep-CSV endpoints (score_w(1)=perf, score_w(0)=usab):
##   col "1" = rank(-perf_mean) = perf_rank ; col "0" = rank(-usab) = usab_rank.
suppressMessages({library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)})
SUM <- "/path/to/multiomeBench/sum_plot"
OUT <- Sys.getenv("WEIGHTING_OUTPUT_DIR", unset = SUM)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

ranks     <- read.csv(file.path(SUM, "weighting_sensitivity_ranks.csv"), check.names = FALSE)
perf_rank <- setNames(ranks[["1"]], ranks$method)
usab_rank <- setNames(ranks[["0"]], ranks$method)

ws   <- seq(0.5, 0.9, by = 0.10)                                     # 50-90%, 10% grid
long <- do.call(rbind, lapply(ws, function(w)
  data.frame(method = ranks$method, w = w,
             wrank = w * perf_rank + (1 - w) * usab_rank, row.names = NULL)))

primary      <- long %>% filter(abs(w - 0.8) < 1e-9)
right_labels <- long %>% filter(abs(w - 0.9) < 1e-9)

method_colors <- c('scglue' = "#1f77b4", 'scglue(multiome)' = "#aec7e8", 'Portal' = "#ffbb78",
                   'scBridge' = "#2ca02c", 'simba' = "#98df8a", 'Seurat(CCA)' = "#d62728",
                   'scJoint' = "#ff9896", 'scVI' = "#9467bd", 'BindSC' = "#e377c2",
                   'MaxFuse' = "#8c564b", 'MIDAS' = "#17becf", 'scButterfly' = "#bcbd22")

p <- ggplot(long, aes(x = w, y = wrank, group = method, color = method)) +
  geom_line(linewidth = 1.1) +
  geom_point(data = long, size = 2.6, stroke = 0) +          # marker at every 10% grid weight (50..90)
  geom_point(data = primary, size = 4.2, stroke = 0) +       # emphasize the 80% primary setting
  geom_point(data = primary, size = 4.2, shape = 21, stroke = 0.8, color = "grey25", fill = NA) +
  geom_vline(xintercept = 0.8, linetype = "dashed", linewidth = 0.7, color = "grey45") +
  annotate("text", x = 0.8, y = 1.1, label = "Primary analysis (80% performance / 20% usability)",
           hjust = 0.5, size = 5.8, color = "grey30") +
  geom_text_repel(data = right_labels, aes(label = method),
                  hjust = 0, direction = "y", size = 6.2, xlim = c(0.915, 1.07),
                  segment.size = 0.3, segment.color = "grey70", box.padding = 0.14, seed = 1) +
  scale_color_manual(values = method_colors) +
  scale_y_reverse(breaks = 1:12, limits = c(12.3, 0.9)) +
  scale_x_continuous(breaks = seq(0.5, 0.9, 0.1), labels = paste0(seq(50, 90, 10), "%"),
                     limits = c(0.5, 1.07), expand = expansion(mult = c(0.03, 0))) +
  labs(x = "Weight on integration performance", y = "Weighted rank (lower = better)") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey88", linewidth = 0.4),
        axis.line.x = element_line(color = "grey25", linewidth = 0.6),
        axis.line.y = element_line(color = "grey25", linewidth = 0.6),
        axis.title.x = element_text(size = 19, margin = margin(t = 8)),
        axis.title.y = element_text(size = 19, margin = margin(r = 8)),
        axis.text = element_text(size = 16, color = "grey20"),
        plot.margin = margin(12, 14, 12, 12)) +
  coord_cartesian(clip = "off")

ggsave(file.path(OUT, "weighting_sensitivity_revised.pdf"), p, width = 11, height = 8, bg = "white")
ggsave(file.path(OUT, "weighting_sensitivity_revised.png"), p, width = 11, height = 8, dpi = 200, bg = "white")

## final integer ranking table (for the text / response document)
fr <- sapply(ws, function(w) rank(round(w * perf_rank + (1 - w) * usab_rank, 6), ties.method = "min"))
colnames(fr) <- paste0(ws * 100, "%"); rownames(fr) <- ranks$method
cat("weight-by-rank final ranking (1 = best):\n"); print(fr)
cat("\ntop method per weight:\n")
print(apply(fr, 2, function(col) paste(rownames(fr)[col == min(col)], collapse = ", ")))
cat("\nWrote weighting_sensitivity_revised.pdf / .png (weight-by-rank, 50-90%, square)\n")
