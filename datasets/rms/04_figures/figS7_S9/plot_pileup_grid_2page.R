#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# FigS6 as TWO self-contained pages, split by method (Fig5B peak-distance rank). Each page = header (CRC-SE +
# gene bar + celltype titles) + MYOD1 | FOXO1 grid (6 cols), Annotation on top + 6 methods + Random on the
# bottom (8 rows) -- mirrors main Fig5C (5-method, 2-region), extended to all 12 methods over two pages.
#   Page 1 (top-6):    Annotation, scglue(multiome), scglue, simba, BindSC, Portal, MIDAS, Random
#   Page 2 (bottom-6): Annotation, scVI, MaxFuse, Seurat(CCA), scBridge, scButterfly, scJoint, Random
# Emits BOTH the grid-only pages (figS6_page*_grid.*, ready to drop under the reused Fig5C header) and the
# auto-assembled full pages with header (figS6_page*_full.*). Reuses plot_figS6.R label recipe.  Rscript plot_pileup_grid_2page.R
suppressMessages({library(dplyr); library(ggplot2); library(tidyr); library(ggrastr); library(cowplot)})
set.seed(123)
mypalette <- c("#1F77B4", "#FF7F0E", "#2CA02C")   # Mesoderm, Myoblast, Myocyte
CTS <- c("Mesoderm", "Myoblast", "Myocyte")

G   <- "/path/to/multiomeBench"
KNN <- file.path(G, "results/rms/knn_test")
NEW <- file.path(G, "RMS/benchmark/fig5b")
S4  <- file.path(G, "RMS/benchmark/figS4a")
OUT <- file.path(G, "RMS/benchmark/figS6")
NEW_METHODS <- c("MaxFuse", "MIDAS", "scButterfly")
RENDER_FULL <- FALSE   # grid-only this run (user assembles the header separately); flip to TRUE for full pages

remap <- function(x) dplyr::case_when(
  x %in% c("Seurat.CCA.", "Seurat(CCA)")           ~ "Seurat \n (CCA)",
  x %in% c("Seurat.WNN.", "Seurat(WNN)")           ~ "Seurat(WNN)",
  x %in% c("scglue.multiome.", "scglue(multiome)") ~ "scglue \n (multiome)",
  x == "Bindsc" ~ "BindSC", x == "annotation" ~ "Annotation", x == "random" ~ "Random", TRUE ~ x)

PAGES <- list(
  page1 = c("scglue \n (multiome)", "scglue", "simba", "BindSC", "Portal", "MIDAS"),
  page2 = c("scVI", "MaxFuse", "Seurat \n (CCA)", "scBridge", "scButterfly", "scJoint"))

tfs <- list(
  myod1 = list(value_old = file.path(KNN, "predicted_ataclabel_myod1value.csv"),
               value_new = file.path(NEW, "new_methods_myod1value.csv"),
               ymax = 10, textx = 17733000, texty = 8.6, xlab = "chr11 position (bp)",
               rects = list(c(17728863, 17736766)), breaks = seq(17725000, 17745000, 15000),
               gene = list(start = 17741115, end = 17743678, name = "MYOD1", col = "#E69F00")),
  foxo1 = list(value_old = file.path(KNN, "predicted_ataclabel_foxo1value.csv"),
               value_new = file.path(NEW, "new_methods_foxo1value.csv"),
               ymax = 8, textx = 41200000, texty = 6.9, xlab = "chr13 position (bp)",
               rects = list(c(41145735, 41171933), c(41224650, 41256947)),
               breaks = seq(41100000, 41300000, 150000),
               gene = list(start = 41129804, end = 41240734, name = "FOXO1", col = "#0072B2")))

## ---------- accuracy + n labels (identical recipe to plot_figS6.R) ----------
KNN.label <- read.csv(file.path(KNN, "knn_k10_pred_label.csv"), row.names = 1)
rownames(KNN.label) <- gsub("_atac", "", rownames(KNN.label)); bc <- rownames(KNN.label)
atac.annot <- read.csv(file.path(KNN, "..", "label.csv"), row.names = 1)
colnames(atac.annot)[colnames(atac.annot) == "cell_type"] <- "annotation"
rownames(atac.annot) <- gsub("_atac", "", rownames(atac.annot)); atac.annot <- atac.annot[bc, ]
atac.annot$random <- sample(unique(atac.annot$annotation), nrow(atac.annot), replace = TRUE)
atac.annot <- atac.annot[bc, c("annotation", "random")]; type <- unique(atac.annot$annotation)
accu_random <- do.call(rbind, lapply(type, function(cat) {
  a <- atac.annot$annotation == cat; p <- atac.annot$random == cat
  data.frame(pipeline = "Random", type = cat, accuracy = sum(a & p)/sum(a), n = sum(p))}))
count_old <- bind_rows(lapply(colnames(KNN.label), function(col) {
  d <- as.data.frame(table(KNN.label[[col]])); colnames(d) <- c("type", "n"); d$pipeline <- col; d}))
count_old$pipeline <- remap(as.character(count_old$pipeline)); count_old$type <- as.character(count_old$type)
accu_old <- read.csv(file.path(KNN, "knn_k10_pred_accu.csv"), row.names = 1)
accu_old$pipeline <- remap(accu_old$pipeline)
accu_old <- accu_old %>% filter(type != "overall") %>% merge(count_old, by = c("pipeline", "type"))
accu_new <- read.csv(file.path(S4, "knn_pred_accu.csv")) %>%
  filter(method %in% NEW_METHODS, type != "overall") %>% rename(pipeline = method)
count_new <- bind_rows(lapply(NEW_METHODS, function(m) {
  lab <- read.csv(file.path(S4, paste0("knn_pred_label__", m, ".csv")), row.names = 1)
  d <- as.data.frame(table(lab[[1]])); colnames(d) <- c("type", "n"); d$pipeline <- m; d$type <- as.character(d$type); d}))
accu_new <- merge(accu_new, count_new, by = c("pipeline", "type"))
accu_annot <- atac.annot %>% group_by(annotation) %>% summarise(n = n(), .groups = "drop") %>%
  transmute(pipeline = "Annotation", type = annotation, accuracy = 1, n = n)
KNN.accu <- bind_rows(accu_old[, c("pipeline","type","accuracy","n")], accu_new[, c("pipeline","type","accuracy","n")],
                      accu_random[, c("pipeline","type","accuracy","n")], accu_annot) %>%
  mutate(label = paste0("accuracy:", round(accuracy, 2), ",n=", n)) %>% drop_na()
KNN.accu$label <- ifelse(KNN.accu$pipeline == "Annotation", paste0("n=", KNN.accu$n), KNN.accu$label)

## ---------- read each region's tracks ONCE, cache to RDS (foxo1 CSV = 88MB) ----------
rds <- file.path(OUT, "region_long.rds")
if (file.exists(rds)) {
  cat("loading cached region tracks (region_long.rds)\n"); region_long <- readRDS(rds)
} else {
  cat("reading region tracks (foxo1 is 88MB, ~1-2 min)...\n")
  region_long <- lapply(tfs, function(cfg) {
    old <- read.csv(cfg$value_old); new <- read.csv(cfg$value_new) %>% filter(pipeline %in% NEW_METHODS)
    bind_rows(old, new) %>% pivot_longer(cols = all_of(CTS), names_to = "type") %>%
      mutate(pipeline = remap(pipeline)) %>% drop_na()
  })
  saveRDS(region_long, rds)
}

# downsample positions for rendering: ~4000 points/track is visually lossless over a ~4in panel (foxo1 has
# ~100k) and cuts the geom_bar cost ~20x so the cowplot-aligned full pages render in reasonable time.
region_long <- setNames(lapply(names(region_long), function(tf) {
  d <- region_long[[tf]]; up <- sort(unique(d$position)); N <- max(1L, floor(length(up) / 4000))
  d[d$position %in% up[seq(1, length(up), by = N)], ]
}), names(region_long))
cat(sprintf("downsampled: myod1=%d, foxo1=%d unique positions\n",
            length(unique(region_long$myod1$position)), length(unique(region_long$foxo1$position))))

## ---------- one region grid, filtered to a page's methods ----------
build_region <- function(tf, methods_page, show_y, show_methodlab) {
  cfg <- tfs[[tf]]; lv <- c("Annotation", methods_page, "Random")
  po <- region_long[[tf]] %>% filter(pipeline %in% lv) %>% mutate(pipeline = factor(pipeline, levels = lv))
  ka <- KNN.accu %>% filter(pipeline %in% lv) %>% mutate(pipeline = factor(pipeline, levels = lv))
  ka$xpos <- cfg$textx; ka$ypos <- cfg$texty
  p <- ggplot(po, aes(position, value, fill = type, color = type)) +
    rasterise(geom_bar(stat = "identity"), dpi = 300) +
    facet_grid(pipeline ~ type) +   # keep the 2-line names so vertical labels fit the row height
    geom_text(data = ka, aes(x = xpos, y = ypos, label = label), size = 3.7) +
    coord_cartesian(ylim = c(0, cfg$ymax)) + xlab(cfg$xlab) + theme_classic() +
    scale_color_manual(values = mypalette) + scale_fill_manual(values = mypalette) +
    theme(legend.position = "none", panel.spacing.x = unit(0.05, "lines"),
          panel.grid = element_blank(), strip.background = element_blank(), strip.text.x = element_blank(),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          text = element_text(size = 13, color = "black"),
          axis.text.x = element_text(size = 11, color = "black"))
  for (r in cfg$rects) p <- p + annotate("rect", xmin = r[1], xmax = r[2], ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "blue")
  p <- p + scale_x_continuous(breaks = cfg$breaks)
  p <- p + if (show_methodlab) theme(strip.text.y = element_text(size = 12, angle = -90)) else theme(strip.text.y = element_blank())
  if (show_y) p <- p + ylab("value") + theme(axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 11))
  else p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(), axis.line.y = element_blank())
  p
}

## ---------- region header: CRC-SE bar + gene model + celltype titles (same facet structure -> aligns) ----------
build_header <- function(tf, show_y, show_methodlab) {
  cfg <- tfs[[tf]]; g <- cfg$gene; xr <- range(region_long[[tf]]$position)
  se_y <- 1.55; se_h <- 0.28; gn_y <- 0.85; gn_h <- 0.20; nm_y <- 2.25   # CRC-SE row, gene row, gene-name label
  # every layer carries type+pipeline so each facet panel has explicit data (mixing has/hasn't breaks the gtable)
  se  <- do.call(rbind, lapply(cfg$rects, function(r) data.frame(xmin = r[1], xmax = r[2], type = CTS))); se$pipeline <- ""
  gn  <- data.frame(xmin = g$start, xmax = g$end, type = CTS); gn$pipeline <- ""
  sg  <- data.frame(x = g$start, xend = g$end, y = gn_y, yend = gn_y, type = CTS); sg$pipeline <- ""
  lbl <- data.frame(type = CTS[1], pipeline = "")
  p <- ggplot() +
    geom_rect(data = se, aes(xmin = xmin, xmax = xmax, ymin = se_y - se_h, ymax = se_y + se_h), fill = "black") +
    geom_segment(data = sg, aes(x = x, xend = xend, y = y, yend = yend), color = g$col, linewidth = 0.5) +
    geom_rect(data = gn, aes(xmin = xmin, xmax = xmax, ymin = gn_y - gn_h, ymax = gn_y + gn_h),
              fill = g$col, color = g$col) +
    geom_text(data = lbl, aes(x = xr[1], y = nm_y, label = g$name), hjust = 0, vjust = 1,
              size = 3.4, fontface = "bold", color = g$col) +   # region gene name, top-left
    facet_grid(pipeline ~ type) +
    coord_cartesian(ylim = c(0.4, 2.6)) + scale_x_continuous(breaks = cfg$breaks, limits = xr) +
    theme_classic() +
    theme(legend.position = "none", panel.spacing.x = unit(0.05, "lines"), panel.grid = element_blank(),
          strip.background = element_blank(), strip.text.x = element_text(size = 12, face = "bold"),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          text = element_text(size = 11, color = "black"), plot.margin = margin(2, 2, 0, 2),
          axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), axis.line.x = element_blank())
  p <- p + if (show_methodlab) theme(strip.text.y = element_text(size = 10, angle = 0, color = "transparent"))
           else theme(strip.text.y = element_blank())
  # "CRC-SE" / "Gene" row labels on the far-left axis (myod1 only) apply to both regions (tracks share rows)
  if (show_y) p <- p + scale_y_continuous(breaks = c(gn_y, se_y), labels = c("Gene", "CRC-SE")) +
                       ylab(NULL) + theme(axis.text.y = element_text(size = 9, color = "black"), axis.ticks.y = element_blank())
  else p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(), axis.line.y = element_blank())
  p
}

## ---------- render both pages ----------
for (pg in names(PAGES)) {
  p_my <- build_region("myod1", PAGES[[pg]], show_y = TRUE,  show_methodlab = FALSE)
  p_fo <- build_region("foxo1", PAGES[[pg]], show_y = FALSE, show_methodlab = TRUE)
  grid_row <- cowplot::plot_grid(p_my, p_fo, ncol = 2, rel_widths = c(1, 1), align = "h", axis = "tb")
  ggsave(file.path(OUT, paste0("figS6_", pg, "_grid.pdf")), grid_row, width = 13, height = 9.0)
  ggsave(file.path(OUT, paste0("figS6_", pg, "_grid.png")), grid_row, width = 13, height = 9.0, dpi = 150)

  if (RENDER_FULL) {
    h_my <- build_header("myod1", show_y = TRUE,  show_methodlab = FALSE)
    h_fo <- build_header("foxo1", show_y = FALSE, show_methodlab = TRUE)
    full <- cowplot::plot_grid(h_my, h_fo, p_my, p_fo, ncol = 2, rel_heights = c(1.0, 8),
                               align = "hv", axis = "tblr")
    ggsave(file.path(OUT, paste0("figS6_", pg, "_full.pdf")), full, width = 13, height = 9.2)
    ggsave(file.path(OUT, paste0("figS6_", pg, "_full.png")), full, width = 13, height = 9.2, dpi = 150)
  }
  cat(sprintf("wrote %s (grid%s)\n", pg, if (RENDER_FULL) " + full" else ""))
}
cat("done\n")
