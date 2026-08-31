#!/usr/bin/env Rscript
## plot_fig3a.R -- Fig3A: per-cell-type metric bubble plot for pbmc3k.
##
## Clean, repo-relative extraction of the Fig3A section (bubble plot) of
## benchmark/plot_cellttype_stat2.R. Four facets -- KS intra-omics cell-type distance, Omics ASW,
## KS inter-omics cell-type distance, adjusted ATAC label-prediction accuracy -- each a dot per
## (method, cell type); bubble size = #cells of that type, colour = cell type.
##
## Differences from the original:
##   * reads the CLEAN fig2b/ metric CSVs (so the 3 added methods MIDAS/MaxFuse/scButterfly are in)
##     instead of the old /path/to/home/... benchmark_matrix CSVs
##   * method order is taken from fig2b_matrix.csv Rank (restricted to the 14 Fig3A methods) rather
##     than a hand-maintained list, so it stays in sync with the Fig2B ranking
##   * drops the Fig2B formattable table that was appended in the original (that lives in fig2b/plot_fig2b.R)
##
## Usage:  Rscript plot_fig3a.R [metrics_dir=../fig2b] [outdir=.]
suppressMessages({library(tidyr); library(dplyr); library(ggplot2)})

args <- commandArgs(trailingOnly = TRUE)
MET  <- ifelse(length(args) >= 1, args[1], "../fig2b")   # holds celltype_metrics.csv, adj_atac_predaccu.csv, label.csv, fig2b_matrix.csv
OUT  <- ifelse(length(args) >= 2, args[2], ".")

tab10_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

celltype_metrics <- read.csv(file.path(MET, "celltype_metrics.csv"))              # method, celltype, ks_celltype, ks_inter_celltype, celltype_omics_ASW
atac_accu        <- read.csv(file.path(MET, "adj_atac_predaccu.csv"))             # method, <7 celltype accu cols>, average_accu
celltype_label   <- read.csv(file.path(MET, "label.csv"))                         # Barcode, cell_type, modality, ...

## bubble size = #cells per (test) cell type
celltype_label %>%
  filter(modality != "train multiomics") %>%
  group_by(cell_type) %>%
  summarise(count = n(), .groups = "drop") -> celltype_count

## method order = Fig2B rank, restricted to the 14 Fig3A methods (original 11 + MIDAS/MaxFuse/scButterfly)
KEEP <- c("scglue(multiome)", "scDART", "scglue", "scVI", "Seurat(CCA)", "Portal", "scJoint",
          "simba", "BindSC", "scBridge", "Cobolt", "MIDAS", "MaxFuse", "scButterfly")
rankdf <- read.csv(file.path(MET, "fig2b_matrix.csv"), check.names = FALSE)
rankdf$method <- gsub('^"|"$', "", rankdf$method)
method.rank <- rankdf %>% filter(method %in% KEEP) %>% arrange(Rank) %>% pull(method)
stopifnot(length(method.rank) == length(KEEP))

## assemble the long table: value per (method, cell type, statistic)
atac_accu %>%
  pivot_longer(cols = `CD14.Monocytes`:`Naive.CD8.T.cells`,
               names_to = "celltype", values_to = "knn pred accu") %>%
  mutate(celltype = case_when(celltype == "CD14.Monocytes" ~ "CD14 Monocytes",
                              celltype == "CD16.Monocytes" ~ "CD16 Monocytes",
                              celltype == "B.cells" ~ "B cells",
                              celltype == "Memory.T.cells" ~ "Memory T cells",
                              celltype == "NK.Effector.T.cells" ~ "NK/Effector T cells",
                              celltype == "Naive.CD4.T.cells" ~ "Naive CD4 T cells",
                              celltype == "Naive.CD8.T.cells" ~ "Naive CD8 T cells",
                              TRUE ~ celltype)) %>%
  merge(celltype_metrics, by = c("method", "celltype")) %>%
  select("method", "celltype", "ks_celltype", "celltype_omics_ASW", "ks_inter_celltype", "knn pred accu") %>%
  merge(celltype_count, by.x = "celltype", by.y = "cell_type") %>%
  pivot_longer(cols = c("ks_celltype", "celltype_omics_ASW", "ks_inter_celltype", "knn pred accu"),
               names_to = "statistics", values_to = "value") %>%
  mutate(statistics = case_when(statistics == "ks_celltype" ~ "KS statistics of intra-omics \n cell type distances",
                                statistics == "celltype_omics_ASW" ~ "Omics ASW",
                                statistics == "ks_inter_celltype" ~ "KS statistics of inter-omics \n cell type distance",
                                statistics == "knn pred accu" ~ "Adjusted ATAC cell label \n prediction accuracy")) %>%
  filter(method %in% method.rank) -> celltype

celltype$statistics <- factor(celltype$statistics,
  levels = c("KS statistics of intra-omics \n cell type distances", "Omics ASW",
             "KS statistics of inter-omics \n cell type distance", "Adjusted ATAC cell label \n prediction accuracy"))
celltype$Methods <- factor(celltype$method, levels = rev(method.rank))            # rank 1 at top
celltype$`Cell type` <- factor(celltype$celltype,
  levels = c("CD14 Monocytes", "B cells", "Memory T cells", "NK/Effector T cells",
             "CD16 Monocytes", "Naive CD4 T cells", "Naive CD8 T cells"))

p <- ggplot(celltype, aes(x = value, y = Methods, size = count, color = `Cell type`)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~statistics, ncol = 4) +
  scale_color_manual(values = tab10_colors) +
  theme_bw() +
  theme(panel.grid.major = element_line(linetype = "dashed", color = "gray", size = 0.1),
        axis.text.x = element_text(angle = 45, hjust = 0.1, vjust = 0),
        panel.background = element_blank()) +
  theme(text = element_text(size = 14, color = "black"),
        plot.title = element_text(size = 14, color = "black", face = "bold"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = NA, colour = NA)) +
  xlab("") + ylab("") +
  theme(legend.position = "right", legend.text = element_text(size = 14),
        legend.box.spacing = unit(11, "pt"),    # clear gap so the legend doesn't touch the last panel
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.y = unit(1, "pt")) +
  guides(size = guide_legend(ncol = 1),
         colour = guide_legend(ncol = 1, override.aes = list(size = 5)))   # legends stacked on the right

ggsave(file.path(OUT, "celltype_stat.pdf"), p, width = 15.5, height = 4.0)
ggsave(file.path(OUT, "celltype_stat.png"), p, width = 15.5, height = 4.0, dpi = 170)
cat("wrote", file.path(OUT, "celltype_stat.pdf"), "+ .png  --", length(method.rank), "methods:\n  ",
    paste(rev(levels(celltype$Methods)), collapse = ", "), "\n")
