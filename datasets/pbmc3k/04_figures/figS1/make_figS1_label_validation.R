#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# Supplementary "Fig. S1. Validation of cell-type labels" — the figure for R4 / R1.6.
# Supersedes make_figSX_label_validation.R (which built A+B only); this adds panel C.
#
#   A = pbmc3k  — canonical markers by 10X-annotated label   (CD14/CD16 + naive CD4/CD8 separability)
#   B = RMS     — myogenesis markers by annotated sub-stage    (ordered myogenic axis)
#   C = pbmc3k  — 10X-annotated vs independent RNA-only SingleR concordance (90.5%); disagreement sits in
#                 the close pairs, which is R1.6's own point made visible.
#
# NOT a panel here: the rank-stability result (ρ). It answers R4 suggestion (b) directly, but re-scoring the
# 23 method latents against the corrected SingleR labels needs the cluster (recompute_sub.sh) — the latents
# are not mounted locally. It lives in the response TEXT for now; add as panel D after that recompute.
#
# Panels are taken from their source scripts (sys.source → grab the ggplot object) so this figure can never
# drift from the standalone panels. Side effect: sourcing re-renders those standalone outputs too (in sync).
# patchwork 1.1.2 is broken with ggplot2 3.5.1 here → cowplot.
#   Rscript make_figS1_label_validation.R          # ~5 min (loads two h5 files)
suppressMessages({library(ggplot2); library(cowplot)})

ROOT <- if (dir.exists("/path/to/project")) "/path/to/multiomeBench" else
                                          "/path/to/multiomeBench"
A_SCRIPT <- file.path(ROOT, "pbmc/pbmc3k/benchmark/label_sensitivity/marker_plot_label.R")
B_SCRIPT <- file.path(ROOT, "RMS/benchmark/label_sensitivity/marker_plot_label.R")
C_SCRIPT <- file.path(ROOT, "pbmc/pbmc3k/benchmark/label_sensitivity/remap_labels.R")   # exposes p1 = confusion
OUT      <- file.path(ROOT, "pbmc/pbmc3k/benchmark/label_sensitivity")

cat("=== panel A: pbmc3k markers ===\n"); eA <- new.env(); sys.source(A_SCRIPT, envir = eA)
cat("=== panel B: RMS markers ===\n");    eB <- new.env(); sys.source(B_SCRIPT, envir = eB)
cat("=== panel C: 10X vs SingleR concordance ===\n"); eC <- new.env(); sys.source(C_SCRIPT, envir = eC)
stopifnot(inherits(eA$p, "ggplot"), inherits(eB$p, "ggplot"), inherits(eC$p1, "ggplot"))

# A over B (both full-width dot plots); C centred below at reduced width to keep its aspect ratio.
ab   <- plot_grid(eA$p, eB$p, ncol = 1, labels = c("A", "B"), label_size = 18, rel_heights = c(4.4, 3.6))
crow <- plot_grid(NULL, plot_grid(eC$p1, labels = "C", label_size = 18), NULL,
                  ncol = 3, rel_widths = c(0.17, 0.66, 0.17))
fig  <- plot_grid(ab, crow, ncol = 1, rel_heights = c(8.0, 5.2))

f <- file.path(OUT, "figS1_label_validation.pdf")
tryCatch({
  ggsave(f, fig, width = 12.2, height = 13.4, limitsize = FALSE)
  ggsave(sub("\\.pdf$", ".png", f), fig, width = 12.2, height = 13.4, dpi = 140, limitsize = FALSE)
  cat("\nwrote", f, "(+ .png)\n")
}, error = function(e) message("!! could not write ", basename(f), " — open in Preview? (", conditionMessage(e), ")"))
