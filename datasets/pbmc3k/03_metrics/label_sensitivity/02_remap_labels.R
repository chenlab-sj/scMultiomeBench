#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# Re-derive the SingleR 9-type labels from the CACHED per-cell fine calls, using the canonical mapping in
# fine_label_map.R. No SingleR / celldex / ExperimentHub / internet needed -- 00_singleR_annotate.R already
# wrote every fine call to singleR_per_cell.csv, and only the fine->coarse MAP changed. Runs on the Mac.
#
# Regenerates:  label_singleR.csv, singleR_per_cell.csv (coarse col), celltype_map.csv,
#               confusion_10x_vs_singleR.csv, confusion_heatmap.pdf/png, fine_label_map_heatmap.pdf/png
#   Rscript 02_remap_labels.R
suppressMessages({library(ggplot2); library(dplyr)})

ROOT <- if (dir.exists("/path/to/project")) "/path/to/multiomeBench" else
                                          "/path/to/multiomeBench"
OUT        <- file.path(ROOT, "pbmc/pbmc3k/benchmark/label_sensitivity")
ORIG_LABEL <- file.path(ROOT, "pbmc/pbmc3k/benchmark/fig2b/label.csv")   # 10X R&D ground-truth labels
source(file.path(OUT, "fine_label_map.R"))                                # -> MAP

ORDER <- c("CD14 Monocytes","CD16 Monocytes","Myeloid DC","Plasmacytoid DC","B cells",
           "Naive CD4 T cells","Naive CD8 T cells","Memory T cells","NK/Effector T cells")

# These CSVs get opened in Excel, which takes an exclusive lock on this SMB mount ("Resource temporarily
# unavailable"). Don't let a locked side-output halt the run -- warn and carry on; re-run to fill it in.
FAILED <- character(0)
safe_write <- function(x, f, ...) invisible(tryCatch({ write.csv(x, f, ...); TRUE },
  error = function(e) { FAILED <<- c(FAILED, basename(f))
    message("  !! could not write ", basename(f), " -- open in Excel? (skipped)"); FALSE }))
safe_save <- function(f, p, ...) invisible(tryCatch({ ggsave(f, p, ...); TRUE },
  error = function(e) { FAILED <<- c(FAILED, basename(f))
    message("  !! could not write ", basename(f), " -- open in Preview? (skipped)"); FALSE }))

# ---- 1) cached fine calls -> coarse via the canonical MAP ----
per  <- read.csv(file.path(OUT, "singleR_per_cell.csv"), stringsAsFactors = FALSE)
fine <- per$singleR_fine; names(fine) <- per$barcode
unmapped <- setdiff(unique(fine), names(MAP))
if (length(unmapped)) cat("WARNING unmapped fine labels (-> 'other'):\n  ", paste(unmapped, collapse = ", "), "\n")
mapped <- MAP[fine]; mapped[is.na(mapped)] <- "other"; names(mapped) <- names(fine)
cat("SingleR mapped-label distribution:\n"); print(table(mapped))

per$singleR_celltype <- unname(mapped[per$barcode])
safe_write(per, file.path(OUT, "singleR_per_cell.csv"), row.names = FALSE)
safe_write(data.frame(fine = names(MAP), benchmark = unname(MAP)),
           file.path(OUT, "celltype_map.csv"), row.names = FALSE)

# ---- 2) new label.csv: test rows get the SingleR label; train rows keep the 10X labels (train isn't scored) ----
lab     <- read.csv(ORIG_LABEL, stringsAsFactors = FALSE, check.names = FALSE)
cell_bc <- sub("_(rna|atac)$", "", gsub('"', "", lab$Barcode))
is_test <- lab$modality != "train multiomics"
lab_new <- lab
hit     <- is_test & cell_bc %in% names(mapped)
lab_new$cell_type[hit] <- mapped[cell_bc[hit]]
miss <- sum(is_test & !(cell_bc %in% names(mapped)))
if (miss) cat(sprintf("%d test rows had no SingleR call (kept original 10X label)\n", miss))
safe_write(lab_new, file.path(OUT, "label_singleR.csv"), row.names = FALSE)

# ---- 3) confusion + agreement on the RNA-side test cells ----
# Restrict to the cell types actually scored in the benchmark (= the types present in the 10X test labels).
# This drops the handful of cells SingleR calls Myeloid DC / "other" -- not benchmark cell types, and with no
# counterpart on the 10X (row) axis, so they only added an empty column.
rna_rows    <- lab$modality == "test scRNA"
truth       <- lab$cell_type[rna_rows]; names(truth) <- cell_bc[rna_rows]
bench_types <- sort(unique(truth))                       # benchmark-scored cell types (test split has no DC)
common      <- intersect(names(truth), names(mapped))
n_all       <- length(common)
common      <- common[truth[common] %in% bench_types & mapped[common] %in% bench_types]
n_bench     <- length(common)
# Additionally keep only SingleR-CONFIDENT cells (pruned.labels non-NA), per author request. Confidence is in
# singleR_confidence.csv (Rscript 01_singleR_confidence.R re-ran SingleR to get it; the cached max-corr score
# tops out at ~0.54 so an absolute >0.8 cutoff is not meaningful, and pruned.labels is SingleR's own
# high-confidence flag). If the file is absent, all benchmark-type cells are used.
conf_file <- file.path(OUT, "singleR_confidence.csv")
CONFIDENT_ONLY <- file.exists(conf_file)
if (CONFIDENT_ONLY) {
  ct        <- read.csv(conf_file, stringsAsFactors = FALSE)
  confident <- ct$barcode[ct$is_confident %in% c(TRUE, "TRUE", "True")]
  common    <- intersect(common, confident)
}
message(sprintf("confusion: %d test cells -> %d benchmark-type%s",
                n_all, n_bench,
                if (CONFIDENT_ONLY) sprintf(" -> %d SingleR-confident (pruned.labels non-NA)", length(common)) else ""))
conf <- table(TenX = truth[common], SingleR = mapped[common])
safe_write(as.data.frame.matrix(conf), file.path(OUT, "confusion_10x_vs_singleR.csv"))
agree <- mean(truth[common] == mapped[common])
cat(sprintf("\n=== agreement (10X-annotated ground truth vs SingleR, %d test cells): %.1f%% ===\n",
            length(common), 100*agree))
per_t <- sapply(sort(unique(truth[common])), function(t) {
  idx <- names(truth[common])[truth[common] == t]; 100*mean(mapped[idx] == t) })
for (t in names(per_t)) cat(sprintf("  %-22s %5.1f%% (n=%d)\n", t, per_t[t], sum(truth[common] == t)))

# ---- 4) HEATMAP 1 (supp fig): 10X x SingleR confusion, row-normalised % ----
lv  <- ORDER[ORDER %in% union(rownames(conf), colnames(conf))]
cdf <- as.data.frame(conf, stringsAsFactors = FALSE) %>%
  group_by(TenX) %>% mutate(pct = 100*Freq/sum(Freq)) %>% ungroup() %>%
  mutate(TenX = factor(TenX, levels = rev(lv)),
         SingleR = factor(SingleR, levels = c(lv, setdiff(unique(SingleR), lv))))
p1 <- ggplot(cdf, aes(SingleR, TenX, fill = pct)) +
  geom_tile(colour = "grey90", linewidth = 0.3) +
  geom_text(aes(label = ifelse(Freq > 0, Freq, "")), size = 3,
            colour = ifelse(cdf$pct > 55, "white", "grey20")) +
  scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0,100), name = "% of\n10X row") +
  labs(x = "SingleR label (independent, RNA-only)", y = "10X-annotated label",
       title = "pbmc3k label concordance") +   # %-concordance + n now reported in the caption, not the title
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title  = element_text(size = 14),
        plot.title  = element_text(size = 15),
        panel.grid  = element_blank())
safe_save(file.path(OUT, "confusion_heatmap.pdf"), p1, width = 8.4, height = 6.1)
safe_save(file.path(OUT, "confusion_heatmap.png"), p1, width = 8.4, height = 6.1, dpi = 150)

# ---- 5) HEATMAP 2 (mapping diagnostic): SingleR FINE x 10X, row-normalised %; row label shows the
#         bucket each fine label was assigned to -- this is what justifies the gd-T / effector-CD4 calls ----
fdf <- as.data.frame(table(fine = fine[common], TenX = truth[common]), stringsAsFactors = FALSE) %>%
  group_by(fine) %>% mutate(pct = 100*Freq/sum(Freq), n = sum(Freq)) %>% ungroup() %>%
  filter(n > 0) %>%
  mutate(bucket = ifelse(is.na(MAP[fine]), "other", MAP[fine]),
         rowlab = sprintf("%s  [%s, n=%d]", fine, bucket, n),
         TenX = factor(TenX, levels = lv))
ord <- fdf %>% distinct(rowlab, bucket, n) %>% arrange(bucket, -n) %>% pull(rowlab)
fdf$rowlab <- factor(fdf$rowlab, levels = rev(ord))
p2 <- ggplot(fdf, aes(TenX, rowlab, fill = pct)) +
  geom_tile(colour = "grey90", linewidth = 0.3) +
  geom_text(aes(label = ifelse(Freq > 0, Freq, "")), size = 2.7,
            colour = ifelse(fdf$pct > 55, "white", "grey20")) +
  scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0,100), name = "% of\nfine row") +
  labs(x = "10X-annotated (ground-truth) label", y = "SingleR fine label  [assigned bucket, n]",
       title = "Fine-label harmonization: where each MonacoImmune label lands") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
        plot.title = element_text(size = 11))
safe_save(file.path(OUT, "fine_label_map_heatmap.pdf"), p2, width = 9.6, height = 7.4)
safe_save(file.path(OUT, "fine_label_map_heatmap.png"), p2, width = 9.6, height = 7.4, dpi = 150)

if (length(FAILED)) {
  cat("\n!! SKIPPED (locked -- close in Excel/Preview and re-run):\n   ",
      paste(unique(FAILED), collapse = ", "), "\n")
} else cat("\nwrote all outputs: label_singleR.csv, singleR_per_cell.csv, celltype_map.csv,",
           "confusion_10x_vs_singleR.csv, confusion_heatmap.pdf/png, fine_label_map_heatmap.pdf/png\n")
