#!/usr/bin/env Rscript
## plot_metrics_matrix.R -- the Fig2B summary matrix (methods x metrics), clean rewrite.
## Faithful to plot_cellttype_stat2.R's formattable section: same merges, same composite
## score, same Rank, same column names + color groups. Differences: reads the clean
## metrics/ CSVs (00_compute_metrics.py + 01_adjust_accuracy.py + peak step), drops the bubble
## plot / dead code, and SAVES the matrix (the old script only rendered to the viewer).
## Usage:  Rscript plot_metrics_matrix.R [indir=.] [outdir=.]
suppressMessages({library(dplyr); library(formattable)})

args   <- commandArgs(trailingOnly = TRUE)
indir  <- ifelse(length(args) >= 1, args[1], ".")
outdir <- ifelse(length(args) >= 2, args[2], indir)
p <- function(f) file.path(indir, f)

## Single-file mode (default): load the shipped final matrix and plot directly.
## Set REBUILD_MATRIX=1 to rebuild it from the per-metric tables (+ the SPLICE_PUBLISHED splice).
matrix_csv <- p("metrics_matrix.csv")
if (Sys.getenv("REBUILD_MATRIX", "0") != "1" && file.exists(matrix_csv)) {
  mat <- read.csv(matrix_csv, check.names = FALSE)  # mirrors write.csv(mat, ..., row.names = FALSE)
} else {

## assemble the 9 per-method Fig2B columns from the 4 source CSVs (sum / celltype / accu / peak)
assemble9 <- function(sum_f, ct_f, accu_f, peak_f, peak_col = "method", remap = FALSE) {
  sm <- read.csv(sum_f)                                             # method, asw, omics_asw, ami, ari, ks.statistic
  ct <- read.csv(ct_f) %>% group_by(method) %>%                     # per-celltype -> means
    summarise(ks_celltype_mean = mean(ks_celltype),
              ks_inter_celltype_mean = mean(ks_inter_celltype), .groups = "drop")
  ac <- read.csv(accu_f) %>% rename(method = 1) %>%                 # method + average_accu
    select(method, average_accu) %>% filter(method != "random_atac")
  pk <- read.csv(peak_f); pk$method <- pk[[peak_col]]               # method + peakdist_adj
  if (remap) pk$method <- recode(pk$method, "Seurat.CCA." = "Seurat(CCA)",
                                 "Seurat.WNN." = "Seurat(WNN)", "scglue.multiome." = "scglue(multiome)")
  pk <- pk %>% filter(method != "random") %>% select(method, peakdist_adj)
  sm %>% merge(ct, "method") %>% merge(ac, "method") %>% merge(pk, "method")
}

## this run's values (the new pipeline) for every method in indir/
new9 <- assemble9(p("sum_metrics.csv"), p("celltype_metrics.csv"),
                  p("adj_atac_predaccu.csv"), p("peakdist.csv"))

## SPLICE (default ON): keep your PUBLISHED values for the original methods + add ONLY the new
## methods from this run, then re-rank. Avoids reviewer questions about slightly-changed numbers.
## (Conos is in the published set, so its value comes from the publication -- no need to recompute.)
## Set SPLICE_PUBLISHED=0 to instead score every method from this run.
if (Sys.getenv("SPLICE_PUBLISHED", "0") == "1") {  # published figure used 0; splice branch is dead code kept for provenance
  bm <- file.path(indir, "..")                                      # benchmark/ (parent of metrics/)
  pub9 <- assemble9(file.path(bm, "metrics/sum_metrics.csv"),
                    file.path(bm, "old/pbmc3k/benchmark_matrix/celltype_metrics.csv"),
                    file.path(bm, "knn_test/adj_atac_predaccu.csv"),
                    file.path(bm, "old/peak_similarity/pbmc3k/peakdist_adj_random.csv"),
                    peak_col = "Method.x", remap = TRUE)
  added   <- new9 %>% filter(!method %in% pub9$method)
  combined <- bind_rows(pub9, added)
  message("SPLICE: ", nrow(pub9), " published + ", nrow(added), " new (",
          paste(added$method, collapse = ", "), ") = ", nrow(combined), " methods")
} else {
  combined <- new9
}

## user: drop scDART + Cobolt from the pbmc_parse (Parse PBMC scRNA + 10X PBMC scATAC) metrics plot
combined <- combined %>% filter(!method %in% c("scDART", "Cobolt"))

## NOTE (Parse cross-platform = UNPAIRED): ks.statistic / ari / ami measure SAME-CELL (paired)
## consistency, which doesn't exist when RNA (Parse) and ATAC (pbmc3k) are different cells. They are
## therefore DROPPED ENTIRELY (not shown) from the matrix + table -- this is real unpaired data, so
## the columns would just be blank. The composite score rescales each group to its valid members:
##   omics gap reduction = omics_asw                                  (was (ks.statistic + omics_asw)/2)
##   alignment accuracy  = (ks_inter_celltype_mean + (average_accu + peakdist_adj)/2)/2  (was .../3 w/ ari,ami)
## bio-conservation is unchanged (both members valid).
mat <- combined %>%
  mutate(score = (ks_celltype_mean + asw) / 2 +                       ## bio-conservation
                 omics_asw +                                          ## omics gap reduction (ks.statistic dropped)
                 (ks_inter_celltype_mean +
                    (average_accu + peakdist_adj) / 2) / 2) %>%       ## alignment accuracy (ari/ami dropped)
  mutate(score = score / 3,
         Rank  = rank(-score, ties.method = "first")) %>%
  select(method, ks_celltype_mean, asw, omics_asw,
         ks_inter_celltype_mean, average_accu, peakdist_adj, score, Rank) %>%
  mutate_at(vars(-c(Rank, method)), round, digits = 2) %>%
  arrange(Rank)

}  ## end REBUILD_MATRIX else-branch (assembly)

## grouped scores (sum_metrics_clean.csv) -- kept for the cross-dataset summary
mat %>%
  transmute(method,
            `Omics gap reduction` = omics_asw,                        ## ks.statistic N/A (unpaired)
            `Bio-conservation`    = (ks_celltype_mean + asw) / 2,
            `Alignment accuracy`  = (ks_inter_celltype_mean +
                                       (average_accu + peakdist_adj) / 2) / 2) %>%  ## (ari+ami) N/A
  mutate(`Overall Performance` = (`Omics gap reduction` + `Bio-conservation` + `Alignment accuracy`) / 3,
         Data = "parse_pbmc3k") %>%
  write.csv(file.path(outdir, "sum_metrics_clean.csv"), row.names = FALSE)

## SAVE the matrix that backs Fig2B (the old script never wrote this)
write.csv(mat, file.path(outdir, "metrics_matrix.csv"), row.names = FALSE)

colnames(mat) <- c("Methods", "Same omics cell type distance statistics", "Cell type ASW",
                   "Omics ASW",
                   "Inter-omics cell type distance statistics",
                   "Adjusted ATAC cell type prediction accuracy",
                   "Adjusted predicted peak score", "Integration performance score", "Rank")

ident <- function(x) x
bar   <- function(color) color_bar(color = color, fun = ident)
ft <- formattable(mat, align = rep("l", ncol(mat)), list(
  area(col = "Same omics cell type distance statistics")    ~ bar("#FC8D62"),
  area(col = "Cell type ASW")                               ~ bar("#FC8D62"),
  ## ks.statistic / ARI / AMI (same-cell metrics) are dropped entirely -- unpaired cross-platform data.
  area(col = "Omics ASW")                                   ~ bar("#8DA0CB"),
  area(col = "Inter-omics cell type distance statistics")   ~ bar("#66C2A5"),
  area(col = "Adjusted ATAC cell type prediction accuracy") ~ bar("#66C2A5"),
  area(col = "Adjusted predicted peak score")               ~ bar("#66C2A5")))

## export the colored table to HTML, then to static PNG + PDF (+ JPG).
## formattable is an HTML widget, so the static image is rendered by a headless browser:
## prefer webshot2 (headless Chrome -> just needs a chrome/chromium on PATH), fall back to
## webshot (needs phantomjs once: webshot::install_phantomjs()). Extension picks the format.
html <- file.path(outdir, "metrics_table.html")
tryCatch({
  htmlwidgets::saveWidget(formattable::as.htmlwidget(ft), html, selfcontained = FALSE)  # FALSE -> no pandoc (Mac)
  shot <- if (requireNamespace("webshot2", quietly = TRUE)) webshot2::webshot
          else if (requireNamespace("webshot", quietly = TRUE)) webshot::webshot
          else NULL
  if (is.null(shot)) {
    message("no webshot2/webshot -> only HTML written. ",
            "install.packages('webshot2') (Chrome) or webshot + webshot::install_phantomjs()")
  } else {
    for (ext in c("png", "pdf", "jpeg"))
      try(shot(html, file.path(outdir, paste0("metrics_table.", ext)), zoom = 2), silent = TRUE)
  }
}, error = function(e) message("table render skipped: ", conditionMessage(e)))

cat("wrote: metrics_matrix.csv, sum_metrics_clean.csv, metrics_table.html (+ png/pdf/jpeg if webshot2/webshot present)\n")
print(ft)
