#!/usr/bin/env Rscript
## RMS Fig5A composite metrics table (PAIRED multiome -> full ARI/AMI/ks.statistic terms, same recipe as
## the BRCA figS4a plot_metrics_matrix.R). The OLD methods keep their PUBLISHED Fig5A values (old/Mast607A/); only
## the 3 NEW methods (MaxFuse, MIDAS, scButterfly) are added from this run (this dir). scDART + Cobolt are
## EXCLUDED, matching the published Fig5A (sum_metrics_clean.csv has 9 methods). Everything is re-ranked.
## The published raw peakdist wasn't saved, so old-method peakdist is reconstructed from the published
## Alignment group (validated to 0.0017); old-method score is pinned to the exact published Overall.
## Usage:  Rscript plot_metrics_matrix.R [indir=.] [outdir=.]
suppressMessages({library(dplyr); library(formattable)})

args   <- commandArgs(trailingOnly = TRUE)
indir  <- ifelse(length(args) >= 1, args[1], ".")
outdir <- ifelse(length(args) >= 2, args[2], indir)
p <- function(f) file.path(indir, f)

## Single-file mode (default): load the shipped final matrix and plot directly.
## Set REBUILD_MATRIX=1 to rebuild it from the per-metric tables (+ the published_reference splice).
matrix_csv <- p("metrics_matrix.csv")
if (Sys.getenv("REBUILD_MATRIX", "0") != "1" && file.exists(matrix_csv)) {
  mat <- read.csv(matrix_csv, check.names = FALSE)   # mirrors write.csv(..., row.names = FALSE)
} else {
  assemble9 <- function(sum_f, ct_f, accu_f, peak_f, peak_col = "method", remap = FALSE) {
    sm <- read.csv(sum_f)
    ct <- read.csv(ct_f) %>% group_by(method) %>%
      summarise(ks_celltype_mean = mean(ks_celltype),
                ks_inter_celltype_mean = mean(ks_inter_celltype), .groups = "drop")
    ac <- read.csv(accu_f) %>% rename(method = 1) %>% select(method, average_accu) %>%
      filter(!method %in% c("random_atac", "random_celltype"))
    pk <- read.csv(peak_f); pk$method <- pk[[peak_col]]
    if (remap) pk$method <- recode(pk$method, "Seurat.CCA." = "Seurat(CCA)", "scglue.multiome." = "scglue(multiome)")
    pk <- pk %>% filter(method != "random") %>% select(method, peakdist_adj)
    sm %>% merge(ct, "method") %>% merge(ac, "method") %>% merge(pk, "method")
  }

  ## ---- NEW methods (this run): MaxFuse / MIDAS / scButterfly (RMS has 3 major cell types, no rare-type issue) ----
  new9 <- assemble9(p("sum_metrics.csv"), p("celltype_metrics.csv"), p("adj_atac_predaccu.csv"), p("peakdist.csv"))
  new9 <- new9 %>% filter(method %in% c("MaxFuse", "MIDAS", "scButterfly"))

  ## ---- OLD methods = published Fig5A (old/Mast607A/); 9 methods, scDART + Cobolt excluded ----
  OLD <- c("BindSC", "Portal", "Seurat(CCA)", "scBridge", "scJoint", "scVI", "scglue", "scglue(multiome)", "simba")
  bm  <- file.path(indir, "..", "published_reference", "Mast607A")  # shipped baselines (results/rms/published_reference/)
  sm  <- read.csv(file.path(bm, "benchmark_matrix/sum_metrics.csv")); sm$method <- recode(sm$method, "Bindsc" = "BindSC")
  ct  <- read.csv(file.path(bm, "benchmark_matrix/celltype_metrics.csv")); ct$method <- recode(ct$method, "Bindsc" = "BindSC")
  ctm <- ct %>% group_by(method) %>%
    summarise(ks_celltype_mean = mean(ks_celltype), ks_inter_celltype_mean = mean(ks_inter_celltype), .groups = "drop")
  ac  <- read.csv(file.path(bm, "knn_test/adj_atac_predaccu.csv")); ac$method <- recode(ac$method, "Bindsc" = "BindSC")
  ac  <- ac %>% select(method, average_accu)
  clean <- read.csv(file.path(bm, "benchmark_matrix/sum_metrics_clean.csv"))
  pub9 <- sm %>% merge(ctm, "method") %>% merge(ac, "method") %>%
    merge(clean[, c("method", "Alignment.accuracy", "Overall.Performance")], "method") %>%
    ## reconstruct the (unsaved) published peakdist from the Alignment group:
    ## Alignment = ((ari+ami)/2 + ks_inter + (accu+peak)/2)/3  =>  peak = 2*(3*Align - (ari+ami)/2 - ks_inter) - accu
    mutate(peakdist_adj = 2 * (3 * Alignment.accuracy - (ari + ami) / 2 - ks_inter_celltype_mean) - average_accu) %>%
    filter(method %in% OLD)

  combined <- bind_rows(pub9 %>% select(-Alignment.accuracy), new9 %>% mutate(Overall.Performance = NA_real_))
  message("SPLICE: ", nrow(pub9), " published (Fig5A) + ", nrow(new9), " new (",
          paste(new9$method, collapse = ", "), ") = ", nrow(combined), " methods")

  ## ---- paired composite score (BRCA figS4a formula); OLD methods pinned to published Overall ----
  ## NOTE: peakdist_adj keeps its real value (scJoint = -1.15, a degenerate worse-than-random peak; the RMS
  ## peak script has no MIN_GROUP guard). The negative is only a DISPLAY issue for the color bar -- handled
  ## below by clipping the BAR WIDTH at 0 (barneg) while the cell TEXT still shows -1.15.
  mat <- combined %>%
    mutate(score = (ks_celltype_mean + asw) / 2 +                       ## bio-conservation
                   (ks.statistic + omics_asw) / 2 +                     ## omics gap reduction
                   ((ari + ami) / 2 + ks_inter_celltype_mean +
                      (average_accu + peakdist_adj) / 2) / 3) %>%       ## alignment accuracy
    mutate(score = score / 3) %>%
    mutate(score = ifelse(is.na(Overall.Performance), score, Overall.Performance),
           Rank  = rank(-score, ties.method = "first")) %>%
    select(method, ks_celltype_mean, asw, ks.statistic, omics_asw, ari, ami,
           ks_inter_celltype_mean, average_accu, peakdist_adj, score, Rank) %>%
    mutate_at(vars(-c(Rank, method)), round, digits = 2) %>%
    arrange(Rank)
}

## grouped scores for the cross-dataset summary
mat %>%
  transmute(method,
            `Omics gap reduction` = (ks.statistic + omics_asw) / 2,
            `Bio-conservation`    = (ks_celltype_mean + asw) / 2,
            `Alignment accuracy`  = ((ari + ami) / 2 + ks_inter_celltype_mean +
                                       (average_accu + peakdist_adj) / 2) / 3) %>%
  mutate(`Overall Performance` = (`Omics gap reduction` + `Bio-conservation` + `Alignment accuracy`) / 3,
         Data = "RMS_Mast607A") %>%
  write.csv(file.path(outdir, "sum_metrics_clean.csv"), row.names = FALSE)

write.csv(mat, file.path(outdir, "metrics_matrix.csv"), row.names = FALSE)

colnames(mat) <- c("Methods", "Same omics cell type distance statistics", "Cell type ASW",
                   "Omics distance statistics", "Omics ASW", "ARI", "AMI",
                   "Inter-omics cell type distance statistics",
                   "Adjusted ATAC cell type prediction accuracy",
                   "Adjusted predicted peak score", "Integration performance score", "Rank")

ident <- function(x) x
bar   <- function(color) color_bar(color = color, fun = ident)
## peak column can go negative (scJoint = -1.15). color_bar sets width = percent(fun(x)); a negative
## width distorts / is invalid. Clip the BAR WIDTH at 0 (text is unaffected -> the cell still shows -1.15).
barneg <- function(color) color_bar(color = color, fun = function(x) pmax(as.numeric(x), 0))
ft <- formattable(mat, align = rep("l", ncol(mat)), list(
  area(col = "Same omics cell type distance statistics")    ~ bar("#FC8D62"),
  area(col = "Cell type ASW")                               ~ bar("#FC8D62"),
  area(col = "Omics distance statistics")                   ~ bar("#8DA0CB"),
  area(col = "Omics ASW")                                   ~ bar("#8DA0CB"),
  area(col = "Inter-omics cell type distance statistics")   ~ bar("#66C2A5"),
  area(col = "ARI")                                         ~ bar("#66C2A5"),
  area(col = "AMI")                                         ~ bar("#66C2A5"),
  area(col = "Adjusted ATAC cell type prediction accuracy") ~ bar("#66C2A5"),
  area(col = "Adjusted predicted peak score")               ~ barneg("#66C2A5")))

html <- file.path(outdir, "metrics_table.html")
tryCatch({
  htmlwidgets::saveWidget(formattable::as.htmlwidget(ft), html, selfcontained = TRUE)
  shot <- if (requireNamespace("webshot2", quietly = TRUE)) webshot2::webshot
          else if (requireNamespace("webshot", quietly = TRUE)) webshot::webshot
          else NULL
  if (is.null(shot)) message("no webshot2/webshot -> only HTML written.")
  else for (ext in c("png", "pdf")) try(shot(html, file.path(outdir, paste0("metrics_table.", ext)), zoom = 2), silent = TRUE)
}, error = function(e) message("table render skipped: ", conditionMessage(e)))

cat("wrote: metrics_matrix.csv, sum_metrics_clean.csv, metrics_table.html\n")
print(ft)
