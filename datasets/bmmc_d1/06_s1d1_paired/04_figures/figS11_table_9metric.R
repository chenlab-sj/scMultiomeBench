#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# Complete-metric (9-metric) SUPPLEMENTARY table for the BMMC s1d1 PAIRED multiome baseline.
# The s1d1_paired run is genuinely paired (same cells), so the three same-cell metrics (ks.statistic,
# ARI, AMI) are valid here and are INCLUDED -- unlike the Fig 6B scatter, which drops them to compare
# fairly against the unpaired cross-site. Single site -> no sample/batch facet, so this uses the
# standard 3-term composite (bio + omics + alignment)/3, identical to the pbmc3k fig2b table.
# Major cell-type basis (test-ATAC > 100), matching the scatter: metrics from major/, peak from peakdist.csv.
#   Rscript figS11_table_9metric.R  ->  s1d1_paired_metrics_matrix_9metric.csv + s1d1_paired_metrics_table_9metric.html
suppressMessages({library(dplyr); library(formattable); library(htmlwidgets)})

D   <- "/path/to/multiomeBench/BMMC_d1/benchmark/s1d1_paired"
clean <- function(x) trimws(gsub('"', '', x))

## assemble the 9 per-method metric columns (major basis for sum/celltype/accu; top-level peakdist)
sm <- read.csv(file.path(D, "major/sum_metrics.csv"), check.names = FALSE); sm$method <- clean(sm$method)
sm <- sm %>% select(method, asw, omics_asw, ami, ari, ks.statistic)
ct <- read.csv(file.path(D, "major/celltype_metrics.csv")); ct$method <- clean(ct$method)
ct <- ct %>% group_by(method) %>% summarise(ks_celltype_mean = mean(ks_celltype),
            ks_inter_celltype_mean = mean(ks_inter_celltype), .groups = "drop")
ac <- read.csv(file.path(D, "major/adj_atac_predaccu.csv"), check.names = FALSE); ac$method <- clean(ac[[1]])
ac <- ac %>% select(method, average_accu) %>% filter(method != "random_atac")
pk <- read.csv(file.path(D, "peakdist.csv")); pk$method <- clean(pk$method); pk <- pk %>% select(method, peakdist_adj)
combined <- sm %>% merge(ct, "method") %>% merge(ac, "method") %>% merge(pk, "method")

## 3-term composite (no batch facet) + rank -- raw peakdist (incl. negatives) feeds the score
mat <- combined %>%
  mutate(score = (ks_celltype_mean + asw) / 2 +                          # biological diversity preservation
                 (ks.statistic + omics_asw) / 2 +                        # omics gaps reduction (same-cell metric kept)
                 ((ari + ami) / 2 + ks_inter_celltype_mean +
                    (average_accu + peakdist_adj) / 2) / 3) %>%          # optimal alignment (ARI/AMI kept)
  mutate(Rank = rank(-score, ties.method = "first"), score = score / 3) %>%
  select(method, ks_celltype_mean, asw, ks.statistic, omics_asw, ami, ari,
         ks_inter_celltype_mean, average_accu, peakdist_adj, score, Rank) %>%
  mutate_at(vars(-c(Rank, method)), round, digits = 2) %>%
  arrange(Rank)

write.csv(mat, file.path(D, "s1d1_paired_metrics_matrix_9metric.csv"), row.names = FALSE)

colnames(mat) <- c("Methods", "Same omics cell type distance statistics", "Cell type ASW",
                   "Omics distance statistics", "Omics ASW", "AMI", "ARI",
                   "Inter-omics cell type distance statistics",
                   "Adjusted ATAC cell type prediction accuracy", "Adjusted predicted peak score",
                   "Integration performance score", "Rank")

ident  <- function(x) x
clampf <- function(x) ifelse(x <= 0, 0, x)                 # accu/peak bars clamp negatives to width 0
bar1 <- function(color) color_bar(color = color, fun = ident)
bar2 <- function(color) color_bar(color = color, fun = clampf)
ft <- formattable(mat, align = rep("l", ncol(mat)), list(
  area(col = "Same omics cell type distance statistics")     ~ bar1("#FC8D62"),
  area(col = "Cell type ASW")                                ~ bar1("#FC8D62"),
  area(col = "Omics distance statistics")                    ~ bar1("#8DA0CB"),
  area(col = "Omics ASW")                                    ~ bar1("#8DA0CB"),
  area(col = "AMI")                                          ~ bar1("#66C2A5"),
  area(col = "ARI")                                          ~ bar1("#66C2A5"),
  area(col = "Inter-omics cell type distance statistics")    ~ bar1("#66C2A5"),
  area(col = "Adjusted ATAC cell type prediction accuracy")  ~ bar2("#66C2A5"),
  area(col = "Adjusted predicted peak score")                ~ bar2("#66C2A5")))

saveWidget(as.htmlwidget(ft), file.path(D, "s1d1_paired_metrics_table_9metric.html"), selfcontained = FALSE)
cat("wrote s1d1_paired_metrics_matrix_9metric.csv + s1d1_paired_metrics_table_9metric.html\n\n")
print(mat)
