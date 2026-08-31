#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# BMMC_d1 SUPPLEMENTARY metrics table = composite "integration performance score" (methods x 11 metrics + score
# + rank), 4-term score (adds the sample-batch term vs pbmc's 3-term). Faithful to the formattable section of
# old/BMMC_d1/plot_cellttype_stat2.R. Adds MaxFuse / MIDAS / MIDAS(batch) / scButterfly.
#
# APPROACH = recompute ALL rows from the new pipeline (consistent with the updated Fig6). Old rows' 10 non-peak
# columns match published exactly; the peak column + composite/rank are recomputed on the new scale (the
# published peak values are not in the repo). Portal is SPLICED from the published raw CSVs (its new-run metrics
# are pending run_portal_sub.sh); its peak uses the published figure value 0.63 -> flagged.
#   Rscript plot_metrics_matrix.R   ->  writes metrics_matrix.csv + metrics_table.html (screenshot to PNG separately)
suppressMessages({library(dplyr); library(formattable); library(htmlwidgets)})

G   <- "/path/to/multiomeBench/BMMC_d1/benchmark"
PUB <- file.path(G, "old/BMMC_d1/benchmark_matrix")   # published raw sources (for Portal splice)
OUT <- file.path(G, "metrics")
dir.create(OUT, showWarnings = FALSE)

# published rows + BindSC/BindSC(batch) + the 4 new methods (exclude Cobolt + scJoint(batch): Cobolt absent from
# the new run; scJoint(batch) not requested)
TARGET <- c("scglue(multiome)", "scglue(multiome,batch)", "scglue", "scglue(batch)", "scVI", "scVI(batch)",
            "scBridge", "simba", "Seurat(CCA)", "Portal", "scJoint",
            "BindSC", "BindSC(batch)",
            "MaxFuse", "MIDAS", "MIDAS(batch)", "scButterfly")  # scDART removed per user

clean <- function(x) trimws(gsub('"', '', x))

## assemble the 11 per-method metric columns from a pipeline dir's 4 CSVs
assemble11 <- function(sum_f, ct_f, accu_f, peak_f) {
  sm <- read.csv(sum_f, check.names = FALSE); sm$method <- clean(sm$method)
  sm <- sm %>% select(method, asw, omics_asw, sample_asw, ami, ari, ks.statistic)
  ct <- read.csv(ct_f); ct$method <- clean(ct$method)
  ct <- ct %>% group_by(method) %>% summarise(ks_celltype_mean = mean(ks_celltype),
              ks_sample_mean = mean(ks_sample), ks_inter_celltype_mean = mean(ks_inter_celltype), .groups = "drop")
  ac <- read.csv(accu_f, check.names = FALSE); ac$method <- clean(ac[[1]])
  ac <- ac %>% select(method, average_accu) %>% filter(method != "random_atac")
  pk <- read.csv(peak_f); pk$method <- clean(pk$method); pk <- pk %>% select(method, peakdist_adj)
  sm %>% merge(ct, "method") %>% merge(ac, "method") %>% merge(pk, "method")
}

## recomputed rows (new pipeline) -- everything except Portal (dropped from the new run)
new11 <- assemble11(file.path(G, "sum_metrics.csv"), file.path(G, "celltype_metrics.csv"),
                    file.path(G, "adj_atac_predaccu.csv"), file.path(G, "peakdist.csv")) %>%
  filter(method %in% TARGET)

## Portal -- splice from published raw CSVs; peak = published figure value (run_portal_sub.sh pending)
portal <- {
  sm <- read.csv(file.path(PUB, "sum_metrics_all.csv"), check.names = FALSE); sm$method <- clean(sm$method)
  sm <- sm %>% filter(method == "Portal") %>% select(method, asw, omics_asw, sample_asw, ami, ari, ks.statistic)
  ct <- read.csv(file.path(PUB, "celltype_metrics2.csv")); ct$method <- clean(ct$method)
  ct <- ct %>% filter(method == "Portal") %>% summarise(method = "Portal",
              ks_celltype_mean = mean(ks_celltype), ks_sample_mean = mean(ks_sample),
              ks_inter_celltype_mean = mean(ks_inter_celltype))
  ac <- read.csv(file.path(G, "old/BMMC_d1/adj_atac_predaccu.csv"), check.names = FALSE); ac$method <- clean(ac[[1]])
  ac <- ac %>% filter(method == "Portal") %>% select(method, average_accu)
  sm %>% merge(ct, "method") %>% merge(ac, "method") %>% mutate(peakdist_adj = 0.63)
}

combined <- bind_rows(new11, portal)
stopifnot(setequal(combined$method, TARGET))

## 4-term composite + rank (raw peakdist, incl. negatives, feeds the score)
mat <- combined %>%
  mutate(score = (ks_celltype_mean + asw) / 2 +                         # biological diversity preservation
                 (ks_sample_mean + sample_asw) / 2 +                    # batch effects correction
                 (ks.statistic + omics_asw) / 2 +                       # omics gaps reduction
                 ((ari + ami) / 2 + ks_inter_celltype_mean +
                    (average_accu + peakdist_adj) / 2) / 3) %>%          # optimal alignment
  mutate(Rank = rank(-score, ties.method = "first"), score = score / 4) %>%
  select(method, ks_celltype_mean, asw, ks_sample_mean, sample_asw, ks.statistic, omics_asw,
         ami, ari, ks_inter_celltype_mean, average_accu, peakdist_adj, score, Rank) %>%
  mutate_at(vars(-c(Rank, method)), round, digits = 2) %>%
  arrange(Rank)

## batch variants -> published "&" display (scglue(multiome,batch)->scglue(multiome)&, scVI(batch)->scVI&, ...)
disp <- function(m) { m <- gsub("\\(multiome,batch\\)", "(multiome)&", m); gsub("\\(batch\\)", "&", m) }
mat$method <- disp(mat$method)

write.csv(mat, file.path(OUT, "metrics_matrix.csv"), row.names = FALSE)

colnames(mat) <- c("Methods", "Same omics cell type distance statistics", "Cell type ASW",
                   "Inter-sample cell type distance statistics", "Sample ASW",
                   "Omics distance statistics", "Omics ASW", "AMI", "ARI",
                   "Inter-omics cell type distance statistics",
                   "Adjusted ATAC cell type prediction accuracy", "Adjusted predicted peak score",
                   "Integration performance score", "Rank")

ident <- function(x) x
clampf <- function(x) ifelse(x <= 0, 0, x)                  # accu/peak bars clamp negatives to width 0
bar1 <- function(color) color_bar(color = color, fun = ident)
bar2 <- function(color) color_bar(color = color, fun = clampf)
ft <- formattable(mat, align = rep("l", ncol(mat)), list(
  area(col = "Same omics cell type distance statistics")     ~ bar1("#FC8D62"),
  area(col = "Cell type ASW")                                ~ bar1("#FC8D62"),
  area(col = "Inter-sample cell type distance statistics")   ~ bar1("#e78ac3"),
  area(col = "Sample ASW")                                   ~ bar1("#e78ac3"),
  area(col = "Omics distance statistics")                    ~ bar1("#8DA0CB"),
  area(col = "Omics ASW")                                    ~ bar1("#8DA0CB"),
  area(col = "Inter-omics cell type distance statistics")    ~ bar1("#66C2A5"),
  area(col = "AMI")                                          ~ bar1("#66C2A5"),
  area(col = "ARI")                                          ~ bar1("#66C2A5"),
  area(col = "Adjusted ATAC cell type prediction accuracy")  ~ bar2("#66C2A5"),
  area(col = "Adjusted predicted peak score")                ~ bar2("#66C2A5")))

saveWidget(as.htmlwidget(ft), file.path(OUT, "metrics_table.html"), selfcontained = FALSE)  # FALSE -> no pandoc (Mac)
cat("wrote metrics_matrix.csv + metrics_table.html\n\n")
print(mat)
