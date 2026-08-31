#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
## Table S4: effect of the sample-batch-correction option on integration performance (BMMC_d1).
## Paired Wilcoxon signed-rank test over the 11 benchmark metrics, batch vs no-batch, per method
## (alternative = "greater": low p = batch correction significantly improves the metrics).
## Reproduces the original 4 methods and ADDS MIDAS. Source: BMMC_d1/benchmark/ (all batch variants).
suppressMessages(library(dplyr))
R <- "/path/to/multiomeBench/BMMC_d1/benchmark"
sm <- read.csv(file.path(R,"sum_metrics.csv"), check.names=FALSE)
ct <- read.csv(file.path(R,"celltype_metrics.csv"), check.names=FALSE)
aa <- read.csv(file.path(R,"adj_atac_predaccu.csv"), check.names=FALSE)
pd <- read.csv(file.path(R,"peakdist.csv"), check.names=FALSE)
names(pd) <- gsub('"','',names(pd)); pd$method <- gsub('"','',pd$method)
ctm <- ct %>% group_by(method) %>%
  summarise(ks_celltype_mean=mean(ks_celltype), ks_sample_mean=mean(ks_sample),
            ks_inter_celltype_mean=mean(ks_inter_celltype), .groups="drop")
M <- sm %>% select(method, asw, omics_asw, sample_asw, ami, ari, ks.statistic) %>%
  left_join(ctm, by="method") %>%
  left_join(aa %>% select(method, average_accu), by="method") %>%
  left_join(pd %>% select(method, peakdist_adj), by="method")
M[ , c("ks_celltype_mean","asw","ks_sample_mean","sample_asw","ks.statistic","omics_asw","ami","ari","ks_inter_celltype_mean","average_accu","peakdist_adj")] <- round(M[ , c("ks_celltype_mean","asw","ks_sample_mean","sample_asw","ks.statistic","omics_asw","ami","ari","ks_inter_celltype_mean","average_accu","peakdist_adj")], 2)
metrics <- c("ks_celltype_mean","asw","ks_sample_mean","sample_asw","ks.statistic",
             "omics_asw","ami","ari","ks_inter_celltype_mean","average_accu","peakdist_adj")
pairs <- list("scglue"="scglue(batch)","scglue(multiome)"="scglue(multiome,batch)",
              "scVI"="scVI(batch)","BindSC"="BindSC(batch)","MIDAS"="MIDAS(batch)")
res <- data.frame(Methods=character(), p_value=numeric(), stringsAsFactors=FALSE)
for (base in names(pairs)) {
  nb <- M[M$method==base, metrics]; ba <- M[M$method==pairs[[base]], metrics]
  if (nrow(nb)==0 || nrow(ba)==0 || anyNA(nb) || anyNA(ba)) { cat("PROBLEM:", base, "\n"); next }
  w <- wilcox.test(unlist(ba), unlist(nb), paired=TRUE, alternative="greater", exact=TRUE)
  res <- rbind(res, data.frame(Methods=base, p_value=round(w$p.value,3)))
}
print(res)
write.csv(res, file.path(R,"tableS4_batch_correction.csv"), row.names=FALSE)
cat("\nwrote tableS4_batch_correction.csv\n")
