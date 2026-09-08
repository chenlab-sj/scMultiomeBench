#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README.
## Weighting-sensitivity analysis for the final ranking (reviewer R1.4/R4-E: "80/20 weighting not justified").
## Reproduces the exact performance + usability pipeline of sum_perfom_plot.R, then sweeps the performance
## weight w from 0 to 1 and asks: is the top-ranked method stable, and how correlated is each ranking to the
## published 80/20 (w=0.8) ordering?  The published score sums 4 per-dataset performance values + 1 usability
## value (=> 4:1 = 80:20), which equals ordering by  score(w) = w*mean_perf + (1-w)*usability  at w=0.8.
suppressMessages({library(dplyr); library(tidyr); library(ggplot2)})
SUM <- "/path/to/multiomeBench/sum_plot"
CL  <- "/path/to/multiomeBench/common"

top10_methods <- c('scglue(multiome)','scVI','scglue','scJoint','Seurat(CCA)','Portal','simba','BindSC',
                   'scBridge','MaxFuse','MIDAS','scButterfly')

## ---- performance: 4 datasets' Overall Performance (old methods; pbmc3k unfiltered as in the original) ----
pbmc3k     <- read.csv(file.path(CL,'pbmc3k/benchmark_matrix/sum_metrics_clean.csv'),  row.names = 1)
RMS        <- read.csv(file.path(CL,'Mast607A/benchmark_matrix/sum_metrics_clean.csv'), row.names = 1)
RMS        <- RMS[RMS$method %in% top10_methods,]
BRCA_1samp <- read.csv(file.path(CL,'HT243B1-S1H4/benchmark_matrix/sum_metrics_clean.csv'), row.names = 1)
BRCA_1samp <- BRCA_1samp[BRCA_1samp$method %in% top10_methods,]
BMMC       <- read.csv(file.path(CL,'BMMC_d1/sum_metrics_clean.csv'), row.names = 1)
BMMC       <- BMMC[BMMC$method %in% top10_methods,]

## ---- usability (comp_res), exactly as sum_perfom_plot.R ----
rev_norm <- function(x) ifelse(is.na(x), NA, 1 - ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))))
comp_res <- read.csv(file.path(SUM,'runtime_memory_sel_plus.csv'))
comp_sum <- comp_res %>% group_by(methods) %>%
  summarise(avg_running_time = mean(running.time, na.rm = TRUE),
            avg_memory_usage  = mean(memory.usage, na.rm = TRUE), .groups = "drop")
comp_sum$avg_running_time <- ifelse(comp_sum$avg_running_time > 200000, 200000, comp_sum$avg_running_time)
comp_res <- comp_res[,c('methods','GPU.requirment','Data.requirement')] %>% merge(comp_sum, by = 'methods')
comp_res$Data <- "Usability"
comp_res$`Running time` <- rev_norm(comp_res$avg_running_time)
comp_res$`Memory usage` <- rev_norm(comp_res$avg_memory_usage)
comp_res <- comp_res %>% mutate(`Overall usability` = (`GPU.requirment` + `Running time` + `Memory usage`)/6 + (`Data.requirement`/2))
colnames(comp_res) <- c('method','GPU requirment','Data requirement','avg_running_time','avg_memory_usage','Data','Running time','Memory usage','Overall usability')

## ---- combine + splice new methods, then GLOBAL min-max normalize the performance (as in the original) ----
combined_res <- bind_rows(list(pbmc3k, RMS, BRCA_1samp, BMMC, comp_res))
newg <- read.csv(file.path(SUM,'new_methods_grouped.csv'), check.names = TRUE)
combined_res <- bind_rows(combined_res, newg)
min_max <- function(x) ifelse(is.na(x), NA, (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
combined_res$perf_norm <- min_max(combined_res$Overall.Performance)

## ---- per-method mean normalized performance + usability (top 12 methods) ----
perf <- combined_res %>% filter(method %in% top10_methods, !is.na(perf_norm)) %>%
  group_by(method) %>% summarise(perf_mean = mean(perf_norm), n_ds = n(), .groups = "drop")
usab <- combined_res %>% filter(method %in% top10_methods, !is.na(`Overall usability`)) %>%
  group_by(method) %>% summarise(usab = mean(`Overall usability`), .groups = "drop")
tab <- merge(perf, usab, by = 'method')
cat("per-method inputs (perf_mean over", unique(perf$n_ds), "datasets; usability):\n")
print(tab %>% arrange(desc(perf_mean)), row.names = FALSE)

## ---- sweep w = 0..1 ----
ws <- seq(0, 1, by = 0.05)
score_w <- function(w) w * tab$perf_mean + (1 - w) * tab$usab
rankmat <- sapply(ws, function(w) rank(-score_w(w), ties.method = "min"))
rownames(rankmat) <- tab$method; colnames(rankmat) <- ws

ref <- rank(-score_w(0.8), ties.method = "min")                       # published 80/20 ordering
rho <- sapply(ws, function(w) suppressWarnings(cor(rank(-score_w(w)), ref, method = "spearman")))
top_by_w <- tab$method[apply(rankmat, 2, which.min)]

cat("\n=== sweep: top method + scglue(multiome) rank + Spearman vs 80/20 ===\n")
summ <- data.frame(w = ws, top_method = top_by_w,
                   `scglue(multiome)_rank` = as.integer(rankmat['scglue(multiome)', ]),
                   `scglue_rank` = as.integer(rankmat['scglue', ]),
                   spearman_vs_80_20 = round(rho, 3), check.names = FALSE)
print(summ, row.names = FALSE)

sg_top_w  <- ws[top_by_w == 'scglue(multiome)']
cat(sprintf("\nscglue(multiome) is #1 for performance weight w in [%.2f, %.2f]  (published 80/20 = w=0.80)\n",
            min(sg_top_w), max(sg_top_w)))
cat(sprintf("Spearman rho vs the 80/20 ordering stays >= %.2f for w in [0.5,1.0]; overall min rho = %.2f\n",
            min(rho[ws >= 0.5]), min(rho)))
write.csv(cbind(method = rownames(rankmat), as.data.frame(rankmat)),
          file.path(SUM,'weighting_sensitivity_ranks.csv'), row.names = FALSE)

## ---- plot: rank vs weight (bump) -- SAME method palette as the manuscript Fig 7B (sum_perfom_plot.R) ----
method_colors <- c('scglue' = "#1f77b4", 'scglue(multiome)' = "#aec7e8", 'Portal' = "#ffbb78",
                   'scBridge' = "#2ca02c", 'simba' = "#98df8a", 'Seurat(CCA)' = "#d62728",
                   'scJoint' = "#ff9896", 'scVI' = "#9467bd", 'MinNet' = "#c5b0d5", 'BindSC' = "#e377c2",
                   'MaxFuse' = "#8c564b", 'MIDAS' = "#17becf", 'scButterfly' = "#bcbd22")
long <- as.data.frame(rankmat) %>% mutate(method = rownames(rankmat)) %>%
  pivot_longer(-method, names_to = "w", values_to = "rank") %>% mutate(w = as.numeric(w))
p <- ggplot(long, aes(x = w, y = rank, group = method, color = method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "grey40") +
  annotate("text", x = 0.8, y = 0.4, label = "published (80/20)", size = 3.2, color = "grey30") +
  scale_color_manual(values = method_colors) +
  scale_y_reverse(breaks = 1:nrow(tab)) +
  scale_x_continuous(breaks = seq(0,1,0.2), labels = paste0(seq(0,100,20),"%")) +
  labs(x = "Weight on integration performance (usability = 1 - w)", y = "Rank",
       title = "Final-ranking sensitivity to the performance-vs-usability weighting") +
  theme_classic() + theme(legend.position = "right", legend.title = element_blank(),
                          text = element_text(size = 12))
ggsave(file.path(SUM,'weighting_sensitivity.pdf'), p, width = 9, height = 5)
ggsave(file.path(SUM,'weighting_sensitivity.png'), p, width = 9, height = 5, dpi = 150)
cat("\nwrote weighting_sensitivity.png/pdf + weighting_sensitivity_ranks.csv\n")
