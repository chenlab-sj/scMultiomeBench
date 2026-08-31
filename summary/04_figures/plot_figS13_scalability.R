#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README.
## Standalone scalability figure (Fig S8): running time + memory usage vs cell number (test RNA + ATAC), per
## method, across all datasets. Reads runtime_memory_sel_plus.csv (one row per dataset x method; cols include
## running.time, memory.usage, Test.cell.counts). Uses the manuscript method palette (Fig 7B).
##   Rscript plot_scalability.R
suppressMessages({library(dplyr); library(reshape2); library(ggplot2)})
SUM <- "/path/to/multiomeBench/sum_plot"

## methods shown in the main figure (12 = 9 published + 3 new: MaxFuse/MIDAS/scButterfly)
methods_sel <- c('scglue(multiome)','scVI','scglue','scJoint','Seurat(CCA)','Portal','simba','BindSC',
                 'scBridge','MaxFuse','MIDAS','scButterfly')
method_colors <- c('scglue' = "#1f77b4", 'scglue(multiome)' = "#aec7e8", 'Portal' = "#ffbb78",
                   'scBridge' = "#2ca02c", 'simba' = "#98df8a", 'Seurat(CCA)' = "#d62728",
                   'scJoint' = "#ff9896", 'scVI' = "#9467bd", 'MinNet' = "#c5b0d5", 'BindSC' = "#e377c2",
                   'MaxFuse' = "#8c564b", 'MIDAS' = "#17becf", 'scButterfly' = "#bcbd22")

comp <- read.csv(file.path(SUM, 'runtime_memory_sel_plus.csv'), check.names = TRUE)
comp <- comp %>% filter(methods %in% methods_sel)
comp <- comp %>% filter(Dataset %in% c("pbmc3k", "BRCA", "Mast607A", "BMMC"))   # 4 main benchmark datasets (one per tissue; drops HT263 + pbmc10k)

## report per-method dataset coverage so missing points are visible
cov <- comp %>% group_by(methods) %>% summarise(n_datasets = n(),
        datasets = paste(sort(unique(Dataset)), collapse = ","), .groups = "drop")
cat("dataset coverage per method:\n"); print(as.data.frame(cov), row.names = FALSE)

long <- melt(comp, id.vars = c("Dataset","methods","Test.cell.counts"),
             measure.vars = c("running.time","memory.usage"),
             variable.name = "Metric", value.name = "Value")
long$Metric <- recode(long$Metric, "running.time" = "Running time (s)", "memory.usage" = "Memory usage (MB)")

p <- ggplot(long, aes(x = Test.cell.counts, y = Value, color = methods)) +
  geom_line(aes(group = methods), linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 3) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = method_colors) +
  labs(x = "Cell number (test scRNA + scATAC)", y = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        text = element_text(size = 13, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 13),
        strip.background = element_rect(fill = NA, colour = NA))

ggsave(file.path(SUM, 'scalability.pdf'), p, width = 9, height = 4.6)
ggsave(file.path(SUM, 'scalability.png'), p, width = 9, height = 4.6, dpi = 150)
cat("\nwrote scalability.pdf/png\n")
