#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# =============================================================================
# 02_compare_donor_composition.R
# After 01_azimuth_annotate_parse.R has been run for each donor, rank the Parse donors
# by how closely their cell-type composition matches pbmc3k -> choose the donor to use.
# Lightweight (reads CSVs + plots); no Seurat/Azimuth -> run interactively or as a tiny job.
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

ROOT    <- "/path/to/multiomeBench"
ANN_DIR <- file.path(ROOT, "pbmc_parse", "azimuth")                              # parse_<DONOR>_composition.csv live here
PBMC3K  <- file.path(ROOT, "pbmc", "pbmc3k", "Data", "pbmc3k_celltype_annotation.csv")
DONORS  <- c("Donor_1", "Donor_2", "Donor_3", "Donor_4")                         # <-- EDIT to your actual sample names
OUT     <- file.path(ANN_DIR, "donor_selection"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

pbmc_parse <- readRDS(file.path(ROOT,  "pbmc_parse/data/e705ee54-5278-41a4-a1f9-0da3ea459cec_processed_matrix.rds"))  # for metadata
pbmc_parse_Mini <- readRDS(file.path(ROOT,  "pbmc_parse/data/c6ef3f1e-2288-4404-b089-50eb701cef06_processed_matrix.rds"))  # for metadata

# the 9 pbmc3k target types (Other/blank ignored in the comparison)
TYPES <- c("CD14 Monocytes","CD16 Monocytes","B cells","Naive CD4 T cells",
           "Naive CD8 T cells","Memory T cells","NK/Effector T cells",
           "Myeloid DC","Plasmacytoid DC")

prop_vec <- function(counts) {                       # counts: named vector -> proportions over TYPES
  v <- setNames(rep(0, length(TYPES)), TYPES)
  common <- intersect(names(counts), TYPES); v[common] <- counts[common]
  v / sum(v)
}

## pbmc3k reference composition
ref_raw <- read.csv(PBMC3K, stringsAsFactors = FALSE)
ref_tab <- table(trimws(ref_raw[[ncol(ref_raw)]]))
ref_p   <- prop_vec(ref_tab)
ref_n   <- sum(ref_tab[intersect(names(ref_tab), TYPES)])

## each donor: distance + coverage
rows <- list(); comp_long <- list()
for (d in DONORS) {
  f <- file.path(ANN_DIR, sprintf("parse_%s_composition.csv", d))
  if (!file.exists(f)) { message("skip (missing): ", f); next }
  cc     <- read.csv(f, stringsAsFactors = FALSE)
  counts <- setNames(cc$n, cc$celltype)
  p      <- prop_vec(counts)
  n      <- sum(counts[intersect(names(counts), TYPES)])
  L1     <- sum(abs(p - ref_p))                                          # total proportion mismatch (0..2)
  m      <- 0.5 * (p + ref_p)                                            # Jensen-Shannon divergence
  JSD    <- 0.5*sum(ifelse(p>0, p*log2(p/m), 0)) + 0.5*sum(ifelse(ref_p>0, ref_p*log2(ref_p/m), 0))
  missing <- TYPES[p == 0]                                              # pbmc3k types absent in this donor
  rows[[d]] <- data.frame(donor = d, n_cells = n, L1 = round(L1, 3), JSD = round(JSD, 3),
                          n_missing = length(missing), missing = paste(missing, collapse = "; "),
                          stringsAsFactors = FALSE)
  comp_long[[d]] <- data.frame(dataset = d, celltype = TYPES, frac = as.numeric(p))
}

## rank: fewest missing types first, then closest composition (L1)
summary_tab <- bind_rows(rows) %>% arrange(n_missing, L1)
write.csv(summary_tab, file.path(OUT, "donor_selection_summary.csv"), row.names = FALSE)
cat(sprintf("\npbmc3k reference: %d cells, %d types present\n", ref_n, sum(ref_p > 0)))
print(summary_tab)
cat(sprintf("\nRECOMMEND -> %s  (missing %d pbmc3k types; L1=%.3f; %d cells)\n",
            summary_tab$donor[1], summary_tab$n_missing[1], summary_tab$L1[1], summary_tab$n_cells[1]))
cat("NB: cell-number match is secondary (you can downsample); coverage + composition matter most.\n")

## visual: pbmc3k vs each donor
plot_df <- bind_rows(data.frame(dataset = "pbmc3k", celltype = TYPES, frac = as.numeric(ref_p)),
                     bind_rows(comp_long))
plot_df$dataset  <- factor(plot_df$dataset, levels = c("pbmc3k", DONORS))
plot_df$celltype <- factor(plot_df$celltype, levels = rev(TYPES))
p <- ggplot(plot_df, aes(celltype, frac, fill = dataset)) +
  geom_col(position = "dodge") + coord_flip() +
  labs(title = "Cell-type composition: pbmc3k vs Parse donors", y = "proportion", x = NULL) +
  theme_bw()
ggsave(file.path(OUT, "donor_composition_comparison.png"), p, width = 9, height = 6, dpi = 200)
message("wrote ", OUT)
