#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# R1.6 / R4 label-validity check: canonical marker expression grouped by the GROUND-TRUTH cell-type label
# used to score the benchmark -- i.e. the labels themselves are under scrutiny, so they (not any annotation
# method's output) define the rows. Shows the labels are biologically real and that the "close pairs" R1.6
# flags (CD14/CD16 mono, naive CD4/CD8 T, Memory/NK-effector T) are separable by canonical markers.
#
# LABEL PROVENANCE: pbmc3k/pbmc10k cell-type labels are from the **10X Genomics R&D team** (Methods: "Cell-
# type labels from the 10X Genomics R&D team were used as surrogate ground truths") -- a THIRD-PARTY
# annotation, NOT derived here with Seurat WNN or any benchmarked method. Don't call these "WNN labels".
#
# Genes are faceted by WHICH CELL TYPE THEY MARK, so the expected signal reads as a diagonal of blocks.
#
#   Rscript plot_marker_label_validation.R                  # 10X ground-truth labels -> marker_plot_label.pdf
#   LABELS=singleR Rscript plot_marker_label_validation.R   # cross-check -> marker_plot_singleR.pdf
suppressMessages({library(Seurat); library(ggplot2)})

ROOT <- if (dir.exists("/path/to/project")) "/path/to/multiomeBench" else
                                          "/path/to/multiomeBench"
TEST_H5 <- file.path(ROOT, "pbmc/pbmc3k/Data/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5")
OUT     <- file.path(ROOT, "pbmc/pbmc3k/benchmark/label_sensitivity")
LABKIND <- Sys.getenv("LABELS", "label")     # "label" = 10X ground truth; "singleR" = cross-check
LABEL   <- if (LABKIND == "singleR") file.path(OUT, "label_singleR.csv") else
                                     file.path(ROOT, "pbmc/pbmc3k/benchmark/fig2b/label.csv")
TITLE   <- if (LABKIND == "singleR") "SingleR (independent, RNA-only) cell-type label" else
                                     "10X-annotated cell-type label"

# rows: cell types, ordered so each type's own markers fall on the diagonal
ORDER <- c("CD14 Monocytes","CD16 Monocytes","Myeloid DC","Plasmacytoid DC","B cells",
           "Naive CD4 T cells","Naive CD8 T cells","Memory T cells","NK/Effector T cells")

# columns: marker panel, grouped by the cell type each gene marks (-> facet strips).
# "Pan-T" is a lineage block (not a scored type) and is always kept.
MARKERS <- list(
  "Pan-T"       = c("CD3D","CD3E"),
  "CD14 Mono"   = c("CD14","LYZ","S100A8"),
  "CD16 Mono"   = c("FCGR3A","MS4A7"),
  "Myeloid DC"  = c("FCER1A","CLEC10A"),
  "pDC"         = c("LILRA4","IRF7"),
  "B"           = c("MS4A1","CD79A"),
  "Naive CD4 T" = c("CCR7","SELL","CD4"),
  "Naive CD8 T" = c("CD8A","CD8B"),
  "Memory T"    = c("IL7R","S100A4"),
  "NK/Effector" = c("GNLY","NKG7","KLRD1","GZMB")
)
# which scored cell type each marker block belongs to (NA = lineage block, always shown)
BLOCK_TYPE <- c("Pan-T" = NA, "CD14 Mono" = "CD14 Monocytes", "CD16 Mono" = "CD16 Monocytes",
                "Myeloid DC" = "Myeloid DC", "pDC" = "Plasmacytoid DC", "B" = "B cells",
                "Naive CD4 T" = "Naive CD4 T cells", "Naive CD8 T" = "Naive CD8 T cells",
                "Memory T" = "Memory T cells", "NK/Effector" = "NK/Effector T cells")

# ---- data + labels ----
x   <- Read10X_h5(TEST_H5)
rna <- if (is.list(x)) x[["Gene Expression"]] else x
seu <- NormalizeData(CreateSeuratObject(rna), verbose = FALSE)

lab     <- read.csv(LABEL, stringsAsFactors = FALSE, check.names = FALSE)
rna_lab <- lab[lab$modality == "test scRNA", c("Barcode", "cell_type")]
rownames(rna_lab) <- sub("_rna$", "", gsub('"', "", rna_lab$Barcode))
common  <- intersect(colnames(seu), rownames(rna_lab))
seu     <- seu[, common]
seu$cell_type <- rna_lab[common, "cell_type"]

present <- ORDER[ORDER %in% unique(seu$cell_type)]
Idents(seu) <- factor(seu$cell_type, levels = present)

# ---- drop marker blocks whose cell type has no cells here (else empty columns with no matching row).
#      pbmc3k's TEST split contains no DC, so the Myeloid DC / pDC blocks drop out automatically. ----
keep  <- names(MARKERS)[is.na(BLOCK_TYPE[names(MARKERS)]) | BLOCK_TYPE[names(MARKERS)] %in% present]
drop  <- setdiff(names(MARKERS), keep)
if (length(drop)) message("marker blocks dropped (no such cells in this split): ", paste(drop, collapse = ", "))
MARKERS <- MARKERS[keep]

GRP   <- setNames(rep(names(MARKERS), lengths(MARKERS)), unlist(MARKERS))
feats <- unlist(MARKERS, use.names = FALSE)
miss  <- setdiff(feats, rownames(seu))
if (length(miss)) message("markers not in data (skipped): ", paste(miss, collapse = ", "))
feats <- feats[feats %in% rownames(seu)]

# ---- dot plot, faceted by the cell type each gene marks ----
p <- DotPlot(seu, features = feats)
p$data$marker_of <- factor(GRP[as.character(p$data$features.plot)], levels = names(MARKERS))
p <- p +
  facet_grid(~ marker_of, scales = "free_x", space = "free_x", switch = "x",
             labeller = labeller(marker_of = c(                       # short strip labels so they fit big
               "Pan-T" = "Pan-T", "CD14 Mono" = "CD14", "CD16 Mono" = "CD16", "Myeloid DC" = "mDC",
               "pDC" = "pDC", "B" = "B", "Naive CD4 T" = "CD4 T", "Naive CD8 T" = "CD8 T",
               "Memory T" = "Mem T", "NK/Effector" = "NK"))) +
  labs(x = "Canonical marker genes (grouped by the cell type they mark)", y = "Cell-type label",
       title = sprintf("pbmc3k marker expression by %s (n=%d test cells)", TITLE, length(common))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 13),
        axis.text.y = element_text(size = 14),
        axis.title  = element_text(size = 15),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text.x = element_text(size = 13, face = "bold", margin = margin(3, 2, 3, 2)),
        strip.placement = "outside",
        panel.spacing.x = unit(0.12, "lines"),
        plot.title = element_text(size = 13))

f <- file.path(OUT, sprintf("marker_plot_%s.pdf", LABKIND))
tryCatch({
  ggsave(f, p, width = 12.0, height = 5.9)
  ggsave(sub("\\.pdf$", ".png", f), p, width = 12.0, height = 5.9, dpi = 150)
  cat("wrote", f, "(+ .png)\n")
}, error = function(e) message("!! could not write ", basename(f), " -- open in Preview? (", conditionMessage(e), ")"))
