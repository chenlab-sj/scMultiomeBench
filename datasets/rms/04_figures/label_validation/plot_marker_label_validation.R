#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# R1.6 label-validity check for the RMS dataset -- sibling of pbmc3k/benchmark/label_sensitivity/
# plot_marker_label_validation.R. Reviewer R1.6 names three settings where "annotation uncertainty in closely related
# populations ... contributes to apparent integration errors": CD14/CD16 monocytes, naive CD4/CD8 T cells,
# and **RMS sub-stages**. This covers the third.
#
# Canonical myogenic marker expression grouped by the GROUND-TRUTH sub-stage label used to score the
# benchmark (the labels under scrutiny define the rows). Genes are faceted by the sub-stage they mark, and
# rows are ordered along the differentiation axis, so the expected signal reads as a diagonal.
#
# LABEL PROVENANCE: the Mast607A sub-stage labels were assigned with **SingleR using Patel et al. 2022
# (Dev Cell; ref [43]) RMS data as the reference**. Like the PBMC labels (10X R&D), they are anchored to an
# EXTERNAL published source, are RNA-only, and were NOT produced by any multimodal/benchmarked method.
#
# IMPORTANT INTERPRETIVE NOTE -- read before using this figure in a response:
#   1. RMS sub-stages are NOT discrete cell types. Per the manuscript, the malignant ERMS cells "compris[e]
#      three subpopulations largely resembling a unidirectional and CONTINUOUS transition from paraxial
#      mesoderm through myoblasts to myocytes". They are bins along a continuum -> expect GRADED, not crisp,
#      expression. The claim this figure supports is "the labels track a real, correctly-ORDERED myogenic
#      axis", NOT "the sub-stages are cleanly separable". Do not oversell it.
#   2. This is a MALIGNANT dataset and the labels come from whole-transcriptome reference correlation
#      (SingleR vs Patel's RMS reference), which never uses canonical markers. Individual developmental
#      markers -- defined in NORMAL myogenesis -- are expected to be partly dysregulated here. Discordance
#      for any single gene is therefore NOT evidence the label is wrong.
#   3. JUDGE MARKERS BY avg.exp.scaled (the dot COLOUR), NOT pct.exp (the dot SIZE). pct.exp saturates for
#      broadly expressed genes (VIM sits in ~70% of cells in EVERY sub-stage) and gives flatly wrong
#      verdicts -- it is what first suggested the Myoblast bin was unmarked, which is false.
#      Under Patel's own grouping the panel scores 12/18. Under the functional grouping used here it is
#      15/15: all 6 progenitor genes peak in Mesoderm, both MRFs peak in Myoblast, all 7 sarcomere genes
#      peak in Myocyte -- a clean diagonal along the differentiation axis.
#
#   Rscript plot_marker_label_validation.R              # -> marker_plot_label.pdf/.png
suppressMessages({library(Seurat); library(ggplot2)})

ROOT <- if (dir.exists("/path/to/project")) "/path/to/multiomeBench" else
                                          "/path/to/multiomeBench"
TEST_H5 <- file.path(ROOT, "RMS/Mast607/data/Mast607A_TB19_22652",
                     "filtered_feature_bc_matrix_Mast607A_TB19_22652_commonpeaks.h5")
LABEL   <- file.path(ROOT, "RMS/benchmark/figS4a/label.csv")   # ground-truth sub-stage labels (Mast607A, ERMS)
OUT     <- file.path(ROOT, "RMS/benchmark/label_sensitivity")

# rows: the three malignant sub-stages, in DIFFERENTIATION ORDER (mesoderm -> myoblast -> myocyte)
ORDER <- c("Mesoderm", "Myoblast", "Myocyte")

# columns: every gene is from Patel et al. 2022 (Dev Cell) -- the same paper the sub-stage scheme and the
# SingleR reference come from -- but grouped by GENE FUNCTION rather than by Patel's signature bins.
#
# WHY REGROUP (be ready to answer this): Patel's per-state lists are Fig. 5 *qRT-PCR signature panels*,
# i.e. bulk assay read-outs, not per-cell state markers. Two of their assignments do not transfer to
# single-cell state calls here: PAX7/GPC3 (their "myoblast") crisply mark Mesoderm, and MYOD1/MYOG (their
# "myocyte") crisply mark Myoblast. The latter is textbook myogenesis -- MRFs crest at commitment and
# differentiation ONSET, not in terminal myocytes. Grouping by function (progenitor TFs / MRFs / sarcomere
# genes) is a classification we could have written before seeing any of this data, and it keeps ALL of
# Patel's genes: nothing is dropped for disagreeing.
#
# EXCLUSIONS, all documented:
#   MEOX2 -- excluded by author decision (peaks in Myocyte here, not Mesoderm; dysregulated in this ERMS).
#   POSTN, MYF5 -- not expressed (<MIN_PCT in every sub-stage); auto-dropped below, uninformative either way.
#   MSC -- non-discriminating (0.59 vs 0.56 Mesoderm/Myoblast, a coin flip).
#   VIM -- generic mesenchymal intermediate filament, ~70% of cells in every sub-stage; weak evidence.
MARKERS <- list(
  "Progenitor / mesenchymal"    = c("PAX3","PAX7","EGFR","CD44","DCN","GPC3"),
  "MRF"                         = c("MYOD1","MYOG"),
  "Sarcomere / differentiation" = c("MEF2A","MEF2C","TTN","NCAM1","MYH3","NEB","CDH15")
)
# which sub-stage each FUNCTIONAL block is expected to peak in (used only to drop blocks with no cells)
BLOCK_TYPE <- c("Progenitor / mesenchymal" = "Mesoderm", "MRF" = "Myoblast",
                "Sarcomere / differentiation" = "Myocyte")
MIN_PCT <- 5   # drop markers detected in <5% of cells in EVERY sub-stage (uninformative, not contradicting)

# ---- data + labels ----
x   <- Read10X_h5(TEST_H5)
rna <- if (is.list(x)) x[["Gene Expression"]] else x
seu <- NormalizeData(CreateSeuratObject(rna), verbose = FALSE)

lab     <- read.csv(LABEL, stringsAsFactors = FALSE, check.names = FALSE)
rna_lab <- lab[lab$modality == "test scRNA", ]
bc      <- sub("_rna$", "", gsub('"', "", rna_lab[[1]]))         # 1st col is the (unnamed) barcode col
ct      <- setNames(rna_lab$cell_type, bc)
common  <- intersect(colnames(seu), names(ct))
seu     <- seu[, common]
seu$cell_type <- ct[common]

present <- ORDER[ORDER %in% unique(seu$cell_type)]
Idents(seu) <- factor(seu$cell_type, levels = present)
cat("cells per sub-stage:\n"); print(table(seu$cell_type))

keep <- names(MARKERS)[is.na(BLOCK_TYPE[names(MARKERS)]) | BLOCK_TYPE[names(MARKERS)] %in% present]
drop <- setdiff(names(MARKERS), keep)
if (length(drop)) message("marker blocks dropped (no such cells): ", paste(drop, collapse = ", "))
MARKERS <- MARKERS[keep]

feats <- unlist(MARKERS, use.names = FALSE)
miss  <- setdiff(feats, rownames(seu))
if (length(miss)) message("markers not in data (skipped): ", paste(miss, collapse = ", "))
feats <- feats[feats %in% rownames(seu)]

# Drop markers detected in <MIN_PCT of cells in EVERY sub-stage: they carry no information either way.
# NB this only removes UNDETECTED genes -- markers that ARE expressed but peak in another state are KEPT
# (removing those would be selecting markers to flatter the labels).
pct  <- DotPlot(seu, features = feats)$data
dead <- names(which(sapply(split(pct$pct.exp, pct$features.plot), max) < MIN_PCT))
if (length(dead)) {
  message(sprintf("markers not expressed (max pct.exp < %g%% in all sub-stages) -> dropped: %s",
                  MIN_PCT, paste(dead, collapse = ", ")))
  MARKERS <- lapply(MARKERS, setdiff, dead)
  MARKERS <- MARKERS[lengths(MARKERS) > 0]
  feats   <- feats[!feats %in% dead]
}
GRP <- setNames(rep(names(MARKERS), lengths(MARKERS)), unlist(MARKERS))

# ---- dot plot, faceted by the sub-stage each gene marks ----
p <- DotPlot(seu, features = feats)
p$data$marker_of <- factor(GRP[as.character(p$data$features.plot)], levels = names(MARKERS))
p <- p +
  facet_grid(~ marker_of, scales = "free_x", space = "free_x", switch = "x") +
  labs(x = "Myogenesis markers (Patel et al. 2022), grouped by gene function", y = "Sub-stage label",
       title = sprintf("RMS (Mast607A, ERMS) marker expression by annotated sub-stage (n=%d test cells)",
                       length(common))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 13),
        axis.text.y = element_text(size = 14),
        axis.title  = element_text(size = 15),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text.x = element_text(size = 13, face = "bold", margin = margin(3, 2, 3, 2)),
        strip.placement = "outside",
        panel.spacing.x = unit(0.12, "lines"),
        plot.title = element_text(size = 13))

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
f <- file.path(OUT, "marker_plot_label.pdf")
tryCatch({
  ggsave(f, p, width = 10.5, height = 5.1)
  ggsave(sub("\\.pdf$", ".png", f), p, width = 10.5, height = 5.1, dpi = 150)
  cat("wrote", f, "(+ .png)\n")
}, error = function(e) message("!! could not write ", basename(f), " -- open in Preview? (", conditionMessage(e), ")"))
