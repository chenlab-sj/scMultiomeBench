#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
## 02_peak_similarity.R -- adjusted peak score (peakdist_adj) for Fig2B, clean + reproducible.
##
## Faithful to old/peak_similarity/pbmc3k/peak_similarity.R (Signac group-coverage bigWig tracks
## per predicted cell type, then track similarity vs the true-annotation track), with:
##   * reads the per-method knn_pred_label__<m>.csv written by 00_compute_metrics.py (combines them)
##   * random baseline = label.csv's FIXED random_atac column (the old script re-sampled with no
##     set.seed -> non-reproducible); this makes the baseline stable across runs
##   * CACHES bigWigs: a group whose <out>/<group>/*.bw already exist is skipped, so the
##     method-independent annotation/random tracks (and any prior methods) are exported once
##   * drops the unused spearman/cosine adjustments (Fig2B uses only the scaled-euclidean one)
##
## peakdist_adj per method = mean over cell types of (random - method)/random  on the
## window-count-scaled euclidean distance of (predicted-group track) vs (true-annotation track).
## 0 == as far from truth as random; 1 == identical to truth; higher == better.
##
## Run on a node with R + Signac/Seurat + rtracklayer/GenomicRanges/AnnotationHub/proxy/lsa.
##   Rscript 02_peak_similarity.R          # paths via env vars below (defaults = the pbmc3k data)
## Output: <PEAK_OUT>/peak_similarity.csv (raw) and peakdist.csv (method, peakdist_adj) in FIG2B_DIR.

EXPORT_BWG <- Sys.getenv("EXPORT_BWG", "/path/to/multiomeBench/common/peak_similarity/export_groupbwg.R")
source(EXPORT_BWG)                          # provides ExportGroupBW()
suppressMessages({library(Signac); library(Seurat)})

## RMS Mast607A uses the hg38 cellranger output (chr-prefixed peaks, hg38 gaps) -- matches the original
## peak_similarity (genome stays "hg38" below). FRAG/H5 now env-overridable; 02_peak_similarity_sub.sh points them at
## the hg38 Mast607A files.
DATA  <- Sys.getenv("PEAK_DATA", "/path/to/data/pbmc3k")
FRAG  <- Sys.getenv("FRAG", file.path(DATA, "pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz"))
H5    <- Sys.getenv("H5",   file.path(DATA, "pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5"))
FIG2B <- Sys.getenv("FIG2B_DIR", ".")       # holds knn_pred_label__*.csv + label.csv
OUT   <- Sys.getenv("PEAK_OUT", file.path(FIG2B, "peak"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- combine the per-method KNN predicted labels (check.names=FALSE keeps exact method names) ----
pred_files <- list.files(FIG2B, pattern = "^knn_pred_label__.*\\.csv$", full.names = TRUE)
stopifnot(length(pred_files) > 0)
KNN.label <- NULL
for (f in pred_files) {
  d <- read.csv(f, row.names = 1, check.names = FALSE)   # one column, named like the method
  if (is.null(KNN.label)) { KNN.label <- d } else {
    m <- merge(KNN.label, d, by = "row.names", all = TRUE)
    rownames(m) <- m$Row.names; m$Row.names <- NULL; KNN.label <- m
  }
}
rownames(KNN.label) <- gsub("_atac", "", rownames(KNN.label))
bc <- rownames(KNN.label)

## ---- annotation (truth) + FIXED random baseline from label.csv ----
annot <- read.csv(file.path(FIG2B, "label.csv"), row.names = 1)
atac.annot <- annot[annot$modality == "test scATAC", c("cell_type", "random_celltype")]  # RMS label uses random_celltype (not random_atac)
colnames(atac.annot) <- c("annotation", "random")
rownames(atac.annot) <- gsub("_atac", "", rownames(atac.annot))
atac.annot <- atac.annot[bc, ]
atac.meta <- merge(atac.annot, KNN.label, by = "row.names", all = TRUE)
rownames(atac.meta) <- atac.meta$Row.names; atac.meta$Row.names <- NULL
atac.meta[] <- lapply(atac.meta, function(x) gsub("[ /]", "_", x))   # spaces/slash -> underscore

## ---- Seurat ChromatinAssay ----
x10 <- Read10X_h5(H5)
atac <- CreateSeuratObject(
  counts = CreateChromatinAssay(counts = x10$Peaks, sep = c(":", "-"), genome = "hg38", fragments = FRAG),
  assay = "peaks", meta.data = atac.meta)

## ---- export group bigWigs (CACHE: skip groups whose .bw already exist) ----
for (col in colnames(atac.meta)) {
  mdir <- gsub("\\.", "_", file.path(OUT, col))
  if (dir.exists(mdir) && length(list.files(mdir, pattern = "\\.bw$")) > 0) {
    message("cached, skip: ", col); next
  }
  dir.create(mdir, showWarnings = FALSE)
  message("export bigWig: ", col)
  ExportGroupBW(atac, assay = "peaks", group.by = col, idents = NULL, normMethod = "RC",
                tileSize = 100000, minCells = 0, cutoff = NULL,
                chromosome = paste0("chr", 1:22), outdir = mdir, verbose = TRUE)
  file.remove(list.files(mdir, pattern = "\\.bed$", full.names = TRUE))   # free disk
}

## ---- gap removal + similarity ----
suppressMessages({library(rtracklayer); library(GenomicRanges); library(AnnotationHub); library(proxy); library(lsa)})
ah <- AnnotationHub()
gaps <- AnnotationHub::query(ah, c("gap", "Homo sapiens", "hg38"))
gaps <- reduce(do.call(c, list(gaps[["AH107355"]], gaps[["AH107356"]], gaps[["AH107357"]],
                               gaps[["AH107358"]], gaps[["AH107359"]])))
gap.ranges <- GRanges(seqnames = gaps@seqnames, ranges = gaps@ranges)

gaps_rm <- function(bw.path, bw.file, gap.ranges) {
  bw <- import(file.path(bw.path, bw.file))
  bw[-queryHits(findOverlaps(bw, gap.ranges))]
}
get_similarity <- function(g1, g2) {        # window-count-scaled euclidean of two coverage tracks
  colnames(mcols(g1)) <- "score1"; colnames(mcols(g2)) <- "score2"
  mg <- merge(g1, g2)
  euc <- dist(rbind(mg$score1, mg$score2))
  as.numeric(euc) / sqrt(length(mg))        # euclidean_distance.scale
}

annot.dir <- gsub("\\.", "_", file.path(OUT, "annotation"))
bw.files  <- list.files(annot.dir)
annot_bw  <- setNames(lapply(bw.files, function(f) gaps_rm(annot.dir, f, gap.ranges)),
                      sub("-.*", "", bw.files))

rows <- list()
for (col in colnames(atac.meta)) {
  if (col == "annotation") next            # annotation is the reference
  mdir <- gsub("\\.", "_", file.path(OUT, col))
  for (f in bw.files) {
    ct <- sub("-.*", "", f)
    if (file.exists(file.path(mdir, f))) {
      d <- get_similarity(gaps_rm(mdir, f, gap.ranges), annot_bw[[ct]])
      rows[[paste(col, ct, sep = "_")]] <- data.frame(Method = col, CellType = ct,
                                                      euclidean_distance.scale = d)
    }
  }
}
sim <- do.call(rbind, rows)
write.csv(sim, file.path(OUT, "peak_similarity.csv"), row.names = FALSE)

## ---- adjust to the random baseline -> peakdist_adj per method ----
suppressMessages(library(dplyr))
rand <- sim %>% filter(Method == "random") %>% select(CellType, rand = euclidean_distance.scale)
peakdist <- sim %>%
  filter(!Method %in% c("random", "annotation")) %>%
  merge(rand, by = "CellType") %>%
  mutate(eucdist_adj = (rand - euclidean_distance.scale) / rand) %>%
  group_by(method = Method) %>%
  summarise(peakdist_adj = mean(eucdist_adj), .groups = "drop")

write.csv(peakdist, file.path(FIG2B, "peakdist.csv"), row.names = FALSE)
cat("wrote", file.path(FIG2B, "peakdist.csv"), "-", nrow(peakdist), "methods\n")
print(peakdist)
