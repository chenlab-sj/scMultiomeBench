#!/usr/bin/env Rscript
# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
## BMMC Fig2b peak step -- adjusted peak score (peakdist_adj), 3 test batches (s2d1/s4d1/s1d1).
## Clean rewrite of peak_similarity/BMMC_d1/peak_similarity.R: reads the per-method
## knn_pred_label__<m>.csv (from 03_knn_bmmc.py), builds one ChromatinAssay per batch from the
## cellranger-arc fragments, exports per-cell-type coverage bigWigs grouped by each method's
## predicted labels (+ truth 'annotation' + 'random' baseline), then
##   peakdist_adj = mean over cell types of (random - method)/random   on the window-count-scaled
##   euclidean distance of (predicted-group track) vs (true-annotation track).
## Fixes the original's atac3-uses-atac2-annot bug; caches bigWigs (a group dir with *.bw is skipped);
## uses old/BMMC_d1/label.csv (same _atacN barcodes as the KNN files) so nothing mis-aligns.
## Output: <PEAK_OUT>/peak_similarity.csv (raw) + peakdist.csv (method, peakdist_adj) in FIG2B_DIR.
EXPORT_BWG <- Sys.getenv("EXPORT_BWG", "/path/to/multiomeBench/common/peak_similarity/export_groupbwg.R")
source(EXPORT_BWG)                                  # ExportGroupBW()
suppressMessages({library(Signac); library(Seurat)})

FIG2B <- Sys.getenv("FIG2B_DIR", ".")              # holds knn_pred_label__*.csv
OUT   <- Sys.getenv("PEAK_OUT", file.path(FIG2B, "peak"))
LABEL <- Sys.getenv("LABEL", file.path(FIG2B, "old/BMMC_d1/label.csv"))
SRA   <- Sys.getenv("BMMC_SRA", "/path/to/data/BMMC/NCBI_sra")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
## _atacN suffix -> batch dir + modality  (verified: atac1=s2d1, atac2=s4d1, atac3=s1d1)
batches <- list(list(suf="atac1", name="s2d1", modality="s2d1 scATAC"),
                list(suf="atac2", name="s4d1", modality="s4d1 scATAC"),
                list(suf="atac3", name="s1d1", modality="s1d1 scATAC"))

## ---- combine the per-method KNN predicted labels (index = <bc>-1_atacN) ----
pred_files <- sort(list.files(FIG2B, pattern = "^knn_pred_label__.*\\.csv$", full.names = TRUE))
stopifnot(length(pred_files) > 0)
## PEAK_LAST (comma-sep method names) -> force these to the END of the column order so the safe tokens
## M1.. for the already-computed methods stay put and their cached peak/M*/ bigWigs remain valid; only
## the appended method (e.g. Portal) gets a new token + fresh export. Keeps a resubmit cheap + correct.
peak_last <- strsplit(Sys.getenv("PEAK_LAST", ""), ",")[[1]]
peak_last <- trimws(peak_last[peak_last != ""])
if (length(peak_last)) {
  is_last <- sapply(pred_files, function(f) any(sapply(peak_last,
                    function(m) grepl(paste0("__", m, "\\.csv$"), f, fixed = FALSE))))
  pred_files <- c(pred_files[!is_last], pred_files[is_last])
  message("PEAK_LAST -> appended last: ", paste(basename(pred_files[is_last]), collapse = ", "))
}
KNN.label <- NULL
for (f in pred_files) {
  d <- read.csv(f, row.names = 1, check.names = FALSE)
  KNN.label <- if (is.null(KNN.label)) d else {
    m <- merge(KNN.label, d, by = "row.names", all = TRUE)
    rownames(m) <- m$Row.names; m$Row.names <- NULL; m }
}
## methods -> safe tokens M1.. (avoid Seurat mangling of "scglue(multiome)" etc.); map back at the end
methods_orig <- colnames(KNN.label)
safe <- paste0("M", seq_along(methods_orig)); names(methods_orig) <- safe
colnames(KNN.label) <- safe

annot <- read.csv(LABEL, row.names = 1)

## ---- per-batch Seurat ChromatinAssay + meta (annotation, random, M1..) ----
objs <- list()
for (b in batches) {
  bc <- rownames(KNN.label)[grep(paste0(b$suf, "$"), rownames(KNN.label))]
  an <- annot[annot$modality == b$modality, c("cell_type", "random_atac")]
  colnames(an) <- c("annotation", "random"); an <- an[bc, , drop = FALSE]
  meta <- merge(an, KNN.label[bc, , drop = FALSE], by = "row.names", all = TRUE)
  rownames(meta) <- gsub(paste0("_", b$suf), "", meta$Row.names); meta$Row.names <- NULL
  meta[] <- lapply(meta, function(x) gsub("[ /+]", "_", x))       # spaces/slash/+ -> underscore
  x10 <- Read10X_h5(file.path(SRA, b$name, "outs", "filtered_feature_bc_matrix.h5"))
  ca  <- CreateChromatinAssay(counts = x10$Peaks, sep = c(":", "-"), genome = "hg38",
                              fragments = file.path(SRA, b$name, "outs", "atac_fragments.tsv.gz"))
  objs[[b$name]] <- CreateSeuratObject(counts = ca, assay = "peaks", meta.data = meta)
}
atac <- merge(objs[[1]], objs[-1])
celltypes <- unique(objs[[1]]$annotation); celltypes <- celltypes[!is.na(celltypes)]
atac <- subset(atac, subset = annotation %in% celltypes)

## ---- export coverage bigWigs per group column (annotation + random + M1..), CACHED ----
group_cols <- c("annotation", "random", safe)
for (col in group_cols) {
  gdir <- file.path(OUT, col)
  if (dir.exists(gdir) && length(list.files(gdir, pattern = "\\.bw$")) > 0) {
    message("cache hit, skip export: ", col); next }
  dir.create(gdir, showWarnings = FALSE, recursive = TRUE)
  message("export bigWigs: ", col)
  ExportGroupBW(atac, assay = "peaks", group.by = col, idents = NULL, normMethod = "RC",
                tileSize = 100000, minCells = 0, cutoff = NULL,
                chromosome = paste0("chr", 1:22), outdir = gdir, verbose = TRUE)
  file.remove(list.files(gdir, pattern = "\\.bed$", full.names = TRUE))
}

## ---- gap mask (hg38) + scaled-euclidean similarity of each track vs the truth track ----
suppressMessages({library(rtracklayer); library(GenomicRanges); library(AnnotationHub); library(lsa)})
ah <- AnnotationHub(); ahData <- AnnotationHub::query(ah, c("gap", "Homo sapiens", "hg38"))
gaps <- reduce(do.call(c, list(ahData[["AH107355"]], ahData[["AH107356"]], ahData[["AH107357"]],
                               ahData[["AH107358"]], ahData[["AH107359"]])))
gap.ranges <- GRanges(seqnames = gaps@seqnames, ranges = gaps@ranges)
gaps_rm <- function(dir, file) { bw <- import(file.path(dir, file))
  o <- findOverlaps(bw, gap.ranges); if (length(o)) bw[-queryHits(o)] else bw }
scaled_euc <- function(g1, g2) { colnames(mcols(g1)) <- "s1"; colnames(mcols(g2)) <- "s2"
  mg <- merge(g1, g2); as.numeric(dist(rbind(mg$s1, mg$s2))) / sqrt(length(mg)) }

annot_bw <- list()
for (file in list.files(file.path(OUT, "annotation")))
  annot_bw[[sub("-.*", "", file)]] <- gaps_rm(file.path(OUT, "annotation"), file)

rows <- list()
for (col in c("random", safe)) {
  gdir <- file.path(OUT, col)
  for (file in list.files(gdir)) {
    ct <- sub("-.*", "", file); if (is.null(annot_bw[[ct]])) next
    rows[[paste(col, ct)]] <- data.frame(
      Method = col, CellType = ct,
      euclidean_distance.scale = scaled_euc(gaps_rm(gdir, file), annot_bw[[ct]]))
  }
}
sim <- do.call(rbind, rows); rownames(sim) <- NULL
write.csv(sim, file.path(OUT, "peak_similarity.csv"), row.names = FALSE)

## ---- peakdist_adj = mean over cell types of (random - method)/random ; map safe tokens -> methods ----
suppressMessages(library(dplyr))
rand <- sim %>% filter(Method == "random") %>% select(CellType, rand = euclidean_distance.scale)
sim %>% filter(Method != "random") %>% merge(rand, by = "CellType") %>%
  mutate(eucdist_adj = (rand - euclidean_distance.scale) / rand) %>%
  group_by(Method) %>% summarise(peakdist_adj = mean(eucdist_adj), .groups = "drop") %>%
  mutate(method = methods_orig[Method]) %>% select(method, peakdist_adj) %>%
  write.csv(file.path(FIG2B, "peakdist.csv"), row.names = FALSE)
cat("wrote", file.path(OUT, "peak_similarity.csv"), "and", file.path(FIG2B, "peakdist.csv"), "\n")
