# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# ArchR gene scores for pbmc3k -> the ALTERNATIVE gene-activity matrix for the R2.1 sensitivity test.
# Adapted from 00_ArchR_forseacell.R, stripped to ONLY what's needed: build Arrow files from the pbmc3k
# ATAC fragments (the SAME fragments the Signac GeneActivity baseline uses) and export the STANDARD ArchR
# GeneScoreMatrix (genes x cells). Removed vs the example: the LSI variable-feature blacklist (we want the
# default full gene-score model, not restricted to variable tiles), clustering, peak calling, peak export.
library(ArchR)
library(Matrix)
set.seed(1)
addArchRThreads(threads = 8)
addArchRGenome("hg38")                     # pbmc3k multiome is GRCh38; matches EnsDb.Hsapiens.v86 (Signac baseline)

OUT        <- "/path/to/multiomeBench/pbmc/pbmc3k/benchmark/geneactivity_archr/archr"
inputFiles <- c(pbmc3k = "/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz")  # same as fig2b (needs .tbi)
proj_name  <- "ArchR_pbmc3k"
setwd(OUT)

# Arrow files WITH the default ArchR GeneScoreMatrix (no blacklist). Permissive QC so the benchmark's
# test+train ATAC cells are retained (each method intersects barcodes downstream).
ArrowFiles <- createArrowFiles(
  inputFiles      = inputFiles,
  sampleNames     = "pbmc3k",
  minTSS          = 1,
  minFrags        = 1000,
  addTileMat      = TRUE,
  addGeneScoreMat = TRUE,                  # <-- the gene-activity matrix we need
  excludeChr      = c("chrM"),
  force           = TRUE
)

proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = proj_name, copyArrows = FALSE)

# ---- export GeneScoreMatrix (genes x cells), sparse. Strip the "pbmc3k#" barcode prefix ArchR adds,
#      so barcodes match the benchmark / fig2b (<bc>-1). ----
gsm <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
m   <- assays(gsm)[["GeneScoreMatrix"]]
rownames(m) <- rowData(gsm)$name
colnames(m) <- gsub("^pbmc3k#", "", colnames(m))

exp <- file.path(OUT, "export"); dir.create(exp, showWarnings = FALSE, recursive = TRUE)
Matrix::writeMM(m, file.path(exp, "archr_gene_scores.mtx"))
write.csv(data.frame(gene = rownames(m)),    file.path(exp, "genes.csv"),    row.names = FALSE, quote = FALSE)
write.csv(data.frame(barcode = colnames(m)), file.path(exp, "barcodes.csv"), row.names = FALSE, quote = FALSE)
saveArchRProject(ArchRProj = proj)

cat(sprintf("ArchR gene scores: %d genes x %d cells -> export/archr_gene_scores.mtx\n", nrow(m), ncol(m)))
