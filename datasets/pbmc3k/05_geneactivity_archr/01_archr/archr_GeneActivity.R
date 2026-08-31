# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
# Drop-in replacement for Signac::GeneActivity(object, features=...) that returns the precomputed
# ArchR GeneScoreMatrix instead of computing Signac gene activity -- the only change each method's
# prep needs for the R2.1 gene-activity sensitivity test: swap  GeneActivity(  ->  archr_GeneActivity(.
#
# Returns a genes x cells sparse matrix aligned to colnames(object): ArchR barcodes are matched to the
# object's cells after stripping a trailing "_atac" (the benchmark Renames ATAC cells <bc>_atac).
# Cells absent from the ArchR run (the ~2.6% dropped by minFrags) are filled with 0.
suppressMessages(library(Matrix))

archr_GeneActivity <- function(object, features = NULL,
    archr_dir = "/path/to/multiomeBench/pbmc/pbmc3k/benchmark/geneactivity_archr/archr/export") {
  m <- Matrix::readMM(file.path(archr_dir, "archr_gene_scores.mtx"))
  rownames(m) <- read.csv(file.path(archr_dir, "genes.csv"),    stringsAsFactors = FALSE)$gene
  colnames(m) <- read.csv(file.path(archr_dir, "barcodes.csv"), stringsAsFactors = FALSE)$barcode

  cells <- colnames(object)
  stems <- sub("_atac$", "", cells)                              # <bc>_atac (or <bc>) -> <bc>
  genes <- if (is.null(features)) rownames(m) else intersect(as.character(features), rownames(m))

  out <- Matrix::Matrix(0, nrow = length(genes), ncol = length(cells),
                        dimnames = list(genes, cells), sparse = TRUE)
  shared <- intersect(stems, colnames(m))                        # cells present in the ArchR run
  if (length(shared)) out[, match(shared, stems)] <- m[genes, shared, drop = FALSE]

  message(sprintf("archr_GeneActivity: %d genes x %d cells (%d/%d cells covered by ArchR; %d filled 0)",
                  length(genes), length(cells), length(shared), length(cells), length(cells) - length(shared)))
  out
}
