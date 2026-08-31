# multiomeBench path configuration for R scripts.
# Every path is read via Sys.getenv() with a placeholder default, following the
# pattern already used by the original 02_peak_similarity.R.
DATA_ROOT        <- Sys.getenv("DATA_ROOT",        "/path/to/data")
PROJECT_ROOT     <- Sys.getenv("PROJECT_ROOT",     normalizePath(".."))
BENCHMARK_FUN_DIR<- Sys.getenv("BENCHMARK_FUN_DIR", file.path(PROJECT_ROOT, "common"))
PUBLISHED_REF    <- Sys.getenv("PUBLISHED_REF",    file.path(PROJECT_ROOT, "published_reference"))
TOOLS_ROOT       <- Sys.getenv("TOOLS_ROOT",       "/path/to/tools")
REF_ROOT         <- Sys.getenv("REF_ROOT",         "/path/to/references")
