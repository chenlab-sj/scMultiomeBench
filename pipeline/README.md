# `pipeline/` — execution order

Five shell drivers that record **the order in which the analysis was run**. They are
documentation-by-execution-order: none of them contains analysis logic. Every line points at a
script under `datasets/` or `summary/` that does the real work.

| driver | stage |
|---|---|
| `00_download.sh` | obtain the raw datasets |
| `01_preprocess.sh` | per-dataset preprocessing and per-method input preparation |
| `02_run_methods.sh` | the 14 integration methods (23 on PBMC 3k), plus replicate seeds |
| `03_compute_metrics.sh` | metric chains, in numeric-prefix order |
| `04_make_figures.sh` | per-dataset figure scripts, then the cross-dataset summary figures |

The script lists were generated from the repository census, not typed by hand, so they cover every
executable file that survives into the public tree.

## Dry run is the default

Every driver defaults to `DRY_RUN=1`, which **prints** the commands instead of running them.
Nothing is downloaded, computed, written, or submitted.

```
./02_run_methods.sh                 # dry run, every dataset
./02_run_methods.sh rms             # dry run, one dataset
DRY_RUN=0 ./02_run_methods.sh rms   # actually execute
```

Datasets accepted as the first argument: `pbmc3k pbmc10k pbmc_parse bmmc_d1 brca rms`, plus
`summary` in stages 03 and 04. With no argument, all of them.

With `DRY_RUN=0` a driver refuses to start while `DATA_ROOT` is still the shipped placeholder
`/path/to/data`. Stage 02 in particular is a multi-week compute job; the dry-run default exists so
that it cannot be launched by reflex.

### Markers in the printed output

| marker | meaning |
|---|---|
| `[dry-run]` | a command that would be executed |
| `[lsf]` | an original LSF submission script — printed, never run (see below) |
| `[config]` | a file you edit rather than execute (e.g. scJoint's `config.py`) |
| `[lib]` | an imported library, not an entry point |
| `[manual]` | stage 00 only: a step with no script; do it by hand from the stated accession |
| `#` | a note |

## Dependency graph

```
                     $DATA_ROOT (raw matrices, fragments)
                              |
  00_download.sh  ------------+
                              |
                              v
  01_preprocess.sh   per-dataset RNA/ATAC objects, common-peak matrices,
                     ground-truth cell-type labels, and the method-specific
                     inputs (gene-activity matrices, region2gene links,
                     scJoint/scBridge h5)
                              |
                              v
  02_run_methods.sh  one latent embedding per method x dataset (x replicate seed)
                              |
                              v
  03_compute_metrics.sh  per-dataset metric tables
                     (KNN label transfer k=10, cosine, distance-weighted;
                      accuracy adjustment; peak-similarity)
                              |
                              +--------------------------+
                              v                          v
  04_make_figures.sh   per-dataset panels        summary/ cross-dataset
                       (Fig2-Fig6, FigS1-S12)    figures (Fig7, FigS2A, FigS13)
                              ^
                              |
                     results/<ds>/published_reference/  (two panels only: pbmc3k Fig2B,
                                            BRCA FigS6)
```

What each stage hands to the next:

| stage | consumes | produces |
|---|---|---|
| 00 | dbGaP / GEO / vendor accessions | raw files under `$DATA_ROOT` |
| 01 | `$DATA_ROOT` | per-dataset objects, common peaks, ground-truth labels, per-method prepared inputs |
| 02 | stage 01 outputs | latent embeddings, one per method per dataset per seed |
| 03 | stage 02 embeddings + stage 01 labels | metric tables (composite score inputs) |
| 04 | stage 03 tables (+ `results/<dataset>/published_reference/` for two panels) | manuscript figures |

Stage 04 is cheap: given the stage-03 metric tables it runs in minutes and needs no method
environment. `results/` holds the collected copies of those tables; the plotting scripts themselves
read each dataset's own `03_metrics/` working directory. Stages 02 and 03 are the expensive ones,
and only stage 02 needs a GPU.

## Script counts per driver

`run` = a step the driver would execute. `not run` = `[lsf]` submitters, `[config]` files and
`[lib]` libraries that are listed for completeness only.

| dataset | 01 preprocess | 02 methods | 03 metrics | 04 figures |
|---|---|---|---|---|
| pbmc3k | 14 run | 52 run + 74 not run | 11 run + 12 not run | 8 run + 6 not run |
| pbmc10k | 1 run | 27 run + 21 not run | 5 run + 4 not run | 1 run + 1 not run |
| pbmc_parse | 5 run | 18 run + 20 not run | 5 run + 4 not run | 3 run + 2 not run |
| bmmc_d1 | 2 run | 50 run + 64 not run | 15 run + 12 not run | 7 run + 2 not run |
| brca | 7 run + 4 not run | 40 run + 55 not run | 8 run + 8 not run | 5 run + 1 not run |
| rms | 5 run + 1 not run | 13 run + 38 not run | 8 run + 7 not run | 8 run + 1 not run |
| summary | — | — | 1 run | 4 run + 1 not run |

BMMC's stage 02 count includes the six `(batch)` variants — `BindSC_batch`, `MIDAS_batch`,
`scJoint_batch`, `scVI_batch`, `scglue_batch`, `scglue_multiome_batch` — which feed Fig 6A and
Table S4, plus the `05_crosssite/` and `06_s1d1_paired/` sub-analyses. BRCA's and PBMC 3k's counts
include the `reps/` replicate-seed folders behind the Fig 4A/B reproducibility panels.

## Method runs need separate conda environments

**There is no single environment that runs all 14 methods.** They carry mutually incompatible pins
(torch, scvi-tools, numpy, networkx), and several are only installable at all on an old Python. The
drivers call a bare `python` / `Rscript`; you are expected to activate the correct environment
yourself, one method at a time, before letting stage 02 run that method's block.

Concretely, the environments split at least three ways:

| kind | examples |
|---|---|
| legacy Python + old torch | MIRA, scButterfly, BABEL-style translation methods, MaxFuse |
| modern Python + recent torch / scvi-tools, GPU | MIDAS, Multigrate, scVI/MultiVI, scGLUE, Cobolt, Portal, simba, scDART, scJoint |
| R (Seurat / Signac stack) | BindSC, Seurat CCA, Seurat WNN, Conos, LIGER, and every `*_peak_similarity.R` |

See `envs/` for the recipes. Which runs needed a GPU is recorded in the corresponding `*_sub.sh`
(for example the `*_subgpu.sh` submitters): the deep-learning methods in the second row above. The R
methods and all of stages 03 and 04 are CPU-only.

## `*_sub.sh` files are the original LSF submitters

Every `*_sub.sh` (and `run_all.sh` / `run_methods.sh`) inside the dataset folders is an **original
LSF submission script from the site where this benchmark was run** (St. Jude). They carry `bsub`
lines, queue names, GPU reservations, walltimes and absolute cluster paths.

They are shipped because they record how each job was actually resourced — queue, walltime, memory
and whether a GPU was reserved. They are **site-specific and will not run anywhere else**, so the
drivers print them under `[lsf]` and never execute them. Port them to your own scheduler, or run the underlying `.py` / `.R` script directly.

## Paths and configuration

`config/config.sh` defines `DATA_ROOT`, `PROJECT_ROOT`, `BENCHMARK_FUN_DIR` and `PUBLISHED_REF`,
with placeholder defaults. Copy `config/config.local.sh.example` to `config/config.local.sh`
(gitignored) and set your real paths. R scripts read the same variables via `Sys.getenv()`, matching
the pattern `peak_similarity.R` already used. Each driver sources `config/config.sh` on startup.

Absolute cluster paths have been scrubbed from every script. Shell scripts pick the config up
automatically; `.py` and `.R` scripts instead carry literal `/path/to/...` placeholders at the top of
the file, which you edit (or override via the environment variables the script reads with
`Sys.getenv()`) before running one by hand. Do not assume a script is portable because the driver
that calls it sources the config.

## Known gaps

- **Large per-cell pileup inputs are not in the repository.** `predicted_ataclabel_foxo1value.csv`
  (88.7 MB) and `predicted_ataclabel_myod1value.csv` (18.1 MB) are gitignored for size and are
  required by the RMS Fig S8/S9 (and Fig 5C) pileup panels; the Fig S7 strips additionally need the
  SJRHB013758_X2-side value tables, which are likewise not committed. All are available from the
  authors on request.
- **Fig S7 assembly is manual.** The three per-locus strip scripts under
  `datasets/rms/04_figures/figS7_multiome_vs_annotation/` draw the panels; combining them into the
  published page was done by hand, so no single script emits the composed figure.
- **BRCA raw data is controlled access.** dbGaP `phs002371.v3.p1`; stages 01–04 for `brca` cannot be
  run without an approved data access request. Every other dataset is open.
