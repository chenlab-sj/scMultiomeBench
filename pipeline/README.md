# `pipeline/` — execution order

Five shell drivers that record **the order in which the analysis was run**. None contains analysis
logic; every line points at a script under `datasets/` or `summary/` that does the real work.

| driver | stage |
|---|---|
| `00_download.sh` | obtain the raw datasets |
| `01_preprocess.sh` | per-dataset preprocessing and per-method input preparation |
| `02_run_methods.sh` | the integration method runs, plus replicate seeds |
| `03_compute_metrics.sh` | metric chains, in numeric-prefix order |
| `04_make_figures.sh` | per-dataset figure scripts, then the cross-dataset summary figures |

## Dry run is the default

Every driver defaults to `DRY_RUN=1`, which **prints** the commands instead of running them.

```bash
./02_run_methods.sh                 # dry run, every dataset
./02_run_methods.sh rms             # dry run, one dataset
DRY_RUN=0 ./02_run_methods.sh rms   # actually execute
```

Datasets accepted as the first argument: `pbmc3k pbmc10k pbmc_parse bmmc_d1 brca rms` (plus
`summary` in stages 03 and 04); with no argument, all of them. With `DRY_RUN=0` a driver refuses to
start while `DATA_ROOT` is still the shipped placeholder — stage 02 is a multi-week compute job.

Markers in the printed output: `[dry-run]` a command that would run; `[lsf]` an original LSF
submission script (printed, never executed — see below); `[config]` a file you edit rather than
execute; `[lib]` an imported library; `[manual]` a by-hand step (stage 00); `#` a note.

## Dependency graph

```
                     $DATA_ROOT (raw matrices, fragments)
                              |
  00_download.sh -> 01_preprocess.sh   per-dataset objects, common peaks, labels,
                              |        per-method prepared inputs
                              v
  02_run_methods.sh  one latent embedding per method x dataset (x replicate seed)
                              |
                              v
  03_compute_metrics.sh  per-dataset metric tables
                              |
                              v
  04_make_figures.sh   per-dataset panels + summary/ cross-dataset figures
```

Stage 04 is cheap: the figures regenerate in minutes from the tables shipped under `results/` and
need no method environment. Stages 02 and 03 are the expensive ones; only stage 02 needs a GPU.

## Environments

There is no single environment that runs all methods — they carry mutually incompatible pins. The
drivers call a bare `python` / `Rscript`; activate the right environment per method before running
its block. Roughly: legacy-Python methods (MIRA, scButterfly, MaxFuse), modern-Python GPU methods
(MIDAS, Multigrate, scVI/MultiVI, scGLUE, Cobolt, Portal, simba, scDART, scJoint), and the R stack
(BindSC, Seurat CCA/WNN, Conos, LIGER, peak similarity). See `envs/` and README section 6.

## `*_sub.sh` files are the original LSF submitters

They record how each job was actually resourced (queue, walltime, memory, GPU) at the originating
cluster. They are site-specific and will not run elsewhere — the drivers print them under `[lsf]`
and never execute them; port them to your own scheduler or run the underlying `.py` / `.R` directly.

## Notes

- The two large RMS pileup tables (`predicted_ataclabel_{foxo1,myod1}value.csv`) are gitignored for
  size; Figs 5B/5C/S7-S9 need them — available from the authors on request.
- The published Fig S7 page was assembled manually from the three per-locus strips.
- BRCA raw data is controlled access (dbGaP `phs002371.v3.p1`); every other dataset is open.
