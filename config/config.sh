#!/bin/bash
# multiomeBench path configuration (shell + R scripts read these).
# Copy config.local.sh.example -> config.local.sh and set your real paths there.
# config.local.sh is gitignored and must never be committed.

# Raw and intermediate data (h5ad/h5/fragments). See DATA.md for how to obtain it.
: "${DATA_ROOT:=/path/to/data}"

# This repository's own root.
: "${PROJECT_ROOT:=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

# Vendored helper libraries (benchmark_fun.py, export_groupbwg.R, plot_regionpeak_fun.R)
: "${BENCHMARK_FUN_DIR:=${PROJECT_ROOT}/common}"

# Published reference values spliced into some figures (see published_reference/README.md)
: "${PUBLISHED_REF:=${PROJECT_ROOT}/published_reference}"

# Cloned method repositories (scJoint, scBridge, scMoMaT were run from git clones,
# not installed packages). See README for the commits used.
: "${TOOLS_ROOT:=/path/to/tools}"

# Genome annotation / reference files. NOTE: RMS uses hg19 (GRCh37); the other
# datasets use hg38 (GRCh38-2020-A).
: "${REF_ROOT:=/path/to/references}"

[ -f "$(dirname "${BASH_SOURCE[0]}")/config.local.sh" ] && . "$(dirname "${BASH_SOURCE[0]}")/config.local.sh"
export DATA_ROOT PROJECT_ROOT BENCHMARK_FUN_DIR PUBLISHED_REF TOOLS_ROOT REF_ROOT
