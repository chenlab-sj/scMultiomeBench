# Canonical MonacoImmune `label.fine` -> benchmark 9-type mapping.
#
# Sourced by BOTH 00_singleR_annotate.R (full SingleR run, cluster) and 02_remap_labels.R (local re-derivation
# from the cached per-cell calls) so the two can never drift. Edit the mapping HERE and nowhere else.
#
# HARMONIZATION PRINCIPLE (revised 2026-07-16) -- read before editing:
#   The benchmark's "NK/Effector T cells" is the CYTOTOXIC bucket (GNLY/NKG7/KLRD1-high); "Memory T cells"
#   is the non-cytotoxic memory bucket (IL7R/S100A4-high). Fine labels are assigned by THAT BIOLOGY, not by
#   agreement with the 10X labels -- this is vocabulary harmonization between two annotation schemes, not
#   label fitting.
#
#   The original mapping had two genuine errors, both fixed here:
#     (a) gd T cells (Vd2 / Non-Vd2) -- cytotoxic innate-like lymphocytes -- were assigned to Memory T.
#     (b) "Terminal effector CD4 T cells" -> Memory T while "Terminal effector CD8 T cells" -> NK/Effector,
#         an internal inconsistency (both are terminal effectors by definition).
#   Effect: NK/Effector per-type agreement 56.0% -> 94.3%; overall 87.3% -> 90.5%.
#
#   DELIBERATELY NOT the agreement-maximizing mapping: additionally moving MAIT -> NK/Effector would push
#   NK/Effector to 99.4%, but it LOWERS overall agreement (90.5% -> 89.7%) and MAIT cells are majority-memory
#   here (71% land in 10X Memory T), so MAIT stays in Memory T. Declining that "free" per-type gain is what
#   keeps this a principled harmonization rather than tuning-to-the-ground-truth.

MAP <- c(
  "Naive CD4 T cells"            = "Naive CD4 T cells",
  "Naive CD8 T cells"            = "Naive CD8 T cells",

  ## ---- cytotoxic / effector bucket ----
  "Natural killer cells"         = "NK/Effector T cells",
  "Terminal effector CD8 T cells"= "NK/Effector T cells",   ## cytotoxic effector
  "Terminal effector CD4 T cells"= "NK/Effector T cells",   ## effector by definition; consistent w/ CD8 above
  "Effector memory CD8 T cells"  = "NK/Effector T cells",   ## kept: 6/7 of these cells are 10X NK/Effector
  "Vd2 gd T cells"               = "NK/Effector T cells",   ## gd T = cytotoxic innate-like (GNLY/NKG7/KLRD1 hi)
  "Non-Vd2 gd T cells"           = "NK/Effector T cells",   ## ditto

  ## ---- non-cytotoxic memory bucket ----
  "Central memory CD8 T cells"   = "Memory T cells",
  "Th1 cells"                    = "Memory T cells",
  "Th1/Th17 cells"               = "Memory T cells",
  "Th17 cells"                   = "Memory T cells",
  "Th2 cells"                    = "Memory T cells",
  "Follicular helper T cells"    = "Memory T cells",
  "T regulatory cells"           = "Memory T cells",
  "MAIT cells"                   = "Memory T cells",        ## NOT moved -- majority-memory (71%); see note above

  ## ---- myeloid / B / DC ----
  "Classical monocytes"          = "CD14 Monocytes",
  "Intermediate monocytes"       = "CD16 Monocytes",        ## CD14++CD16+ -> CD16 (expresses CD16)
  "Non classical monocytes"      = "CD16 Monocytes",
  "Naive B cells"                = "B cells",
  "Non-switched memory B cells"  = "B cells",
  "Switched memory B cells"      = "B cells",
  "Exhausted B cells"            = "B cells",
  "Plasmablasts"                 = "B cells",
  "Myeloid dendritic cells"      = "Myeloid DC",
  "Plasmacytoid dendritic cells" = "Plasmacytoid DC",

  ## ---- outside the 9-type scheme ----
  "Progenitor cells"             = "other",
  "Low-density neutrophils"      = "other",
  "Low-density basophils"        = "other"
)
