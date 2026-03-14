# ================================================================
# PrepareShinyData.R
# Run once after each model update to generate ShinyApp/app_data.RDS
# Does NOT require the fitted model — only PrepareData
# ================================================================

rm(list = ls())

library(here)
library(tidyverse)

source(here("Code", "MPVA_Functions1_4.R"))
suppressMessages(source(here("Code", "PrepareData_MPVA1_2.R")))

# ── Population-level core index ───────────────────────────────────────────────
pop_core_vec <- CoreArea_Num[match(seq_len(nPop), LocalPops_Num)]

# ── Display names from PatchWB_In ("BLTPatchID.WaterbodyName") ────────────────
wb_raw      <- sapply(PatchWB_In, function(x) strsplit(x, "\\.")[[1]][2])
wb_readable <- trimws(gsub("([a-z0-9'])([A-Z])", "\\1 \\2", wb_raw))

# Population display: "Boulder Creek (09040001BLD)"
pop_display <- sapply(seq_len(nPop), function(p) {
  first <- which(LocalPops_Num == p)[1]
  paste0(wb_readable[first], " (", U_LocalPops_In[p], ")")
})

# Site display: "Boulder Creek (09040001BLD.BoulderCreek)"
site_id      <- sapply(PatchWB_In, function(x) strsplit(x, "\\.")[[1]][1])
site_display <- paste0(wb_readable, " (", site_id, ")")

# ── Save ──────────────────────────────────────────────────────────────────────
dir.create(here("ShinyApp"), showWarnings = FALSE)

saveRDS(
  list(
    n_in           = n_in,
    DatYears       = DatYears,
    PatchWB_In     = PatchWB_In,
    site_display   = site_display,
    LocalPops_Num  = LocalPops_Num,
    U_LocalPops_In = U_LocalPops_In,
    pop_display    = pop_display,
    pop_core_vec   = pop_core_vec,
    U_CoreAreas_In = U_CoreAreas_In
  ),
  here("ShinyApp", "app_data.RDS")
)

cat("Saved ShinyApp/app_data.RDS\n")
cat("  Populations:", nPop, "\n")
cat("  Core areas: ", nCore, "\n")
cat("  Sites:      ", nrow(n_in), "\n")
cat("  Years:      ", min(DatYears), "-", max(DatYears), "\n")
