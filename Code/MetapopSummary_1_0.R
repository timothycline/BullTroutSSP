# ================================================================
# Metapopulation Summary v1.0
# Generates Results/MetapopSummary_D{data}_M{model}.RDS
#
# Computes core-area level viability metrics from:
#   - fit_rep ADREPORT logN (current abundance)
#   - fit_rep ADREPORT b0r_core (core-level growth rate)
#   - logN_proj from Qext RDS (metapopulation quasi-extinction)
#   - prioritization data RDS (per-pop λ status)
# ================================================================

rm(list = ls())

library(TMB)
library(here)
library(tidyverse)

# ── Configuration ─────────────────────────────────────────────────
DataVersion     <- "1_8"
ModelVersion    <- "1_22"
RunVersion      <- "1_21"
FunctionVersion <- "1_13"

SaveName <- paste0("MPVA_Fit_D", DataVersion, "_M", ModelVersion, "_R", RunVersion)
OutName  <- paste0("MetapopSummary_D", DataVersion, "_M", ModelVersion)

# ── Load fit ──────────────────────────────────────────────────────
cat("Loading fit:", SaveName, "\n")
fit <- readRDS(here("Results", paste0(SaveName, ".RDS")))
list2env(fit, envir = .GlobalEnv)

source(here("Code", paste0("MPVA_Functions", FunctionVersion, ".R")))
suppressMessages(
  source(here("Code", paste0("PrepareData_MPVA", DataVersion, ".R")))
)

nPop   <- data_list$nPop
nCore  <- data_list$nCore
nYears <- data_list$nYears
pop_core_vec <- data_list$pop_core + 1L   # 1-indexed

cat("nPop:", nPop, "| nCore:", nCore, "| nYears:", nYears, "\n\n")

# ── Load Qext RDS ────────────────────────────────────────────────
cat("Loading Qext RDS...\n")
qext     <- readRDS(here("Results", paste0(SaveName, "_Qext.RDS")))
logN_proj <- qext$logN_proj    # nsim x nPop x nYears_proj
nsim      <- qext$nsim
nYears_proj <- qext$nYears_proj
proj_years  <- qext$proj_years
cat("logN_proj:", paste(dim(logN_proj), collapse=" x "), "\n\n")

# ── Load prioritization data (per-pop λ) ──────────────────────────
pdat <- readRDS(here("Results",
                     paste0("PrioritizationData_D", DataVersion, "_M", ModelVersion, ".RDS")))
pop_lambda  <- pdat$status_table$lambda_med    # median λ per pop
pop_qext20  <- pdat$status_table$qext_p20

# ── Population display names ──────────────────────────────────────
wb_raw      <- sapply(PatchWB_In, function(x) strsplit(x, "\\.")[[1]][2])
wb_readable <- trimws(gsub("([a-z0-9'])([A-Z])", "\\1 \\2", wb_raw))
pop_display <- sapply(seq_len(nPop), function(p) {
  wb_readable[which(LocalPops_Num == p)[1]]
})

# ── Extract ADREPORT logN and b0r_core ────────────────────────────
logN_vals   <- fit_rep$value[names(fit_rep$value) == "logN"]
logN_se     <- fit_rep$sd[names(fit_rep$value)    == "logN"]
logN_mat    <- matrix(logN_vals, nrow = nPop, ncol = nYears)   # nPop x nYears
logN_se_mat <- matrix(logN_se,   nrow = nPop, ncol = nYears)

b0r_core_vals <- fit_rep$value[names(fit_rep$value) == "b0r_core"]  # nCore
b0r_core_se   <- fit_rep$sd[names(fit_rep$value)    == "b0r_core"]

cat("ADREPORT extracted.\n")
cat("  logN last-year range:", round(range(logN_mat[, nYears]), 2), "\n")
cat("  b0r_core range:", round(range(b0r_core_vals), 3), "\n\n")

# ── Per-population current abundance (last data year) ─────────────
N_last_med <- exp(logN_mat[, nYears])
N_last_lo  <- exp(logN_mat[, nYears] - 1.645 * logN_se_mat[, nYears])
N_last_hi  <- exp(logN_mat[, nYears] + 1.645 * logN_se_mat[, nYears])

# ── Core-area summaries ───────────────────────────────────────────
cat("Computing core-area summaries...\n")

# Metapop Qext: P(total core N < threshold at any point in projection)
# Using joint simulation draws — accounts for correlated population dynamics
N_proj_arr <- exp(logN_proj)   # nsim x nPop x nYears_proj

metapop_qext <- sapply(seq_len(nCore), function(c) {
  pops  <- which(pop_core_vec == c)
  if (length(pops) == 0) return(NA_real_)
  # Total N across all populations in core, per simulation, per year
  N_core <- apply(N_proj_arr[, pops, , drop = FALSE], c(1, 3), sum)  # nsim x nYears_proj
  # Threshold: 20% of current median core N (population-size-scaled threshold)
  # Also compute at fixed thresholds for comparability
  mean(apply(N_core < 20 * length(pops), 1, any))   # Qext at 20 * n_pops threshold
})

# More useful: P(majority of pops < 20) at any projection year
metapop_qext_majority <- sapply(seq_len(nCore), function(c) {
  pops  <- which(pop_core_vec == c)
  if (length(pops) == 0) return(NA_real_)
  n_pops <- length(pops)
  # For each sim × year, count how many pops < 20
  N_mat <- N_proj_arr[, pops, , drop = FALSE]  # nsim x n_pops x nYears_proj
  n_below <- apply(N_mat < 20, c(1, 3), sum)   # nsim x nYears_proj
  # P(majority < 20 at any point)
  mean(apply(n_below > n_pops / 2, 1, any))
})

core_summary <- data.frame(
  core          = U_CoreAreas_In,
  n_pops        = tabulate(pop_core_vec, nbins = nCore),
  # Connected habitat
  hab_km        = sapply(seq_len(nCore), function(c)
                    sum(data_list$extent[pop_core_vec == c])),
  # Current abundance
  N_med         = sapply(seq_len(nCore), function(c)
                    sum(N_last_med[pop_core_vec == c])),
  N_lo90        = sapply(seq_len(nCore), function(c)
                    sum(N_last_lo[pop_core_vec == c])),
  N_hi90        = sapply(seq_len(nCore), function(c)
                    sum(N_last_hi[pop_core_vec == c])),
  # Core-level intrinsic growth rate
  lambda_core   = round(exp(b0r_core_vals), 3),
  lambda_lo90   = round(exp(b0r_core_vals - 1.645 * b0r_core_se), 3),
  lambda_hi90   = round(exp(b0r_core_vals + 1.645 * b0r_core_se), 3),
  # Population status breakdown
  n_growing     = sapply(seq_len(nCore), function(c)
                    sum(pop_lambda[pop_core_vec == c] >= 1, na.rm = TRUE)),
  n_declining   = sapply(seq_len(nCore), function(c)
                    sum(pop_lambda[pop_core_vec == c] < 1, na.rm = TRUE)),
  # Metapopulation quasi-extinction
  metapop_qext  = round(metapop_qext, 3),
  metapop_qext_maj = round(metapop_qext_majority, 3)
) %>%
  mutate(
    density_N_km  = round(N_med  / hab_km, 2),
    density_lo_km = round(N_lo90 / hab_km, 2),
    density_hi_km = round(N_hi90 / hab_km, 2),
    pct_declining = round(100 * n_declining / n_pops, 0)
  )

cat("Core summary assembled.\n")
print(select(core_summary, core, n_pops, density_N_km, lambda_core, metapop_qext))

# ── Per-population table (for pop-level breakdown plot) ───────────
pop_table <- data.frame(
  pop_id    = U_LocalPops_In,
  pop_name  = pop_display,
  core      = U_CoreAreas_In[pop_core_vec],
  hab_km    = data_list$extent,
  N_med     = round(N_last_med, 1),
  N_lo90    = round(N_last_lo,  1),
  N_hi90    = round(N_last_hi,  1),
  lambda    = round(pop_lambda, 3),
  qext_p20  = round(pop_qext20, 3)
) %>%
  mutate(
    density_N_km = round(N_med / hab_km, 2),
    status       = ifelse(lambda >= 1, "Growing", "Declining")
  )

# ── Save ──────────────────────────────────────────────────────────
out <- list(
  core_summary  = core_summary,
  pop_table     = pop_table,
  meta = list(
    DataVersion  = DataVersion,
    ModelVersion = ModelVersion,
    RunVersion   = RunVersion,
    DatYears     = DatYears,
    nYears       = nYears,
    nPop         = nPop,
    nCore        = nCore,
    nYears_proj  = nYears_proj,
    proj_years   = proj_years,
    run_date     = Sys.Date()
  )
)

out_path <- here("Results", paste0(OutName, ".RDS"))
saveRDS(out, out_path)
cat("\nSaved:", out_path, "\n")
