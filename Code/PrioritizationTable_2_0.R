# ================================================================
# Prioritization Data Prep v2.0
# Generates:
#   Results/PrioritizationData_D{data}_M{model}.RDS  (Quarto page)
#   ShinyApp/scenario_app/prioritization_data.RDS    (Shiny app)
#
# Requires TMB (for summary.sdreport dispatch).
# Run once after each new model fit to update both products.
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

nsamp <- 1000   # marginal posterior draws for scenario engine
seed  <- 42
set.seed(seed)

SaveName <- paste0("MPVA_Fit_D", DataVersion, "_M", ModelVersion, "_R", RunVersion)
OutName  <- paste0("PrioritizationData_D", DataVersion, "_M", ModelVersion)

# ── Load fit ──────────────────────────────────────────────────────
cat("Loading fit:", SaveName, "\n")
fit <- readRDS(here("Results", paste0(SaveName, ".RDS")))
list2env(fit, envir = .GlobalEnv)

source(here("Code", paste0("MPVA_Functions", FunctionVersion, ".R")))
suppressMessages(
  source(here("Code", paste0("PrepareData_MPVA", DataVersion, ".R")))
)

nPop  <- data_list$nPop
nCore <- data_list$nCore
nYears <- data_list$nYears
pop_core_vec <- data_list$pop_core + 1L
has_LKT_pop  <- data_list$has_LKT[pop_core_vec]   # 0/1 vector, length nPop

cat("nPop:", nPop, "| nCore:", nCore, "| nYears:", nYears, "\n")

# ── Extract MAP estimates + SEs from sdreport ──────────────────────
re_summ <- summary(fit_rep, select = "random")

# Helper: extract random effect vector (Estimate + SE)
get_re <- function(nm) {
  idx <- rownames(re_summ) == nm
  list(est = re_summ[idx, "Estimate"],
       se  = re_summ[idx, "Std. Error"])
}

# Helper: extract fixed effect scalar
get_fe <- function(nm) {
  list(est = fit_rep$value[names(fit_rep$value) == nm],
       se  = fit_rep$sd[names(fit_rep$value)    == nm])
}

b0r      <- get_re("b0r_pop")
bFlow    <- get_re("bFlow_pop")
bTemp    <- get_re("bTemp_pop")
bINV     <- get_re("bINV_pop")
bLKT     <- get_fe("bLKT_M")
log_bT2  <- get_fe("log_bTemp2")   # bTemp2 = -exp(log_bTemp2); may be absent if mapped out

cat("Parameters extracted from sdreport.\n")
cat("  b0r_pop  : n =", length(b0r$est),  " | median =", round(median(b0r$est), 3), "\n")
cat("  bINV_pop : n =", length(bINV$est), " | median =", round(median(bINV$est), 3), "\n")
cat("  bLKT_M   : est =", round(bLKT$est, 3), "\n\n")

# ── Draw marginal posterior samples ──────────────────────────────
# Independent Normal(mu, sigma) draws per parameter.
# Loses parameter cross-correlation but is fast and TMB-free.
# Adequate for scenario comparison and display of uncertainty.

# Each column p gets n independent draws from N(est[p], se[p]).
# Replaces NA SEs with 0 (MAP estimate only, no uncertainty for that parameter).
draw_mat <- function(est, se, n) {
  se_safe <- ifelse(is.na(se) | se <= 0, 0, se)
  sapply(seq_along(est), function(p) rnorm(n, est[p], se_safe[p]))
}

b0r_samp    <- draw_mat(b0r$est,   b0r$se,   nsamp)   # nsamp x nPop
bFlow_samp  <- draw_mat(bFlow$est, bFlow$se, nsamp)   # nsamp x nPop
bTemp_samp  <- draw_mat(bTemp$est, bTemp$se, nsamp)   # nsamp x nPop
bINV_samp   <- draw_mat(bINV$est,  bINV$se,  nsamp)   # nsamp x nPop
bLKT_samp   <- rnorm(nsamp, bLKT$est, bLKT$se)        # nsamp (global scalar)

# bTemp2 drawn on log scale then transformed to preserve sign constraint.
# Falls back to 0 if log_bTemp2 was mapped out or absent in this model version.
if (length(log_bT2$est) > 0 && !all(is.na(log_bT2$se))) {
  se_bT2       <- ifelse(is.na(log_bT2$se), 0, log_bT2$se)
  log_bT2_samp <- rnorm(nsamp, log_bT2$est, se_bT2)
  bTemp2_samp  <- -exp(log_bT2_samp)                  # nsamp, always negative
} else {
  cat("  log_bTemp2 not found in fit_rep — bTemp2 fixed at 0.\n")
  bTemp2_samp  <- rep(0, nsamp)
}

cat("Drew", nsamp, "marginal posterior samples.\n\n")

# ── Covariate current-year and historical values ──────────────────
# NA covariate values (e.g., missing LKT data for non-lake-connected cores,
# or missing flow/temp in a given year) are replaced with 0 = historical mean.
# This prevents NA propagation in lambda_quants while remaining conservative.
safe0 <- function(x) ifelse(is.na(x), 0, x)

INV_last       <- safe0(data_list$INV[,  nYears])      # nPop, z-scored
Temp_last      <- safe0(data_list$Temp[, nYears])      # nPop, z-scored
Flow_last      <- safe0(data_list$Flow[, nYears])      # nPop, z-scored
LKT_last_core  <- data_list$LKT[,  nYears]             # nCore, z-scored (may have NAs)
LKT_last_pop   <- safe0(LKT_last_core[pop_core_vec])   # nPop

# Historical range (for slider bounds in Shiny)
INV_range  <- round(range(data_list$INV,  na.rm = TRUE), 1)
Temp_range <- round(range(data_list$Temp, na.rm = TRUE), 1)
Flow_range <- round(range(data_list$Flow, na.rm = TRUE), 1)
LKT_range  <- round(range(data_list$LKT,  na.rm = TRUE), 1)

# ── Core lambda computation helper ────────────────────────────────
# Compute posterior λ quantiles given per-pop z-scored covariate vectors.
# inv_z, temp_z, flow_z: length-nPop vectors; lkt_z: length-nPop vector.
lambda_quants <- function(inv_z, temp_z, flow_z, lkt_z,
                          probs = c(0.05, 0.5, 0.95)) {
  nP   <- nPop
  nsmp <- nsamp

  log_lam <- b0r_samp

  if (covar_config$use_Flow)
    log_lam <- log_lam + bFlow_samp * matrix(flow_z, nsmp, nP, byrow = TRUE)

  if (covar_config$use_Temp)
    log_lam <- log_lam +
      bTemp_samp * matrix(temp_z, nsmp, nP, byrow = TRUE) +
      matrix(bTemp2_samp, nsmp, nP) * matrix(temp_z^2, nsmp, nP, byrow = TRUE)

  if (covar_config$use_INV)
    log_lam <- log_lam + bINV_samp * matrix(inv_z, nsmp, nP, byrow = TRUE)

  if (covar_config$use_LKT)
    log_lam <- log_lam +
      matrix(bLKT_samp, nsmp, nP) *
      matrix(has_LKT_pop * lkt_z, nsmp, nP, byrow = TRUE)

  # Return matrix: nPop x length(probs)
  t(apply(exp(log_lam), 2, quantile, probs = probs, na.rm = TRUE))
}

# ── Compute status table quantities ───────────────────────────────
cat("Computing status table...\n")

# Current conditions (last observed year)
lam_current    <- lambda_quants(INV_last, Temp_last, Flow_last, LKT_last_pop)

# Sensitivity scenarios — one lever at a time, all others at current
# INV at historical mean (z=0): simulates average invasive pressure
lam_inv_mean   <- lambda_quants(rep(0, nPop), Temp_last, Flow_last, LKT_last_pop)
# INV at historical min: most favorable invasive condition
lam_inv_min    <- lambda_quants(rep(INV_range[1], nPop), Temp_last, Flow_last, LKT_last_pop)
# Temp 1 SD cooler: beneficial
lam_temp_1sd   <- lambda_quants(INV_last, Temp_last - 1, Flow_last, LKT_last_pop)
# Flow 1 SD higher: beneficial
lam_flow_1sd   <- lambda_quants(INV_last, Temp_last, Flow_last + 1, LKT_last_pop)
# LKT at historical min: most favorable (least lake trout)
lam_lkt_min    <- lambda_quants(INV_last, Temp_last, Flow_last,
                                rep(LKT_range[1], nPop))

cat("  Done.\n\n")

# ── Load Qext probabilities (optional) ───────────────────────────
qext_file <- here("Results", paste0(SaveName, "_Qext.RDS"))
if (file.exists(qext_file)) {
  qext <- readRDS(qext_file)
  qext_p20  <- qext$qext_results$thr20$p_ever
  qext_p10  <- qext$qext_results$thr10$p_ever
  cat("Loaded Qext RDS — Pr(N < 20) and Pr(N < 10) available.\n\n")
} else {
  qext_p20 <- rep(NA_real_, nPop)
  qext_p10 <- rep(NA_real_, nPop)
  cat("No Qext RDS found — Qext columns will be NA.\n\n")
}

# ── Population display names ──────────────────────────────────────
wb_raw      <- sapply(PatchWB_In, function(x) strsplit(x, "\\.")[[1]][2])
wb_readable <- trimws(gsub("([a-z0-9'])([A-Z])", "\\1 \\2", wb_raw))
pop_display <- sapply(seq_len(nPop), function(p) {
  first_site <- which(LocalPops_Num == p)[1]
  wb_readable[first_site]
})

# ── Assemble status table ─────────────────────────────────────────
status_table <- data.frame(
  pop_id          = U_LocalPops_In,
  pop_name        = pop_display,
  core            = U_CoreAreas_In[pop_core_vec],
  has_lkt         = has_LKT_pop,
  # Current λ (last year conditions)
  lambda_med      = round(lam_current[, 2], 3),
  lambda_lo90     = round(lam_current[, 1], 3),
  lambda_hi90     = round(lam_current[, 3], 3),
  # Quasi-extinction
  qext_p20        = round(qext_p20, 3),
  qext_p10        = round(qext_p10, 3),
  # Sensitivity: λ under best-case for each lever
  lam_inv_mean    = round(lam_inv_mean[,  2], 3),   # INV at historical avg
  lam_inv_best    = round(lam_inv_min[,   2], 3),   # INV at historical min
  lam_temp_1sd    = round(lam_temp_1sd[,  2], 3),   # Temp -1 SD
  lam_flow_1sd    = round(lam_flow_1sd[,  2], 3),   # Flow +1 SD
  lam_lkt_best    = ifelse(has_LKT_pop == 1,
                           round(lam_lkt_min[, 2], 3), NA_real_)
)

cat("Status table: nrow =", nrow(status_table), "\n")
cat("  λ range:", round(range(status_table$lambda_med, na.rm = TRUE), 3), "\n")
cat("  λ NAs  :", sum(is.na(status_table$lambda_med)), "\n")
cat("  Qext > 0.1:", sum(status_table$qext_p20 > 0.1, na.rm = TRUE), "pops\n\n")

# ── Assemble output ───────────────────────────────────────────────
out <- list(

  # Pre-computed status table for reactable
  status_table = status_table,

  # Posterior samples for live scenario computation (pure R matrix ops in app)
  samples = list(
    b0r_pop   = b0r_samp,    # nsamp x nPop
    bFlow_pop = bFlow_samp,  # nsamp x nPop
    bTemp_pop = bTemp_samp,  # nsamp x nPop
    bINV_pop  = bINV_samp,   # nsamp x nPop
    bLKT_M    = bLKT_samp,   # nsamp (global scalar)
    bTemp2    = bTemp2_samp  # nsamp (global scalar, always < 0)
  ),

  # Covariate current values and historical ranges
  cov_info = list(
    INV_last    = INV_last,         # nPop, z-scored last year
    Temp_last   = Temp_last,
    Flow_last   = Flow_last,
    LKT_last    = LKT_last_pop,     # nPop (mapped from core)
    has_LKT_pop = has_LKT_pop,
    INV_range   = INV_range,
    Temp_range  = Temp_range,
    Flow_range  = Flow_range,
    LKT_range   = LKT_range
  ),

  # Model metadata
  meta = list(
    DataVersion  = DataVersion,
    ModelVersion = ModelVersion,
    RunVersion   = RunVersion,
    DatYears     = DatYears,
    nYears       = nYears,
    nPop         = nPop,
    nCore        = nCore,
    nsamp        = nsamp,
    pop_names    = U_LocalPops_In,
    pop_display  = pop_display,
    core_names   = U_CoreAreas_In,
    pop_core_vec = pop_core_vec,
    has_LKT_pop  = has_LKT_pop,
    covar_config = covar_config,
    run_date     = Sys.Date()
  )
)

# ── Save to Results/ (for prioritization.qmd) ────────────────────
out_results <- here("Results", paste0(OutName, ".RDS"))
saveRDS(out, out_results)
cat("Saved:", out_results, "\n")
cat("File size:", round(file.size(out_results) / 1e6, 1), "MB\n\n")

# ── Copy to Shiny app directory ───────────────────────────────────
shiny_dir <- here("ShinyApp", "scenario_app")
if (!dir.exists(shiny_dir)) dir.create(shiny_dir, recursive = TRUE)
out_shiny <- file.path(shiny_dir, "prioritization_data.RDS")
file.copy(out_results, out_shiny, overwrite = TRUE)
cat("Copied to Shiny app dir:", out_shiny, "\n")
