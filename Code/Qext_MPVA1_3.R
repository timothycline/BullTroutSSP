# ================================================================
# MPVA TMB — Quasi-Extinction Probability  (simplified v1.3)
# Qext_MPVA1_3.R
#
# Simplified from v1.2: no TMB rebuild, no joint precision matrix,
# no sparse Cholesky.  Uses only the saved fit RDS.
#
# Sampling strategy (two-tier):
#   Fixed outer params  (log_sigmaN, logit_rho, bLKT_M, bTemp2,
#                        sigma hyperparameters, …):
#     Draw from N(par.fixed, cov.fixed) via small Cholesky.
#     cov.fixed is the marginal covariance of these ~15 params —
#     always well-conditioned.
#   Random effects (b0r_pop, log_b0phi_pop, bFlow_pop, bTemp_pop,
#                   bINV_pop):
#     Draw independently from N(MAP, marginal_SE^2).
#     Captures individual-population uncertainty; ignores
#     cross-parameter correlations (b0r × b0phi etc.), which are
#     a second-order effect for 7-yr projections.
#
# Projection starts from MAP logN at the last data year.
# Historical bands use marginal SEs on logN from sdreport
# (lognormal approximation: CI = exp(logN_map ± z * logN_sd)).
#
# Outputs match v1.2 so RiskComposite_1_0.R works unchanged.
# ================================================================

rm(list = ls())

library(here)
library(tidyverse)

# ================================================================
# CONFIGURATION
# ================================================================
DataVersion  <- '1_8'
ModelVersion <- '1_22'
RunVersion   <- '1_21'

SaveName <- paste0('MPVA_Fit_D', DataVersion, '_M', ModelVersion, '_R', RunVersion)

nsim        <- 5000
nYears_proj <- 7
thresholds  <- c(5, 10, 20, 50)
seed        <- 42

future_cov_assumption <- 'mean'   # 'mean' (z-score 0) or 'last'

# ================================================================
# LOAD FIT
# ================================================================
cat("Loading:", SaveName, "\n")
fit <- readRDS(here('Results', paste0(SaveName, '.RDS')))
list2env(fit, envir = .GlobalEnv)

nPop   <- data_list$nPop
nCore  <- data_list$nCore
nYears <- data_list$nYears
pop_core_vec <- data_list$pop_core + 1L
has_LKT_pop  <- data_list$has_LKT[pop_core_vec]
cfg          <- covar_config

cat("nPop:", nPop, "| nCore:", nCore, "| nYears:", nYears, "\n")
cat("Projecting", nYears_proj, "years beyond", DatYears[nYears], "\n\n")

# ================================================================
# EXTRACT MAP ESTIMATES AND STANDARD ERRORS FROM fit_rep
# ================================================================
# fit_rep$value / $sd  — ADREPORT quantities (logN, b0r_pop if reported)
# fit_rep$par.fixed    — outer fixed parameters
# summary(fit_rep, select="random") — random effect estimates + marginal SEs
#   (always available regardless of what was ADREPORT-ed in the C++ code)
# ================================================================
get_val <- function(nm) fit_rep$value[names(fit_rep$value) == nm]
get_sd  <- function(nm) fit_rep$sd[   names(fit_rep$value) == nm]

# Random effects: use summary() so we always get both estimate and SE
re_summ  <- summary(fit_rep, select = "random")   # matrix: Estimate | Std. Error
re_nm    <- rownames(re_summ)
get_re   <- function(nm) unname(re_summ[re_nm == nm, "Estimate"])
get_re_sd <- function(nm) unname(re_summ[re_nm == nm, "Std. Error"])

# logN and its SEs come from ADREPORT (used for historical bands)
logN_map <- matrix(get_val("logN"), nrow = nPop, ncol = nYears)
logN_sd  <- matrix(get_sd( "logN"), nrow = nPop, ncol = nYears)

# Random effects — values and marginal SEs
epsN_map   <- matrix(get_re("epsN"),          nrow = nPop)  # nPop x (nYears-1)
b0r_map    <- get_re("b0r_pop");      b0r_sd    <- get_re_sd("b0r_pop")
logphi_map <- get_re("log_b0phi_pop"); logphi_sd <- get_re_sd("log_b0phi_pop")
bFlow_map  <- get_re("bFlow_pop");    bFlow_sd  <- get_re_sd("bFlow_pop")
bTemp_map  <- get_re("bTemp_pop");    bTemp_sd  <- get_re_sd("bTemp_pop")
bINV_map   <- get_re("bINV_pop");     bINV_sd   <- get_re_sd("bINV_pop")

# Sanity checks — fail fast with a clear message
stopifnot("logN not in fit_rep$value — check ADREPORT in C++" =
            length(logN_map) == nPop * nYears)
stopifnot("b0r_pop not found in random effects" = length(b0r_map)    == nPop)
stopifnot("log_b0phi_pop not found"             = length(logphi_map) == nPop)

cat("Random effect SEs: b0r_pop mean SD =",   round(mean(b0r_sd),    3), "\n")
cat("                   log_b0phi_pop mean SD =", round(mean(logphi_sd), 3), "\n\n")

# bTemp2 is a fixed outer param (negative to enforce concave response)
bTemp2_map <- if (cfg$use_Temp && any(names(fit_rep$par.fixed) == "log_bTemp2"))
                -exp(fit_rep$par.fixed[["log_bTemp2"]]) else 0

# ================================================================
# FIXED OUTER PARAMETERS — small Cholesky of cov.fixed
# ================================================================
mu_fix  <- fit_rep$par.fixed
nm_fix  <- names(mu_fix)
cov_fix <- fit_rep$cov.fixed

# Tiny diagonal regularisation as safety (negligible effect)
eps_reg <- 1e-10 * mean(diag(cov_fix))
L_fix   <- chol(cov_fix + diag(eps_reg, nrow(cov_fix)))

cat("Fixed outer params:", length(mu_fix), "\n")
cat("  Names:", paste(nm_fix, collapse = ", "), "\n\n")

# ================================================================
# FUTURE COVARIATES
# ================================================================
if (future_cov_assumption == 'mean') {
  fut <- list(
    Flow = matrix(0, nPop,  nYears_proj),
    Temp = matrix(0, nPop,  nYears_proj),
    INV  = matrix(0, nPop,  nYears_proj),
    LKT  = matrix(0, nCore, nYears_proj)
  )
} else if (future_cov_assumption == 'last') {
  fut <- list(
    Flow = matrix(data_list$Flow[, nYears], nPop,  nYears_proj),
    Temp = matrix(data_list$Temp[, nYears], nPop,  nYears_proj),
    INV  = matrix(data_list$INV[,  nYears], nPop,  nYears_proj),
    LKT  = matrix(data_list$LKT[,  nYears], nCore, nYears_proj)
  )
} else {
  stop("future_cov_assumption must be 'mean' or 'last'")
}

# ================================================================
# FORWARD PROJECTION LOOP
# ================================================================
logN_start <- logN_map[, nYears]            # MAP N at last data year
eps_start  <- epsN_map[, ncol(epsN_map)]    # MAP last process error (AR(1) seed)

set.seed(seed)
proj_years <- DatYears[nYears] + seq_len(nYears_proj)
logN_proj  <- array(NA_real_, c(nsim, nPop, nYears_proj))

cat("Running", nsim, "projections...\n")
pb <- txtProgressBar(min = 0, max = nsim, style = 3)

for (i in seq_len(nsim)) {

  # --- Draw fixed outer params (correlated via cov.fixed Cholesky) ---
  th_fix   <- mu_fix + drop(t(L_fix) %*% rnorm(length(mu_fix)))
  sigmaN_i <- exp(th_fix[nm_fix == "log_sigmaN"][1])
  rho_i    <- tanh(th_fix[nm_fix == "logit_rho"][1])
  bLKT_i   <- if (cfg$use_LKT  && any(nm_fix == "bLKT_M"))
                th_fix[nm_fix == "bLKT_M"][1]   else 0
  bTemp2_i <- if (cfg$use_Temp && any(nm_fix == "log_bTemp2"))
                -exp(th_fix[nm_fix == "log_bTemp2"][1]) else bTemp2_map

  # --- Draw random effects from marginal normals ---
  b0r_i    <- rnorm(nPop, b0r_map,    b0r_sd)
  logphi_i <- rnorm(nPop, logphi_map, logphi_sd)
  bFlow_i  <- if (cfg$use_Flow && length(bFlow_map))
                rnorm(nPop, bFlow_map, bFlow_sd) else rep(0, nPop)
  bTemp_i  <- if (cfg$use_Temp && length(bTemp_map))
                rnorm(nPop, bTemp_map, bTemp_sd) else rep(0, nPop)
  bINV_i   <- if (cfg$use_INV  && length(bINV_map))
                rnorm(nPop, bINV_map,  bINV_sd)  else rep(0, nPop)

  phi_i        <- exp(logphi_i)
  sigma_innov  <- sigmaN_i * sqrt(max(1 - rho_i^2, 1e-8))

  # --- Project nYears_proj steps ---
  logN_t <- logN_start
  eps_t  <- eps_start

  for (t in seq_len(nYears_proj)) {
    R <- b0r_i
    if (cfg$use_Flow) R <- R + bFlow_i * fut$Flow[, t]
    if (cfg$use_Temp) R <- R + bTemp_i * fut$Temp[, t] + bTemp2_i * fut$Temp[, t]^2
    if (cfg$use_INV)  R <- R + bINV_i  * fut$INV[, t]
    if (cfg$use_LKT)  R <- R + bLKT_i  * has_LKT_pop * fut$LKT[pop_core_vec, t]

    logdens <- pmin(logN_t - log(data_list$extent), 10)
    R       <- R - phi_i * exp(logdens) / 1000

    eps_t  <- rho_i * eps_t + rnorm(nPop, 0, sigma_innov)
    logN_t <- pmax(pmin(logN_t + R + eps_t, 30), -30)
    logN_proj[i, , t] <- logN_t
  }
  setTxtProgressBar(pb, i)
}
close(pb)
cat("\nProjection complete.\n\n")

# ================================================================
# QUASI-EXTINCTION PROBABILITIES
# ================================================================
qext_results <- lapply(thresholds, function(thr) {
  N_arr  <- exp(logN_proj)
  p_yr   <- apply(N_arr < thr, c(2, 3), mean)
  p_ever <- apply(apply(N_arr < thr, c(1, 2), any), 2, mean)
  N_q    <- apply(N_arr, c(2, 3), quantile,
                  probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
  list(threshold = thr, p_yr = p_yr, p_ever = p_ever,
       N_p05 = N_q[1,,], N_p25 = N_q[2,,], N_p50 = N_q[3,,],
       N_p75 = N_q[4,,], N_p95 = N_q[5,,])
})
names(qext_results) <- paste0("thr", thresholds)

# ================================================================
# SUMMARY TABLE
# ================================================================
qext_table <- lapply(qext_results, function(r) {
  data.frame(
    pop_name    = U_LocalPops_In,
    core_name   = U_CoreAreas_In[pop_core_vec],
    threshold   = r$threshold,
    p_qext      = round(r$p_ever, 3),
    p_below_end = round(r$p_yr[, nYears_proj], 3),
    N_med_end   = round(r$N_p50[, nYears_proj], 1),
    N_lo_end    = round(r$N_p05[, nYears_proj], 1),
    N_hi_end    = round(r$N_p95[, nYears_proj], 1)
  )
}) |> bind_rows()

cat("--- Quasi-Extinction Summary (", nYears_proj, "-year horizon) ---\n", sep = "")
for (thr in thresholds) {
  df_t <- filter(qext_table, threshold == thr) |> arrange(desc(p_qext))
  cat(sprintf("\nThreshold N < %d:\n", thr))
  risk_pops <- filter(df_t, p_qext > 0.05)
  if (nrow(risk_pops)) {
    print(select(risk_pops, pop_name, core_name, p_qext, p_below_end, N_med_end),
          row.names = FALSE)
  } else {
    cat("  No populations with P(qext) > 5%\n")
  }
}

# ================================================================
# SAVE
# ================================================================
QextSave <- paste0(SaveName, '_Qext')
saveRDS(
  list(
    qext_results          = qext_results,
    qext_table            = qext_table,
    logN_proj             = logN_proj,
    thresholds            = thresholds,
    nYears_proj           = nYears_proj,
    proj_years            = proj_years,
    future_cov_assumption = future_cov_assumption,
    SaveName              = SaveName,
    pop_core_vec          = pop_core_vec,
    nsim                  = nsim
  ),
  file = here('Results', paste0(QextSave, '.RDS'))
)
cat("\nSaved to:", here('Results', paste0(QextSave, '.RDS')), "\n")

# ================================================================
# PLOTS
# ================================================================
pdf_path <- here('Results', paste0(QextSave, '_Plots.pdf'))
pdf(pdf_path, width = 11, height = 8.5)

# Historical CI from sdreport marginal SEs — lognormal approximation
# N ~ lognormal(logN_map, logN_sd^2), so CI = exp(logN_map ± z * logN_sd)
z90 <- qnorm(0.95)   # 1.645
z50 <- qnorm(0.75)   # 0.674

df_hist <- expand.grid(pop = seq_len(nPop), y_idx = seq_len(nYears)) |>
  mutate(
    pop_name  = U_LocalPops_In[pop],
    core_name = U_CoreAreas_In[pop_core_vec[pop]],
    year      = DatYears[y_idx],
    N_p50     = exp(logN_map[cbind(pop, y_idx)]),
    N_p05     = exp(logN_map[cbind(pop, y_idx)] - z90 * logN_sd[cbind(pop, y_idx)]),
    N_p25     = exp(logN_map[cbind(pop, y_idx)] - z50 * logN_sd[cbind(pop, y_idx)]),
    N_p75     = exp(logN_map[cbind(pop, y_idx)] + z50 * logN_sd[cbind(pop, y_idx)]),
    N_p95     = exp(logN_map[cbind(pop, y_idx)] + z90 * logN_sd[cbind(pop, y_idx)])
  )

r1 <- qext_results[[1]]
df_proj <- expand.grid(pop = seq_len(nPop), h = seq_len(nYears_proj)) |>
  mutate(
    pop_name  = U_LocalPops_In[pop],
    core_name = U_CoreAreas_In[pop_core_vec[pop]],
    year      = proj_years[h],
    N_p05     = r1$N_p05[cbind(pop, h)],
    N_p25     = r1$N_p25[cbind(pop, h)],
    N_p50     = r1$N_p50[cbind(pop, h)],
    N_p75     = r1$N_p75[cbind(pop, h)],
    N_p95     = r1$N_p95[cbind(pop, h)]
  )

# Plot 1: Trajectory ribbons by core area
# Page layout: compute facet grid dimensions and size the plot to fit naturally
# (avoids single-row panels stretching to fill the full 8.5-in page height)
traj_ncol <- 4   # max columns per page

for (ca in unique(df_hist$core_name)) {
  d_h <- filter(df_hist, core_name == ca)
  d_p <- filter(df_proj, core_name == ca)

  n_pops   <- n_distinct(d_h$pop_name)
  n_cols   <- min(traj_ncol, ceiling(sqrt(n_pops)))
  n_rows   <- ceiling(n_pops / n_cols)

  # ~2.2 in per panel row + 1.2 in for title/subtitle/axis labels
  plot_h_in  <- n_rows * 2.2 + 1.2
  height_frac <- min(1, plot_h_in / 8.5)

  p_traj <- ggplot() +
    geom_ribbon(data = d_h, aes(x = year, ymin = N_p05, ymax = N_p95),
                fill = "steelblue", alpha = 0.15) +
    geom_ribbon(data = d_h, aes(x = year, ymin = N_p25, ymax = N_p75),
                fill = "steelblue", alpha = 0.30) +
    geom_line(data = d_h, aes(x = year, y = N_p50),
              color = "steelblue", linewidth = 0.6) +
    geom_vline(xintercept = DatYears[nYears] + 0.5,
               linetype = 2, color = "grey50") +
    geom_ribbon(data = d_p, aes(x = year, ymin = N_p05, ymax = N_p95),
                fill = "tomato", alpha = 0.15) +
    geom_ribbon(data = d_p, aes(x = year, ymin = N_p25, ymax = N_p75),
                fill = "tomato", alpha = 0.30) +
    geom_line(data = d_p, aes(x = year, y = N_p50),
              color = "tomato", linewidth = 0.7) +
    geom_hline(yintercept = 20, linetype = 3,
               color = "firebrick", linewidth = 0.5) +
    facet_wrap(~pop_name, scales = "free_y", ncol = n_cols) +
    labs(x = "Year", y = "Adult Abundance (N)",
         title = paste0(ca, " - Historical (blue) + ",
                        nYears_proj, "-yr projection (red) | 50% and 90% CI"),
         subtitle = paste0("Future covariates: ", future_cov_assumption,
                           " | nsim = ", nsim,
                           " | dotted line: N = 20")) +
    theme_bw(base_size = 8)

  grid::grid.newpage()
  print(p_traj,
        vp = grid::viewport(width = 1, height = height_frac,
                            just = c("left", "top"), x = 0, y = 1))
}

# Plot 2: Quasi-extinction bar charts
for (thr in thresholds) {
  df_t <- filter(qext_table, threshold == thr) |>
    arrange(desc(p_qext)) |>
    mutate(pop_name = factor(pop_name, levels = pop_name))

  p_bar <- ggplot(df_t, aes(x = pop_name, y = p_qext)) +
    geom_col(aes(fill = p_qext > 0.05), alpha = 0.8, show.legend = FALSE) +
    scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "tomato")) +
    geom_hline(yintercept = 0.05, linetype = 2, color = "grey30") +
    annotate("text", x = nPop * 0.9, y = 0.07,
             label = "5% risk", size = 3, color = "grey30") +
    coord_flip() +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    labs(x = NULL, y = "Quasi-extinction probability",
         title = sprintf("P(N < %d at any point in %d-year horizon)",
                         thr, nYears_proj)) +
    theme_bw(base_size = 8)
  print(p_bar)
}

# Plot 3: P(below threshold) through time — per core area
df_time <- lapply(qext_results, function(r) {
  as.data.frame(r$p_yr) |>
    setNames(as.character(proj_years)) |>
    mutate(pop      = seq_len(nPop),
           pop_name  = U_LocalPops_In,
           core_name = U_CoreAreas_In[pop_core_vec],
           threshold = r$threshold) |>
    pivot_longer(as.character(proj_years), names_to = "year",
                 values_to = "p_below") |>
    mutate(year = as.integer(year))
}) |> bind_rows()

for (ca in unique(df_time$core_name)) {
  p_time <- df_time |>
    filter(core_name == ca) |>
    ggplot(aes(x = year, y = p_below,
               color = factor(threshold), group = factor(threshold))) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.5) +
    facet_wrap(~pop_name) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, NA)) +
    scale_color_brewer(palette = "OrRd", name = "Threshold\n(N adults)") +
    labs(x = "Year", y = "P(N < threshold)",
         title = paste0(ca, " - Extinction risk through time")) +
    theme_bw(base_size = 8)
  print(p_time)
}

dev.off()
cat("Plots saved to:", pdf_path, "\n")
