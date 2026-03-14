# ================================================================
# MPVA TMB — Helper Functions v1_1
# Source this at the top of Run and Diagnostics scripts
# Changes from v1_0: hierarchical covariate slope parameters
# ================================================================

# Required packages
required_pkgs <- c("tidyverse", "TMB")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

# ================================================================
# 1. INITIALIZATION
# ================================================================

#' Back-calculate initial log-abundances from redd counts
init_pars_from_redds <- function(data_list,
                                 N_floor  = 1,
                                 Ns_cap   = 1e7,
                                 use_year = c("earliest", "maxdata", "medianyear")) {
  use_year   <- match.arg(use_year)
  ReddCounts <- data_list$ReddCounts
  pS         <- data_list$pS
  site_pop   <- data_list$site_pop
  extent     <- data_list$extent
  nPop       <- data_list$nPop
  nYears     <- data_list$nYears
  nSites     <- data_list$nSites
  nCore      <- data_list$nCore

  stopifnot(nrow(ReddCounts) == nSites, ncol(ReddCounts) == nYears)
  stopifnot(length(pS) == nSites, length(site_pop) == nSites, length(extent) == nPop)

  Nhat_py <- matrix(NA_real_, nPop, nYears)
  for (p in 0:(nPop - 1)) {
    sites_p <- which(site_pop == p)
    if (!length(sites_p)) next
    for (y in seq_len(nYears)) {
      rc  <- ReddCounts[sites_p, y]
      ok  <- is.finite(rc) & !is.na(rc) & is.finite(pS[sites_p]) &
        !is.na(pS[sites_p]) & (pS[sites_p] > 0)
      if (!any(ok)) next
      Ns  <- 2 * rc[ok] / pS[sites_p][ok]
      Ns  <- Ns[is.finite(Ns) & Ns > 0]
      if (!length(Ns)) next
      Nhat_py[p + 1, y] <- median(pmin(Ns, Ns_cap))
    }
  }

  pick_year <- function(v) {
    idx <- which(is.finite(v))
    if (!length(idx)) return(NA_integer_)
    switch(use_year,
           earliest   = idx[1],
           maxdata    = idx[which.max(v[idx])],
           medianyear = idx[ceiling(length(idx) / 2)]
    )
  }

  logN0 <- rep(log(10), nPop)
  for (p in seq_len(nPop)) {
    y0 <- pick_year(Nhat_py[p, ])
    logN0[p] <- if (!is.na(y0)) {
      log(max(Nhat_py[p, y0], N_floor))
    } else {
      log(max(N_floor, 0.01 * extent[p]))
    }
  }

  list(
    log_sigmaN           = log(0.2),
    log_sigma_b0r_pop    = log(0.3),
    log_sigma_b0r_core   = log(0.3),
    log_sigma_b0phi_pop  = log(0.3),
    log_sigma_b0phi_core = log(1e-4),
    b0r_M                = 0,
    log_b0phi_M          = log(2),
    b0r_pop              = rep(0, nPop),
    b0r_core             = rep(0, nCore),
    log_b0phi_pop        = rep(log(2), nPop),
    log_b0phi_core       = rep(log(2), nCore),
    logN0                = logN0,
    epsN                 = matrix(0, nPop, nYears - 1),
    bAdf       = 0,
    bVB        = 0,
    bFlood     = 0,
    bRD        = 0,
    # Flow: hierarchical
    bFlow_M        = 0,
    bFlow_pop      = rep(0, nPop),
    log_sigma_bFlow = log(0.3),
    # Temp linear: hierarchical; quadratic global
    bTemp_M        = 0,
    bTemp_pop      = rep(0, nPop),
    log_sigma_bTemp = log(0.3),
    log_bTemp2     = log(0.05),
    # Invasives: hierarchical
    bINV_M         = 0,
    bINV_pop       = rep(0, nPop),
    log_sigma_bINV  = log(0.3),
    # Lake trout: core-level hierarchical
    bLKT_M         = 0,
    bLKT_core      = rep(0, nCore),
    log_sigma_bLKT  = log(0.3)
  )
}

# ================================================================
# 2. OPTIMIZATION SETUP
# ================================================================

#' Build named parameter bounds for nlminb
make_bounds <- function(obj) {
  pn    <- names(obj$par)
  lower <- rep(-Inf, length(pn)); names(lower) <- pn
  upper <- rep( Inf, length(pn)); names(upper) <- pn

  lower["log_sigmaN"] <- log(1e-4); upper["log_sigmaN"] <- log(5)

  if ("b0r_M" %in% pn) { lower["b0r_M"] <- -2; upper["b0r_M"] <- 2 }
  idx <- grep("^b0r_pop", pn)
  if (length(idx)) { lower[idx] <- -2; upper[idx] <- 2 }

  lower["log_b0phi_M"] <- log(1e-4); upper["log_b0phi_M"] <- log(10)
  idx <- grep("^log_b0phi_pop", pn)
  if (length(idx)) { lower[idx] <- log(1e-4); upper[idx] <- log(10) }

  idx <- grep("^logN0", pn)
  if (length(idx)) { lower[idx] <- log(1e-3); upper[idx] <- log(1e10) }

  if ("log_bTemp2" %in% pn) {
    lower["log_bTemp2"] <- log(1e-4)
    upper["log_bTemp2"] <- log(10)
  }

  # Covariate slope SD bounds
  for (nm in c("log_sigma_bFlow", "log_sigma_bTemp",
               "log_sigma_bINV",  "log_sigma_bLKT")) {
    if (nm %in% pn) { lower[nm] <- log(1e-4); upper[nm] <- log(5) }
  }

  list(lower = lower, upper = upper)
}

#' Build map list to fix unused covariate parameters
#' LKT uses nCore-length factor; all others use nPop-length factor
make_covar_map <- function(config, nPop, nCore) {
  covar_list <- list(
    list(flag = "use_Adfluvial", scalars = "bAdf",
         pop_res = NULL, core_res = NULL),
    list(flag = "use_VB",        scalars = "bVB",
         pop_res = NULL, core_res = NULL),
    list(flag = "use_Flood",     scalars = "bFlood",
         pop_res = NULL, core_res = NULL),
    list(flag = "use_RD",        scalars = "bRD",
         pop_res = NULL, core_res = NULL),
    list(flag = "use_Flow",      scalars = c("bFlow_M", "log_sigma_bFlow"),
         pop_res = "bFlow_pop", core_res = NULL),
    list(flag = "use_Temp",      scalars = c("bTemp_M", "log_bTemp2", "log_sigma_bTemp"),
         pop_res = "bTemp_pop", core_res = NULL),
    list(flag = "use_INV",       scalars = c("bINV_M", "log_sigma_bINV"),
         pop_res = "bINV_pop",  core_res = NULL),
    list(flag = "use_LKT",       scalars = c("bLKT_M", "log_sigma_bLKT"),
         pop_res = NULL, core_res = "bLKT_core")
  )

  map <- list()
  for (cov in covar_list) {
    if (config[[cov$flag]] == 0) {
      for (b in cov$scalars) map[[b]] <- factor(NA)
      if (!is.null(cov$pop_res))
        map[[cov$pop_res]]  <- factor(rep(NA, nPop))
      if (!is.null(cov$core_res))
        map[[cov$core_res]] <- factor(rep(NA, nCore))
    }
  }
  map
}

# ================================================================
# 3. FIT DIAGNOSTICS
# ================================================================

#' Convergence diagnostics — requires live TMB object
check_fit <- function(obj, opt) {
  cat("\n--- Optimizer ---\n")
  cat("Convergence:", opt$convergence, "\n")
  if (!is.null(opt$message)) cat("Message:", opt$message, "\n")
  cat("Objective:", opt$objective, "\n")

  g <- obj$gr(opt$par)
  cat("\n--- Gradient ---\n")
  cat("max|grad|:", max(abs(g)), "\n")
  cat("median|grad|:", median(abs(g)), "\n")

  grad_df <- data.frame(par = names(opt$par), grad = g) %>%
    arrange(desc(abs(grad))) %>% head(10)
  cat("\nLargest gradients:\n"); print(grad_df)

  cat("\n--- sdreport ---\n")
  rep <- try(TMB::sdreport(obj), silent = TRUE)
  if (inherits(rep, "try-error")) {
    cat("sdreport FAILED\n")
    return(invisible(list(grad = g, sdreport = NULL)))
  }
  s <- summary(rep)
  cat("Parameters reported:", nrow(s), "\n")
  cat("Any non-finite SEs:", any(!is.finite(s[, "Std. Error"])), "\n")

  invisible(list(grad = g, sdreport = rep))
}

#' Full fit summary — works from saved RDS without TMB object
summarize_fit <- function(fit_rep, optB, data_list, covar_config,
                          U_LocalPops_In, U_CoreAreas_In, pop_core_vec) {

  cat("================================================================\n")
  cat("Model Fit Summary\n")
  cat("================================================================\n")

  cat("\n--- Optimizer ---\n")
  cat("Convergence:", optB$convergence, "\n")
  if (!is.null(optB$message)) cat("Message:", optB$message, "\n")
  cat("Objective:", optB$objective, "\n")

  cat("\n--- Gradient ---\n")
  if (!is.null(fit_rep$gradient.fixed)) {
    cat("max|grad|:", max(abs(fit_rep$gradient.fixed)), "\n")
    cat("median|grad|:", median(abs(fit_rep$gradient.fixed)), "\n")
    grad_df <- data.frame(
      par  = names(fit_rep$gradient.fixed),
      grad = fit_rep$gradient.fixed
    ) %>% arrange(desc(abs(grad))) %>% head(10)
    cat("\nLargest gradients:\n"); print(grad_df)
  } else {
    cat("Gradient not stored in sdreport.\n")
  }

  cat("\n--- Standard Errors ---\n")
  s <- summary(fit_rep)
  cat("Parameters reported:", nrow(s), "\n")
  cat("Any non-finite SEs:", any(!is.finite(s[, "Std. Error"])), "\n")

  cat("\n--- Variance Components ---\n")
  sigma_names <- c("sigmaN", "sigma_b0r_pop", "sigma_b0r_core",
                   "sigma_b0phi_pop", "sigma_b0phi_core",
                   "sigma_bFlow", "sigma_bTemp", "sigma_bINV", "sigma_bLKT")
  for (nm in sigma_names) {
    val <- fit_rep$value[names(fit_rep$value) == nm]
    se  <- fit_rep$sd[names(fit_rep$value) == nm]
    if (length(val)) cat(sprintf("  %-25s %.3f (SE: %.3f)\n", nm, val[1], se[1]))
  }

  cat("\n--- Global Growth Rate ---\n")
  b0r_M    <- fit_rep$value[names(fit_rep$value) == "b0r_M"]
  b0r_M_se <- fit_rep$sd[names(fit_rep$value) == "b0r_M"]
  cat(sprintf("  b0r_M: %.3f (SE: %.3f)\n", b0r_M, b0r_M_se))

  cat("\n--- Core Area Growth Rates ---\n")
  b0r_core    <- fit_rep$value[names(fit_rep$value) == "b0r_core"]
  b0r_core_se <- fit_rep$sd[names(fit_rep$value) == "b0r_core"]
  data.frame(
    core_name = U_CoreAreas_In,
    b0r_core  = round(b0r_core, 3),
    se        = round(b0r_core_se, 3),
    lower     = round(b0r_core - 1.96 * b0r_core_se, 3),
    upper     = round(b0r_core + 1.96 * b0r_core_se, 3)
  ) %>% print()

  cat("\n--- Density Dependence (fish/km) ---\n")
  cat("b0phi_M:   ", exp(fit_rep$value[names(fit_rep$value) == "log_b0phi_M"])   / 1000, "\n")
  cat("b0phi_core:", exp(fit_rep$value[names(fit_rep$value) == "log_b0phi_core"]) / 1000, "\n")
  cat("b0phi_pop: ", exp(fit_rep$value[names(fit_rep$value) == "log_b0phi_pop"])  / 1000, "\n")

  covar_summary(fit_rep, covar_config)

  b0r_pop_vals <- fit_rep$value[names(fit_rep$value) == "b0r_pop"]
  b0r_pop_se   <- fit_rep$sd[names(fit_rep$value) == "b0r_pop"]

  b0r_df <- data.frame(
    pop       = seq_len(data_list$nPop),
    pop_name  = U_LocalPops_In,
    core_name = U_CoreAreas_In[pop_core_vec],
    b0r       = b0r_pop_vals,
    b0r_se    = b0r_pop_se
  )

  print(
    ggplot(b0r_df, aes(x = core_name, y = b0r)) +
      geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
      geom_boxplot(fill = "steelblue", alpha = 0.5, outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
      geom_errorbar(aes(ymin = b0r - 1.96 * b0r_se,
                        ymax = b0r + 1.96 * b0r_se),
                    width = 0.1, alpha = 0.4) +
      labs(x = "Core Area", y = "b0r (pop level)",
           title = "Population Growth Rates by Core Area") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )

  invisible(b0r_df)
}

# ================================================================
# 4. OBSERVATION MODEL HELPERS
# ================================================================

#' Extract latent abundance trajectories into tidy dataframe
latent_df <- function(data_list, rep) {
  logN_vals <- rep$value[names(rep$value) == "logN"]
  logN_mat  <- matrix(logN_vals, nrow = data_list$nPop, ncol = data_list$nYears)
  tibble(
    pop  = rep(seq_len(data_list$nPop),  each  = data_list$nYears),
    year = rep(seq_len(data_list$nYears), times = data_list$nPop),
    logN = as.vector(t(logN_mat)),
    N    = exp(logN)
  )
}

#' Build observed vs fitted dataframe at site level
obs_fit_df <- function(data_list, rep) {
  logN_vals <- rep$value[names(rep$value) == "logN"]
  logN_mat  <- matrix(logN_vals, nrow = data_list$nPop, ncol = data_list$nYears)
  expand.grid(site = seq_len(data_list$nSites),
              year = seq_len(data_list$nYears)) %>%
    as_tibble() %>%
    mutate(
      pop    = data_list$site_pop[site] + 1L,
      y_obs  = as.vector(data_list$ReddCounts),
      N_pred = exp(logN_mat[cbind(pop, year)]),
      mu     = 0.5 * N_pred * data_list$pS[site]
    ) %>%
    filter(is.finite(y_obs), !is.na(y_obs))
}

# ================================================================
# 5. PLOTTING
# ================================================================

#' Standard diagnostic plots
plot_fit <- function(data_list, rep) {
  dfN      <- latent_df(data_list, rep)
  df_fit   <- obs_fit_df(data_list, rep)
  df_resid <- df_fit %>%
    mutate(
      pearson = (y_obs - mu) / sqrt(pmax(mu, 1e-12)),
      pop     = factor(pop),
      site    = factor(site)
    )

  print(ggplot(dfN, aes(x = year, y = N, group = pop)) +
          geom_line() + facet_wrap(~pop, scales = "free_y") +
          labs(x = "Year", y = "Latent N",
               title = "Latent abundance by population") +
          theme_bw())

  print(ggplot(df_fit, aes(x = year)) +
          geom_point(aes(y = y_obs), alpha = 0.7) +
          geom_line(aes(y = mu), linewidth = 0.7) +
          facet_wrap(~pop, scales = "free_y") +
          labs(x = "Year", y = "Redd counts",
               title = "Observed vs fitted by population") +
          theme_bw())

  print(ggplot(df_fit, aes(x = year)) +
          geom_point(aes(y = y_obs), alpha = 0.7) +
          geom_line(aes(y = mu), linewidth = 0.7) +
          facet_wrap(~site, scales = "free_y") +
          labs(x = "Year", y = "Redd counts",
               title = "Observed vs fitted by site") +
          theme_bw())

  print(ggplot(df_resid, aes(x = year, y = pearson)) +
          geom_hline(yintercept = 0, linetype = 2) +
          geom_point(alpha = 0.6) + facet_wrap(~pop) +
          labs(x = "Year", y = "Pearson residual",
               title = "Residuals by population") +
          theme_bw())

  print(ggplot(df_resid, aes(x = pearson)) +
          geom_histogram(bins = 40) +
          labs(x = "Pearson residual",
               title = "Residual distribution") +
          theme_bw())

  bad <- df_resid %>%
    group_by(pop, site) %>%
    summarize(n = n(), rmse = sqrt(mean((y_obs - mu)^2)),
              mean_pearson = mean(pearson), sd_pearson = sd(pearson),
              .groups = "drop") %>%
    arrange(desc(rmse))
  cat("\nWorst-fitting sites:\n"); print(head(bad, 15))

  invisible(list(latent = dfN, fit = df_fit, resid = df_resid, bad_sites = bad))
}

#' Covariate parameter estimates — hyper-means and SDs
covar_summary <- function(rep, config) {
  covar_list <- list(
    list(flag = "use_Adfluvial", betas = "bAdf",   labels = "Adfluvial",
         sd_par = NULL, res_par = NULL),
    list(flag = "use_VB",        betas = "bVB",    labels = "Valley bottom",
         sd_par = NULL, res_par = NULL),
    list(flag = "use_Flood",
         betas = "bFlood", labels = "W95 flood frequency",
         sd_par = NULL, res_par = NULL),
    list(flag = "use_RD",
         betas = "bRD",    labels = "Road density",
         sd_par = NULL, res_par = NULL),
    list(flag = "use_Flow",
         betas = "bFlow_M",      labels = "Flow (hyper-mean)",
         sd_par = "sigma_bFlow", res_par = "bFlow_pop"),
    list(flag = "use_Temp",
         betas = c("bTemp_M", "bTemp2"), labels = c("Temp (hyper-mean)", "Temp (quadratic)"),
         sd_par = "sigma_bTemp", res_par = "bTemp_pop"),
    list(flag = "use_INV",
         betas = "bINV_M",       labels = "INV (hyper-mean)",
         sd_par = "sigma_bINV",  res_par = "bINV_pop"),
    list(flag = "use_LKT",
         betas = "bLKT_M",       labels = "LKT (hyper-mean)",
         sd_par = "sigma_bLKT",  res_par = "bLKT_core")
  )

  rows <- list()
  for (cov in covar_list) {
    if (config[[cov$flag]] == 1) {
      for (i in seq_along(cov$betas)) {
        b   <- cov$betas[i]
        idx <- which(names(rep$value) == b)
        if (!length(idx)) next
        rows[[length(rows) + 1]] <- data.frame(
          covariate = cov$labels[i],
          parameter = b,
          estimate  = rep$value[idx[1]],
          se        = rep$sd[idx[1]]
        )
      }
    }
  }

  if (!length(rows)) { cat("No covariates active.\n"); return(invisible(NULL)) }

  df <- bind_rows(rows) %>%
    mutate(
      z     = estimate / se,
      lower = estimate - 1.96 * se,
      upper = estimate + 1.96 * se
    )

  cat("\n--- Covariate Estimates ---\n")
  print(df)
  invisible(df)
}

#' Plot implied thermal performance curve using hyper-mean temperature parameters
plot_temp_curve <- function(rep, temp_mean = 0, temp_sd = 1) {
  b1 <- rep$value[names(rep$value) == "bTemp_M"]
  b2 <- rep$value[names(rep$value) == "bTemp2"]

  if (!length(b1) || !length(b2)) {
    cat("Temperature parameters not found in rep.\n")
    return(invisible(NULL))
  }

  temp_z   <- seq(-3, 3, length.out = 200)
  temp_raw <- temp_z * temp_sd + temp_mean
  effect   <- b1 * temp_z + b2 * temp_z^2
  opt_z    <- -b1 / (2 * b2)
  opt_raw  <- opt_z * temp_sd + temp_mean

  cat("Thermal optimum (z-scored):", round(opt_z, 3), "\n")
  cat("Thermal optimum (raw scale):", round(opt_raw, 2), "\n")

  print(
    data.frame(temp_raw = temp_raw, effect = effect) %>%
      ggplot(aes(x = temp_raw, y = effect)) +
      geom_line(linewidth = 1, color = "steelblue") +
      geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
      geom_vline(xintercept = opt_raw, linetype = 2, color = "tomato") +
      annotate("text", x = opt_raw, y = max(effect),
               label = paste0("optimum = ", round(opt_raw, 1)),
               hjust = -0.1, color = "tomato", size = 3.5) +
      labs(x = "Temperature", y = "Effect on log growth rate",
           title = "Thermal performance curve (global hyper-mean)") +
      theme_bw()
  )

  invisible(list(b1 = b1, b2 = b2, opt_raw = opt_raw))
}

# ================================================================
# 6. PROJECTION
# ================================================================

#' Forward project population trajectories from final fitted state
project_populations <- function(fit_rep, data_list,
                                nProj     = 5,
                                nSim      = 10000,
                                threshold = 50,
                                seed      = 42) {
  set.seed(seed)

  logN_vals  <- fit_rep$value[names(fit_rep$value) == "logN"]
  logN_mat   <- matrix(logN_vals, nrow = data_list$nPop, ncol = data_list$nYears)
  logN_final <- logN_mat[, data_list$nYears]

  r_pop     <- fit_rep$value[names(fit_rep$value) == "b0r_pop"]
  b0phi_pop <- -exp(fit_rep$value[names(fit_rep$value) == "log_b0phi_pop"])
  sigmaN    <- fit_rep$value[names(fit_rep$value) == "sigmaN"]
  extent_km <- data_list$extent / 1000

  proj_array <- array(NA_real_, dim = c(data_list$nPop, nProj, nSim))

  for (sim in seq_len(nSim)) {
    logN_curr <- logN_final
    for (t in seq_len(nProj)) {
      dens      <- exp(logN_curr) / extent_km
      R         <- r_pop + b0phi_pop * dens
      eps       <- rnorm(data_list$nPop, 0, sigmaN)
      logN_curr <- logN_curr + R + eps
      proj_array[, t, sim] <- exp(logN_curr)
    }
  }

  qe_prob <- apply(proj_array, 1, function(pop_sims) {
    mean(apply(pop_sims, 2, function(traj) any(traj < threshold)))
  })

  proj_summary <- apply(proj_array, c(1, 2), function(x) {
    c(median = median(x),
      lo90   = quantile(x, 0.05),
      hi90   = quantile(x, 0.95))
  })

  cat("\n--- Quasi-Extinction Probabilities ---\n")
  cat("Threshold:", threshold, "adults |",
      "Projection:", nProj, "years |",
      "Simulations:", nSim, "\n\n")

  qe_df <- data.frame(
    pop     = seq_len(data_list$nPop),
    qe_prob = round(qe_prob, 3)
  ) %>% arrange(desc(qe_prob))
  print(qe_df)

  invisible(list(
    proj_array   = proj_array,
    proj_summary = proj_summary,
    qe_prob      = qe_prob,
    qe_df        = qe_df
  ))
}
