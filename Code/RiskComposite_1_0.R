# ================================================================
# Bull Trout MPVA — Composite Risk Plot
# RiskComposite_1_0.R
#
# Combines two orthogonal dimensions of population risk:
#   X-axis: b0r_pop — intrinsic log growth rate at low density
#           (from TMB model posterior; negative = declining, positive = growing)
#   Y-axis: P(N < threshold at any point in projection horizon)
#           (from Qext_MPVA1_2 simulation)
#   Point size: current abundance (MAP N at last data year)
#   Color: core area
#
# Quadrant interpretation:
#   Q1 (declining + high Qext): immediate crisis
#   Q2 (growing  + high Qext): small but recovering — stochastic risk
#   Q3 (declining + low Qext): large but declining — long-run threat
#   Q4 (growing  + low Qext): healthy
#
# Requires:
#   - Results/<SaveName>.RDS               (main model fit)
#   - Results/<SaveName>_Qext.RDS          (quasi-extinction output)
# ================================================================

rm(list = ls())

library(here)
library(tidyverse)
library(ggrepel)

# ================================================================
# CONFIGURATION
# ================================================================
DataVersion  <- '1_8'
ModelVersion <- '1_22'
RunVersion   <- '1_21'

SaveName <- paste0('MPVA_Fit_D', DataVersion, '_M', ModelVersion, '_R', RunVersion)
QextName <- paste0(SaveName, '_Qext')

# Which Qext threshold to use for the Y-axis
plot_threshold <- 20   # P(N < plot_threshold)

# Risk cutoffs for quadrant lines
lambda_cutoff <- 0      # b0r = 0 separates declining from growing
qext_cutoff   <- 0.10   # 10% quasi-extinction probability

# Label populations above this Qext probability OR with b0r below this
label_qext_above  <- 0.15
label_lambda_below <- -0.15

# ================================================================
# LOAD DATA
# ================================================================
cat("Loading fit:", SaveName, "\n")
fit <- readRDS(here('Results', paste0(SaveName, '.RDS')))

cat("Loading Qext:", QextName, "\n")
qext <- readRDS(here('Results', paste0(QextName, '.RDS')))

# Unpack fit
fit_rep      <- fit$fit_rep
data_list    <- fit$data_list
DatYears     <- fit$DatYears
U_LocalPops  <- fit$U_LocalPops_In
U_CoreAreas  <- fit$U_CoreAreas_In

nPop         <- data_list$nPop
nYears       <- data_list$nYears
pop_core_vec <- data_list$pop_core + 1L

# ================================================================
# EXTRACT LAMBDA (b0r_pop) FROM FIT
# ================================================================
b0r_vals <- fit_rep$value[names(fit_rep$value) == "b0r_pop"]
b0r_se   <- fit_rep$sd[names(fit_rep$value)    == "b0r_pop"]

stopifnot(length(b0r_vals) == nPop)

# ================================================================
# EXTRACT CURRENT ABUNDANCE (MAP N at final year)
# ================================================================
logN_map <- matrix(fit_rep$value[names(fit_rep$value) == "logN"],
                   nrow = nPop, ncol = nYears)
N_current <- exp(logN_map[, nYears])

# ================================================================
# EXTRACT QEXT PROBABILITIES
# ================================================================
thr_key <- paste0("thr", plot_threshold)
if (!thr_key %in% names(qext$qext_results)) {
  available <- paste(names(qext$qext_results), collapse = ", ")
  stop("Threshold ", plot_threshold, " not in Qext results. Available: ", available)
}
p_qext_vec <- qext$qext_results[[thr_key]]$p_ever   # length nPop

# ================================================================
# BUILD COMPOSITE DATA FRAME
# ================================================================
df <- data.frame(
  pop_name  = U_LocalPops,
  core_name = U_CoreAreas[pop_core_vec],
  b0r       = b0r_vals,
  b0r_se    = b0r_se,
  b0r_lo    = b0r_vals - 1.96 * b0r_se,
  b0r_hi    = b0r_vals + 1.96 * b0r_se,
  p_qext    = p_qext_vec,
  N_current = N_current
) |>
  mutate(
    # Risk quadrant label
    quadrant = case_when(
      b0r < lambda_cutoff & p_qext >= qext_cutoff ~ "Declining + High risk",
      b0r >= lambda_cutoff & p_qext >= qext_cutoff ~ "Growing + High risk",
      b0r < lambda_cutoff & p_qext < qext_cutoff  ~ "Declining + Low risk",
      TRUE                                         ~ "Stable/Growing + Low risk"
    ),
    # Which populations to label
    label_flag = p_qext > label_qext_above | b0r < label_lambda_below
  )

# ================================================================
# RISK SUMMARY
# ================================================================
cat("\n=== Composite Risk Summary ===\n")
cat(sprintf("Threshold: N < %d | Projection: %d years | nsim: %d\n",
            plot_threshold, qext$nYears_proj, qext$nsim))
cat(sprintf("Lambda cutoff: %g | Qext cutoff: %g\n\n",
            lambda_cutoff, qext_cutoff))

quadrant_summary <- df |>
  count(quadrant) |>
  mutate(pct = round(100 * n / nPop, 1)) |>
  arrange(desc(n))
print(quadrant_summary)

cat("\nHigh-risk populations (declining AND P(qext) >", qext_cutoff, "):\n")
df |>
  filter(b0r < lambda_cutoff, p_qext >= qext_cutoff) |>
  arrange(b0r) |>
  select(pop_name, core_name, b0r, b0r_se, p_qext, N_current) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  print(row.names = FALSE)

# ================================================================
# COMPOSITE SCATTER PLOT
# ================================================================
quadrant_colors <- c(
  "Declining + High risk"      = "#d73027",   # red
  "Growing + High risk"        = "#fc8d59",   # orange
  "Declining + Low risk"       = "#4575b4",   # blue
  "Stable/Growing + Low risk"  = "#74add1"    # light blue
)

# Clamp N_current for size scale (avoids huge outliers dominating)
size_cap   <- quantile(df$N_current, 0.95, na.rm = TRUE)
df$N_sized <- pmin(df$N_current, size_cap)

p_composite <- ggplot(df, aes(x = b0r, y = p_qext)) +

  # Quadrant shading
  annotate("rect", xmin = -Inf, xmax = lambda_cutoff,
           ymin = qext_cutoff, ymax = Inf,
           fill = "#d73027", alpha = 0.06) +
  annotate("rect", xmin = lambda_cutoff, xmax = Inf,
           ymin = qext_cutoff, ymax = Inf,
           fill = "#fc8d59", alpha = 0.06) +
  annotate("rect", xmin = -Inf, xmax = lambda_cutoff,
           ymin = -Inf, ymax = qext_cutoff,
           fill = "#4575b4", alpha = 0.06) +
  annotate("rect", xmin = lambda_cutoff, xmax = Inf,
           ymin = -Inf, ymax = qext_cutoff,
           fill = "#74add1", alpha = 0.06) +

  # Reference lines
  geom_vline(xintercept = lambda_cutoff, linetype = 2,
             color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = qext_cutoff, linetype = 2,
             color = "grey40", linewidth = 0.5) +

  # Error bars for b0r uncertainty
  geom_errorbarh(aes(xmin = b0r_lo, xmax = b0r_hi, color = quadrant),
                 height = 0, alpha = 0.4, linewidth = 0.4) +

  # Points
  geom_point(aes(size = N_sized, color = quadrant, fill = quadrant),
             alpha = 0.75, shape = 21, stroke = 0.4) +

  # Labels for high-risk or steeply declining populations
  geom_label_repel(
    data = filter(df, label_flag),
    aes(label = pop_name, color = quadrant),
    size       = 2.5,
    label.size = 0.2,
    max.overlaps = 20,
    show.legend  = FALSE
  ) +

  scale_color_manual(values = quadrant_colors, name = "Risk quadrant") +
  scale_fill_manual( values = quadrant_colors, name = "Risk quadrant") +
  scale_size_area(max_size = 10, name = paste0("Current N\n(MAP, ", DatYears[nYears], ")")) +
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(0, min(1, max(df$p_qext) * 1.15))) +

  # Quadrant annotation labels
  annotate("text", x = lambda_cutoff - 0.01, y = 1,
           label = "Declining", hjust = 1, vjust = 1,
           size = 3, color = "grey40", fontface = "italic") +
  annotate("text", x = lambda_cutoff + 0.01, y = 1,
           label = "Growing", hjust = 0, vjust = 1,
           size = 3, color = "grey40", fontface = "italic") +
  annotate("text", x = min(df$b0r) * 0.95, y = qext_cutoff + 0.01,
           label = paste0("High risk (>", 100 * qext_cutoff, "%)"),
           hjust = 0, vjust = 0, size = 2.8, color = "grey40", fontface = "italic") +

  labs(
    x        = expression(b[0*r]*" - Intrinsic log growth rate (per year)"),
    y        = paste0("P(N < ", plot_threshold, " at any point in ",
                      qext$nYears_proj, " years)"),
    title    = "Bull Trout Local Population Composite Risk",
    subtitle = paste0("Model v", ModelVersion, " / Run v", RunVersion,
                      " | Data v", DataVersion,
                      " | Future covariates: ", qext$future_cov_assumption,
                      " | nsim = ", qext$nsim)
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

print(p_composite)

# ================================================================
# PER-CORE VERSION (faceted)
# ================================================================
p_core <- p_composite +
  facet_wrap(~core_name, scales = "free") +
  theme(strip.text = element_text(size = 7)) +
  labs(title = "Bull Trout Composite Risk - by Core Area")

print(p_core)

# ================================================================
# RISK TABLE: all populations sorted by composite threat
# ================================================================
# Composite score: simple sum of standardised b0r (inverted) and p_qext
df_table <- df |>
  mutate(
    b0r_rank   = rank(-b0r),        # high rank = more declining
    qext_rank  = rank(-p_qext),     # high rank = higher Qext risk
    N_rank     = rank(N_current),   # high rank = smaller current N
    # Equal-weight composite of three axes
    risk_score = (b0r_rank + qext_rank + N_rank) / 3
  ) |>
  arrange(desc(risk_score)) |>
  select(pop_name, core_name, b0r, b0r_se, p_qext, N_current, quadrant, risk_score) |>
  mutate(across(c(b0r, b0r_se, p_qext, risk_score), \(x) round(x, 3)),
         N_current = round(N_current, 1))

cat("\n=== Full Risk Table (sorted by composite score) ===\n")
print(df_table, row.names = FALSE)

# ================================================================
# SAVE OUTPUTS
# ================================================================
pdf_path   <- here('Results', paste0(SaveName, '_RiskComposite.pdf'))
csv_path   <- here('Results', paste0(SaveName, '_RiskTable.csv'))

pdf(pdf_path, width = 11, height = 8.5)
print(p_composite)
print(p_core)
dev.off()

write.csv(df_table, csv_path, row.names = FALSE)

cat("\nComposite plot saved to:", pdf_path, "\n")
cat("Risk table saved to:",     csv_path, "\n")
