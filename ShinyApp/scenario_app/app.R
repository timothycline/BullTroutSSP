# ================================================================
# Bull Trout Conservation Prioritization Tool
# Shiny Scenario Explorer
#
# Requires prioritization_data.RDS in the same directory.
# Generate it by running: Code/PrioritizationTable_2_0.R
#
# Deploy to shinyapps.io:
#   rsconnect::deployApp("ShinyApp/scenario_app")
# ================================================================

library(shiny)
library(reactable)
library(tidyverse)

# ── Load pre-computed data ────────────────────────────────────────
pdat <- readRDS("prioritization_data.RDS")
st   <- pdat$status_table
samp <- pdat$samples
cov  <- pdat$cov_info
meta <- pdat$meta

nPop        <- meta$nPop
pop_display <- meta$pop_display
core_names  <- meta$core_names[meta$pop_core_vec]   # per-pop core label
has_lkt     <- cov$has_LKT_pop
cfg         <- meta$covar_config
nsamp       <- meta$nsamp

# ── Scenario computation ─────────────────────────────────────────
# Pure R matrix ops on pre-stored posterior samples — no TMB needed.
# inv_z, temp_z, flow_z, lkt_z: per-pop z-scored covariate vectors.
compute_lambda <- function(inv_z, temp_z, flow_z, lkt_z) {
  nP    <- nPop
  nsmp  <- nsamp

  log_lam <- samp$b0r_pop   # nsamp x nPop, start with intrinsic rate

  if (cfg$use_Flow)
    log_lam <- log_lam + samp$bFlow_pop * matrix(flow_z, nsmp, nP, byrow = TRUE)

  if (cfg$use_Temp)
    log_lam <- log_lam +
      samp$bTemp_pop * matrix(temp_z, nsmp, nP, byrow = TRUE) +
      matrix(samp$bTemp2, nsmp, nP) * matrix(temp_z^2, nsmp, nP, byrow = TRUE)

  if (cfg$use_INV)
    log_lam <- log_lam + samp$bINV_pop * matrix(inv_z, nsmp, nP, byrow = TRUE)

  if (cfg$use_LKT)
    log_lam <- log_lam +
      matrix(samp$bLKT_M, nsmp, nP) *
      matrix(has_lkt * lkt_z, nsmp, nP, byrow = TRUE)

  # Return nPop x 3 matrix (10th, 50th, 90th percentile)
  t(apply(exp(log_lam), 2, quantile, probs = c(0.10, 0.50, 0.90), na.rm = TRUE))
}

# Pre-compute current-conditions λ for the overlay
current_lam <- compute_lambda(cov$INV_last, cov$Temp_last, cov$Flow_last, cov$LKT_last)

core_choices <- c("All core areas", sort(unique(core_names)))

# ── UI ───────────────────────────────────────────────────────────
ui <- navbarPage(
  title = "Bull Trout Prioritization Tool",
  id    = "nav",

  # ── Tab 1: Population Status ──────────────────────────────────
  tabPanel("Population Status",
    fluidRow(
      column(12,
        tags$h4("Current Population Status — All Local Populations"),
        tags$p(
          "Model-estimated intrinsic growth rate (λ) and quasi-extinction risk",
          "under current conditions (last observed year of each covariate).",
          tags$br(),
          tags$b("λ > 1:"), " population growing under current conditions.",
          tags$b(" λ < 1:"), " declining.",
          tags$br(),
          tags$i("Sensitivity columns show λ if one threat lever is improved while others remain at current levels.")
        ),
        tags$hr(),
        reactableOutput("status_tbl", height = "680px")
      )
    )
  ),

  # ── Tab 2: Scenario Explorer ──────────────────────────────────
  tabPanel("Scenario Explorer",
    sidebarLayout(
      sidebarPanel(width = 3,
        tags$h5("Set Scenario Conditions"),
        tags$p(
          "Sliders set covariate values in", tags$b("standard deviation units"),
          "relative to each population's historical mean.",
          "0 = historical average; −1 = one SD below mean (less pressure).",
          style = "font-size: 12px; color: #555;"
        ),
        tags$hr(),

        sliderInput("inv_z",
                    label     = "Invasive species pressure (SD)",
                    min       = max(-3, floor(cov$INV_range[1])),
                    max       = min( 3, ceiling(cov$INV_range[2])),
                    value     = 0,
                    step      = 0.25),

        sliderInput("temp_z",
                    label     = "Summer temperature (SD)",
                    min       = max(-3, floor(cov$Temp_range[1])),
                    max       = min( 5, ceiling(cov$Temp_range[2])),
                    value     = 0,
                    step      = 0.25),

        sliderInput("flow_z",
                    label     = "Summer streamflow (SD)",
                    min       = max(-3, floor(cov$Flow_range[1])),
                    max       = min( 3, ceiling(cov$Flow_range[2])),
                    value     = 0,
                    step      = 0.25),

        conditionalPanel(
          condition = "true",   # always show; LKT slider only affects LKT pops
          sliderInput("lkt_z",
                      label   = "Lake trout pressure (SD; affects lake-connected pops only)",
                      min     = max(-3, floor(cov$LKT_range[1])),
                      max     = min( 3, ceiling(cov$LKT_range[2])),
                      value   = 0,
                      step    = 0.25)
        ),

        tags$hr(),

        checkboxInput("show_current",
                      "Show current conditions (open circles)",
                      value = TRUE),

        selectInput("core_filter",
                    "Focus on core area:",
                    choices  = core_choices,
                    selected = "All core areas"),

        tags$hr(),
        actionButton("reset_btn", "Reset to historical mean", class = "btn-sm btn-outline-secondary"),

        tags$hr(),
        tags$p(
          tags$b("Filled dots:"), " scenario λ with 80% CI.", tags$br(),
          tags$b("Open circles:"), " current conditions.", tags$br(),
          tags$span("Blue = λ ≥ 1 (growing)  |  Red = λ < 1 (declining)",
                    style = "font-size: 11px;")
        )
      ),

      mainPanel(width = 9,
        tags$h5(textOutput("scenario_subtitle")),
        plotOutput("scenario_plot", height = "680px")
      )
    )
  ),

  # ── Tab 3: About ──────────────────────────────────────────────
  tabPanel("About",
    fluidRow(
      column(8, offset = 1,
        tags$h4("About This Tool"),
        tags$p(
          "This tool supports prioritization of bull trout conservation actions",
          "in the Upper Columbia River Basin (Montana), as part of the Bull Trout",
          "Interagency Recovery Team (BIRT) conservation framework."
        ),
        tags$p(
          "Population-level estimates are derived from a hierarchical state-space",
          "Ricker model (TMB) fitted to annual redd count time series across",
          meta$nPop, "local populations nested within",
          length(meta$core_names), "core areas.",
          tags$br(),
          sprintf("Model v%s | Data v%s | Run date: %s",
                  meta$ModelVersion, meta$DataVersion, meta$run_date)
        ),
        tags$h5("Interpreting the Scenario Explorer"),
        tags$p(
          "The Scenario Explorer shows", tags$b("intrinsic growth rate (λ)"),
          "— the rate of population change at each population's own density-dependence",
          "equilibrium, given the specified covariate conditions.",
          "Covariate sliders are calibrated in standard deviation (SD) units",
          "relative to each population's historical mean, so 0 = average historical",
          "pressure and −1 = one SD below average (more favorable).",
          "Uncertainty bars reflect the marginal posterior distribution of",
          "population-level covariate effects."
        ),
        tags$h5("Sensitivity Columns (Population Status tab)"),
        tags$ul(
          tags$li(tags$b("λ (INV avg):"), " INV set to historical mean (z = 0); other covariates at current."),
          tags$li(tags$b("λ (INV best):"), " INV set to historical minimum; other covariates at current."),
          tags$li(tags$b("λ (Temp −1SD):"), " Temperature 1 SD below current."),
          tags$li(tags$b("λ (Flow +1SD):"), " Flow 1 SD above current."),
          tags$li(tags$b("λ (LKT best):"), " Lake trout at historical minimum (lake-connected pops only).")
        ),
        tags$hr(),
        tags$p(
          tags$b("Contacts:"),
          "Timothy Cline (Montana State University)",
          "· Clint Muhlfeld / Robert Al-Chokhachy (U.S. Geological Survey)",
          "· Ryan Kovach / David Schmetterling (Montana FWP)"
        )
      )
    )
  )
)

# ── Server ──────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Reset sliders to 0
  observeEvent(input$reset_btn, {
    updateSliderInput(session, "inv_z",  value = 0)
    updateSliderInput(session, "temp_z", value = 0)
    updateSliderInput(session, "flow_z", value = 0)
    updateSliderInput(session, "lkt_z",  value = 0)
  })

  # Reactive scenario λ
  scenario_lam <- reactive({
    compute_lambda(
      inv_z  = rep(input$inv_z,  nPop),
      temp_z = rep(input$temp_z, nPop),
      flow_z = rep(input$flow_z, nPop),
      lkt_z  = rep(input$lkt_z,  nPop)
    )
  })

  # Subtitle showing active scenario
  output$scenario_subtitle <- renderText({
    sprintf("INV = %.2f SD | Temp = %.2f SD | Flow = %.2f SD | LKT = %.2f SD",
            input$inv_z, input$temp_z, input$flow_z, input$lkt_z)
  })

  # ── Status table ───────────────────────────────────────────────
  output$status_tbl <- renderReactable({

    df <- st |>
      select(core, pop_name, pop_id,
             lambda_med, lambda_lo90, lambda_hi90,
             qext_p20, qext_p10,
             lam_inv_mean, lam_inv_best,
             lam_temp_1sd, lam_flow_1sd, lam_lkt_best)

    # Colour helper: white → blue for above-1 λ, red → white for below-1
    lambda_style <- function(value) {
      if (is.na(value)) return(list())
      if (value >= 1) {
        intensity <- min(1, (value - 1) / 0.5)
        bg <- rgb(1 - intensity * 0.6, 1 - intensity * 0.4, 1)
      } else {
        intensity <- min(1, (1 - value) / 0.5)
        bg <- rgb(1, 1 - intensity * 0.7, 1 - intensity * 0.7)
      }
      list(background = bg, fontWeight = "bold")
    }

    qext_style <- function(value) {
      if (is.na(value)) return(list())
      intensity <- min(1, value)
      bg <- rgb(1, 1 - intensity * 0.7, 1 - intensity * 0.7)
      list(background = bg)
    }

    reactable(
      df,
      filterable      = TRUE,
      searchable      = TRUE,
      sortable        = TRUE,
      resizable       = TRUE,
      defaultPageSize = 30,
      highlight       = TRUE,
      striped         = FALSE,
      columns = list(
        core        = colDef(name = "Core Area",   minWidth = 130, sticky = "left"),
        pop_name    = colDef(name = "Population",  minWidth = 130),
        pop_id      = colDef(name = "Pop ID",      maxWidth = 90),
        lambda_med  = colDef(name = "λ (median)",  format = colFormat(digits = 2),
                             style = lambda_style,  minWidth = 85),
        lambda_lo90 = colDef(name = "λ (5%)",      format = colFormat(digits = 2), maxWidth = 70),
        lambda_hi90 = colDef(name = "λ (95%)",     format = colFormat(digits = 2), maxWidth = 70),
        qext_p20    = colDef(name = "Qext (N<20)", format = colFormat(percent = TRUE, digits = 1),
                             style = qext_style,   maxWidth = 90),
        qext_p10    = colDef(name = "Qext (N<10)", format = colFormat(percent = TRUE, digits = 1),
                             style = qext_style,   maxWidth = 90),
        lam_inv_mean = colDef(name = "λ (INV avg)",  format = colFormat(digits = 2),
                              style = lambda_style,  maxWidth = 90),
        lam_inv_best = colDef(name = "λ (INV best)", format = colFormat(digits = 2),
                              style = lambda_style,  maxWidth = 90),
        lam_temp_1sd = colDef(name = "λ (Temp −1SD)", format = colFormat(digits = 2),
                              style = lambda_style,  maxWidth = 95),
        lam_flow_1sd = colDef(name = "λ (Flow +1SD)", format = colFormat(digits = 2),
                              style = lambda_style,  maxWidth = 95),
        lam_lkt_best = colDef(name = "λ (LKT best)", format = colFormat(digits = 2),
                              style = lambda_style,  na = "—", maxWidth = 90)
      ),
      columnGroups = list(
        colGroup(name = "Current Status",      columns = c("lambda_med", "lambda_lo90", "lambda_hi90")),
        colGroup(name = "Quasi-Extinction",    columns = c("qext_p20", "qext_p10")),
        colGroup(name = "Sensitivity (best-case λ)", columns = c("lam_inv_mean", "lam_inv_best",
                                                                   "lam_temp_1sd", "lam_flow_1sd",
                                                                   "lam_lkt_best"))
      )
    )
  })

  # ── Scenario plot ─────────────────────────────────────────────
  output$scenario_plot <- renderPlot({
    lam <- scenario_lam()

    df <- data.frame(
      pop_id    = meta$pop_names,
      pop_name  = pop_display,
      core      = core_names,
      lam_lo    = lam[, 1],
      lam_med   = lam[, 2],
      lam_hi    = lam[, 3],
      lam_cur   = current_lam[, 2]
    )

    if (input$core_filter != "All core areas") {
      df <- filter(df, core == input$core_filter)
    }

    df <- df |>
      arrange(lam_med) |>
      mutate(
        pop_label = paste0(pop_name, " (", pop_id, ")"),
        pop_label = factor(pop_label, levels = unique(pop_label)),
        direction = ifelse(lam_med >= 1, "Growing (λ \u2265 1)", "Declining (λ < 1)")
      )

    p <- ggplot(df, aes(y = pop_label)) +
      geom_vline(xintercept = 1, linetype = 2, color = "grey40", linewidth = 0.8)

    if (input$show_current) {
      p <- p + geom_point(aes(x = lam_cur), color = "grey55", size = 2.5, shape = 1, stroke = 1)
    }

    p <- p +
      geom_errorbarh(aes(xmin = lam_lo, xmax = lam_hi, color = direction),
                     height = 0.35, linewidth = 0.5, alpha = 0.6) +
      geom_point(aes(x = lam_med, color = direction), size = 3) +
      scale_color_manual(
        values = c("Growing (λ \u2265 1)" = "#1565C0", "Declining (λ < 1)" = "#C62828"),
        name   = NULL
      ) +
      scale_x_continuous(labels = function(x) sprintf("%.2f", x),
                         expand = expansion(mult = 0.07)) +
      facet_grid(core ~ ., scales = "free_y", space = "free_y") +
      labs(
        x        = "Intrinsic growth rate (\u03bb)",
        y        = NULL,
        subtitle = "Filled = scenario | Open circle = current conditions | Bars = 80% CI"
      ) +
      theme_bw(base_size = 11) +
      theme(
        legend.position  = "bottom",
        strip.text.y     = element_text(angle = 0, size = 8.5, face = "bold"),
        panel.grid.minor = element_blank(),
        axis.text.y      = element_text(size = 7.5),
        plot.subtitle    = element_text(color = "grey40", size = 10)
      )

    print(p)
  })
}

shinyApp(ui, server)
