library(shiny)
library(tidyverse)
library(DT)
library(blastula)

# ── Load pre-computed data ────────────────────────────────────────────────────
app_data       <- readRDS("app_data.RDS")
n_in           <- app_data$n_in
DatYears       <- app_data$DatYears
PatchWB_In     <- app_data$PatchWB_In
site_display   <- app_data$site_display
LocalPops_Num  <- app_data$LocalPops_Num
U_LocalPops_In <- app_data$U_LocalPops_In
pop_display    <- app_data$pop_display
pop_core_vec   <- app_data$pop_core_vec
U_CoreAreas_In <- app_data$U_CoreAreas_In

# ── Email config (from environment variables) ─────────────────────────────────
SMTP_USER  <- Sys.getenv("BLT_SMTP_USER",  unset = "")
SMTP_PASS  <- Sys.getenv("BLT_SMTP_PASSWORD", unset = "")
EMAIL_TO   <- Sys.getenv("BLT_EMAIL_TO",   unset = "")
email_configured <- nchar(SMTP_USER) > 0 && nchar(SMTP_PASS) > 0 && nchar(EMAIL_TO) > 0

# ── Helpers ───────────────────────────────────────────────────────────────────
send_submission_email <- function(df) {
  rows_html <- apply(df, 1, function(r) {
    sprintf(
      "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
      r["SiteName"], r["Year"], r["OldValue"], r["NewValue"],
      r["Submitter"], r["Agency"], r["CoreArea"], r["Notes"]
    )
  })
  body_html <- paste0(
    "<h3>Bull Trout Redd Count Submission</h3>",
    "<table border='1' cellpadding='5' cellspacing='0' style='border-collapse:collapse;font-family:sans-serif;font-size:13px;'>",
    "<thead><tr style='background:#dce6f1;'>",
    "<th>Site</th><th>Year</th><th>Old Value</th><th>New Value</th>",
    "<th>Submitter</th><th>Agency</th><th>Core Area</th><th>Notes</th>",
    "</tr></thead><tbody>",
    paste(rows_html, collapse = ""),
    "</tbody></table>",
    "<p style='color:#666;font-size:12px;margin-top:16px;'>",
    "Submitted: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " UTC</p>"
  )

  email <- compose_email(body = md(paste0(
    "## Bull Trout Redd Count Submission\n\n",
    "**Submitter:** ", df$Submitter[1], " (", df$Agency[1], ")  \n",
    "**Core area:** ", df$CoreArea[1], "  \n",
    "**Year:** ", df$Year[1], "  \n",
    "**Records submitted:** ", nrow(df), "\n\n",
    "See HTML version for full table."
  )))

  smtp_send(
    email,
    from    = SMTP_USER,
    to      = EMAIL_TO,
    subject = paste0("BLT Redd Submission \u2014 ", df$CoreArea[1], " ", df$Year[1]),
    credentials = creds_anonymous(
      host     = "smtp.gmail.com",
      port     = 587,
      use_ssl  = FALSE,
      user     = SMTP_USER,
      password = SMTP_PASS
    )
  )
}

log_to_csv <- function(df) {
  log_file <- "submissions_log.csv"
  write.table(df, log_file,
              sep = ",", row.names = FALSE,
              col.names = !file.exists(log_file),
              append = file.exists(log_file))
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: 'Helvetica Neue', Arial, sans-serif; }
    .sidebar-panel { background: #f8f9fa; padding: 16px; border-radius: 6px; }
    .site-row { padding: 8px 0; border-bottom: 1px solid #eee; }
    .current-val { color: #888; font-size: 12px; }
    .status-ok  { color: #2e7d32; font-weight: bold; margin-top: 8px; }
    .status-err { color: #c62828; font-weight: bold; margin-top: 8px; }
    h4 { color: #2c5282; margin-top: 18px; }
  "))),

  titlePanel(
    div(
      h3("Bull Trout Redd Count Data Entry",
         style = "margin-bottom: 2px;"),
      p("Bull Trout Interagency Recovery Team — Annual Data Verification",
        style = "color:#666; font-size:13px; margin-top:0;")
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      class = "sidebar-panel",

      h4("Your information"),
      textInput("submitter", "Full name", placeholder = "First Last"),
      textInput("agency",    "Agency",    placeholder = "USFWS / MFWP / USFS / MSU"),

      hr(),
      h4("Select population"),
      selectInput("core_area",  "Core area",  choices = sort(U_CoreAreas_In)),
      selectInput("population", "Population", choices = NULL),

      hr(),
      h4("Survey year"),
      numericInput("year", NULL,
                   value = as.integer(format(Sys.Date(), "%Y")),
                   min   = 1980,
                   max   = as.integer(format(Sys.Date(), "%Y")) + 1,
                   step  = 1),

      hr(),
      textAreaInput("notes", "Notes / reason for correction", rows = 3,
                    placeholder = "New survey data, correcting entry error, etc."),

      hr(),
      actionButton("submit", "Submit", class = "btn-primary", width = "100%"),
      uiOutput("submit_status")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        # ── Tab 1: Data entry ──────────────────────────────────────────────
        tabPanel(
          "Enter / Verify Data",
          br(),
          wellPanel(
            style = "background:#eef4fb; border:none;",
            p(strong("Instructions:"), "Enter redd counts for each survey site
              in the selected population for the chosen year.",
              tags$ul(
                tags$li("Leave a field", strong("blank"), "if the site was not surveyed this year."),
                tags$li("Enter", strong("0"), "if the site was surveyed and no redds were observed."),
                tags$li("The", em("current value"), "shown is from the model database — verify or correct it.")
              )
            )
          ),
          h4(textOutput("pop_title")),
          uiOutput("site_inputs"),
          br()
        ),

        # ── Tab 2: Historical data ─────────────────────────────────────────
        tabPanel(
          "Historical Data",
          br(),
          h4(textOutput("pop_title2")),
          plotOutput("ts_plot", height = "320px"),
          br(),
          h5("All survey records"),
          p(style = "color:#666; font-size:12px;",
            "Rows = survey sites, columns = years. Blank = not surveyed."),
          DTOutput("historical_table")
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Population choices update when core area changes ───────────────────────
  observe({
    ca_idx   <- which(U_CoreAreas_In == input$core_area)
    pop_idxs <- which(pop_core_vec == ca_idx)
    choices  <- setNames(as.character(pop_idxs), pop_display[pop_idxs])
    updateSelectInput(session, "population", choices = choices)
  })

  pop_idx <- reactive({
    req(input$population)
    as.integer(input$population)
  })

  site_idxs <- reactive({
    which(LocalPops_Num == pop_idx())
  })

  # ── Population title (shared across tabs) ──────────────────────────────────
  output$pop_title  <- renderText({ pop_display[pop_idx()] })
  output$pop_title2 <- renderText({ pop_display[pop_idx()] })

  # ── Dynamic site entry inputs ───────────────────────────────────────────────
  output$site_inputs <- renderUI({
    si     <- site_idxs()
    yr     <- input$year
    yr_col <- which(DatYears == yr)

    header <- fluidRow(
      column(6, tags$b("Survey site")),
      column(3, tags$b(paste("Redd count \u2014", yr))),
      column(3, tags$b("Current value in database"))
    )

    rows <- lapply(si, function(s) {
      cur <- if (length(yr_col) == 1) n_in[s, yr_col] else NA
      cur_label <- if (is.na(cur)) "not surveyed" else as.character(as.integer(cur))

      fluidRow(
        class = "site-row",
        column(6, tags$span(site_display[s])),
        column(3,
          numericInput(
            inputId = paste0("site_", s),
            label   = NULL,
            value   = NA,
            min     = 0,
            step    = 1
          )
        ),
        column(3,
          tags$span(cur_label, class = "current-val",
                    style = if (is.na(cur)) "color:#aaa;" else "color:#555;")
        )
      )
    })

    tagList(header, hr(), rows)
  })

  # ── Time series plot ────────────────────────────────────────────────────────
  output$ts_plot <- renderPlot({
    si <- site_idxs()

    df_sites <- as.data.frame(n_in[si, , drop = FALSE]) |>
      setNames(as.character(DatYears)) |>
      mutate(site = site_display[si]) |>
      pivot_longer(-site, names_to = "year", values_to = "redds") |>
      mutate(year = as.integer(year))

    df_total <- df_sites |>
      group_by(year) |>
      summarise(total   = sum(redds, na.rm = TRUE),
                any_obs = any(!is.na(redds)),
                .groups = "drop") |>
      mutate(total = if_else(any_obs, total, NA_real_))

    ggplot() +
      geom_col(data = filter(df_total, !is.na(total)),
               aes(x = year, y = total),
               fill = "steelblue", alpha = 0.2, width = 0.8) +
      geom_line(data = filter(df_total, !is.na(total)),
                aes(x = year, y = total),
                color = "steelblue", linewidth = 1) +
      geom_point(data = filter(df_sites, !is.na(redds)),
                 aes(x = year, y = redds, color = site),
                 size = 2.5, alpha = 0.8) +
      scale_color_brewer(palette = "Set2", name = NULL) +
      labs(x = "Year", y = "Redd count",
           caption = "Blue bars/line = population total  |  coloured points = individual survey sites") +
      theme_bw(base_size = 13) +
      theme(legend.position = "bottom",
            legend.text     = element_text(size = 9))
  })

  # ── Historical table ────────────────────────────────────────────────────────
  output$historical_table <- renderDT({
    si  <- site_idxs()
    mat <- n_in[si, , drop = FALSE]
    df  <- as.data.frame(mat)
    colnames(df) <- DatYears
    rownames(df) <- site_display[si]
    df[is.na(df)] <- ""

    datatable(df,
              rownames  = TRUE,
              options   = list(scrollX    = TRUE,
                               pageLength = 25,
                               dom        = "t"))
  })

  # ── Submit ──────────────────────────────────────────────────────────────────
  observeEvent(input$submit, {

    # Validation
    if (trimws(input$submitter) == "") {
      output$submit_status <- renderUI(
        tags$p("Please enter your name before submitting.", class = "status-err")
      )
      return()
    }

    si     <- site_idxs()
    yr     <- input$year
    yr_col <- which(DatYears == yr)

    rows <- lapply(si, function(s) {
      new_val <- input[[paste0("site_", s)]]
      if (is.null(new_val) || is.na(new_val)) return(NULL)
      old_val <- if (length(yr_col) == 1) n_in[s, yr_col] else NA
      data.frame(
        Timestamp  = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        Submitter  = trimws(input$submitter),
        Agency     = trimws(input$agency),
        CoreArea   = input$core_area,
        Population = U_LocalPops_In[pop_idx()],
        SiteID     = PatchWB_In[s],
        SiteName   = site_display[s],
        Year       = yr,
        OldValue   = if (is.na(old_val)) "not surveyed" else as.integer(old_val),
        NewValue   = as.integer(new_val),
        Notes      = trimws(input$notes),
        stringsAsFactors = FALSE
      )
    })

    rows <- Filter(Negate(is.null), rows)

    if (length(rows) == 0) {
      output$submit_status <- renderUI(
        tags$p("No values entered \u2014 nothing to submit.", class = "status-err")
      )
      return()
    }

    df_submit <- bind_rows(rows)

    success <- tryCatch({
      if (email_configured) {
        send_submission_email(df_submit)
      } else {
        log_to_csv(df_submit)
      }
      TRUE
    }, error = function(e) {
      message("Submission error: ", e$message)
      FALSE
    })

    if (success) {
      n   <- nrow(df_submit)
      how <- if (email_configured) "emailed to the data team" else "logged locally"
      output$submit_status <- renderUI(
        tags$p(
          sprintf("\u2713 %d record%s submitted \u2014 %s. Thank you, %s!",
                  n, if (n > 1) "s" else "", how, df_submit$Submitter[1]),
          class = "status-ok"
        )
      )
    } else {
      output$submit_status <- renderUI(
        tags$p(
          "Submission failed. Please contact the model team.",
          class = "status-err"
        )
      )
    }
  })
}

shinyApp(ui, server)
