# Bull Trout Redd Count Data Entry App

Shiny app for Bull Trout Interagency Recovery Team biologists to verify,
correct, and add new annual redd count data. Submissions are emailed to the
designated data manager and incorporated before the next annual model run.

---

## First-time setup

### 1. Generate the app data file

From the project root in R:

```r
source("Code/PrepareShinyData.R")
```

This creates `ShinyApp/app_data.RDS`. Re-run this every year after the model
is updated with new data.

### 2. Set up a Gmail account for sending submissions

1. Create (or designate) a Gmail account for the project, e.g. `bulltroutssp@gmail.com`
2. Enable 2-Factor Authentication on that account
3. Go to **Google Account → Security → 2-Step Verification → App passwords**
4. Create an app password named "BullTroutShiny" — copy the 16-character password

### 3. Configure environment variables

**For local testing** — add to `~/.Renviron` (one time):

```
BLT_SMTP_USER=bulltroutssp@gmail.com
BLT_SMTP_PASSWORD=xxxx xxxx xxxx xxxx
BLT_EMAIL_TO=your.email@msu.edu
```

Restart R after editing `.Renviron`.

**For shinyapps.io** — set via the dashboard:
1. Go to [shinyapps.io](https://www.shinyapps.io) → your account → the app
2. **Settings → Environment Variables**
3. Add the same three variables: `BLT_SMTP_USER`, `BLT_SMTP_PASSWORD`, `BLT_EMAIL_TO`

> **Note:** If env vars are not set, submissions are written to
> `ShinyApp/submissions_log.csv` instead (useful for local development).

---

## Running locally

```r
shiny::runApp("ShinyApp")
```

---

## Deploying to shinyapps.io

Install `rsconnect` if needed:

```r
install.packages("rsconnect")
```

Authenticate (one time):

```r
rsconnect::setAccountInfo(
  name   = "your-shinyapps-username",
  token  = "YOUR_TOKEN",
  secret = "YOUR_SECRET"
)
```

Deploy:

```r
rsconnect::deployApp(
  appDir  = "ShinyApp",
  appName = "BullTroutReddEntry"
)
```

The app will be live at:
`https://your-username.shinyapps.io/BullTroutReddEntry/`

Share this URL with recovery team biologists.

---

## Annual update workflow

Each year after new redd count data is incorporated into the model:

1. Update the raw data CSV (`Data/Raw/USGS_Harmonized_BLT_Redd_Data_*.csv`)
2. Re-run `source("Code/PrepareShinyData.R")` to regenerate `app_data.RDS`
3. Re-deploy: `rsconnect::deployApp("ShinyApp", appName = "BullTroutReddEntry")`
4. Notify biologists that the app is open for the new survey year

---

## How submissions work

1. Biologist opens the app URL in their browser
2. Selects their core area → population → survey year
3. Enters redd counts for each survey site (or corrects existing values)
4. Clicks **Submit**
5. An email is sent to `BLT_EMAIL_TO` with a table of all submitted values
6. Data manager reviews the email and applies corrections to the database
   before the next model run

---

## Required R packages

```r
install.packages(c("shiny", "tidyverse", "DT", "blastula"))
```
