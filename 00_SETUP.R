### 00. SETUP #### 
suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2); library(glmmTMB); library(purrr); library(broom.mixed)
})

# Enter the date for this analysis
analysis_date <- "2026.06.09" # update each time 
# Alternative:
# analysis_date <- format(Sys.Date(), "%Y.%m.%d")

### directories ####

data_raw_dir   <- "data_raw"
data_clean_dir <- "data_clean"
docs_dir       <- "docs"

dir.create(data_raw_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(data_clean_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(docs_dir,       showWarnings = FALSE, recursive = TRUE)

### output structure ####

output_dir <- paste0("Analysis_", analysis_date)

fits_dir   <- file.path(output_dir, "fits")
plots_dir  <- file.path(output_dir, "plots")
stats_dir  <- file.path(output_dir, "stats")
summ_dir   <- file.path(output_dir, "summaries")
eda_dir    <- file.path(output_dir, "eda")
tables_dir <- file.path(output_dir, "tables")

purrr::walk(
  c(output_dir, fits_dir, plots_dir, stats_dir, summ_dir, eda_dir, tables_dir),
  ~ dir.create(.x, showWarnings = FALSE, recursive = TRUE)
)


### custom theme and colour palettes ####  
theme_clean <- theme_minimal(base_family = "Arial") +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    plot.title = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# colour palettes
fg_cols <- c(
  "Grazer" = "#66c2a4",
  "Invertivore" = "#41b6c4",
  "Mesopredator" = "#2c7fb8",
  "HTLP" = "#253494"
) 

reef_cols <- c(
  "High Density" = "#253494",  # deep ocean blue
  "Low Density"  = "#66c2a4"   # seafoam
)


### helpers ####
format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3))
}

save_obj <- function(x, filename, dir = summ_dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(x, file.path(dir, filename))
  invisible(TRUE)
}

model_export <- function(model, model_name, output_dir, sigfigs = 3) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  library(broom.mixed)
  library(dplyr)
  
  fx <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  
  fx_out <- fx %>%
    transmute(
      Effect   = term,
      Estimate = signif(estimate, sigfigs),
      SE       = signif(std.error, sigfigs),
      CI       = paste0("[",
                        signif(conf.low,  sigfigs), ", ",
                        signif(conf.high, sigfigs), "]"),
      p_value  = ifelse(p.value < 0.001,
                        "<0.001",
                        formatC(p.value, format = "f", digits = 3))
    )
  
  write.csv(
    fx_out,
    file = file.path(output_dir, paste0(model_name, "_summarytable.csv")),
    row.names = FALSE
  )
  
  invisible(fx_out)
}

# only using the three grazer taxa 
functional_taxa <- tibble::tribble(
  ~Species,            ~Trophic_Group, ~Genus,                ~Species_epithet, ~sci_name,
  "Parrotfish",        "Grazer",          "Scarus",              "spp.",           "Scarus spp.",
  "Rabbitfish",        "Grazer",          "Siganus",             "spp.",           "Siganus spp.",
  "Butterflyfish",     "Grazer",          "Chaetodon",           "spp.",          "Chaetodon spp."
)


