### 00. SETUP #### 
suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2); library(glmmTMB); library(purrr); library(broom.mixed); library(janitor); library(DHARMa)
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
  "High Density" = "#FF9683",  # orange siesta
  "Control"  = "#95B971" # Atlantis Green
)



substrate_cols_pal <- c(
  "HC"  = "#ff8c94",  # coral pink
  "SC"  = "#E2B5F8",  # soft coral / pale peach
  "AB"  = "#d8c39f",  # abiotic / sand
  "TUR" = "#95B971",  # turf / seafoam green
  "MAC" = "#238b45",  # macroalgae / deep green
  "SP"  = "#6372EB",  # sponge / purple
  "OB"  = "#8c7a6b",  # other biotic / reef brown
  "UKN" = "#bdbdbd"   # unknown / grey
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

### lookup tables #### 

# abundance data: coarse grazer groups
functional_taxa <- tibble::tribble(
  ~family,          ~trophic_group, ~taxon_family,        ~genus,       ~species_epithet, ~sci_name,
  "Parrotfish",     "Grazer",       "Labridae: Scarinae", "Scarus",     "spp.",           "Scarus spp.",
  "Rabbitfish",     "Grazer",       "Siganidae",          "Siganus",    "spp.",           "Siganus spp.",
  "Butterflyfish",  "Grazer",       "Chaetodontidae",     "Chaetodon",  "spp.",           "Chaetodon spp."
)

# fish bites data: species-level observations
fish_bite_species_lookup <- tibble::tribble(
  ~family,          ~species,       ~trophic_group, ~taxon_family,        ~genus,       ~species_epithet, ~sci_name,                    ~taxon_notes,
  "Butterflyfish",  "8-banded",     "Grazer",       "Chaetodontidae",     "Chaetodon",  "octofasciatus",  "Chaetodon octofasciatus",  "",
  "Butterflyfish",  "Weibel",       "Grazer",       "Chaetodontidae",     "Chaetodon",  "wiebeli",        "Chaetodon wiebeli",        "",
  "Parrotfish",     "Blue-barred",  "Grazer",       "Labridae: Scarinae", "Scarus",     "ghobban",        "Scarus ghobban",           "",
  "Parrotfish",     "Purple",       "Grazer",       "Labridae: Scarinae", "Scarus",     "globiceps",      "Scarus globiceps",         "check ID",
  "Parrotfish",     "Surf",         "Grazer",       "Labridae: Scarinae", "Scarus",     "rivulatus",      "Scarus rivulatus",         "check ID",
  "Rabbitfish",     "Virgate",      "Grazer",       "Siganidae",          "Siganus",    "virgatus",       "Siganus virgatus",         "",
  "Rabbitfish",     "Yellow",       "Grazer",       "Siganidae",          "Siganus",    "guttatus",       "Siganus guttatus",         "check ID"
)


