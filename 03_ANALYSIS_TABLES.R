#### 03_ANALYSIS_TABLES.R ####
### Build manuscript analysis tables for CoTS, benthos, fish abundance, and feeding models

source("00_SETUP.R")

#### 01. LOAD CLEAN DATA ####

fish_bites <- read_csv(file.path(data_clean_dir, "fish_bites_clean.csv"), show_col_types = FALSE)
fish_abund <- read_csv(file.path(data_clean_dir, "fish_abund_clean.csv"), show_col_types = FALSE)
cots_survey <- read_csv(file.path(data_clean_dir, "cots_survey_clean.csv"), show_col_types = FALSE)
substrate_wide <- read_csv(file.path(data_clean_dir, "substrate_transect_cover_wide_clean.csv"), show_col_types = FALSE)
substrate_transect <- read_csv(file.path(data_clean_dir, "substrate_transect_cover_clean.csv"), show_col_types = FALSE)

fish_bites <- fish_bites %>% mutate(date = as.Date(date))
fish_abund <- fish_abund %>% mutate(date = as.Date(date))
cots_survey <- cots_survey %>% mutate(date = as.Date(date))

#### 02. HELPERS ####

rescale_model_vars <- function(dat) {
  dat %>%
    mutate(
      cots_mean_ha_sc = as.numeric(scale(cots_mean_ha)),
      HC_sc = as.numeric(scale(HC)),
      MAC_sc = as.numeric(scale(MAC)),
      TUR_sc = as.numeric(scale(TUR))
    )
}

#### 03. COTS SITE CONTEXT ####
# One row per site. Mean CoTS density is the primary continuous predictor.
# The 15 CoTS ha-1 threshold is retained for sensitivity analyses.

cots_site_context <- cots_survey %>%
  group_by(site_code, site, site_type) %>%
  summarise(
    cots_mean_ha = mean(cots_density_ha, na.rm = TRUE),
    cots_median_ha = median(cots_density_ha, na.rm = TRUE),
    cots_max_ha = max(cots_density_ha, na.rm = TRUE),
    n_cots_surveys = n(),
    .groups = "drop"
  ) %>%
  mutate(
    cots_15_mean = factor(if_else(cots_mean_ha >= 15, "Above_15", "Below_15"),
                          levels = c("Below_15", "Above_15")),
    cots_15_max = factor(if_else(cots_max_ha >= 15, "Above_15", "Below_15"),
                         levels = c("Below_15", "Above_15")),
    site_type = factor(site_type, levels = c("Control", "High Density"))
  )

write_csv(cots_site_context, file.path(tables_dir, "cots_site_context.csv"))

#### 04. SUBSTRATE SITE CONTEXT ####
# Site-level substrate availability for feeding models.

substrate_site_context <- substrate_wide %>%
  group_by(site_code, site) %>%
  summarise(
    HC = mean(HC, na.rm = TRUE),
    TUR = mean(TUR, na.rm = TRUE),
    MAC = mean(MAC, na.rm = TRUE),
    AB = mean(AB, na.rm = TRUE),
    OB = mean(OB, na.rm = TRUE),
    n_substrate_transects = n(),
    .groups = "drop"
  ) %>%
  mutate(
    algal_cover = TUR + MAC,
    benthic_state = algal_cover - HC
  )

write_csv(substrate_site_context, file.path(tables_dir, "substrate_site_context.csv"))

#### 05. BENTHOS MODEL TABLE ####
# Transect-level percent cover with site-level CoTS context.

benthos_mod <- substrate_wide %>%
  left_join(
    cots_site_context %>%
      select(site_code, site_type, cots_mean_ha, cots_median_ha, cots_max_ha, cots_15_mean, cots_15_max),
    by = "site_code"
  ) %>%
  mutate(
    algal_cover = TUR + MAC,
    benthic_state = algal_cover - HC,
    site_type = factor(site_type, levels = c("Control", "High Density"))
  ) %>%
  rescale_model_vars()

write_csv(benthos_mod, file.path(tables_dir, "benthos_mod.csv"))

#### 06. FISH ABUNDANCE MODEL TABLE ####
# Raw family-level counts are retained. Density is descriptive only.

fish_abund_mod <- fish_abund %>%
  left_join(
    cots_site_context %>%
      select(site_code, site_type, cots_mean_ha, cots_median_ha, cots_max_ha, cots_15_mean, cots_15_max),
    by = "site_code"
  ) %>%
  left_join(
    substrate_site_context %>%
      select(site_code, HC, TUR, MAC, algal_cover, benthic_state, n_substrate_transects),
    by = "site_code"
  ) %>%
  mutate(
    site_type = factor(site_type, levels = c("Control", "High Density")),
    family = fct_relevel(factor(family), "Parrotfish", "Rabbitfish", "Butterflyfish"),
    sci_name = fct_relevel(factor(sci_name), "Scarus spp.", "Siganus spp.", "Chaetodon spp.")
  ) %>%
  rescale_model_vars()

write_csv(fish_abund_mod, file.path(tables_dir, "fish_abund_mod.csv"))

#### 07. FISH BITE MODEL TABLE ####
# Observation-level feeding table for total bite-rate and feeding preference models.

fish_bites_mod <- fish_bites %>%
  left_join(
    cots_site_context %>%
      select(site_code, site_type, cots_mean_ha, cots_median_ha, cots_max_ha, cots_15_mean, cots_15_max),
    by = "site_code"
  ) %>%
  left_join(
    substrate_site_context %>%
      select(site_code, HC, TUR, MAC, algal_cover, benthic_state, n_substrate_transects),
    by = "site_code"
  ) %>%
  mutate(
    site_type = factor(site_type, levels = c("Control", "High Density")),
    family = fct_relevel(factor(family), "Parrotfish", "Rabbitfish", "Butterflyfish"),
    obs_id = row_number()
  ) %>%
  rescale_model_vars()

write_csv(fish_bites_mod, file.path(tables_dir, "fish_bites_mod.csv"))

#### 08. FEEDING PREFERENCE MODEL TABLE ####
# One row per observation x focal substrate. Focal substrates are TUR, MAC, and HC.

feeding_pref_mod <- fish_bites_mod %>%
  select(
    obs_id, site_code, site, site_type, date, buddy,
    family, species, sci_name, bites_total, minutes_obs,
    live_coral_bites, macro_algal_bites, turf_filamentous_algal_bites,
    HC, MAC, TUR, cots_mean_ha, cots_mean_ha_sc
  ) %>%
  pivot_longer(
    cols = c(live_coral_bites, macro_algal_bites, turf_filamentous_algal_bites),
    names_to = "bite_substrate",
    values_to = "bites_sub"
  ) %>%
  mutate(
    substrate = case_when(
      bite_substrate == "live_coral_bites" ~ "HC",
      bite_substrate == "macro_algal_bites" ~ "MAC",
      bite_substrate == "turf_filamentous_algal_bites" ~ "TUR",
      TRUE ~ NA_character_
    ),
    avail_pct = case_when(
      substrate == "HC" ~ HC,
      substrate == "MAC" ~ MAC,
      substrate == "TUR" ~ TUR,
      TRUE ~ NA_real_
    ),
    bites_other = bites_total - bites_sub,
    substrate = factor(substrate, levels = c("TUR", "MAC", "HC")),
    family = fct_relevel(factor(family), "Parrotfish", "Rabbitfish", "Butterflyfish"),
    obs_id = factor(obs_id),
    avail_pct_sc = as.numeric(scale(avail_pct))
  ) %>%
  filter(!is.na(bites_sub), !is.na(bites_other), bites_other >= 0)

write_csv(feeding_pref_mod, file.path(tables_dir, "feeding_pref_mod.csv"))

#### 09. ANALYSIS TABLE CHECKS ####

table_dims <- tibble(
  table = c("cots_site_context", "substrate_site_context", "benthos_mod", "fish_abund_mod", "fish_bites_mod", "feeding_pref_mod"),
  rows = c(nrow(cots_site_context), nrow(substrate_site_context), nrow(benthos_mod), nrow(fish_abund_mod), nrow(fish_bites_mod), nrow(feeding_pref_mod)),
  cols = c(ncol(cots_site_context), ncol(substrate_site_context), ncol(benthos_mod), ncol(fish_abund_mod), ncol(fish_bites_mod), ncol(feeding_pref_mod))
)

write_csv(table_dims, file.path(tables_dir, "analysis_table_dimensions.csv"))
print(table_dims)

join_checks <- list(
  fish_bites_mod = fish_bites_mod %>% summarise(across(c(site_type, cots_mean_ha, HC, TUR, MAC), ~ sum(is.na(.)))),
  fish_abund_mod = fish_abund_mod %>% summarise(across(c(site_type, cots_mean_ha, HC, TUR, MAC), ~ sum(is.na(.)))),
  benthos_mod = benthos_mod %>% summarise(across(c(site_type, cots_mean_ha, HC, TUR, MAC), ~ sum(is.na(.))))
)

saveRDS(join_checks, file.path(stats_dir, "analysis_table_join_checks.rds"))
