### 03_MODELS.R ### 

# helpers ####

rescale_model_vars <- function(dat) {
  dat %>%
    mutate(
      cots_mean_ha_sc = as.numeric(scale(cots_mean_ha)),
      HC_sc = as.numeric(scale(HC)),
      MAC_sc = as.numeric(scale(MAC)),
      TUR_sc = as.numeric(scale(TUR))
    )
}and

sig_effects <- function(model, alpha = 0.05) {
  out <- broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term != "(Intercept)", !is.na(p.value), p.value < alpha) %>%
    mutate(sig = paste0(term, " (p = ", format_p(p.value), ")")) %>%
    pull(sig)
  
  if (length(out) == 0) return("none")
  paste(out, collapse = "; ")
}

export_model_list <- function(models, stats_out = stats_dir, fits_out = fits_dir) {
  purrr::iwalk(models, ~ {
    model_export(.x, .y, stats_out)
    save_obj(.x, paste0(.y, ".rds"), dir = fits_out)
  })
  invisible(models)
}

make_model_overview <- function(models, questions, responses, model_families, file_name, stats_out = stats_dir) {
  overview <- tibble(
    model = names(models),
    question = questions,
    response = responses,
    family = model_families,
    n_obs = purrr::map_int(models, nobs),
    AIC = purrr::map_dbl(models, AIC),
    significant_fixed_effects = purrr::map_chr(models, sig_effects)
  )
  
  write_csv(overview, file.path(stats_out, file_name))
  print(overview)
  invisible(overview)
}

run_dharma <- function(model, model_name, out_dir = plots_dir, n = 1000, save_plot = FALSE) {
  
  sim <- DHARMa::simulateResiduals(model, n = n)
  
  tests <- tibble(
    model = model_name,
    uniformity_p = DHARMa::testUniformity(sim, plot = FALSE)$p.value,
    dispersion_p = DHARMa::testDispersion(sim, plot = FALSE)$p.value,
    outlier_p = DHARMa::testOutliers(sim, plot = FALSE)$p.value,
    zero_inflation_p = DHARMa::testZeroInflation(sim, plot = FALSE)$p.value
  )
  
  if (isTRUE(save_plot)) {
    png(
      filename = file.path(out_dir, paste0(model_name, "_DHARMa.png")),
      width = 1800,
      height = 1400,
      res = 200
    )
    plot(sim)
    dev.off()
  }
  
  invisible(list(sim = sim, tests = tests))
}

run_dharma_tests <- function(model, model_name, n = 1000) {
  run_dharma(model, model_name, n = n, save_plot = FALSE)$tests
}

# A. BUILD ANALYSIS READY TABLES #####

### 01. COTS SITE CONTEXT ####
# Scale: site-level summary.
# Input: cots_survey, where each row is one CoTS survey.
# Output: one row per site_code.
# CoTS density is summarised across all available surveys per site.
# This creates site-level CoTS context for joining onto benthos, fish abundance, and fish bite tables.

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
    cots_15_mean = if_else(cots_mean_ha >= 15, "Above_15", "Below_15"),
    cots_15_max = if_else(cots_max_ha >= 15, "Above_15", "Below_15"),
    cots_15_mean = factor(cots_15_mean, levels = c("Below_15", "Above_15")),
    cots_15_max = factor(cots_15_max, levels = c("Below_15", "Above_15"))
  )

write_csv(cots_site_context, file.path(tables_dir, "cots_site_context.csv"))

### 02. SUBSTRATE SITE CONTEXT ####
# Scale: site-level summary.
# Input: substrate_wide, where each row is one substrate transect.
# Output: one row per site_code.
# Substrate percent cover is averaged across transects within each site.
# This creates site-level benthic context for joining onto fish abundance and fish bite tables.

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

### 03. BENTHOS MODEL TABLE ####
# Scale: transect-level model table.
# Input: substrate_wide, where each row is one substrate transect.
# Output: one row per substrate transect.
# Substrate variables remain at transect level.
# Site-level CoTS context is joined onto each transect by site_code.
# This table tests whether benthic condition differs among CoTS density contexts.

benthos_mod <- substrate_wide %>%
  left_join(
    cots_site_context %>%
      select(site_code, site_type, cots_mean_ha, cots_max_ha, cots_15_mean, cots_15_max),
    by = "site_code"
  ) %>%
  mutate(
    algal_cover = TUR + MAC,
    benthic_state = algal_cover - HC,
    site_type = factor(site_type, levels = c("Control", "High Density"))
  ) %>%
  rescale_model_vars()

write_csv(benthos_mod, file.path(tables_dir, "benthos_mod.csv"))

### 04. FISH ABUNDANCE MODEL TABLE ####
# Scale: survey/family-level model table.
# Input: fish_abund, where each row is one survey_id x researcher x family abundance record.
# Output: one row per original fish abundance record.
# Fish counts are not averaged before modelling.
# Site-level CoTS and benthic context are joined by site_code.
# Survey area is retained so fish_count can be modelled with offset(log(area_ha)).

fish_abund_mod <- fish_abund %>%
  left_join(
    cots_site_context %>%
      select(site_code, site_type, cots_mean_ha, cots_max_ha, cots_15_mean, cots_15_max),
    by = "site_code"
  ) %>%
  left_join(
    substrate_site_context %>%
      select(site_code, HC, TUR, MAC, algal_cover, benthic_state, n_substrate_transects),
    by = "site_code"
  ) %>%
  mutate(
    site_type = factor(site_type, levels = c("Control", "High Density")),
    family = fct_relevel(family, "Parrotfish", "Rabbitfish", "Butterflyfish"),
    sci_name = fct_relevel(sci_name, "Scarus spp.", "Siganus spp.", "Chaetodon spp.")
  ) %>%
  rescale_model_vars()

write_csv(fish_abund_mod, file.path(tables_dir, "fish_abund_mod.csv"))

### 05A. FISH BITE MODEL TABLE ####
# Scale: observation-level table used to build the feeding preference model.
# Input: fish_bites, where each row is one fish feeding observation period.
# Output: one row per original fish bite observation.
# Bite observations are not averaged.
# Site-level CoTS and benthic context are joined by site_code.
# Total bite rate is not modelled formally; feeding behaviour is modelled as substrate-specific bite use.

fish_bites_mod <- fish_bites %>%
  left_join(
    cots_site_context %>%
      select(site_code, site_type, cots_mean_ha, cots_max_ha, cots_15_mean, cots_15_max),
    by = "site_code"
  ) %>%
  left_join(
    substrate_site_context %>%
      select(site_code, HC, TUR, MAC, algal_cover, benthic_state, n_substrate_transects),
    by = "site_code"
  ) %>%
  mutate(
    site_type = factor(site_type, levels = c("Control", "High Density")),
    family = fct_relevel(family, "Parrotfish", "Rabbitfish", "Butterflyfish"),
    obs_id = row_number()
  ) %>%
  rescale_model_vars()

write_csv(fish_bites_mod, file.path(tables_dir, "fish_bites_mod.csv"))

### 05B. FEEDING PREFERENCE MODEL TABLE ####
# Scale: substrate-specific feeding table.
# Input: fish_bites_mod, where each row is one fish feeding observation period.
# Output: one row per observation x focal substrate.
# Focal substrates are HC, MAC, and TUR.

feeding_pref_mod <- fish_bites_mod %>%
  select(
    obs_id, site_code, site, site_type, date, buddy,
    family, species, sci_name, bites_total,
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
      bite_substrate == "turf_filamentous_algal_bites" ~ "TUR"
    ),
    avail_pct = case_when(
      substrate == "HC" ~ HC,
      substrate == "MAC" ~ MAC,
      substrate == "TUR" ~ TUR
    ),
    bites_other = bites_total - bites_sub,
    substrate = factor(substrate, levels = c("TUR", "MAC", "HC")),
    family = fct_relevel(family, "Parrotfish", "Rabbitfish", "Butterflyfish"),
    obs_id = factor(obs_id),
    avail_pct_sc = as.numeric(scale(avail_pct))
  ) %>%
  filter(!is.na(bites_sub), !is.na(bites_other), bites_other >= 0)

write_csv(feeding_pref_mod, file.path(tables_dir, "feeding_pref_mod.csv"))



### 06. ANALYSIS TABLE CHECKS ####
# Scale: diagnostic checks only.
# These checks confirm table dimensions, site/type balance, and missing values after joins.
# No data are summarised here for modelling, except as quick printed diagnostics.
table_dims <- tibble(
  table = c(
    "cots_site_context", "substrate_site_context", "benthos_mod",
    "fish_abund_mod", "fish_bites_mod", "feeding_pref_mod"
  ),
  rows = c(
    nrow(cots_site_context), nrow(substrate_site_context), nrow(benthos_mod),
    nrow(fish_abund_mod), nrow(fish_bites_mod), nrow(feeding_pref_mod)
  ),
  cols = c(
    ncol(cots_site_context), ncol(substrate_site_context), ncol(benthos_mod),
    ncol(fish_abund_mod), ncol(fish_bites_mod), ncol(feeding_pref_mod)
  )
)

write_csv(table_dims, file.path(tables_dir, "analysis_table_dimensions.csv"))

cat("\nCHECK 1: fish_bites_mod balance by site and CoTS density context\n")
print(fish_bites_mod %>% count(site_code, site_type))

cat("\nCHECK 2: fish_abund_mod balance by site, CoTS density context, and fish family\n")
print(fish_abund_mod %>% count(site_code, site_type, family))

cat("\nCHECK 3: benthos_mod balance by site and CoTS density context\n")
print(benthos_mod %>% count(site_code, site_type))

cat("\nCHECK 4: missing values in fish_bites_mod after joins\n")
print(
  fish_bites_mod %>%
    summarise(across(c(site_type, cots_mean_ha, HC, algal_cover, benthic_state), ~ sum(is.na(.))))
)

cat("\nCHECK 5: missing values in fish_abund_mod after joins\n")
print(
  fish_abund_mod %>%
    summarise(across(c(site_type, cots_mean_ha, HC, algal_cover, benthic_state), ~ sum(is.na(.))))
)

cat("\nCHECK 6: missing values in benthos_mod after joins\n")
print(
  benthos_mod %>%
    summarise(across(c(site_type, cots_mean_ha, benthic_state), ~ sum(is.na(.))))
)


# B. MODELS #####

### 07. MAIN MODELS ####
# Main model set.
# These models keep CoTS density central by using continuous mean CoTS density.
# Benthos models test CoTS density against HC, MAC, and TUR separately.
# Fish abundance tests whether family-level grazer abundance varies with CoTS density.

main_models <- list(
  m_benthos_hc = glmmTMB(
    HC ~ cots_mean_ha_sc + (1 | site_code),
    family = gaussian(),
    data = benthos_mod
  ),
  
  m_benthos_mac = glmmTMB(
    MAC ~ cots_mean_ha_sc + (1 | site_code),
    family = gaussian(),
    data = benthos_mod
  ),
  
  m_benthos_tur = glmmTMB(
    TUR ~ cots_mean_ha_sc + (1 | site_code),
    family = gaussian(),
    data = benthos_mod
  ),
  
  m_fish_abund = glmmTMB(
    fish_count ~ family * cots_mean_ha_sc + offset(log(area_ha)) +
      (1 | site_code) + (1 | date),
    family = nbinom2(),
    data = fish_abund_mod
  )
)

list2env(main_models, envir = .GlobalEnv)

purrr::iwalk(main_models, ~ {
  cat("\n==============================\n")
  cat(.y, "\n")
  cat("==============================\n")
  print(summary(.x))
})

export_model_list(main_models)

model_overview <- make_model_overview(
  models = main_models,
  questions = c(
    "Hard coral cover differs with variation in CoTS density",
    "Macroalgae cover differs with variation in CoTS density",
    "Turf algae cover differs with variation in CoTS density",
    "Fish abundance differs by family and continuous CoTS density"
  ),
  responses = c("HC", "MAC", "TUR", "fish_count"),
  model_families = c("gaussian", "gaussian", "gaussian", "negative binomial"),
  file_name = "model_overview.csv"
)

### 08. FEEDING PREFERENCE MODEL ####
# Formal feeding behaviour model.
# Question: does substrate-specific feeding change with substrate availability, fish family, and CoTS context?
# Scale: substrate-specific rows nested within each fish feeding observation.
# Response: bites to each focal substrate vs all other bites.
# Focal substrates: HC, MAC, TUR.

feeding_models <- list(
  m_feeding_pref_cots_sub = glmmTMB(
    cbind(bites_sub, bites_other) ~ substrate * family + substrate * cots_mean_ha_sc +
      avail_pct_sc + (1 | obs_id),
    family = betabinomial(link = "logit"),
    data = feeding_pref_mod
  )
)

list2env(feeding_models, envir = .GlobalEnv)

summary(m_feeding_pref_cots_sub)
export_model_list(feeding_models)

feeding_pref_overview <- make_model_overview(
  models = feeding_models,
  questions = "Substrate-specific feeding differs by substrate, family, availability, and CoTS context, with substrate-specific CoTS associations",
  responses = "cbind(bites_sub, bites_other)",
  model_families = "beta-binomial",
  file_name = "feeding_pref_overview.csv"
)

#### FEEDING PREFERENCE: ECOLOGICAL INTERACTION MODEL ####

m_feeding_pref_cots_sub_family <- glmmTMB(
  cbind(bites_sub, bites_other) ~ substrate * family * cots_mean_ha_sc +
    avail_pct_sc +
    (1 | obs_id),
  family = betabinomial(link = "logit"),
  data = feeding_pref_mod
)

saveRDS(
  m_feeding_pref_cots_sub_family,
  file.path(fits_dir, "m_feeding_pref_cots_sub_family.rds")
)


# C. SENSITIVITY ####

### 09. 15 COTS HA-1 THRESHOLD SENSITIVITY MODELS ####
# Sensitivity check only.
# These repeat the main benthos and fish abundance questions using the 15 CoTS ha-1 threshold.

threshold_models <- list(
  m_benthos_hc_15 = glmmTMB(
    HC ~ cots_15_mean + (1 | site_code),
    family = gaussian(),
    data = benthos_mod
  ),
  
  m_benthos_mac_15 = glmmTMB(
    MAC ~ cots_15_mean + (1 | site_code),
    family = gaussian(),
    data = benthos_mod
  ),
  
  m_benthos_tur_15 = glmmTMB(
    TUR ~ cots_15_mean + (1 | site_code),
    family = gaussian(),
    data = benthos_mod
  ),
  
  m_fish_abund_15 = glmmTMB(
    fish_count ~ family * cots_15_mean + offset(log(area_ha)) +
      (1 | site_code) + (1 | date),
    family = nbinom2(),
    data = fish_abund_mod
  )
)

list2env(threshold_models, envir = .GlobalEnv)
export_model_list(threshold_models)

model_overview_15 <- make_model_overview(
  models = threshold_models,
  questions = c(
    "Hard coral cover differs above vs below 15 CoTS ha-1",
    "Macroalgae cover differs above vs below 15 CoTS ha-1",
    "Turf algae cover differs above vs below 15 CoTS ha-1",
    "Fish abundance differs by family and 15 CoTS ha-1 threshold"
  ),
  responses = c("HC", "MAC", "TUR", "fish_count"),
  model_families = c("gaussian", "gaussian", "gaussian", "negative binomial"),
  file_name = "model_overview_15.csv"
)

### 10. AOW LEUK SENSITIVITY MODELS ####
# Sensitivity check only.
# AL is retained in the main analysis, but models are repeated without AL because its substrate context may behave differently.

run_al_sensitivity <- function(
    benthos_data,
    abund_data,
    feeding_pref_data,
    exclude_site = "AL",
    stats_out = stats_dir,
    fits_out = fits_dir
) {
  
  benthos_no_al <- benthos_data %>%
    filter(site_code != exclude_site) %>%
    rescale_model_vars()
  
  abund_no_al <- abund_data %>%
    filter(site_code != exclude_site) %>%
    rescale_model_vars()
  
  feeding_pref_no_al <- feeding_pref_data %>%
    filter(site_code != exclude_site) %>%
    mutate(
      cots_mean_ha_sc = as.numeric(scale(cots_mean_ha)),
      avail_pct_sc = as.numeric(scale(avail_pct)),
      obs_id = factor(obs_id)
    )
  
  models <- list(
    m_benthos_hc_no_al = glmmTMB(
      HC ~ cots_mean_ha_sc + (1 | site_code),
      family = gaussian(),
      data = benthos_no_al
    ),
    
    m_benthos_mac_no_al = glmmTMB(
      MAC ~ cots_mean_ha_sc + (1 | site_code),
      family = gaussian(),
      data = benthos_no_al
    ),
    
    m_benthos_tur_no_al = glmmTMB(
      TUR ~ cots_mean_ha_sc + (1 | site_code),
      family = gaussian(),
      data = benthos_no_al
    ),
    
    m_fish_abund_no_al = glmmTMB(
      fish_count ~ family * cots_mean_ha_sc + offset(log(area_ha)) +
        (1 | site_code) + (1 | date),
      family = nbinom2(),
      data = abund_no_al
    ),
    
    m_feeding_pref_no_al = glmmTMB(
      cbind(bites_sub, bites_other) ~ substrate * family + avail_pct_sc + cots_mean_ha_sc +
        (1 | obs_id),
      family = betabinomial(link = "logit"),
      data = feeding_pref_no_al
    )
  )
  
  export_model_list(models, stats_out = stats_out, fits_out = fits_out)
  models
}

al_models <- run_al_sensitivity(
  benthos_data = benthos_mod,
  abund_data = fish_abund_mod,
  feeding_pref_data = feeding_pref_mod,
  exclude_site = "AL"
)

al_sensitivity <- tibble(
  model_pair = c(
    "HC ~ CoTS",
    "MAC ~ CoTS",
    "TUR ~ CoTS",
    "Fish abundance ~ family * CoTS",
    "Feeding preference ~ substrate * family + availability + CoTS"
  ),
  full_model = c(names(main_models), names(feeding_models)),
  no_al_model = names(al_models),
  full_AIC = c(
    purrr::map_dbl(main_models, AIC),
    purrr::map_dbl(feeding_models, AIC)
  ),
  no_al_AIC = purrr::map_dbl(al_models, AIC),
  full_sig_effects = c(
    purrr::map_chr(main_models, sig_effects),
    purrr::map_chr(feeding_models, sig_effects)
  ),
  no_al_sig_effects = purrr::map_chr(al_models, sig_effects)
)

write_csv(al_sensitivity, file.path(stats_dir, "al_sensitivity_overview.csv"))
print(al_sensitivity)

# D. DIAGNOSTICS #####

diagnostic_tests_all <- bind_rows(
  purrr::imap_dfr(main_models, run_dharma_tests),
  purrr::imap_dfr(feeding_models, run_dharma_tests),
  purrr::imap_dfr(threshold_models, run_dharma_tests),
  purrr::imap_dfr(al_models, run_dharma_tests)
)

write_csv(diagnostic_tests_all, file.path(stats_dir, "dharma_diagnostic_tests.csv"))
print(diagnostic_tests_all)

# E. FINAL MODEL SUMMARIES ####

final_models <- c(
  main_models,
  feeding_models,
  threshold_models,
  al_models
)

final_summary_file <- file.path(stats_dir, "00_final_model_summaries.txt")

sink(final_summary_file)
tryCatch({
  cat("\n\nFINAL MODEL SUMMARIES\n")
  cat("=====================\n")
  
  purrr::iwalk(final_models, ~ {
    cat("\n\n------------------------------\n")
    cat(.y, "\n")
    cat("------------------------------\n")
    print(summary(.x))
  })
}, finally = {
  sink()
})

cat("\nFinal model summaries saved to:\n"); cat(final_summary_file, "\n")

final_model_fixed_effects <- purrr::imap_dfr(
  final_models,
  ~ broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE) %>%
    mutate(
      model = .y,
      p_formatted = if_else(
        p.value < 0.001,
        "<0.001",
        formatC(p.value, format = "f", digits = 3)
      )
    ),
  .id = NULL
) %>%
  select(
    model, term, estimate, std.error, conf.low, conf.high,
    statistic, p.value, p_formatted
  )

write_csv(
  final_model_fixed_effects,
  file.path(stats_dir, "final_model_fixed_effects.csv")
)

print(final_model_fixed_effects)
