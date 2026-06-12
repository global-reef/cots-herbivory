### 04_MODELS.R ####
### Manuscript model fitting for CoTS, benthos, fish abundance, and feeding responses

source("00_SETUP.R")

#### 01. LOAD ANALYSIS TABLES ####

benthos_mod <- read_csv(file.path(tables_dir, "benthos_mod.csv"), show_col_types = FALSE) %>%
  mutate(site_code = factor(site_code), site_type = factor(site_type, levels = c("Control", "High Density")))

fish_abund_mod <- read_csv(file.path(tables_dir, "fish_abund_mod.csv"), show_col_types = FALSE) %>%
  mutate(
    date = as.Date(date),
    site_code = factor(site_code),
    site_type = factor(site_type, levels = c("Control", "High Density")),
    family = fct_relevel(factor(family), "Parrotfish", "Rabbitfish", "Butterflyfish")
  )

fish_bites_mod <- read_csv(file.path(tables_dir, "fish_bites_mod.csv"), show_col_types = FALSE) %>%
  mutate(
    date = as.Date(date),
    site_code = factor(site_code),
    site_type = factor(site_type, levels = c("Control", "High Density")),
    family = fct_relevel(factor(family), "Parrotfish", "Rabbitfish", "Butterflyfish")
  )

feeding_pref_mod <- read_csv(file.path(tables_dir, "feeding_pref_mod.csv"), show_col_types = FALSE) %>%
  mutate(
    date = as.Date(date),
    site_code = factor(site_code),
    site_type = factor(site_type, levels = c("Control", "High Density")),
    substrate = factor(substrate, levels = c("TUR", "MAC", "HC")),
    family = fct_relevel(factor(family), "Parrotfish", "Rabbitfish", "Butterflyfish"),
    obs_id = factor(obs_id)
  )

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

run_dharma_tests <- function(model, model_name, n = 1000) {
  sim <- DHARMa::simulateResiduals(model, n = n)
  tibble(
    model = model_name,
    uniformity_p = DHARMa::testUniformity(sim, plot = FALSE)$p.value,
    dispersion_p = DHARMa::testDispersion(sim, plot = FALSE)$p.value,
    outlier_p = DHARMa::testOutliers(sim, plot = FALSE)$p.value,
    zero_inflation_p = DHARMa::testZeroInflation(sim, plot = FALSE)$p.value
  )
}

#### 03. MAIN BENTHIC AND FISH ABUNDANCE MODELS ####

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
export_model_list(main_models)

make_model_overview(
  models = main_models,
  questions = c(
    "Hard coral cover changes with continuous mean site-level CoTS density",
    "Macroalgae cover changes with continuous mean site-level CoTS density",
    "Turf algae cover changes with continuous mean site-level CoTS density",
    "Fish abundance differs among families across continuous mean site-level CoTS density"
  ),
  responses = c("HC", "MAC", "TUR", "fish_count"),
  model_families = c("gaussian", "gaussian", "gaussian", "negative binomial"),
  file_name = "model_overview_main.csv"
)

#### 04. MAIN SUBSTRATE-SPECIFIC FEEDING MODEL ####
# This is the main feeding allocation model used for Figure 14 and text estimates.

feeding_models <- list(
  m_feeding_pref_cots_sub = glmmTMB(
    cbind(bites_sub, bites_other) ~ substrate * family + substrate * cots_mean_ha_sc +
      avail_pct_sc + (1 | obs_id),
    family = betabinomial(link = "logit"),
    data = feeding_pref_mod
  )
)

list2env(feeding_models, envir = .GlobalEnv)
export_model_list(feeding_models)

make_model_overview(
  models = feeding_models,
  questions = "Substrate-specific bite probability differs by focal substrate, fish family, substrate availability, and CoTS density",
  responses = "cbind(bites_sub, bites_other)",
  model_families = "beta-binomial",
  file_name = "model_overview_feeding_preference.csv"
)

#### 05. FEEDING INTERACTION CHECKS REPORTED IN TEXT ####
# These check whether CoTS-related feeding responses differ among fish families.
# They are retained because the results text reports limited support for family-specific interactions.

feeding_interaction_models <- list(
  m_feeding_pref_cots_family = glmmTMB(
    cbind(bites_sub, bites_other) ~ substrate * family + family * cots_mean_ha_sc +
      avail_pct_sc + (1 | obs_id),
    family = betabinomial(link = "logit"),
    data = feeding_pref_mod
  ),
  m_feeding_pref_cots_sub_family = glmmTMB(
    cbind(bites_sub, bites_other) ~ substrate * family * cots_mean_ha_sc +
      avail_pct_sc + (1 | obs_id),
    family = betabinomial(link = "logit"),
    data = feeding_pref_mod
  )
)

list2env(feeding_interaction_models, envir = .GlobalEnv)
export_model_list(feeding_interaction_models)

feeding_interaction_compare <- AIC(
  m_feeding_pref_cots_sub,
  m_feeding_pref_cots_family,
  m_feeding_pref_cots_sub_family
) %>%
  as_tibble(rownames = "model") %>%
  arrange(AIC) %>%
  mutate(delta_AIC = AIC - min(AIC))

write_csv(feeding_interaction_compare, file.path(stats_dir, "feeding_interaction_model_compare.csv"))
print(feeding_interaction_compare)

#### 06. SECONDARY TOTAL BITE-RATE HURDLE MODELS ####
# Secondary feeding intensity analysis reported in Results.
# Part 1 models feeding occurrence. Part 2 models positive bite counts only.

fish_bites_rate_mod <- fish_bites_mod %>%
  mutate(
    any_feeding = as.integer(bites_total > 0),
    minutes_obs = as.numeric(minutes_obs)
  ) %>%
  filter(!is.na(bites_total), !is.na(minutes_obs), minutes_obs > 0, !is.na(site_type), !is.na(family))

fish_bites_positive <- fish_bites_rate_mod %>%
  filter(bites_total > 0)

bite_rate_models <- list(
  m_bite_any = glmmTMB(
    any_feeding ~ family * site_type + (1 | site_code) + (1 | date),
    family = binomial(link = "logit"),
    data = fish_bites_rate_mod
  ),
  m_bite_positive = glmmTMB(
    bites_total ~ family * site_type + offset(log(minutes_obs)) +
      (1 | site_code) + (1 | date),
    family = nbinom2(),
    data = fish_bites_positive
  )
)

list2env(bite_rate_models, envir = .GlobalEnv)
export_model_list(bite_rate_models)

make_model_overview(
  models = bite_rate_models,
  questions = c(
    "Probability of any feeding differs by family and CoTS site type",
    "Positive bite rate differs by family and CoTS site type"
  ),
  responses = c("any_feeding", "positive bites_total with offset(log(minutes_obs))"),
  model_families = c("binomial", "negative binomial"),
  file_name = "model_overview_bite_rate_hurdle.csv"
)
### sensitivity models ####
#### 07. 15 COTS HA-1 THRESHOLD SENSITIVITY MODELS ####

threshold_models <- list(
  m_benthos_hc_15 = glmmTMB(HC ~ cots_15_mean + (1 | site_code), family = gaussian(), data = benthos_mod),
  m_benthos_mac_15 = glmmTMB(MAC ~ cots_15_mean + (1 | site_code), family = gaussian(), data = benthos_mod),
  m_benthos_tur_15 = glmmTMB(TUR ~ cots_15_mean + (1 | site_code), family = gaussian(), data = benthos_mod),
  m_fish_abund_15 = glmmTMB(
    fish_count ~ family * cots_15_mean + offset(log(area_ha)) +
      (1 | site_code) + (1 | date),
    family = nbinom2(),
    data = fish_abund_mod
  )
)

list2env(threshold_models, envir = .GlobalEnv)
export_model_list(threshold_models)

make_model_overview(
  models = threshold_models,
  questions = c(
    "Hard coral cover differs above vs below 15 CoTS ha-1",
    "Macroalgae cover differs above vs below 15 CoTS ha-1",
    "Turf algae cover differs above vs below 15 CoTS ha-1",
    "Fish abundance differs by family and 15 CoTS ha-1 threshold"
  ),
  responses = c("HC", "MAC", "TUR", "fish_count"),
  model_families = c("gaussian", "gaussian", "gaussian", "negative binomial"),
  file_name = "model_overview_threshold_15.csv"
)

#### 08. AOW LEUK EXCLUSION SENSITIVITY MODELS ####

benthos_no_al <- benthos_mod %>% filter(site_code != "AL") %>% rescale_model_vars()
fish_abund_no_al <- fish_abund_mod %>% filter(site_code != "AL") %>% rescale_model_vars()
feeding_pref_no_al <- feeding_pref_mod %>%
  filter(site_code != "AL") %>%
  mutate(
    cots_mean_ha_sc = as.numeric(scale(cots_mean_ha)),
    avail_pct_sc = as.numeric(scale(avail_pct)),
    obs_id = factor(obs_id)
  )

al_models <- list(
  m_benthos_hc_no_al = glmmTMB(HC ~ cots_mean_ha_sc + (1 | site_code), family = gaussian(), data = benthos_no_al),
  m_benthos_mac_no_al = glmmTMB(MAC ~ cots_mean_ha_sc + (1 | site_code), family = gaussian(), data = benthos_no_al),
  m_benthos_tur_no_al = glmmTMB(TUR ~ cots_mean_ha_sc + (1 | site_code), family = gaussian(), data = benthos_no_al),
  m_fish_abund_no_al = glmmTMB(
    fish_count ~ family * cots_mean_ha_sc + offset(log(area_ha)) +
      (1 | site_code) + (1 | date),
    family = nbinom2(),
    data = fish_abund_no_al
  ),
  m_feeding_pref_cots_sub_no_al = glmmTMB(
    cbind(bites_sub, bites_other) ~ substrate * family + substrate * cots_mean_ha_sc +
      avail_pct_sc + (1 | obs_id),
    family = betabinomial(link = "logit"),
    data = feeding_pref_no_al
  )
)

list2env(al_models, envir = .GlobalEnv)
export_model_list(al_models)

al_sensitivity <- tibble(
  model_pair = c("HC ~ CoTS", "MAC ~ CoTS", "TUR ~ CoTS", "Fish abundance ~ family * CoTS", "Feeding preference ~ substrate * CoTS"),
  full_model = c(names(main_models), names(feeding_models)),
  no_al_model = names(al_models),
  full_AIC = c(purrr::map_dbl(main_models, AIC), purrr::map_dbl(feeding_models, AIC)),
  no_al_AIC = purrr::map_dbl(al_models, AIC),
  full_sig_effects = c(purrr::map_chr(main_models, sig_effects), purrr::map_chr(feeding_models, sig_effects)),
  no_al_sig_effects = purrr::map_chr(al_models, sig_effects)
)

write_csv(al_sensitivity, file.path(stats_dir, "al_sensitivity_overview.csv"))
print(al_sensitivity)

#### 09. DIAGNOSTICS ####

diagnostic_models <- c(
  main_models,
  feeding_models,
  feeding_interaction_models,
  bite_rate_models,
  threshold_models,
  al_models
)

diagnostic_tests_all <- purrr::imap_dfr(diagnostic_models, run_dharma_tests)
write_csv(diagnostic_tests_all, file.path(stats_dir, "dharma_diagnostic_tests.csv"))
print(diagnostic_tests_all)

#### 10. FINAL MODEL SUMMARIES ####

final_summary_file <- file.path(stats_dir, "00_final_model_summaries.txt")

sink(final_summary_file)
tryCatch({
  cat("FINAL MODEL SUMMARIES\n")
  cat("=====================\n")
  purrr::iwalk(diagnostic_models, ~ {
    cat("\n\n------------------------------\n")
    cat(.y, "\n")
    cat("------------------------------\n")
    print(summary(.x))
  })
}, finally = sink())

final_model_fixed_effects <- purrr::imap_dfr(
  diagnostic_models,
  ~ broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE) %>%
    mutate(
      model = .y,
      p_formatted = if_else(p.value < 0.001, "<0.001", formatC(p.value, format = "f", digits = 3))
    ),
  .id = NULL
) %>%
  select(model, term, estimate, std.error, conf.low, conf.high, statistic, p.value, p_formatted)

write_csv(final_model_fixed_effects, file.path(stats_dir, "final_model_fixed_effects.csv"))
print(final_model_fixed_effects)
