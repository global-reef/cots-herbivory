### 05_RESULTS.R ####
### Manuscript result tables, text values, and model-derived predictions

source("00_SETUP.R")

#### 01. LOAD DATA AND MODELS ####

cots_survey <- read_csv(file.path(data_clean_dir, "cots_survey_clean.csv"), show_col_types = FALSE) %>% mutate(date = as.Date(date))
fish_abund <- read_csv(file.path(data_clean_dir, "fish_abund_clean.csv"), show_col_types = FALSE) %>% mutate(date = as.Date(date))
fish_bites <- read_csv(file.path(data_clean_dir, "fish_bites_clean.csv"), show_col_types = FALSE) %>% mutate(date = as.Date(date))
substrate_wide <- read_csv(file.path(data_clean_dir, "substrate_transect_cover_wide_clean.csv"), show_col_types = FALSE)
substrate_point <- read_csv(file.path(data_clean_dir, "substrate_point_clean.csv"), show_col_types = FALSE)

cots_site_context <- read_csv(file.path(tables_dir, "cots_site_context.csv"), show_col_types = FALSE)
benthos_mod <- read_csv(file.path(tables_dir, "benthos_mod.csv"), show_col_types = FALSE)
fish_abund_mod <- read_csv(file.path(tables_dir, "fish_abund_mod.csv"), show_col_types = FALSE) %>% mutate(date = as.Date(date))
fish_bites_mod <- read_csv(file.path(tables_dir, "fish_bites_mod.csv"), show_col_types = FALSE) %>% mutate(date = as.Date(date))
feeding_pref_mod <- read_csv(file.path(tables_dir, "feeding_pref_mod.csv"), show_col_types = FALSE) %>% mutate(date = as.Date(date))

m_benthos_hc <- readRDS(file.path(fits_dir, "m_benthos_hc.rds"))
m_benthos_mac <- readRDS(file.path(fits_dir, "m_benthos_mac.rds"))
m_benthos_tur <- readRDS(file.path(fits_dir, "m_benthos_tur.rds"))
m_fish_abund <- readRDS(file.path(fits_dir, "m_fish_abund.rds"))
m_feeding_pref_cots_sub <- readRDS(file.path(fits_dir, "m_feeding_pref_cots_sub.rds"))
m_feeding_pref_cots_family <- readRDS(file.path(fits_dir, "m_feeding_pref_cots_family.rds"))
m_feeding_pref_cots_sub_family <- readRDS(file.path(fits_dir, "m_feeding_pref_cots_sub_family.rds"))
m_bite_any <- readRDS(file.path(fits_dir, "m_bite_any.rds"))
m_bite_positive <- readRDS(file.path(fits_dir, "m_bite_positive.rds"))

final_model_fixed_effects <- read_csv(file.path(stats_dir, "final_model_fixed_effects.csv"), show_col_types = FALSE)

#### 02. HELPERS ####

site_code_order <- c("GR", "RR", "TW", "AL", "SI", "TB")

pred_link_inv <- function(model, newdata, inv_fun) {
  pr <- predict(model, newdata = newdata, type = "link", se.fit = TRUE, re.form = NA)
  newdata %>%
    mutate(
      fit_link = pr$fit,
      se_link = pr$se.fit,
      lwr_link = fit_link - 1.96 * se_link,
      upr_link = fit_link + 1.96 * se_link,
      fit = inv_fun(fit_link),
      lwr = inv_fun(lwr_link),
      upr = inv_fun(upr_link)
    )
}

inv_logit <- function(x) plogis(x)

#### 03. SURVEY EFFORT VALUES ####

effort_summary <- tibble(
  dataset = c("CoTS abundance", "Fish abundance", "Fish feeding", "Benthic substrate"),
  n_units = c(
    nrow(cots_survey),
    n_distinct(fish_abund$survey_id),
    nrow(fish_bites),
    n_distinct(paste(substrate_wide$site_code, substrate_wide$transect, sep = "_"))
  ),
  total_count = c(
    sum(cots_survey$cots_count, na.rm = TRUE),
    sum(fish_abund$fish_count, na.rm = TRUE),
    sum(fish_bites$bites_total, na.rm = TRUE),
    sum(substrate_wide$total_points, na.rm = TRUE)
  ),
  start_date = c(
    min(cots_survey$date, na.rm = TRUE),
    min(fish_abund$date, na.rm = TRUE),
    min(fish_bites$date, na.rm = TRUE),
    NA
  ),
  end_date = c(
    max(cots_survey$date, na.rm = TRUE),
    max(fish_abund$date, na.rm = TRUE),
    max(fish_bites$date, na.rm = TRUE),
    NA
  )
)

write_csv(effort_summary, file.path(tables_dir, "results_effort_summary.csv"))
print(effort_summary)

substrate_effort_summary <- substrate_point %>%
  summarise(
    n_points = n(),
    n_transects = n_distinct(paste(site_code, transect, sep = "_")),
    n_photoquadrats = n_distinct(paste(site_code, transect, photo_id, sep = "_"))
  )

write_csv(substrate_effort_summary, file.path(tables_dir, "results_substrate_effort_summary.csv"))
print(substrate_effort_summary)

#### 04. TABLE 1 SITE-LEVEL CONTEXT ####

tab_cots <- cots_survey %>%
  group_by(site_code) %>%
  summarise(
    cots_mean = mean(cots_density_ha, na.rm = TRUE),
    cots_sd = sd(cots_density_ha, na.rm = TRUE),
    cots_max = max(cots_density_ha, na.rm = TRUE),
    .groups = "drop"
  )

tab_substrate <- substrate_wide %>%
  group_by(site_code) %>%
  summarise(
    HC_mean = mean(HC, na.rm = TRUE), HC_sd = sd(HC, na.rm = TRUE),
    TUR_mean = mean(TUR, na.rm = TRUE), TUR_sd = sd(TUR, na.rm = TRUE),
    MAC_mean = mean(MAC, na.rm = TRUE), MAC_sd = sd(MAC, na.rm = TRUE),
    .groups = "drop"
  )

table_1_context <- tab_cots %>%
  left_join(tab_substrate, by = "site_code") %>%
  mutate(site_code = factor(site_code, levels = site_code_order)) %>%
  arrange(site_code) %>%
  transmute(
    Site = site_code,
    `CoTS density` = paste0(round(cots_mean, 1), " ± ", round(cots_sd, 1)),
    `Max CoTS density` = round(cots_max, 1),
    `Hard coral (%)` = paste0(round(HC_mean, 1), " ± ", round(HC_sd, 1)),
    `Turf algae (%)` = paste0(round(TUR_mean, 1), " ± ", round(TUR_sd, 1)),
    `Macroalgae (%)` = paste0(round(MAC_mean, 1), " ± ", round(MAC_sd, 1))
  )

write_csv(table_1_context, file.path(tables_dir, "table_1_site_context_compact.csv"))
print(table_1_context)

#### 05. FULL VS 2025-ONWARD COTS CHECK ####

cots_compare_full_vs_2025 <- cots_survey %>%
  group_by(site_code) %>%
  summarise(
    cots_mean_ha_full = mean(cots_density_ha, na.rm = TRUE),
    cots_max_ha_full = max(cots_density_ha, na.rm = TRUE),
    n_cots_surveys_full = n(),
    .groups = "drop"
  ) %>%
  left_join(
    cots_survey %>%
      filter(date >= as.Date("2025-01-01")) %>%
      group_by(site_code) %>%
      summarise(
        cots_mean_ha_2025 = mean(cots_density_ha, na.rm = TRUE),
        cots_median_ha_2025 = median(cots_density_ha, na.rm = TRUE),
        cots_max_ha_2025 = max(cots_density_ha, na.rm = TRUE),
        n_cots_surveys_2025 = n(),
        .groups = "drop"
      ),
    by = "site_code"
  ) %>%
  mutate(
    rank_full = dense_rank(desc(cots_mean_ha_full)),
    rank_2025 = dense_rank(desc(cots_mean_ha_2025)),
    rank_changed = rank_full != rank_2025
  ) %>%
  arrange(rank_full)

write_csv(cots_compare_full_vs_2025, file.path(tables_dir, "cots_compare_full_vs_2025.csv"))
print(cots_compare_full_vs_2025)

#### 06. MAIN MODEL EFFECTS FOR TEXT ####

benthos_results_text <- final_model_fixed_effects %>%
  filter(model %in% c("m_benthos_hc", "m_benthos_mac", "m_benthos_tur"), term == "cots_mean_ha_sc") %>%
  mutate(
    substrate = case_when(
      model == "m_benthos_hc" ~ "Hard coral",
      model == "m_benthos_mac" ~ "Macroalgae",
      model == "m_benthos_tur" ~ "Turf algae"
    ),
    estimate = round(estimate, 2),
    conf.low = round(conf.low, 2),
    conf.high = round(conf.high, 2)
  ) %>%
  select(substrate, estimate, conf.low, conf.high, p_value = p_formatted)

write_csv(benthos_results_text, file.path(tables_dir, "results_text_benthos.csv"))
print(benthos_results_text)

fish_abund_results_text <- final_model_fixed_effects %>%
  filter(model == "m_fish_abund") %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ round(.x, 3))) %>%
  select(term, estimate, conf.low, conf.high, p_value = p_formatted)

write_csv(fish_abund_results_text, file.path(tables_dir, "results_text_fish_abundance.csv"))
print(fish_abund_results_text)

feeding_pref_results_text <- final_model_fixed_effects %>%
  filter(model %in% c("m_feeding_pref_cots_sub", "m_feeding_pref_cots_family", "m_feeding_pref_cots_sub_family")) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ round(.x, 3))) %>%
  select(model, term, estimate, conf.low, conf.high, p_value = p_formatted)

write_csv(feeding_pref_results_text, file.path(tables_dir, "results_text_feeding_preference.csv"))
print(feeding_pref_results_text)

bite_rate_results_text <- final_model_fixed_effects %>%
  filter(model %in% c("m_bite_any", "m_bite_positive")) %>%
  mutate(across(c(estimate, conf.low, conf.high), ~ round(.x, 3))) %>%
  select(model, term, estimate, conf.low, conf.high, p_value = p_formatted)

write_csv(bite_rate_results_text, file.path(tables_dir, "results_text_bite_rate_hurdle.csv"))
print(bite_rate_results_text)

#### 07. FISH ABUNDANCE PREDICTED VALUES FOR TEXT AND FIGURES ####

fish_pred_grid <- expand_grid(
  family = levels(factor(fish_abund_mod$family)),
  cots_mean_ha_sc = seq(min(fish_abund_mod$cots_mean_ha_sc, na.rm = TRUE), max(fish_abund_mod$cots_mean_ha_sc, na.rm = TRUE), length.out = 100)
) %>%
  mutate(
    area_ha = 1,
    site_code = fish_abund_mod$site_code[1],
    date = fish_abund_mod$date[1]
  )

cots_lookup_fish <- fish_abund_mod %>% distinct(cots_mean_ha, cots_mean_ha_sc) %>% arrange(cots_mean_ha_sc)

fish_pred <- pred_link_inv(m_fish_abund, fish_pred_grid, inv_fun = exp) %>%
  mutate(
    cots_mean_ha = approx(cots_lookup_fish$cots_mean_ha_sc, cots_lookup_fish$cots_mean_ha, xout = cots_mean_ha_sc, rule = 2)$y
  ) %>%
  left_join(functional_taxa %>% select(family, sci_name), by = "family")

fish_pred_summary <- fish_pred %>%
  group_by(family, sci_name) %>%
  summarise(
    pred_low_cots = fit[which.min(cots_mean_ha)],
    pred_high_cots = fit[which.max(cots_mean_ha)],
    pct_change = ((pred_high_cots - pred_low_cots) / pred_low_cots) * 100,
    .groups = "drop"
  ) %>%
  mutate(across(c(pred_low_cots, pred_high_cots, pct_change), ~ round(.x, 1)))

write_csv(fish_pred, file.path(tables_dir, "plotdata_fish_abundance_predictions.csv"))
write_csv(fish_pred_summary, file.path(tables_dir, "results_text_fish_abundance_prediction_summary.csv"))
print(fish_pred_summary)

#### 08. FEEDING PREFERENCE PREDICTED VALUES FOR TEXT AND FIGURES ####

feed_pref_cots_grid <- expand_grid(
  substrate = levels(factor(feeding_pref_mod$substrate)),
  family = levels(factor(feeding_pref_mod$family)),
  cots_mean_ha_sc = seq(min(feeding_pref_mod$cots_mean_ha_sc, na.rm = TRUE), max(feeding_pref_mod$cots_mean_ha_sc, na.rm = TRUE), length.out = 100)
) %>%
  mutate(
    avail_pct_sc = 0,
    obs_id = feeding_pref_mod$obs_id[1]
  )

cots_lookup_feed <- feeding_pref_mod %>% distinct(cots_mean_ha, cots_mean_ha_sc) %>% arrange(cots_mean_ha_sc)

feed_pref_cots_pred <- pred_link_inv(m_feeding_pref_cots_sub, feed_pref_cots_grid, inv_fun = inv_logit) %>%
  mutate(
    cots_mean_ha = approx(cots_lookup_feed$cots_mean_ha_sc, cots_lookup_feed$cots_mean_ha, xout = cots_mean_ha_sc, rule = 2)$y
  )

feed_pref_summary <- feed_pref_cots_pred %>%
  group_by(substrate, family) %>%
  summarise(
    pred_low_cots = fit[which.min(cots_mean_ha)] * 100,
    pred_high_cots = fit[which.max(cots_mean_ha)] * 100,
    pct_point_change = pred_high_cots - pred_low_cots,
    .groups = "drop"
  ) %>%
  mutate(across(c(pred_low_cots, pred_high_cots, pct_point_change), ~ round(.x, 1)))

write_csv(feed_pref_cots_pred, file.path(tables_dir, "plotdata_feeding_preference_predictions.csv"))
write_csv(feed_pref_summary, file.path(tables_dir, "results_text_feeding_preference_prediction_summary.csv"))
print(feed_pref_summary)

#### 09. SENSITIVITY TABLES FOR TEXT ####

sens_threshold <- read_csv(file.path(stats_dir, "model_overview_threshold_15.csv"), show_col_types = FALSE)
sens_al <- read_csv(file.path(stats_dir, "al_sensitivity_overview.csv"), show_col_types = FALSE)

aow_leuk_benthos_effects <- final_model_fixed_effects %>%
  filter(model %in% c("m_benthos_hc_no_al", "m_benthos_mac_no_al", "m_benthos_tur_no_al"), term == "cots_mean_ha_sc") %>%
  mutate(
    substrate = case_when(
      model == "m_benthos_hc_no_al" ~ "Hard coral",
      model == "m_benthos_mac_no_al" ~ "Macroalgae",
      model == "m_benthos_tur_no_al" ~ "Turf algae"
    ),
    estimate = round(estimate, 2),
    conf.low = round(conf.low, 2),
    conf.high = round(conf.high, 2)
  ) %>%
  select(substrate, estimate, conf.low, conf.high, p_value = p_formatted)

write_csv(sens_threshold, file.path(tables_dir, "results_sensitivity_threshold_15.csv"))
write_csv(sens_al, file.path(tables_dir, "results_sensitivity_aow_leuk_overview.csv"))
write_csv(aow_leuk_benthos_effects, file.path(tables_dir, "results_sensitivity_aow_leuk_benthos_effects.csv"))
print(aow_leuk_benthos_effects)
