#### RESULTS: effort summary + compact Table 1 ####

site_code_order <- c("GR", "RR", "TW", "AL", "SI", "TB")

#### 01. Overall survey effort text values ####

effort_summary <- tibble(
  dataset = c("CoTS abundance", "Fish abundance", "Fish feeding", "Benthic substrate"),
  n_units = c(
    nrow(cots_survey),
    n_distinct(fish_abund$survey_id),
    nrow(fish_bites),
    n_distinct(substrate_wide$transect)
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

effort_summary
write_csv(effort_summary, file.path(tables_dir, "results_effort_summary.csv"))


substrate_effort_summary <- substrate_point %>%
  summarise(
    n_points = n(),
    n_transects = n_distinct(paste(site_code, transect, sep = "_")),
    n_photoquadrats = n_distinct(paste(site_code, transect, photo_id, sep = "_"))
  )
substrate_effort_summary
#### 02. Compact Table 1: site-level ecological context ####

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
    HC_mean = mean(HC, na.rm = TRUE),
    HC_sd = sd(HC, na.rm = TRUE),
    TUR_mean = mean(TUR, na.rm = TRUE),
    TUR_sd = sd(TUR, na.rm = TRUE),
    MAC_mean = mean(MAC, na.rm = TRUE),
    MAC_sd = sd(MAC, na.rm = TRUE),
    .groups = "drop"
  )

table_1_context <- tab_cots %>%
  left_join(tab_substrate, by = "site_code") %>%
  mutate(
    site_code = factor(site_code, levels = site_code_order)
  ) %>%
  arrange(site_code) %>%
  transmute(
    Site = site_code,
    `CoTS density` = paste0(round(cots_mean, 1), " ± ", round(cots_sd, 1)),
    `Max CoTS` = round(cots_max, 1),
    `Hard coral` = paste0(round(HC_mean, 1), " ± ", round(HC_sd, 1)),
    `Turf algae` = paste0(round(TUR_mean, 1), " ± ", round(TUR_sd, 1)),
    `Macroalgae` = paste0(round(MAC_mean, 1), " ± ", round(MAC_sd, 1))
  )

table_1_context

write_csv(table_1_context, file.path(tables_dir, "table_1_site_context_compact.csv"))



#### 
cots_survey_2025 <- cots_survey %>%
  filter(date >= as.Date("2025-01-01"))

cots_site_2025 <- cots_survey_2025 %>%
  group_by(site_code) %>%
  summarise(
    cots_mean_ha_2025 = mean(cots_density_ha, na.rm = TRUE),
    cots_median_ha_2025 = median(cots_density_ha, na.rm = TRUE),
    cots_max_ha_2025 = max(cots_density_ha, na.rm = TRUE),
    n_cots_surveys_2025 = n(),
    .groups = "drop"
  )
cots_compare_full_vs_2025 <- cots_survey %>%
  group_by(site_code) %>%
  summarise(
    cots_mean_ha_full = mean(cots_density_ha, na.rm = TRUE),
    cots_max_ha_full = max(cots_density_ha, na.rm = TRUE),
    n_cots_surveys_full = n(),
    .groups = "drop"
  ) %>%
  left_join(cots_site_2025, by = "site_code") %>%
  arrange(desc(cots_mean_ha_full))

cots_compare_full_vs_2025

#### BENTHIC MODEL RESULTS FOR TEXT ####

benthos_results_text <- read_csv(file.path(stats_dir, "final_model_fixed_effects.csv"),
                                 show_col_types = FALSE) %>%
  filter(model %in% c("m_benthos_hc", "m_benthos_mac", "m_benthos_tur"),
         term == "cots_mean_ha_sc") %>%
  mutate(
    substrate = case_when(
      model == "m_benthos_hc" ~ "Hard coral",
      model == "m_benthos_mac" ~ "Macroalgae",
      model == "m_benthos_tur" ~ "Turf algae"
    ),
    estimate = round(estimate, 2),
    conf.low = round(conf.low, 2),
    conf.high = round(conf.high, 2),
    p_value = p_formatted
  ) %>%
  select(substrate, estimate, conf.low, conf.high, p_value)

benthos_results_text


#### FISH ABUND MODEL RESULTS FOR TEXT ####
fish_abund_results_text <- read_csv(file.path(stats_dir, "final_model_fixed_effects.csv"),
                                    show_col_types = FALSE) %>%
  filter(model == "m_fish_abund") %>%
  mutate(
    estimate = round(estimate, 3),
    conf.low = round(conf.low, 3),
    conf.high = round(conf.high, 3)
  ) %>%
  select(term, estimate, conf.low, conf.high, p_formatted)

fish_abund_results_text

fish_pred_summary <- fish_pred %>%
  group_by(sci_name) %>%
  summarise(
    pred_low_cots = fit[which.min(cots_mean_ha)],
    pred_high_cots = fit[which.max(cots_mean_ha)],
    pct_change = ((pred_high_cots - pred_low_cots) / pred_low_cots) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    across(c(pred_low_cots, pred_high_cots, pct_change), ~ round(.x, 1))
  )

fish_pred_summary

#### FEEDING MODEL RESULTS FOR TEXT ####
feeding_pref_results_text <- read_csv(file.path(stats_dir, "final_model_fixed_effects.csv"),
                                      show_col_types = FALSE) %>%
  filter(model %in% c("m_feeding_pref", "m_feeding_pref_cots_sub")) %>%
  mutate(
    estimate = round(estimate, 3),
    conf.low = round(conf.low, 3),
    conf.high = round(conf.high, 3)
  ) %>%
  select(model, term, estimate, conf.low, conf.high, p_formatted)

feeding_pref_results_text

feed_pref_summary <- feed_pref_cots_pred %>%
  group_by(substrate, family) %>%
  summarise(
    pred_low_cots = fit[which.min(cots_mean_ha)],
    pred_high_cots = fit[which.max(cots_mean_ha)],
    change = pred_high_cots - pred_low_cots,
    .groups = "drop"
  ) %>%
  mutate(
    across(c(pred_low_cots, pred_high_cots, change), ~ round(.x, 3))
  )

feed_pref_summary

#### FEEDING INTERACTION MODEL RESULTS ####

feeding_interaction_results <- broom.mixed::tidy(
  m_feeding_pref_cots_sub_family,
  effects = "fixed",
  conf.int = TRUE
) %>%
  mutate(
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    conf.low = round(conf.low, 3),
    conf.high = round(conf.high, 3),
    statistic = round(statistic, 3),
    p_formatted = case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  )

feeding_interaction_results

write_csv(
  feeding_interaction_results,
  file.path(stats_dir, "feeding_pref_cots_sub_family_fixed_effects.csv")
)
feeding_cots_terms <- feeding_interaction_results %>%
  filter(str_detect(term, "cots_mean_ha_sc")) %>%
  select(term, estimate, conf.low, conf.high, p_formatted)

feeding_cots_terms

#### 02. 15 CoTS threshold sensitivity ####

sens_threshold <- read_csv(
  file.path(stats_dir, "model_overview_15.csv"),
  show_col_types = FALSE
)

sens_threshold

#### 03. Aow Leuk exclusion sensitivity ####

sens_al <- read_csv(
  file.path(stats_dir, "al_sensitivity_overview.csv"),
  show_col_types = FALSE
)

sens_al

names(sens_threshold)
names(sens_al)

#### AOW LEUK SENSITIVITY: BENTHOS ESTIMATES ####

al_benthos_effects <- bind_rows(
  broom.mixed::tidy(m_benthos_hc, effects = "fixed", conf.int = TRUE) %>%
    mutate(model = "HC full"),
  broom.mixed::tidy(al_models$m_benthos_hc_no_al, effects = "fixed", conf.int = TRUE) %>%
    mutate(model = "HC no AL"),
  
  broom.mixed::tidy(m_benthos_mac, effects = "fixed", conf.int = TRUE) %>%
    mutate(model = "MAC full"),
  broom.mixed::tidy(al_models$m_benthos_mac_no_al, effects = "fixed", conf.int = TRUE) %>%
    mutate(model = "MAC no AL"),
  
  broom.mixed::tidy(m_benthos_tur, effects = "fixed", conf.int = TRUE) %>%
    mutate(model = "TUR full"),
  broom.mixed::tidy(al_models$m_benthos_tur_no_al, effects = "fixed", conf.int = TRUE) %>%
    mutate(model = "TUR no AL")
) %>%
  filter(term == "cots_mean_ha_sc") %>%
  mutate(
    estimate = round(estimate, 2),
    conf.low = round(conf.low, 2),
    conf.high = round(conf.high, 2),
    p_formatted = case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  ) %>%
  select(model, estimate, conf.low, conf.high, p_formatted)

al_benthos_effects