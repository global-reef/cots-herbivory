### 02_EXPLORE.R ###

#### 00. helpers ####
substrate_order <- c("HC", "SC", "AB", "TUR", "MAC", "SP", "OB", "UKN")

substrate_cols <- intersect(substrate_order, names(substrate_wide))
key_substrates <- c("HC", "TUR", "MAC", "AB", "OB")

save_plot <- function(plot, filename, width = 8, height = 5) {
  ggsave(file.path(eda_dir, filename), plot, width = width, height = height, dpi = 300)
}

miss_summary <- function(x) {
  x %>%
    summarise(across(everything(), ~ sum(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
    mutate(pct_missing = n_missing / nrow(x) * 100) %>%
    filter(n_missing > 0) %>%
    arrange(desc(n_missing))
}

response_summary <- function(x, response) {
  x %>%
    summarise(
      n = n(),
      min = min({{ response }}, na.rm = TRUE),
      mean = mean({{ response }}, na.rm = TRUE),
      median = median({{ response }}, na.rm = TRUE),
      max = max({{ response }}, na.rm = TRUE),
      zeros = sum({{ response }} == 0, na.rm = TRUE),
      pct_zero = mean({{ response }} == 0, na.rm = TRUE) * 100
    )
}

iqr_outliers <- function(x, group_var, response) {
  x %>%
    group_by({{ group_var }}) %>%
    mutate(
      q1 = quantile({{ response }}, 0.25, na.rm = TRUE),
      q3 = quantile({{ response }}, 0.75, na.rm = TRUE),
      iqr = q3 - q1,
      outlier = {{ response }} < q1 - 1.5 * iqr | {{ response }} > q3 + 1.5 * iqr
    ) %>%
    filter(outlier) %>%
    ungroup()
}


#### 01. STRUCTURE CHECKS ####

str(fish_bites)
str(fish_abund)
str(cots_survey)
str(substrate_wide)


#### 02. MISSINGNESS ####

missingness <- bind_rows(
  miss_summary(fish_bites) %>% mutate(dataset = "fish_bites"),
  miss_summary(fish_abund) %>% mutate(dataset = "fish_abund"),
  miss_summary(cots_survey) %>% mutate(dataset = "cots_survey"),
  miss_summary(substrate_wide) %>% mutate(dataset = "substrate_wide")
) %>%
  select(dataset, variable, n_missing, pct_missing)

missingness %>% print(n = Inf)
write_csv(missingness, file.path(eda_dir, "missingness_summary.csv"))


#### 03. RANGE / IMPOSSIBLE VALUE CHECKS ####

range_flags_fish_bites <- fish_bites %>%
  filter(time_s <= 0 | bites_total < 0 | bites_min < 0 | date > Sys.Date()) %>%
  select(site, site_code, date, family, species, time_s, bites_total, bites_min)

range_flags_fish_abund <- fish_abund %>%
  filter(fish_count < 0 | fish_density_ha < 0 | area_ha <= 0 | transect_m <= 0 |
           visibility_m < 0 | date > Sys.Date()) %>%
  select(site, site_code, survey_id, date, researcher, family, fish_count, fish_density_ha, area_ha, visibility_m)

range_flags_cots <- cots_survey %>%
  filter(cots_count < 0 | cots_density_ha < 0 | area_ha <= 0 |
           duration_min <= 0 | vis_m < 0 | avg_depth_m < 0 | date > Sys.Date()) %>%
  select(site, site_code, survey_id, date, cots_count, cots_density_ha, duration_min, vis_m, avg_depth_m)

range_flags_fish_bites
range_flags_fish_abund
range_flags_cots


#### 04. SUBSTRATE COLLINEARITY ####

substrate_covars <- substrate_wide %>%
  select(any_of(substrate_cols))

substrate_cor <- cor(substrate_covars, use = "pairwise.complete.obs")

substrate_collinearity <- as.data.frame(as.table(substrate_cor)) %>%
  rename(var1 = Var1, var2 = Var2, r = Freq) %>%
  filter(var1 != var2, abs(r) > 0.5) %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(as.character(var1), as.character(var2))), collapse = " ~ ")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  arrange(desc(abs(r)))

substrate_collinearity
write_csv(substrate_collinearity, file.path(eda_dir, "substrate_collinearity_pairs.csv"))


#### 05. BALANCE / EFFORT CHECKS ####

balance_fish_bites <- fish_bites %>%
  group_by(site, site_code) %>%
  summarise(
    n_obs = n(),
    n_dates = n_distinct(date),
    first_date = min(date, na.rm = TRUE),
    last_date = max(date, na.rm = TRUE),
    .groups = "drop"
  )

balance_fish_abund <- fish_abund %>%
  group_by(site, site_code) %>%
  summarise(
    n_records = n(),
    n_surveys = n_distinct(survey_id),
    n_dates = n_distinct(date),
    n_researchers = n_distinct(researcher),
    total_fish = sum(fish_count, na.rm = TRUE),
    mean_density_ha = mean(fish_density_ha, na.rm = TRUE),
    max_density_ha = max(fish_density_ha, na.rm = TRUE),
    first_date = min(date, na.rm = TRUE),
    last_date = max(date, na.rm = TRUE),
    .groups = "drop"
  )

balance_cots <- cots_survey %>%
  group_by(site, site_code, site_type) %>%
  summarise(
    n_surveys = n(),
    n_dates = n_distinct(date),
    total_cots = sum(cots_count, na.rm = TRUE),
    mean_density_ha = mean(cots_density_ha, na.rm = TRUE),
    max_density_ha = max(cots_density_ha, na.rm = TRUE),
    pct_zero = mean(cots_count == 0, na.rm = TRUE) * 100,
    first_date = min(date, na.rm = TRUE),
    last_date = max(date, na.rm = TRUE),
    .groups = "drop"
  )

balance_substrate_effort <- substrate_transect %>%
  group_by(site, site_code) %>%
  summarise(
    n_transects = n_distinct(transect),
    mean_points = mean(total_points, na.rm = TRUE),
    min_points = min(total_points, na.rm = TRUE),
    max_points = max(total_points, na.rm = TRUE),
    .groups = "drop"
  )

balance_substrate_cover <- substrate_wide %>%
  group_by(site, site_code) %>%
  summarise(
    n_transects = n_distinct(transect),
    mean_HC = mean(HC, na.rm = TRUE),
    mean_TUR = mean(TUR, na.rm = TRUE),
    mean_MAC = mean(MAC, na.rm = TRUE),
    mean_AB = mean(AB, na.rm = TRUE),
    .groups = "drop"
  )

balance_fish_bites
balance_fish_abund
balance_cots
balance_substrate_effort
balance_substrate_cover

write_csv(balance_fish_bites, file.path(eda_dir, "balance_fish_bites.csv"))
write_csv(balance_fish_abund, file.path(eda_dir, "balance_fish_abund.csv"))
write_csv(balance_cots, file.path(eda_dir, "balance_cots.csv"))
write_csv(balance_substrate_effort, file.path(eda_dir, "balance_substrate_effort.csv"))
write_csv(balance_substrate_cover, file.path(eda_dir, "balance_substrate_cover.csv"))


#### 06. JOIN FEASIBILITY / TEMPORAL COVERAGE ####

coverage_site <- fish_bites %>%
  group_by(site_code) %>%
  summarise(
    bite_dates = n_distinct(date),
    first_bite = min(date, na.rm = TRUE),
    last_bite = max(date, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  full_join(
    fish_abund %>%
      group_by(site_code) %>%
      summarise(
        fish_surveys = n_distinct(survey_id),
        first_fish = min(date, na.rm = TRUE),
        last_fish = max(date, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "site_code"
  ) %>%
  full_join(
    cots_survey %>%
      group_by(site_code) %>%
      summarise(
        cots_surveys = n_distinct(survey_id),
        first_cots = min(date, na.rm = TRUE),
        last_cots = max(date, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "site_code"
  ) %>%
  full_join(
    substrate_wide %>%
      group_by(site_code) %>%
      summarise(
        substrate_transects = n_distinct(transect),
        mean_HC = mean(HC, na.rm = TRUE),
        mean_TUR = mean(TUR, na.rm = TRUE),
        mean_MAC = mean(MAC, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "site_code"
  )

coverage_site
write_csv(coverage_site, file.path(eda_dir, "coverage_site.csv"))


#### 07. RESPONSE DISTRIBUTIONS ####

response_distributions <- bind_rows(
  response_summary(fish_bites, bites_total) %>% mutate(dataset = "fish_bites", response = "bites_total"),
  response_summary(fish_abund, fish_count) %>% mutate(dataset = "fish_abund", response = "fish_count"),
  response_summary(cots_survey, cots_count) %>% mutate(dataset = "cots_survey", response = "cots_count")
) %>%
  select(dataset, response, everything())

response_distributions
write_csv(response_distributions, file.path(eda_dir, "response_distributions.csv"))

p_bites_hist <- ggplot(fish_bites, aes(x = bites_total)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ family, scales = "free_y") +
  labs(x = "Total bites", y = "Count") +
  theme_minimal()

p_fish_hist <- ggplot(fish_abund, aes(x = fish_count)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ family, scales = "free_y") +
  labs(x = "Fish count", y = "Count") +
  theme_minimal()

p_cots_hist <- ggplot(cots_survey, aes(x = cots_count)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ site_code, scales = "free_y") +
  labs(x = "CoTS count", y = "Count") +
  theme_minimal()

p_bites_hist
p_fish_hist
p_cots_hist

save_plot(p_bites_hist, "hist_fish_bites_total.png")
save_plot(p_fish_hist, "hist_fish_abund_count.png")
save_plot(p_cots_hist, "hist_cots_count.png")


#### 08. OUTLIER CHECKS ####

outliers_fish_bites <- iqr_outliers(fish_bites, family, bites_total) %>%
  select(site, site_code, date, family, species, time_s, bites_total, bites_min)

outliers_fish_abund <- iqr_outliers(fish_abund, family, fish_count) %>%
  select(site, site_code, survey_id, date, researcher, family, fish_count, fish_density_ha)

outliers_cots <- cots_survey %>%
  mutate(dummy_group = "cots") %>%
  iqr_outliers(dummy_group, cots_count) %>%
  select(site, site_code, site_type, survey_id, date, cots_count, cots_density_ha)

outliers_fish_bites
outliers_fish_abund
outliers_cots

write_csv(outliers_fish_bites, file.path(eda_dir, "outliers_fish_bites.csv"))
write_csv(outliers_fish_abund, file.path(eda_dir, "outliers_fish_abund.csv"))
write_csv(outliers_cots, file.path(eda_dir, "outliers_cots.csv"))


#### 09. MEAN-VARIANCE / OVERDISPERSION CHECKS ####

mean_variance_fish_bites <- fish_bites %>%
  group_by(family) %>%
  summarise(
    n = n(),
    mean_bites = mean(bites_total, na.rm = TRUE),
    var_bites = var(bites_total, na.rm = TRUE),
    var_mean_ratio = var_bites / mean_bites,
    .groups = "drop"
  )

mean_variance_fish_abund <- fish_abund %>%
  group_by(family) %>%
  summarise(
    n = n(),
    mean_count = mean(fish_count, na.rm = TRUE),
    var_count = var(fish_count, na.rm = TRUE),
    var_mean_ratio = var_count / mean_count,
    .groups = "drop"
  )

mean_variance_cots <- cots_survey %>%
  group_by(site_type) %>%
  summarise(
    n = n(),
    mean_count = mean(cots_count, na.rm = TRUE),
    var_count = var(cots_count, na.rm = TRUE),
    var_mean_ratio = var_count / mean_count,
    .groups = "drop"
  )

mean_variance_fish_bites
mean_variance_fish_abund
mean_variance_cots

write_csv(mean_variance_fish_bites, file.path(eda_dir, "mean_variance_fish_bites.csv"))
write_csv(mean_variance_fish_abund, file.path(eda_dir, "mean_variance_fish_abund.csv"))
write_csv(mean_variance_cots, file.path(eda_dir, "mean_variance_cots.csv"))


#### 10. RESPONSE BY GROUP ####

response_group_fish_bites <- fish_bites %>%
  group_by(site, site_code, family) %>%
  summarise(
    n_obs = n(),
    mean_bites_min = mean(bites_min, na.rm = TRUE),
    sd_bites_min = sd(bites_min, na.rm = TRUE),
    mean_total_bites = mean(bites_total, na.rm = TRUE),
    .groups = "drop"
  )

response_group_fish_abund <- fish_abund %>%
  group_by(site, site_code, family) %>%
  summarise(
    n_records = n(),
    total_fish = sum(fish_count, na.rm = TRUE),
    mean_density_ha = mean(fish_density_ha, na.rm = TRUE),
    sd_density_ha = sd(fish_density_ha, na.rm = TRUE),
    .groups = "drop"
  )

response_group_cots <- cots_survey %>%
  group_by(site, site_code, site_type) %>%
  summarise(
    n_surveys = n(),
    mean_cots_density_ha = mean(cots_density_ha, na.rm = TRUE),
    sd_cots_density_ha = sd(cots_density_ha, na.rm = TRUE),
    max_cots_density_ha = max(cots_density_ha, na.rm = TRUE),
    pct_zero = mean(cots_count == 0, na.rm = TRUE) * 100,
    .groups = "drop"
  )

response_group_fish_bites
response_group_fish_abund
response_group_cots

write_csv(response_group_fish_bites, file.path(eda_dir, "response_group_fish_bites.csv"))
write_csv(response_group_fish_abund, file.path(eda_dir, "response_group_fish_abund.csv"))
write_csv(response_group_cots, file.path(eda_dir, "response_group_cots.csv"))


#### 11. DEPENDENCE / REPEATED MEASURES ####

repeat_fish_abund <- fish_abund %>%
  distinct(site_code, survey_id, date, researcher) %>%
  count(site_code, date, name = "n_researcher_surveys") %>%
  arrange(desc(n_researcher_surveys))

repeat_cots <- cots_survey %>%
  count(site_code, date, name = "n_cots_surveys") %>%
  arrange(desc(n_cots_surveys))

repeat_fish_bites <- fish_bites %>%
  count(site_code, date, buddy, family, name = "n_obs") %>%
  arrange(site_code, date, buddy, family)

repeat_fish_abund
repeat_cots
repeat_fish_bites

write_csv(repeat_fish_abund, file.path(eda_dir, "repeat_fish_abund_site_date.csv"))
write_csv(repeat_cots, file.path(eda_dir, "repeat_cots_site_date.csv"))
write_csv(repeat_fish_bites, file.path(eda_dir, "repeat_fish_bites_site_date_buddy_family.csv"))


#### 12. RESPONSE-PREDICTOR SCREENING ####

site_substrate_mean <- substrate_wide %>%
  group_by(site_code) %>%
  summarise(
    mean_HC = mean(HC, na.rm = TRUE),
    mean_TUR = mean(TUR, na.rm = TRUE),
    mean_AB = mean(AB, na.rm = TRUE),
    .groups = "drop"
  )

fish_abund_sub <- fish_abund %>%
  left_join(site_substrate_mean, by = "site_code")

cots_sub <- cots_survey %>%
  left_join(site_substrate_mean, by = "site_code")

p_fish_hc <- ggplot(fish_abund_sub, aes(x = mean_HC, y = fish_density_ha)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(~ family, scales = "free_y") +
  labs(x = "Mean hard coral cover (%)", y = "Fish density (ha-1)") +
  theme_minimal()

p_cots_hc <- ggplot(cots_sub, aes(x = mean_HC, y = cots_density_ha)) +
  geom_point() +
  geom_smooth(method = "loess", se = TRUE) +
  facet_wrap(~ site_type) +
  labs(x = "Mean hard coral cover (%)", y = "CoTS density (ha-1)") +
  theme_minimal()

p_fish_hc
p_cots_hc

save_plot(p_fish_hc, "screen_fish_density_vs_mean_HC.png")
save_plot(p_cots_hc, "screen_cots_density_vs_mean_HC.png")


#### 13. SUBSTRATE COMPOSITION AROUND ISLAND ####

site_substrate <- substrate_transect %>%
  group_by(site_code, site, substrate) %>%
  summarise(
    mean_pct = mean(pct, na.rm = TRUE),
    sd_pct = sd(pct, na.rm = TRUE),
    n_transects = n_distinct(transect),
    .groups = "drop"
  )

p_substrate_site <- ggplot(site_substrate, aes(x = site, y = mean_pct, fill = substrate)) +
  geom_col() +
  scale_fill_manual(
    values = substrate_cols_pal[substrate_order],
    breaks = substrate_order, drop = FALSE) +
  labs(x = NULL, y = "Mean cover (%)", fill = "Substrate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

site_type_lookup <- cots_survey %>%
  distinct(site_code, site_type)

substrate_transect_plot <- substrate_transect %>%
  left_join(site_type_lookup, by = "site_code") %>%
  mutate(
    substrate = factor(substrate, levels = substrate_order),
    site_type = factor(site_type, levels = c("Control", "High Density"))
  )

p_key_substrate_box <- substrate_transect_plot %>%
  filter(substrate %in% key_substrates) %>%
  ggplot(aes(x = site, y = pct, fill = site_type)) +
  geom_boxplot(outlier.alpha = 0.4) +
  facet_wrap(~ substrate, scales = "free_y") +
  scale_fill_manual(values = reef_cols, drop = FALSE) +
  labs(x = NULL, y = "Cover per transect (%)", fill = "Site type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_hc_box <- substrate_transect_plot %>%
  filter(substrate == "HC") %>%
  ggplot(aes(x = site, y = pct, fill = site_type)) +
  geom_boxplot(outlier.alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  scale_fill_manual(values = reef_cols, drop = FALSE) +
  labs(x = NULL, y = "Hard coral cover per transect (%)", fill = "Site type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_substrate_site
p_key_substrate_box
p_hc_box

write_csv(site_substrate, file.path(eda_dir, "site_substrate_summary.csv"))
save_plot(p_substrate_site, "substrate_mean_cover_by_site.png")
save_plot(p_key_substrate_box, "substrate_key_categories_boxplot.png", width = 9, height = 6)
save_plot(p_hc_box, "substrate_hard_coral_boxplot.png")


#### 14. SUBSTRATE FIRST VS LAST TRANSECT SCREEN ####

substrate_first_last <- substrate_wide %>%
  mutate(transect_num = readr::parse_number(as.character(transect))) %>%
  group_by(site_code, site) %>%
  filter(transect_num %in% c(min(transect_num, na.rm = TRUE), max(transect_num, na.rm = TRUE))) %>%
  mutate(period = if_else(transect_num == min(transect_num, na.rm = TRUE), "First", "Last")) %>%
  ungroup()

substrate_first_last_long <- substrate_first_last %>%
  pivot_longer(
    cols = any_of(substrate_cols),
    names_to = "substrate",
    values_to = "pct"
  )

substrate_first_last_diff <- substrate_first_last_long %>%
  select(site_code, site, period, substrate, pct) %>%
  pivot_wider(names_from = period, values_from = pct) %>%
  mutate(diff_last_minus_first = Last - First) %>%
  arrange(site_code, desc(abs(diff_last_minus_first)))

p_first_last_comp <- ggplot(substrate_first_last_long, aes(x = period, y = pct, fill = substrate)) +
  geom_col() +
  scale_fill_manual(
    values = substrate_cols_pal[substrate_order],
    breaks = substrate_order,
    drop = FALSE
  ) +
  facet_wrap(~ site) +
  labs(x = NULL, y = "Substrate cover (%)", fill = "Substrate") +
  theme_minimal()

p_first_last_diff <- ggplot(substrate_first_last_diff, aes(x = substrate, y = diff_last_minus_first)) +
  geom_col() + 
  facet_wrap(~ site) +
  labs(x = "Substrate", y = "Change in cover, last minus first (%)") +
  theme_minimal() 

substrate_first_last_diff
p_first_last_comp
p_first_last_diff

write_csv(substrate_first_last_diff, file.path(eda_dir, "substrate_first_last_diff.csv"))
save_plot(p_first_last_comp, "substrate_first_last_composition.png", width = 9, height = 6)
save_plot(p_first_last_diff, "substrate_first_last_diff.png", width = 9, height = 6)


#### 15. GREEN ROCK SUBSTRATE CHECK ####

gr_compare <- substrate_wide %>%
  filter(site_code == "GR") %>%
  mutate(transect_num = readr::parse_number(as.character(transect))) %>%
  filter(transect == "GR5" | transect_num == max(transect_num, na.rm = TRUE)) %>%
  mutate(period = if_else(transect == "GR5", "GR5", paste0("Most recent: ", transect))) %>%
  select(site_code, site, transect, period, any_of(substrate_cols))

gr_compare_long <- gr_compare %>%
  pivot_longer(
    cols = any_of(substrate_cols),
    names_to = "substrate",
    values_to = "pct"
  )

gr_hc_by_transect <- substrate_wide %>%
  filter(site_code == "GR") %>%
  mutate(transect_num = readr::parse_number(as.character(transect))) %>%
  arrange(transect_num) %>%
  select(site, transect, HC)

gr_hc_summary <- substrate_wide %>%
  filter(site_code == "GR") %>%
  summarise(
    mean_HC = mean(HC, na.rm = TRUE),
    sd_HC = sd(HC, na.rm = TRUE),
    min_HC = min(HC, na.rm = TRUE),
    max_HC = max(HC, na.rm = TRUE),
    n_transects = n()
  )

p_gr_compare <- ggplot(gr_compare_long, aes(x = period, y = pct, fill = substrate)) +
  geom_col() +
  scale_fill_manual(
    values = substrate_cols_pal[substrate_order],
    breaks = substrate_order, drop = FALSE) +
  labs(
    x = NULL,
    y = "Substrate cover (%)",
    fill = "Substrate",
    title = "Green Rock substrate composition: GR5 vs most recent transect"
  ) +
  theme_minimal()

gr_hc_by_transect
gr_hc_summary
p_gr_compare

write_csv(gr_compare_long, file.path(eda_dir, "green_rock_gr5_vs_latest_long.csv"))
write_csv(gr_hc_by_transect, file.path(eda_dir, "green_rock_hc_by_transect.csv"))
write_csv(gr_hc_summary, file.path(eda_dir, "green_rock_hc_summary.csv"))
save_plot(p_gr_compare, "green_rock_gr5_vs_latest.png")

list.files(eda_dir)
cat("\nEDA complete. Outputs saved to:", eda_dir, "\n")


