##### 04_PLOTS.R #### 
### 04_PLOTS.R ####
### Publication figures for CoTS / fish feeding manuscript

source("00_SETUP.R")

#### 01. LOAD MODEL-READY DATA ####

cots_site_context <- read_csv(file.path(tables_dir, "cots_site_context.csv"), show_col_types = FALSE)
substrate_site_context <- read_csv(file.path(tables_dir, "substrate_site_context.csv"), show_col_types = FALSE)
benthos_mod <- read_csv(file.path(tables_dir, "benthos_mod.csv"), show_col_types = FALSE)
fish_abund_mod <- read_csv(file.path(tables_dir, "fish_abund_mod.csv"), show_col_types = FALSE)
fish_bites_mod <- read_csv(file.path(tables_dir, "fish_bites_mod.csv"), show_col_types = FALSE)
feeding_pref_mod <- read_csv(file.path(tables_dir, "feeding_pref_mod.csv"), show_col_types = FALSE)

#### 02. LOAD FINAL MODEL OBJECTS ####

m_benthos_hc  <- readRDS(file.path(fits_dir, "m_benthos_hc.rds"))
m_benthos_mac <- readRDS(file.path(fits_dir, "m_benthos_mac.rds"))
m_benthos_tur <- readRDS(file.path(fits_dir, "m_benthos_tur.rds"))
m_fish_abund  <- readRDS(file.path(fits_dir, "m_fish_abund.rds"))
m_feeding_pref <- readRDS(file.path(fits_dir, "m_feeding_pref.rds"))

#### 03. PLOT HELPERS ####

substrate_order <- c("HC", "SC", "AB", "TUR", "MAC", "SP", "OB", "UKN")
key_substrates <- c("HC", "MAC", "TUR")
site_code_order <- c("GR", "RR", "TW", "AL", "SI", "TB")

save_pub <- function(plot, filename, width = 180, height = 120) {
  ggsave(
    filename = file.path(plots_dir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "mm",
    dpi = 600,
    bg = "white"
  )
}

pred_link <- function(model, newdata) {
  pr <- predict(model, newdata = newdata, type = "link", se.fit = TRUE, re.form = NA)
  
  newdata %>%
    mutate(
      fit_link = pr$fit,
      se_link = pr$se.fit,
      lwr_link = fit_link - 1.96 * se_link,
      upr_link = fit_link + 1.96 * se_link,
      fit = fit_link,
      lwr = lwr_link,
      upr = upr_link
    )
}

pred_response <- function(model, newdata) {
  pr <- predict(model, newdata = newdata, type = "response", se.fit = TRUE, re.form = NA)
  
  newdata %>%
    mutate(
      fit = pr$fit,
      se = pr$se.fit,
      lwr = fit - 1.96 * se,
      upr = fit + 1.96 * se
    )
}

#### 04. SITE ORDER AND SITE LABELS ####
# Order sites from highest to lowest mean CoTS density.
# Use site codes on the x-axis for publication figures.

site_code_order <- c("GR", "RR", "TW", "AL", "SI", "TB")

site_lookup <- cots_site_context %>%
  mutate(
    site_code = factor(site_code, levels = site_code_order),
    site_type = factor(site_type, levels = c("Control", "High Density")),
    site_label = as.character(site_code)
  ) %>%
  arrange(site_code) %>%
  select(site_code, site_label, site_type)

site_label_order <- site_code_order

#### 05. FIGURE 1A: SUBSTRATE COMPOSITION BY SITE ####

site_substrate <- benthos_mod %>%
  select(site_code, transect, all_of(substrate_order)) %>%
  pivot_longer(
    cols = all_of(substrate_order),
    names_to = "substrate",
    values_to = "pct"
  ) %>%
  left_join(site_lookup, by = "site_code") %>%
  group_by(site_code, site_label, site_type, substrate) %>%
  summarise(
    mean_pct = mean(pct, na.rm = TRUE),
    sd_pct = sd(pct, na.rm = TRUE),
    n_transects = n_distinct(transect),
    .groups = "drop"
  ) %>%
  mutate(
    site_label = factor(site_label, levels = site_label_order),
    site_type = factor(site_type, levels = c("Control", "High Density")),
    substrate = factor(substrate, levels = substrate_order)
  )

p_fig1a_substrate <- ggplot(site_substrate, aes(x = site_label, y = mean_pct, fill = substrate)) +
  geom_col(width = 0.8, colour = "white", linewidth = 0.15) +
  scale_fill_manual(
    values = substrate_cols_pal[substrate_order],
    breaks = substrate_order,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Mean cover (%)",
    fill = "Substrate"
  ) +
  theme_clean +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p_fig1a_substrate
save_pub(p_fig1a_substrate, "Fig1A_substrate_composition_by_site.png", width = 190, height = 115)

#### 06. FIGURE 1B: COTS DENSITY BY SITE ####

cots_plot <- cots_site_context %>%
  left_join(site_lookup %>% select(site_code, site_label), by = "site_code") %>%
  mutate(
    site_label = factor(site_label, levels = site_label_order),
    site_type = factor(site_type, levels = c("Control", "High Density"))
  )

p_fig1b_cots <- ggplot(cots_plot, aes(x = site_label, y = cots_mean_ha, fill = site_type)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  geom_point(aes(y = cots_max_ha), shape = 21, size = 2.6, fill = "white", colour = "black") +
  geom_hline(yintercept = 15, linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = reef_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = expression("CoTS density (ha"^-1*")"),
    fill = "Site type"
  ) +
  theme_clean +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

p_fig1b_cots
save_pub(p_fig1b_cots, "Fig1B_cots_density_by_site.png", width = 190, height = 100)

#### 07. OPTIONAL COMBINED FIGURE 1 ####

if (requireNamespace("patchwork", quietly = TRUE)) {
  p_fig1 <- p_fig1a_substrate / p_fig1b_cots +
    patchwork::plot_annotation(tag_levels = "A")
  
  p_fig1
  save_pub(p_fig1, "Fig1_site_context_combined.png", width = 190, height = 210)
}

#### 08. FIGURE 2: COTS DENSITY VS KEY BENTHIC SUBSTRATES ####

benthos_long <- benthos_mod %>%
  select(site_code, transect, cots_mean_ha, cots_mean_ha_sc, HC, MAC, TUR) %>%
  left_join(site_lookup, by = "site_code") %>%
  pivot_longer(
    cols = c(HC, MAC, TUR),
    names_to = "substrate",
    values_to = "pct"
  ) %>%
  mutate(
    site_label = factor(site_label, levels = site_label_order),
    site_type = factor(site_type, levels = c("Control", "High Density")),
    substrate = factor(substrate, levels = c("HC", "MAC", "TUR")),
    al_label = if_else(site_code == "AL", "AL", NA_character_)
  )

make_benthos_pred <- function(model, substrate_name, data = benthos_mod) {
  newdat <- tibble(
    cots_mean_ha_sc = seq(
      min(data$cots_mean_ha_sc, na.rm = TRUE),
      max(data$cots_mean_ha_sc, na.rm = TRUE),
      length.out = 100
    )
  )
  
  cots_lookup <- data %>%
    distinct(cots_mean_ha, cots_mean_ha_sc) %>%
    arrange(cots_mean_ha_sc)
  
  pred_response(model, newdat) %>%
    mutate(
      cots_mean_ha = approx(
        x = cots_lookup$cots_mean_ha_sc,
        y = cots_lookup$cots_mean_ha,
        xout = cots_mean_ha_sc,
        rule = 2
      )$y,
      substrate = substrate_name
    )
}

benthos_pred <- bind_rows(
  make_benthos_pred(m_benthos_hc, "HC"),
  make_benthos_pred(m_benthos_mac, "MAC"),
  make_benthos_pred(m_benthos_tur, "TUR")
) %>%
  mutate(substrate = factor(substrate, levels = c("HC", "MAC", "TUR")))

p_fig2_benthos <- ggplot() +
  geom_ribbon(
    data = benthos_pred,
    aes(x = cots_mean_ha, ymin = lwr, ymax = upr),
    alpha = 0.18
  ) +
  geom_line(
    data = benthos_pred,
    aes(x = cots_mean_ha, y = fit),
    linewidth = 0.8
  ) +
  geom_point(
    data = benthos_long,
    aes(x = cots_mean_ha, y = pct, fill = site_type),
    shape = 21,
    size = 2.5,
    colour = "black",
    alpha = 0.85
  ) +
  # geom_text(
  #   data = benthos_long %>% filter(site_code == "AL"), # this is for AL labels 
  #   aes(x = cots_mean_ha, y = pct, label = al_label),
  #   nudge_y = 2,
  #   size = 3
  # ) +
  facet_wrap(~ substrate, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = reef_cols, drop = FALSE) +
  labs(
    x = expression("Mean CoTS density (ha"^-1*")"),
    y = "Transect cover (%)",
    fill = "Site type"
  ) +
  theme_clean +
  theme(legend.position = "right")

p_fig2_benthos
save_pub(p_fig2_benthos, "Fig2_cots_density_vs_key_substrates.png", width = 190, height = 90)


#### 09. FIGURE 3: FISH ABUNDANCE MODEL PREDICTIONS ####
# Model: fish_count ~ family * cots_mean_ha_sc + offset(log(area_ha)) +
#   (1 | site_code) + (1 | date)
# Prediction uses area_ha = 1, so fitted values are fish ha-1.
pred_link_inv <- function(model, newdata, inv_fun = exp) {
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

fish_pred_grid <- expand_grid(
  family = levels(factor(fish_abund_mod$family)),
  cots_mean_ha_sc = seq(
    min(fish_abund_mod$cots_mean_ha_sc, na.rm = TRUE),
    max(fish_abund_mod$cots_mean_ha_sc, na.rm = TRUE),
    length.out = 100
  )
) %>%
  mutate(
    area_ha = 1,
    site_code = fish_abund_mod$site_code[1],
    date = fish_abund_mod$date[1]
  )

cots_lookup_fish <- fish_abund_mod %>%
  distinct(cots_mean_ha, cots_mean_ha_sc) %>%
  arrange(cots_mean_ha_sc)

sci_name_labs <- c(
  "Scarus spp." = "italic('Scarus')~spp.",
  "Siganus spp." = "italic('Siganus')~spp.",
  "Chaetodon spp." = "italic('Chaetodon')~spp."
)

fish_pred <- pred_link_inv(m_fish_abund, fish_pred_grid, inv_fun = exp) %>%
  mutate(
    cots_mean_ha = approx(
      x = cots_lookup_fish$cots_mean_ha_sc,
      y = cots_lookup_fish$cots_mean_ha,
      xout = cots_mean_ha_sc,
      rule = 2
    )$y
  ) %>%
  left_join(functional_taxa %>% select(family, sci_name), by = "family") %>%
  mutate(
    sci_name = if_else(is.na(sci_name), as.character(family), as.character(sci_name)),
    sci_name = factor(sci_name, levels = c("Scarus spp.", "Siganus spp.", "Chaetodon spp.")),
    sci_name_lab = recode(as.character(sci_name), !!!sci_name_labs)
  )

fish_raw_site <- fish_abund_mod %>%
  group_by(site_code, site, site_type, family, sci_name, cots_mean_ha) %>%
  summarise(
    mean_density_ha = mean(fish_density_ha, na.rm = TRUE),
    se_density_ha = sd(fish_density_ha, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    site_type = factor(site_type, levels = c("Control", "High Density")),
    sci_name = factor(as.character(sci_name), levels = c("Scarus spp.", "Siganus spp.", "Chaetodon spp.")),
    sci_name_lab = recode(as.character(sci_name), !!!sci_name_labs)
  )

fish_raw_obs_trim <- fish_abund_mod %>%
  filter(fish_density_ha <= quantile(fish_density_ha, 0.95, na.rm = TRUE)) %>%
  mutate(
    site_type = factor(site_type, levels = c("Control", "High Density")),
    sci_name = factor(as.character(sci_name), levels = c("Scarus spp.", "Siganus spp.", "Chaetodon spp.")),
    sci_name_lab = recode(as.character(sci_name), !!!sci_name_labs)
  )
# trim top 5 percent of observations for plotting sake 

ymax_fish <- quantile(fish_abund_mod$fish_density_ha, 0.95, na.rm=TRUE)
p_fig3_fish_abund <- ggplot() +
  geom_ribbon(
    data = fish_pred,
    aes(x = cots_mean_ha, ymin = lwr, ymax = upr),
    alpha = 0.18
  ) +
  geom_line(
    data = fish_pred,
    aes(x = cots_mean_ha, y = fit),
    linewidth = 0.9
  ) +
  geom_point(
    data = fish_raw_obs_trim,
    aes(x = cots_mean_ha, y = fish_density_ha, fill = site_type),
    shape = 21,
    size = 1.8,
    colour = "black",
    alpha = 0.45,
    position = position_jitter(width = 0.4, height = 0)
  ) +
  facet_wrap(~ sci_name_lab, scales = "free_y", labeller = label_parsed) +
  scale_fill_manual(values = reef_cols, drop = FALSE) +
  labs(
    x = expression("Mean CoTS density (ha"^-1*")"),
    y = expression("Predicted fish density (ha"^-1*")"),
    fill = "Site type"
  ) +
  theme_clean +
  theme(legend.position = "right")

p_fig3_fish_abund
save_pub(p_fig3_fish_abund, "Fig3_fish_abundance_predictions.png", width = 180, height = 120)


#### 10. FIGURE 4: FEEDING PREFERENCE ACROSS COTS DENSITY ####

family_labs_parsed <- c(
  "Parrotfish" = "italic('Scarus')~spp.",
  "Rabbitfish" = "italic('Siganus')~spp.",
  "Butterflyfish" = "italic('Chaetodon')~spp."
)

feeding_obs_prop <- feeding_pref_mod %>%
  mutate(
    bite_prop = bites_sub / (bites_sub + bites_other),
    substrate = factor(substrate, levels = c("TUR", "MAC", "HC")),
    family = factor(family, levels = c("Parrotfish", "Rabbitfish", "Butterflyfish"))
  ) %>%
  filter(!is.na(bite_prop))

feed_pref_cots_grid <- expand_grid(
  substrate = c("TUR", "MAC", "HC"),
  family = c("Parrotfish", "Rabbitfish", "Butterflyfish"),
  cots_mean_ha_sc = seq(
    min(feeding_pref_mod$cots_mean_ha_sc, na.rm = TRUE),
    max(feeding_pref_mod$cots_mean_ha_sc, na.rm = TRUE),
    length.out = 100
  )
) %>%
  mutate(
    substrate = factor(substrate, levels = c("TUR", "MAC", "HC")),
    family = factor(family, levels = c("Parrotfish", "Rabbitfish", "Butterflyfish")),
    avail_pct_sc = 0,
    obs_id = feeding_pref_mod$obs_id[1]
  )

cots_lookup_feed <- feeding_pref_mod %>%
  distinct(cots_mean_ha, cots_mean_ha_sc) %>%
  arrange(cots_mean_ha_sc)

feed_pref_cots_pred <- pred_link_inv(
  m_feeding_pref_cots_sub_family,
  feed_pref_cots_grid,
  inv_fun = plogis
) %>%
  mutate(
    cots_mean_ha = approx(
      x = cots_lookup_feed$cots_mean_ha_sc,
      y = cots_lookup_feed$cots_mean_ha,
      xout = cots_mean_ha_sc,
      rule = 2
    )$y
  )

p_fig4_feed_cots_sub_family <- ggplot() +
  geom_jitter(
    data = feeding_obs_prop,
    aes(x = cots_mean_ha, y = bite_prop, fill = family),
    width = 0.45,
    height = 0,
    shape = 21,
    size = 1.2,
    alpha = 0.14,
    colour = "black"
  ) +
  geom_ribbon(
    data = feed_pref_cots_pred,
    aes(x = cots_mean_ha, ymin = lwr, ymax = upr, fill = family, group = family),
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(
    data = feed_pref_cots_pred,
    aes(x = cots_mean_ha, y = fit, colour = family, group = family),
    linewidth = 0.95
  ) +
  facet_wrap(
    ~ substrate,
    nrow = 1,
    labeller = as_labeller(c(
      "TUR" = "Turf algae",
      "MAC" = "Macroalgae",
      "HC" = "Hard coral"
    ))
  ) +
  scale_fill_manual(
    values = family_cols,
    labels = parse(text = family_labs_parsed),
    drop = FALSE
  ) +
  scale_colour_manual(
    values = family_cols,
    labels = parse(text = family_labs_parsed),
    drop = FALSE
  ) +
  labs(
    x = expression("Mean CoTS density (ha"^-1*")"),
    y = "Probability of substrate-specific bites",
    colour = NULL
  ) +
  theme_clean +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  ) +
  guides(fill = "none")

p_fig4_feed_cots_sub_family

save_pub(
  p_fig4_feed_cots_sub_family,
  "Fig4_feeding_preference_cots_by_substrate_family_interaction.png",
  width = 190,
  height = 105
)



#### 11. SUPPLEMENTARY FIGURE: RAW COTS SURVEY VARIATION ####

cots_raw_plot <- cots_site_context %>%
  select(site_code, site, site_type) %>%
  right_join(
    read_csv(file.path(data_clean_dir, "cots_survey_clean.csv"), show_col_types = FALSE) %>%
      select(site_code, survey_id, cots_density_ha),
    by = "site_code"
  ) %>%
  mutate(
    site_code = factor(site_code, levels = site_code_order),
    site_type = factor(site_type, levels = c("Control", "High Density"))
  )

p_supp_cots_raw <- ggplot(cots_raw_plot, aes(x = site_code, y = cots_density_ha, fill = site_type)) +
  geom_boxplot(outlier.alpha = 0.35, width = 0.65) +
  geom_jitter(width = 0.15, alpha = 0.45, size = 1.6, shape = 21, colour = "black") +
  geom_hline(yintercept = 15, linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = reef_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = expression("Survey-level CoTS density (ha"^-1*")"),
    fill = "Site type"
  ) +
  theme_clean

p_supp_cots_raw
save_pub(p_supp_cots_raw, "Supp_raw_cots_density_by_site.png", width = 180, height = 110)


#### 12. SUPPLEMENTARY FIGURE: RAW SUBSTRATE TRANSECT VARIATION ####

p_supp_substrate_raw <- benthos_long %>%
  ggplot(aes(x = site_label, y = pct, fill = site_type)) +
  geom_boxplot(outlier.alpha = 0.35, width = 0.65) +
  geom_jitter(width = 0.12, alpha = 0.4, size = 1.4, shape = 21, colour = "black") +
  facet_wrap(~ substrate, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = reef_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "Transect cover (%)",
    fill = "Site type"
  ) +
  theme_clean +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_supp_substrate_raw
save_pub(p_supp_substrate_raw, "Supp_raw_substrate_transect_variation.png", width = 190, height = 95)


#### 13. TABLE 1: SITE-LEVEL CONTEXT TABLE ####
site_code_order <- c("GR", "RR", "TW", "AL", "SI", "TB")

site_name_lookup <- c(
  "GR" = "Green Rock Wall",
  "RR" = "Red Rock Wall",
  "TW" = "Twins Wall",
  "AL" = "Aow Leuk",
  "SI" = "Shark Island",
  "TB" = "Tanote Bay"
)

table1_site_context <- cots_site_context %>%
  mutate(site_code = factor(site_code, levels = site_code_order)) %>%
  select(
    site_code, site_type,
    cots_mean_ha, cots_median_ha, cots_max_ha,
    cots_15_mean, n_cots_surveys
  ) %>%
  left_join(
    substrate_site_context %>%
      mutate(site_code = factor(site_code, levels = site_code_order)) %>%
      select(site_code, HC, MAC, TUR, AB, OB, algal_cover, benthic_state, n_substrate_transects),
    by = "site_code"
  ) %>%
  mutate(
    Site = site_name_lookup[as.character(site_code)],
    across(
      c(cots_mean_ha, cots_median_ha, cots_max_ha, HC, MAC, TUR, AB, OB, algal_cover, benthic_state),
      ~ round(.x, 2)
    )
  ) %>%
  arrange(site_code) %>%
  rename(
    `Site code` = site_code,
    `Site type` = site_type,
    `Mean CoTS density ha-1` = cots_mean_ha,
    `Median CoTS density ha-1` = cots_median_ha,
    `Max CoTS density ha-1` = cots_max_ha,
    `Mean 15 CoTS category` = cots_15_mean,
    `CoTS surveys` = n_cots_surveys,
    `Hard coral (%)` = HC,
    `Macroalgae (%)` = MAC,
    `Turf algae (%)` = TUR,
    `Abiotic (%)` = AB,
    `Other biotic (%)` = OB,
    `Algal cover (%)` = algal_cover,
    `Benthic state` = benthic_state,
    `Substrate transects` = n_substrate_transects
  ) %>%
  relocate(`Site code`, Site)

table1_site_context
write_csv(table1_site_context, file.path(tables_dir, "Table1_site_context.csv"))

table1_context_compact <- table1_site_context %>%
  transmute(
    Site = `Site code`,
    `CoTS density` = `Mean CoTS density ha-1`,
    `Max CoTS` = `Max CoTS density ha-1`,
    `Hard coral` = `Hard coral (%)`,
    `Turf algae` = `Turf algae (%)`,
    Macroalgae = `Macroalgae (%)`
  )

table1_context_compact
write_csv(table1_context_compact, file.path(tables_dir, "table1_context_compact.csv"))
#### 14. TABLE 2: MAIN MODEL FIXED EFFECTS ####

main_model_names <- c(
  "m_benthos_hc",
  "m_benthos_mac",
  "m_benthos_tur",
  "m_fish_abund",
  "m_feeding_pref"
)

table2_main_models <- read_csv(file.path(stats_dir, "final_model_fixed_effects.csv"), show_col_types = FALSE) %>%
  filter(model %in% main_model_names) %>%
  mutate(
    model_label = case_when(
      model == "m_benthos_hc" ~ "Benthos: hard coral",
      model == "m_benthos_mac" ~ "Benthos: macroalgae",
      model == "m_benthos_tur" ~ "Benthos: turf algae",
      model == "m_fish_abund" ~ "Fish abundance",
      model == "m_feeding_pref" ~ "Feeding preference"
    ),
    response = case_when(
      model == "m_benthos_hc" ~ "HC",
      model == "m_benthos_mac" ~ "MAC",
      model == "m_benthos_tur" ~ "TUR",
      model == "m_fish_abund" ~ "fish_count",
      model == "m_feeding_pref" ~ "cbind(bites_sub, bites_other)"
    ),
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    conf.low = round(conf.low, 3),
    conf.high = round(conf.high, 3),
    statistic = round(statistic, 3)
  ) %>%
  select(
    model_label, response, term, estimate, std.error,
    conf.low, conf.high, statistic, p_formatted
  ) %>%
  rename(
    Model = model_label,
    Response = response,
    Term = term,
    Estimate = estimate,
    SE = std.error,
    `CI low` = conf.low,
    `CI high` = conf.high,
    Statistic = statistic,
    `p-value` = p_formatted
  )

table2_main_models
write_csv(table2_main_models, file.path(tables_dir, "Table2_main_model_fixed_effects.csv"))


#### 15. SUPPLEMENTARY TABLE: AL SENSITIVITY OVERVIEW ####

if (file.exists(file.path(stats_dir, "al_sensitivity_overview.csv"))) {
  supp_al_sensitivity <- read_csv(file.path(stats_dir, "al_sensitivity_overview.csv"), show_col_types = FALSE)
  supp_al_sensitivity
  write_csv(supp_al_sensitivity, file.path(tables_dir, "Supp_AL_sensitivity_overview.csv"))
}


#### 16. SUPPLEMENTARY TABLE: 15 COTS THRESHOLD MODEL OVERVIEW ####

if (file.exists(file.path(stats_dir, "model_overview_15.csv"))) {
  supp_threshold_overview <- read_csv(file.path(stats_dir, "model_overview_15.csv"), show_col_types = FALSE)
  supp_threshold_overview
  write_csv(supp_threshold_overview, file.path(tables_dir, "Supp_15COTS_threshold_model_overview.csv"))
}


#### 17. SUPPLEMENTARY TABLE: DHARMA DIAGNOSTICS ####

if (file.exists(file.path(stats_dir, "dharma_diagnostic_tests.csv"))) {
  supp_dharma <- read_csv(file.path(stats_dir, "dharma_diagnostic_tests.csv"), show_col_types = FALSE) %>%
    mutate(
      across(
        c(uniformity_p, dispersion_p, outlier_p, zero_inflation_p),
        ~ round(.x, 3)
      )
    )
  
  supp_dharma
  write_csv(supp_dharma, file.path(tables_dir, "Supp_DHARMa_diagnostic_tests.csv"))
}


#### 18. SAVE PLOT DATA ####

write_csv(site_substrate, file.path(tables_dir, "plotdata_fig1_site_substrate.csv"))
write_csv(cots_plot, file.path(tables_dir, "plotdata_fig1_cots_site.csv"))
write_csv(benthos_long, file.path(tables_dir, "plotdata_fig2_benthos_raw.csv"))
write_csv(benthos_pred, file.path(tables_dir, "plotdata_fig2_benthos_predictions.csv"))
write_csv(fish_pred, file.path(tables_dir, "plotdata_fig3_fish_predictions.csv"))
write_csv(feed_pred, file.path(tables_dir, "plotdata_fig4_feeding_predictions.csv"))

cat("\n04_PLOTS complete.\n")
cat("Figures saved to:", plots_dir, "\n")
cat("Tables saved to:", tables_dir, "\n")
