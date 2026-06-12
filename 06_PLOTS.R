### 06_PLOTS.R ####
### Publication figures for CoTS / fish feeding manuscript

source("00_SETUP.R")

#### 01. LOAD MODEL-READY DATA AND PREDICTIONS ####

cots_site_context <- read_csv(file.path(tables_dir, "cots_site_context.csv"), show_col_types = FALSE)
benthos_mod <- read_csv(file.path(tables_dir, "benthos_mod.csv"), show_col_types = FALSE)
fish_abund_mod <- read_csv(file.path(tables_dir, "fish_abund_mod.csv"), show_col_types = FALSE)
feeding_pref_mod <- read_csv(file.path(tables_dir, "feeding_pref_mod.csv"), show_col_types = FALSE)

fish_pred <- read_csv(file.path(tables_dir, "plotdata_fish_abundance_predictions.csv"), show_col_types = FALSE)
feed_pref_cots_pred <- read_csv(file.path(tables_dir, "plotdata_feeding_preference_predictions.csv"), show_col_types = FALSE)

m_benthos_hc <- readRDS(file.path(fits_dir, "m_benthos_hc.rds"))
m_benthos_mac <- readRDS(file.path(fits_dir, "m_benthos_mac.rds"))
m_benthos_tur <- readRDS(file.path(fits_dir, "m_benthos_tur.rds"))

#### 02. PLOT HELPERS ####

substrate_order <- c("HC", "SC", "AB", "TUR", "MAC", "SP", "OB", "UKN")
site_code_order <- c("GR", "RR", "TW", "AL", "SI", "TB")

save_pub <- function(plot, filename, width = 180, height = 120) {
  ggsave(file.path(plots_dir, filename), plot, width = width, height = height, units = "mm", dpi = 600, bg = "white")
}

pred_response <- function(model, newdata) {
  pr <- predict(model, newdata = newdata, type = "response", se.fit = TRUE, re.form = NA)
  newdata %>% mutate(fit = pr$fit, se = pr$se.fit, lwr = fit - 1.96 * se, upr = fit + 1.96 * se)
}

site_lookup <- cots_site_context %>%
  mutate(
    site_code = factor(site_code, levels = site_code_order),
    site_type = factor(site_type, levels = c("Control", "High Density")),
    site_label = as.character(site_code)
  ) %>%
  arrange(site_code) %>%
  select(site_code, site_label, site_type)

#### 03. FIGURE 11A: SUBSTRATE COMPOSITION BY SITE ####

site_substrate <- benthos_mod %>%
  select(site_code, transect, all_of(substrate_order)) %>%
  pivot_longer(cols = all_of(substrate_order), names_to = "substrate", values_to = "pct") %>%
  left_join(site_lookup, by = "site_code") %>%
  group_by(site_code, site_label, site_type, substrate) %>%
  summarise(mean_pct = mean(pct, na.rm = TRUE), sd_pct = sd(pct, na.rm = TRUE), n_transects = n_distinct(transect), .groups = "drop") %>%
  mutate(
    site_label = factor(site_label, levels = site_code_order),
    site_type = factor(site_type, levels = c("Control", "High Density")),
    substrate = factor(substrate, levels = substrate_order)
  )

p_fig11a_substrate <- ggplot(site_substrate, aes(x = site_label, y = mean_pct, fill = substrate)) +
  geom_col(width = 0.8, colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = substrate_cols_pal[substrate_order], breaks = substrate_order, drop = FALSE) +
  labs(x = NULL, y = "Mean cover (%)", fill = "Substrate") +
  theme_clean +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

save_pub(p_fig11a_substrate, "Fig11A_substrate_composition_by_site.png", width = 190, height = 115)

#### 04. FIGURE 11B: COTS DENSITY BY SITE ####

cots_plot <- cots_site_context %>%
  left_join(site_lookup %>% select(site_code, site_label), by = "site_code") %>%
  mutate(
    site_label = factor(site_label, levels = site_code_order),
    site_type = factor(site_type, levels = c("Control", "High Density"))
  )

p_fig11b_cots <- ggplot(cots_plot, aes(x = site_label, y = cots_mean_ha, fill = site_type)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  geom_point(aes(y = cots_max_ha), shape = 21, size = 2.6, fill = "white", colour = "black") +
  geom_hline(yintercept = 15, linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = reef_cols, drop = FALSE) +
  labs(x = NULL, y = expression("CoTS density (ha"^-1*")"), fill = "Site type") +
  theme_clean +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

save_pub(p_fig11b_cots, "Fig11B_cots_density_by_site.png", width = 190, height = 100)

if (requireNamespace("patchwork", quietly = TRUE)) {
  p_fig11 <- p_fig11a_substrate / p_fig11b_cots + patchwork::plot_annotation(tag_levels = "A")
  save_pub(p_fig11, "Fig11_site_context_combined.png", width = 190, height = 210)
}

#### 05. FIGURE 12: COTS DENSITY VS KEY BENTHIC SUBSTRATES ####

benthos_long <- benthos_mod %>%
  select(site_code, transect, cots_mean_ha, cots_mean_ha_sc, HC, MAC, TUR) %>%
  left_join(site_lookup, by = "site_code") %>%
  pivot_longer(cols = c(HC, MAC, TUR), names_to = "substrate", values_to = "pct") %>%
  mutate(
    site_label = factor(site_label, levels = site_code_order),
    site_type = factor(site_type, levels = c("Control", "High Density")),
    substrate = factor(substrate, levels = c("HC", "MAC", "TUR"))
  )

make_benthos_pred <- function(model, substrate_name, data = benthos_mod) {
  newdat <- tibble(cots_mean_ha_sc = seq(min(data$cots_mean_ha_sc, na.rm = TRUE), max(data$cots_mean_ha_sc, na.rm = TRUE), length.out = 100))
  cots_lookup <- data %>% distinct(cots_mean_ha, cots_mean_ha_sc) %>% arrange(cots_mean_ha_sc)
  pred_response(model, newdat) %>%
    mutate(
      cots_mean_ha = approx(cots_lookup$cots_mean_ha_sc, cots_lookup$cots_mean_ha, xout = cots_mean_ha_sc, rule = 2)$y,
      substrate = substrate_name
    )
}

benthos_pred <- bind_rows(
  make_benthos_pred(m_benthos_hc, "HC"),
  make_benthos_pred(m_benthos_mac, "MAC"),
  make_benthos_pred(m_benthos_tur, "TUR")
) %>% mutate(substrate = factor(substrate, levels = c("HC", "MAC", "TUR")))

p_fig12_benthos <- ggplot() +
  geom_ribbon(data = benthos_pred, aes(x = cots_mean_ha, ymin = lwr, ymax = upr), alpha = 0.18) +
  geom_line(data = benthos_pred, aes(x = cots_mean_ha, y = fit), linewidth = 0.8) +
  geom_point(data = benthos_long, aes(x = cots_mean_ha, y = pct, fill = site_type), shape = 21, size = 2.5, colour = "black", alpha = 0.85) +
  facet_wrap(~ substrate, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = reef_cols, drop = FALSE) +
  labs(x = expression("Mean CoTS density (ha"^-1*")"), y = "Transect cover (%)", fill = "Site type") +
  theme_clean +
  theme(legend.position = "right")

save_pub(p_fig12_benthos, "Fig12_cots_density_vs_key_substrates.png", width = 190, height = 90)

#### 06. FIGURE 13: FISH ABUNDANCE MODEL PREDICTIONS ####

family_labs_parsed <- c(
  "Parrotfish" = "italic('Scarus')~spp.",
  "Rabbitfish" = "italic('Siganus')~spp.",
  "Butterflyfish" = "italic('Chaetodon')~spp."
)

fish_raw_obs_trim <- fish_abund_mod %>%
  filter(fish_density_ha <= quantile(fish_density_ha, 0.95, na.rm = TRUE)) %>%
  mutate(site_type = factor(site_type, levels = c("Control", "High Density")))

p_fig13_fish_abund <- ggplot() +
  geom_ribbon(data = fish_pred, aes(x = cots_mean_ha, ymin = lwr, ymax = upr), alpha = 0.18) +
  geom_line(data = fish_pred, aes(x = cots_mean_ha, y = fit), linewidth = 0.9) +
  geom_point(data = fish_raw_obs_trim, aes(x = cots_mean_ha, y = fish_density_ha, fill = site_type), shape = 21, size = 1.8, colour = "black", alpha = 0.45, position = position_jitter(width = 0.45, height = 0)) +
  facet_wrap(~ family, nrow = 1, labeller = as_labeller(family_labs_parsed, label_parsed)) +
  scale_fill_manual(values = reef_cols, drop = FALSE) +
  labs(x = expression("Mean CoTS density (ha"^-1*")"), y = expression("Fish density (individuals ha"^-1*")"), fill = "Site type") +
  theme_clean +
  theme(legend.position = "right")

save_pub(p_fig13_fish_abund, "Fig13_fish_abundance_predictions.png", width = 190, height = 90)

#### 07. FIGURE 14: FEEDING PREFERENCE ACROSS COTS DENSITY ####

feeding_obs_prop <- feeding_pref_mod %>%
  mutate(
    bite_prop = bites_sub / (bites_sub + bites_other),
    substrate = factor(substrate, levels = c("TUR", "MAC", "HC")),
    family = factor(family, levels = c("Parrotfish", "Rabbitfish", "Butterflyfish"))
  ) %>%
  filter(!is.na(bite_prop))

feed_pref_cots_pred <- feed_pref_cots_pred %>%
  mutate(
    substrate = factor(substrate, levels = c("TUR", "MAC", "HC")),
    family = factor(family, levels = c("Parrotfish", "Rabbitfish", "Butterflyfish"))
  )

p_fig14_feed_cots <- ggplot() +
  geom_jitter(data = feeding_obs_prop, aes(x = cots_mean_ha, y = bite_prop, fill = family), width = 0.45, height = 0, shape = 21, size = 1.2, alpha = 0.14, colour = "black") +
  geom_ribbon(data = feed_pref_cots_pred, aes(x = cots_mean_ha, ymin = lwr, ymax = upr, fill = family, group = family), alpha = 0.12, colour = NA) +
  geom_line(data = feed_pref_cots_pred, aes(x = cots_mean_ha, y = fit, colour = family, group = family), linewidth = 0.95) +
  facet_wrap(~ substrate, nrow = 1, labeller = as_labeller(c("TUR" = "Turf algae", "MAC" = "Macroalgae", "HC" = "Hard coral"))) +
  scale_fill_manual(values = family_cols, labels = parse(text = family_labs_parsed), drop = FALSE) +
  scale_colour_manual(values = family_cols, labels = parse(text = family_labs_parsed), drop = FALSE) +
  labs(x = expression("Mean CoTS density (ha"^-1*")"), y = "Bite probability", fill = "Family", colour = "Family") +
  theme_clean +
  theme(legend.position = "right")

save_pub(p_fig14_feed_cots, "Fig14_feeding_preference_cots_substrate_family.png", width = 190, height = 90)

#### 08. SAVE PLOT DATA USED DIRECTLY IN THIS SCRIPT ####

write_csv(site_substrate, file.path(tables_dir, "plotdata_fig11a_site_substrate.csv"))
write_csv(cots_plot, file.path(tables_dir, "plotdata_fig11b_cots.csv"))
write_csv(benthos_long, file.path(tables_dir, "plotdata_fig12_benthos_observed.csv"))
write_csv(benthos_pred, file.path(tables_dir, "plotdata_fig12_benthos_predicted.csv"))
write_csv(fish_raw_obs_trim, file.path(tables_dir, "plotdata_fig13_fish_observed_trimmed.csv"))
write_csv(feeding_obs_prop, file.path(tables_dir, "plotdata_fig14_feeding_observed.csv"))
