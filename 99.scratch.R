### quick n dirty explorations ####

site_substrate <- substrate_transect %>%
  group_by(site_code, site, substrate) %>%
  summarise(
    mean_pct = mean(pct, na.rm = TRUE),
    sd_pct = sd(pct, na.rm = TRUE),
    n_transects = n(),
    .groups = "drop"
  )

ggplot(site_substrate, aes(x = site, y = mean_pct, fill = substrate)) +
  geom_col() +
  labs(x = NULL, y = "Mean cover (%)", fill = "Substrate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(substrate_transect, aes(x = site, y = pct, fill = site)) +
  geom_boxplot(outlier.alpha = 0.4) +
  facet_wrap(~ substrate, scales = "free_y") +
  labs(x = NULL, y = "Cover per transect (%)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


key_substrates <- c("HC", "TUR", "AB", "MAC", "OB")

substrate_transect %>%
  filter(substrate %in% key_substrates) %>%
  ggplot(aes(x = site, y = pct, fill = site)) +
  geom_boxplot(outlier.alpha = 0.4) +
  facet_wrap(~ substrate, scales = "free_y") +
  labs(x = NULL, y = "Cover per transect (%)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
# AL is bucking the trend so we may want to do a senstivity analysis, removing AL to see if the trends still hold / how much variation that's accounting for 
# but actually, it may not be COTs thats driving the changes, there are just weird reefs 


substrate_transect %>%
  filter(substrate == "HC") %>%
  ggplot(aes(x = site, y = pct)) +
  geom_boxplot(outlier.alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  labs(x = NULL, y = "Hard coral cover per transect (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


site_substrate %>%
  arrange(site_code, desc(mean_pct))


### compare first and last substrate - is there a substantial change ### 
substrate_first_last <- substrate_wide %>%
  mutate(
    transect_num = readr::parse_number(as.character(transect))
  ) %>%
  group_by(site_code, site) %>%
  filter(transect_num %in% c(min(transect_num), max(transect_num))) %>%
  mutate(period = if_else(transect_num == min(transect_num), "First", "Last")) %>%
  ungroup()
substrate_first_last_long <- substrate_first_last %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to = "substrate",
    values_to = "pct"
  ) %>%
  filter(substrate != "transect_num")

ggplot(substrate_first_last_long, aes(x = period, y = pct, fill = substrate)) +
  geom_col() +
  facet_wrap(~ site) +
  labs(x = NULL, y = "Substrate cover (%)", fill = "Substrate") +
  theme_minimal()

substrate_first_last_diff <- substrate_first_last_long %>%
  select(site_code, site, period, substrate, pct) %>%
  pivot_wider(names_from = period, values_from = pct) %>%
  mutate(diff_last_minus_first = Last - First) %>%
  arrange(site_code, desc(abs(diff_last_minus_first)))

substrate_first_last_diff

substrate_first_last_diff %>%
  ggplot(aes(x = substrate, y = diff_last_minus_first)) +
  geom_col() +
  facet_wrap(~ site) +
  labs(x = "Substrate", y = "Change in cover, last minus first (%)") +
  theme_minimal()


# check GR5 against latest 
gr_compare <- substrate_wide %>%
  filter(site_code == "GR") %>%
  mutate(transect_num = readr::parse_number(as.character(transect))) %>%
  filter(transect == "GR5" | transect_num == max(transect_num, na.rm = TRUE)) %>%
  mutate(period = if_else(transect == "GR5", "GR5", paste0("Most recent: ", transect))) %>%
  select(site_code, site, transect, period, AB, HC, MAC, OB, TUR, UKN, O, SP, SC)

gr_compare_long <- gr_compare %>%
  pivot_longer(
    cols = c(AB, HC, MAC, OB, TUR, UKN, O, SP, SC),
    names_to = "substrate",
    values_to = "pct"
  )

ggplot(gr_compare_long, aes(x = period, y = pct, fill = substrate)) +
  geom_col() +
  labs(
    x = NULL,
    y = "Substrate cover (%)",
    fill = "Substrate",
    title = "Green Rock substrate composition: GR5 vs most recent transect"
  ) +
  theme_minimal()


substrate_wide %>%
  filter(site_code == "GR") %>%
  mutate(transect_num = readr::parse_number(as.character(transect))) %>%
  arrange(transect_num) %>%
  select(site, transect, HC) %>%
  print(n = Inf)

substrate_wide %>%
  filter(site_code == "GR") %>%
  summarise(
    mean_HC = mean(HC, na.rm = TRUE),
    sd_HC = sd(HC, na.rm = TRUE),
    n_transects = n()
  )
