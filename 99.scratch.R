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
