#### Feeding preference interaction candidates ####

m_feeding_pref_main <- m_feeding_pref

m_feeding_pref_cots_sub <- glmmTMB(
  cbind(bites_sub, bites_other) ~ substrate * family + substrate * cots_mean_ha_sc +
    avail_pct_sc + (1 | obs_id),
  family = betabinomial(link = "logit"),
  data = feeding_pref_mod
)

m_feeding_pref_cots_family <- glmmTMB(
  cbind(bites_sub, bites_other) ~ substrate * family + family * cots_mean_ha_sc +
    avail_pct_sc + (1 | obs_id),
  family = betabinomial(link = "logit"),
  data = feeding_pref_mod
)

m_feeding_pref_cots_sub_family <- glmmTMB(
  cbind(bites_sub, bites_other) ~ substrate * family + substrate * cots_mean_ha_sc +
    family * cots_mean_ha_sc + avail_pct_sc + (1 | obs_id),
  family = betabinomial(link = "logit"),
  data = feeding_pref_mod
)

feeding_pref_interaction_models <- list(
  m_feeding_pref_main = m_feeding_pref_main,
  m_feeding_pref_cots_sub = m_feeding_pref_cots_sub,
  m_feeding_pref_cots_family = m_feeding_pref_cots_family,
  m_feeding_pref_cots_sub_family = m_feeding_pref_cots_sub_family
)

purrr::iwalk(feeding_pref_interaction_models, ~ {
  cat("\n==============================\n")
  cat(.y, "\n")
  cat("==============================\n")
  print(summary(.x))
})

feeding_pref_interaction_compare <- tibble(
  model = names(feeding_pref_interaction_models),
  n_obs = purrr::map_int(feeding_pref_interaction_models, nobs),
  AIC = purrr::map_dbl(feeding_pref_interaction_models, AIC),
  significant_fixed_effects = purrr::map_chr(feeding_pref_interaction_models, sig_effects)
) %>%
  arrange(AIC)

feeding_pref_interaction_compare
write_csv(
  feeding_pref_interaction_compare,
  file.path(stats_dir, "feeding_pref_interaction_model_compare.csv")
)

feeding_pref_interaction_dharma <- purrr::imap_dfr(
  feeding_pref_interaction_models,
  run_dharma_tests
)

feeding_pref_interaction_dharma
write_csv(
  feeding_pref_interaction_dharma,
  file.path(stats_dir, "feeding_pref_interaction_dharma.csv")
)



#### RELEVEL FEEDING PREF DATA BEFORE INTERACTION MODELS ####

feeding_pref_mod <- feeding_pref_mod %>%
  mutate(
    substrate = factor(substrate, levels = c("TUR", "MAC", "HC")),
    family = factor(family, levels = c("Parrotfish", "Rabbitfish", "Butterflyfish")),
    obs_id = factor(obs_id)
  )

m_feeding_pref_main_refit <- glmmTMB(
  cbind(bites_sub, bites_other) ~ substrate * family + avail_pct_sc + cots_mean_ha_sc +
    (1 | obs_id),
  family = betabinomial(link = "logit"),
  data = feeding_pref_mod
)

m_feeding_pref_cots_sub <- glmmTMB(
  cbind(bites_sub, bites_other) ~ substrate * family + substrate * cots_mean_ha_sc +
    avail_pct_sc + (1 | obs_id),
  family = betabinomial(link = "logit"),
  data = feeding_pref_mod
)

summary(m_feeding_pref_cots_sub)

anova(m_feeding_pref_main_refit, m_feeding_pref_cots_sub)

AIC(m_feeding_pref_main_refit, m_feeding_pref_cots_sub)

run_dharma_tests(m_feeding_pref_cots_sub, "m_feeding_pref_cots_sub")