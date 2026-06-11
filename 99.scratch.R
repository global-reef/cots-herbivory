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



#### QUICK CHECK: feeding preference CoTS-by-family models ####

library(glmmTMB)
library(DHARMa)
library(broom.mixed)
library(performance)

#### 01. Main model for comparison ####

m_feeding_pref_main <- m_feeding_pref

#### 02. Does CoTS effect differ among families? ####

m_feeding_pref_cots_family <- glmmTMB(
  cbind(bites_sub, bites_other) ~ substrate * family +
    family * cots_mean_ha_sc +
    avail_pct_sc +
    (1 | obs_id),
  family = betabinomial(link = "logit"),
  data = feeding_pref_mod
)

#### 03. Does CoTS effect differ by substrate and family? ####

m_feeding_pref_cots_sub_family <- glmmTMB(
  cbind(bites_sub, bites_other) ~ substrate * family * cots_mean_ha_sc +
    avail_pct_sc +
    (1 | obs_id),
  family = betabinomial(link = "logit"),
  data = feeding_pref_mod
)

#### 04. Model comparison ####

AIC(
  m_feeding_pref_main,
  m_feeding_pref_cots_family,
  m_feeding_pref_cots_sub_family
)

anova(
  m_feeding_pref_main,
  m_feeding_pref_cots_family,
  m_feeding_pref_cots_sub_family
)

#### 05. Quick model checks ####

check_convergence(m_feeding_pref_cots_family)
check_singularity(m_feeding_pref_cots_family)

check_convergence(m_feeding_pref_cots_sub_family)
check_singularity(m_feeding_pref_cots_sub_family)

#### 06. DHARMa diagnostics ####

set.seed(123)

sim_cots_family <- simulateResiduals(m_feeding_pref_cots_family, n = 1000)
plot(sim_cots_family)
testUniformity(sim_cots_family)
testDispersion(sim_cots_family)
testOutliers(sim_cots_family)
testZeroInflation(sim_cots_family)

sim_cots_sub_family <- simulateResiduals(m_feeding_pref_cots_sub_family, n = 1000)
plot(sim_cots_sub_family)
testUniformity(sim_cots_sub_family)
testDispersion(sim_cots_sub_family)
testOutliers(sim_cots_sub_family)
testZeroInflation(sim_cots_sub_family)

#### 07. Fixed effects summaries ####

feeding_cots_family_results <- tidy(
  m_feeding_pref_cots_family,
  effects = "fixed",
  conf.int = TRUE
)

feeding_cots_sub_family_results <- tidy(
  m_feeding_pref_cots_sub_family,
  effects = "fixed",
  conf.int = TRUE
)

feeding_cots_family_results
feeding_cots_sub_family_results
#### 08. Simple decision table ####

model_compare_feeding <- tibble(
  model = c("main", "family × CoTS", "substrate × family × CoTS"),
  AIC = AIC(m_feeding_pref_main,
            m_feeding_pref_cots_family,
            m_feeding_pref_cots_sub_family)$AIC
) %>%
  mutate(delta_AIC = AIC - min(AIC))

model_compare_feeding

#### FEEDING RATE MODEL: total bite rate by CoTS site type ####

fish_bites_rate_mod <- fish_bites_mod %>%
  mutate(
    site_type = factor(site_type, levels = c("Control", "High Density")),
    family = fct_relevel(family, "Parrotfish", "Rabbitfish", "Butterflyfish"),
    minutes_obs = as.numeric(minutes_obs)
  ) %>%
  filter(
    !is.na(bites_total),
    !is.na(minutes_obs),
    minutes_obs > 0,
    !is.na(site_type),
    !is.na(family)
  )

m_bite_rate_site_type <- glmmTMB(
  bites_total ~ family * site_type +
    offset(log(minutes_obs)) +
    (1 | site_code) +
    (1 | date),
  family = nbinom2,
  data = fish_bites_rate_mod
)

summary(m_bite_rate_site_type)

bite_rate_site_type_results <- broom.mixed::tidy(
  m_bite_rate_site_type,
  effects = "fixed",
  conf.int = TRUE,
  exponentiate = TRUE
) %>%
  mutate(
    estimate = round(estimate, 3),
    conf.low = round(conf.low, 3),
    conf.high = round(conf.high, 3),
    p_formatted = case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  )

bite_rate_site_type_results

set.seed(123)

sim_bite_rate_site_type <- DHARMa::simulateResiduals(m_bite_rate_site_type, n = 1000)
plot(sim_bite_rate_site_type)

DHARMa::testUniformity(sim_bite_rate_site_type)
DHARMa::testDispersion(sim_bite_rate_site_type)
DHARMa::testOutliers(sim_bite_rate_site_type)
DHARMa::testZeroInflation(sim_bite_rate_site_type)

### try a hurdle model ### 
#### FEEDING INTENSITY PART 1: any feeding vs no feeding ####

fish_bites_rate_mod <- fish_bites_mod %>%
  mutate(
    site_type = factor(site_type, levels = c("Control", "High Density")),
    family = fct_relevel(family, "Parrotfish", "Rabbitfish", "Butterflyfish"),
    minutes_obs = as.numeric(minutes_obs),
    any_bites = as.integer(bites_total > 0)
  ) %>%
  filter(
    !is.na(bites_total),
    !is.na(minutes_obs),
    minutes_obs > 0,
    !is.na(site_type),
    !is.na(family)
  )

m_bite_any <- glmmTMB(
  any_bites ~ family * site_type +
    (1 | site_code) +
    (1 | date),
  family = binomial(link = "logit"),
  data = fish_bites_rate_mod
)

summary(m_bite_any)

set.seed(123)
sim_bite_any <- DHARMa::simulateResiduals(m_bite_any, n = 1000)
plot(sim_bite_any)
DHARMa::testUniformity(sim_bite_any)
DHARMa::testDispersion(sim_bite_any)
DHARMa::testOutliers(sim_bite_any)
DHARMa::testZeroInflation(sim_bite_any)

#### FEEDING INTENSITY PART 2: positive bite counts only ####

fish_bites_positive <- fish_bites_rate_mod %>%
  filter(bites_total > 0)

m_bite_positive <- glmmTMB(
  bites_total ~ family * site_type +
    offset(log(minutes_obs)) +
    (1 | site_code) +
    (1 | date),
  family = nbinom2,
  data = fish_bites_positive
)

summary(m_bite_positive)

set.seed(123)
sim_bite_positive <- DHARMa::simulateResiduals(m_bite_positive, n = 1000)
plot(sim_bite_positive)
DHARMa::testUniformity(sim_bite_positive)
DHARMa::testDispersion(sim_bite_positive)
DHARMa::testOutliers(sim_bite_positive)
DHARMa::testZeroInflation(sim_bite_positive)

m_bite_rate_zinb <- glmmTMB(
  bites_total ~ family * site_type +
    offset(log(minutes_obs)) +
    (1 | site_code) +
    (1 | date),
  ziformula = ~ family,
  family = nbinom2,
  data = fish_bites_rate_mod
)

summary(m_bite_rate_zinb)

set.seed(123)
sim_bite_rate_zinb <- DHARMa::simulateResiduals(m_bite_rate_zinb, n = 1000)
plot(sim_bite_rate_zinb)
DHARMa::testUniformity(sim_bite_rate_zinb)
DHARMa::testDispersion(sim_bite_rate_zinb)
DHARMa::testOutliers(sim_bite_rate_zinb)
DHARMa::testZeroInflation(sim_bite_rate_zinb)
