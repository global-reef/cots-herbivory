### data exploration
#### 0) Setup ####

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(glmmTMB)
  library(lme4)
})

# Assumes these objects already exist in your environment:
# cots            (site, survey_id, cots_count, area_m2, cots_density_ha)
# fish            (survey_id, site, date, researcher, family, fish_count, area_m2, fish_density_ha)
# feed_long       (site, row_id, family, species, date, time_s, bites_total, bites_min, ...)
# substrate_wide  (site, transect, HC, MAC, TUR)


#### 1) Quick structure checks ####

cat("\n## STR() CHECKS\n")
print(str(cots))
print(str(fish))
print(str(feed_long))
print(str(substrate_wide))

cat("\n## BASIC COUNTS\n")
cat("CoTS rows:", nrow(cots), " | sites:", n_distinct(cots$site), "\n")
cat("Fish rows:", nrow(fish), " | sites:", n_distinct(fish$site), " | families:", n_distinct(fish$family), "\n")
cat("Feed rows:", nrow(feed_long), " | sites:", n_distinct(feed_long$site), " | families:", n_distinct(feed_long$family), "\n")
cat("Substrate rows:", nrow(substrate_wide), " | sites:", n_distinct(substrate_wide$site), "\n")


#### 2) Feeding exploration ####

## 2.1 Distribution by family
p_feed_hist <- ggplot(feed_long, aes(bites_min)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ family, scales = "free_x") +
  labs(title = "Feeding rate distribution (bites/min) by family")

## 2.2 Boxplots by site (per family)
p_feed_site <- ggplot(feed_long, aes(site, bites_min)) +
  geom_boxplot() +
  facet_wrap(~ family, scales = "free_y") +
  labs(title = "Feeding rates by site (bites/min)", x = "Site", y = "Bites per min")

## 2.3 Check for zeros / near-zeros
feed_zero_check <- feed_long %>%
  summarise(
    n = n(),
    n_zero = sum(bites_min == 0, na.rm = TRUE),
    n_lt_0_5 = sum(bites_min < 0.5, na.rm = TRUE),
    min_val = min(bites_min, na.rm = TRUE),
    max_val = max(bites_min, na.rm = TRUE)
  )

print(feed_zero_check)
print(p_feed_hist)
print(p_feed_site)


#### 3) Fish exploration ####

## 3.1 Distribution by family
p_fish_hist <- ggplot(fish, aes(fish_density_ha)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ family, scales = "free_x") +
  labs(title = "Fish density distribution (per ha) by family")

## 3.2 Boxplots by site
p_fish_site <- ggplot(fish, aes(site, fish_density_ha)) +
  geom_boxplot() +
  facet_wrap(~ family, scales = "free_y") +
  labs(title = "Fish density by site (per ha)", x = "Site", y = "Fish per ha")

## 3.3 Replicate structure (researchers)
fish_rep_check <- fish %>%
  count(site, date, family) %>%
  summarise(
    n_groups = n(),
    min_reps = min(n),
    max_reps = max(n)
  )

print(fish_rep_check)
print(p_fish_hist)
print(p_fish_site)


#### 4) CoTS exploration ####

## 4.1 Distribution overall
p_cots_hist <- ggplot(cots, aes(cots_density_ha)) +
  geom_histogram(bins = 15) +
  labs(title = "CoTS density distribution (per ha)", x = "CoTS per ha")

## 4.2 Boxplots by site
p_cots_site <- ggplot(cots, aes(site, cots_density_ha)) +
  geom_boxplot() +
  labs(title = "CoTS density by site (per ha)", x = "Site", y = "CoTS per ha")

## 4.3 Within-site spread summary
cots_spread <- cots %>%
  group_by(site) %>%
  summarise(
    n = n(),
    mean = mean(cots_density_ha, na.rm = TRUE),
    sd = sd(cots_density_ha, na.rm = TRUE),
    med = median(cots_density_ha, na.rm = TRUE),
    iqr = IQR(cots_density_ha, na.rm = TRUE),
    .groups = "drop"
  )

print(cots_spread)
print(p_cots_hist)
print(p_cots_site)


#### 5) Substrate exploration ####

## 5.1 Long version for plotting
substrate_long <- substrate_wide %>%
  pivot_longer(cols = c(HC, MAC, TUR), names_to = "substrate", values_to = "pct")

## 5.2 Distributions
p_sub_hist <- ggplot(substrate_long, aes(pct)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ substrate, scales = "free_x") +
  labs(title = "Substrate % cover distributions", x = "% cover")

## 5.3 Boxplots by site
p_sub_site <- ggplot(substrate_long, aes(site, pct)) +
  geom_boxplot() +
  facet_wrap(~ substrate, scales = "free_y") +
  labs(title = "Substrate % cover by site", x = "Site", y = "% cover")

## 5.4 Per-transect sum check (only these 3 categories, so should be <= 100)
sub_sum_check <- substrate_wide %>%
  mutate(sum_3cats = HC + MAC + TUR) %>%
  summarise(
    min_sum = min(sum_3cats, na.rm = TRUE),
    max_sum = max(sum_3cats, na.rm = TRUE)
  )

print(sub_sum_check)
print(p_sub_hist)
print(p_sub_site)


#### 6) Better-than-mean site-level summaries (partial pooling) ####

## 6.1 CoTS: site random-intercept model (BLUP-style site effect)
m_cots_site <- lmer(cots_density_ha ~ 1 + (1 | site), data = cots)

cots_site_blup <- ranef(m_cots_site)$site %>%
  as.data.frame() %>%
  rownames_to_column("site") %>%
  rename(cots_site_effect = `(Intercept)`)

## Also store overall intercept (grand mean) so you can recover site-level means if needed
cots_grand_mean <- fixef(m_cots_site)[["(Intercept)"]]

cots_site_blup <- cots_site_blup %>%
  mutate(cots_site_mean_ha = cots_grand_mean + cots_site_effect)

print(summary(m_cots_site))
print(cots_site_blup)

## 6.2 Substrate: site random-intercept models (MAC + TUR)
m_MAC_site <- lmer(MAC ~ 1 + (1 | site), data = substrate_wide)
m_TUR_site <- lmer(TUR ~ 1 + (1 | site), data = substrate_wide)

MAC_grand <- fixef(m_MAC_site)[["(Intercept)"]]
TUR_grand <- fixef(m_TUR_site)[["(Intercept)"]]

benthic_site_blup <- tibble(
  site = rownames(ranef(m_MAC_site)$site),
  MAC_site_effect = ranef(m_MAC_site)$site[, 1],
  TUR_site_effect = ranef(m_TUR_site)$site[, 1]
) %>%
  mutate(
    MAC_site_mean = MAC_grand + MAC_site_effect,
    TUR_site_mean = TUR_grand + TUR_site_effect
  )

print(summary(m_MAC_site))
print(summary(m_TUR_site))
print(benthic_site_blup)


#### 7) Join-ready analysis table (feeding backbone) ####

## 7.1 Fish: summarise to site-date-family (keeps feeding variance; reduces researcher replicates)
fish_svy <- fish %>%
  group_by(site, date, family) %>%
  summarise(
    fish_density_ha = mean(fish_density_ha, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  )

## 7.2 Build table
dat <- feed_long %>%
  mutate(family = factor(family, levels = c("parrot","rabbit","butterfly"))) %>%
  left_join(fish_svy, by = c("site","date","family")) %>%
  left_join(cots_site_blup %>% select(site, cots_site_mean_ha), by = "site") %>%
  left_join(benthic_site_blup %>% select(site, MAC_site_mean, TUR_site_mean), by = "site")

## 7.3 Missingness check
dat %>%
  summarise(
    n = n(),
    miss_fish = sum(is.na(fish_density_ha)),
    miss_cots = sum(is.na(cots_site_mean_ha)),
    miss_MAC  = sum(is.na(MAC_site_mean)),
    miss_TUR  = sum(is.na(TUR_site_mean))
  ) %>%
  print()

## quick look
dat %>% glimpse()


#### 8) Optional: quick distribution checks on join-ready covariates ####

p_cov <- dat %>%
  select(bites_min, fish_density_ha, cots_site_mean_ha, MAC_site_mean, TUR_site_mean, family) %>%
  pivot_longer(-family, names_to = "var", values_to = "val") %>%
  ggplot(aes(val)) +
  geom_histogram(bins = 20) +
  facet_grid(var ~ family, scales = "free_x") +
  labs(title = "Distributions of response + covariates (join-ready table)")

print(p_cov)


# Based on the EDA, the defensible modelling sequence is:
#   Model class
# Gamma GLMM with log link
# Random intercept for site
# Family as fixed effect
# Hypothesis flow (aligned to your plan)
# Does CoTS relate to benthos?
#   (already partially answered descriptively, can be formalised)
# Does benthos relate to feeding, accounting for family?
#   Does CoTS still matter once benthos is included?
#   (attenuation = support for indirect pathway)
# You are not testing strict mediation in the causal-inference sense. You are testing consistent directional support, which is appropriate here.


