############################################################
#### 04_DIET_PREFERENCE.R
#### Dietary preference model:
#### substrate choice ~ availability + CoTS + family
############################################################

#### 0) Setup ####

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmmTMB)
})

# Assumes in memory:
# feed_long        (site, family, bites_total, mac, tur, hc, ...)
# substrate_wide   (site, transect, MAC, TUR, HC)
# cots_site        (site, cots_site_mean_ha)


############################################################
#### 1) Substrate availability at site level
#### (FIX many-to-many joins)
############################################################

substrate_site <- substrate_wide %>%
  group_by(site) %>%
  summarise(
    MAC = mean(MAC, na.rm = TRUE),
    TUR = mean(TUR, na.rm = TRUE),
    HC  = mean(HC,  na.rm = TRUE),
    .groups = "drop"
  )


############################################################
#### 2) Build diet choice dataset
############################################################

diet_long <- feed_long %>%
  left_join(substrate_site, by = "site") %>%
  left_join(cots_site, by = "site") %>%
  mutate(
    CoTS_sc = as.numeric(scale(cots_site_mean_ha))
  ) %>%
  pivot_longer(
    cols = c(mac, tur, hc),
    names_to = "substrate",
    values_to = "bites_sub"
  ) %>%
  mutate(
    substrate = recode(substrate,
                       mac = "MAC",
                       tur = "TUR",
                       hc  = "HC"),
    bites_other = bites_total - bites_sub,
    avail_pct = case_when(
      substrate == "MAC" ~ MAC,
      substrate == "TUR" ~ TUR,
      substrate == "HC"  ~ HC
    ) / 100
  ) %>%
  filter(
    bites_total > 0,
    bites_other >= 0,
    !is.na(avail_pct)
  ) %>%
  mutate(
    family = factor(family),
    substrate = factor(substrate, levels = c("MAC", "TUR", "HC"))
  )


############################################################
#### 3) Dietary preference model (resource selection)
############################################################

m_diet_pref <- glmmTMB(
  cbind(bites_sub, bites_other) ~
    substrate * family +        # family-specific preference
    avail_pct +                 # availability control
    CoTS_sc +                   # overall CoTS effect
    substrate:CoTS_sc +         # change in preference with CoTS
    (1 | site),                 # site-level structure
  family = binomial(link = "logit"),
  data = diet_long
)

summary(m_diet_pref)


############################################################
#### 4) Optional: simpler nested models (for comparison)
############################################################

## 4.1 No CoTS (pure preference + availability)
m_diet_noCoTS <- glmmTMB(
  cbind(bites_sub, bites_other) ~
    substrate * family +
    avail_pct +
    (1 | site),
  family = binomial(link = "logit"),
  data = diet_long
)

## 4.2 No family (community-level preference only)
m_diet_nofam <- glmmTMB(
  cbind(bites_sub, bites_other) ~
    substrate +
    avail_pct +
    CoTS_sc +
    substrate:CoTS_sc +
    (1 | site),
  family = binomial(link = "logit"),
  data = diet_long
)

AIC(m_diet_pref, m_diet_noCoTS, m_diet_nofam)


############################################################
#### 5) Helper: odds ratios for interpretation
############################################################

get_or <- function(model, term) {
  est <- summary(model)$coefficients$cond[term, "Estimate"]
  se  <- summary(model)$coefficients$cond[term, "Std. Error"]
  
  tibble(
    term = term,
    OR = exp(est),
    OR_low = exp(est - 1.96 * se),
    OR_high = exp(est + 1.96 * se)
  )
}

## Examples:
get_or(m_diet_pref, "substrateTUR")
get_or(m_diet_pref, "substrateHC")
get_or(m_diet_pref, "substrateTUR:CoTS_sc")
summary(m_diet_pref)
while (sink.number() > 0) sink()

cat("\n==============================\n")
cat("DIET PREFERENCE MODEL OUTPUT\n")
cat("==============================\n\n")

print(summary(m_diet_pref))


