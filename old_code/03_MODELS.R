#### 0) Setup ####

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmmTMB)
  library(lme4)
  library(DHARMa)
})

# Assumes in memory:
# cots (site, cots_density_ha, ...)
# substrate_wide (site, MAC, TUR, ...)
# feed_long (site, date, family, bites_min, ...)
# fish (site, date, family, fish_density_ha, ...)


#### 1) Minimal site-level predictors (partial pooling) ####

## 1.1 CoTS: site pressure (BLUP mean per site)
m_cots_site <- lmer(cots_density_ha ~ 1 + (1 | site), data = cots)
cots_grand <- fixef(m_cots_site)[["(Intercept)"]]

cots_site <- ranef(m_cots_site)$site %>%
  as.data.frame() %>%
  rownames_to_column("site") %>%
  transmute(
    site,
    cots_site_mean_ha = cots_grand + `(Intercept)`
  )

## 1.2 Substrate: site benthos (MAC + TUR; BLUP means)
m_MAC_site <- lmer(MAC ~ 1 + (1 | site), data = substrate_wide)
m_TUR_site <- lmer(TUR ~ 1 + (1 | site), data = substrate_wide)
# --- Live coral cover (HC) site-level BLUP mean ---

m_HC_site <- lmer(HC ~ 1 + (1 | site), data = substrate_wide)

HC_grand <- fixef(m_HC_site)[["(Intercept)"]]

HC_site <- tibble(
  site = rownames(ranef(m_HC_site)$site),
  HC_site_mean = HC_grand + ranef(m_HC_site)$site[, 1]
)

m_CoTS_to_HC <- lm(
  HC_site_mean ~ cots_site_mean_ha,
  data = HC_site %>% left_join(cots_site, by = "site")
)

summary(m_CoTS_to_HC)

MAC_grand <- fixef(m_MAC_site)[["(Intercept)"]]
TUR_grand <- fixef(m_TUR_site)[["(Intercept)"]]

benthic_site <- tibble(
  site = rownames(ranef(m_MAC_site)$site),
  MAC_site_mean = MAC_grand + ranef(m_MAC_site)$site[, 1],
  TUR_site_mean = TUR_grand + ranef(m_TUR_site)$site[, 1]
)

## 1.3 Fish: summarise to site-date-family for joining to feeding
fish_svy <- fish %>%
  group_by(site, date, family) %>%
  summarise(
    fish_density_ha = mean(fish_density_ha, na.rm = TRUE),
    .groups = "drop"
  )


#### 2) Backbone analysis table (feeding) ####

dat <- feed_long %>%
  mutate(family = factor(family, levels = c("parrot","rabbit","butterfly"))) %>%
  left_join(fish_svy, by = c("site","date","family")) %>%
  left_join(cots_site, by = "site") %>%
  left_join(benthic_site, by = "site") %>%
  filter(!is.na(bites_min), bites_min > 0) %>%
  mutate(
    MAC_sc  = as.numeric(scale(MAC_site_mean)),
    TUR_sc  = as.numeric(scale(TUR_site_mean)),
    CoTS_sc = as.numeric(scale(cots_site_mean_ha)),
    fish_sc = as.numeric(scale(fish_density_ha))
  )

# optional: keep a no-fish version for max sample size
dat_nofish <- dat %>%
  select(-fish_density_ha, -fish_sc)


#### 3) Minimal model set (TUR + MAC) ####

## 3.1 CoTS -> benthos (site-level link tests)
m_CoTS_to_MAC <- lm(MAC_site_mean ~ cots_site_mean_ha, data = benthic_site %>% left_join(cots_site, by = "site"))
m_CoTS_to_TUR <- lm(TUR_site_mean ~ cots_site_mean_ha, data = benthic_site %>% left_join(cots_site, by = "site"))

summary(m_CoTS_to_MAC)
summary(m_CoTS_to_TUR)

## 3.2 Feeding ~ CoTS (total effect)
m_total <- glmmTMB(
  bites_min ~ CoTS_sc + family + (1 | site),
  family = Gamma(link = "log"),
  data = dat_nofish
)

## 3.3 Feeding ~ benthos (MAC + TUR)
m_benthic <- glmmTMB(
  bites_min ~ MAC_sc + TUR_sc + family + (1 | site),
  family = Gamma(link = "log"),
  data = dat_nofish
)

## 3.4 Feeding ~ CoTS + benthos (attenuation check)
m_full <- glmmTMB(
  bites_min ~ CoTS_sc + MAC_sc + TUR_sc + family + (1 | site),
  family = Gamma(link = "log"),
  data = dat_nofish
)

summary(m_total)
summary(m_benthic)
summary(m_full) # Higher macroalgal cover is associated with lower feeding rates, after accounting for family and CoTS. 

AIC(m_total, m_benthic, m_full)


#### 4) Optional single add-on: fish density (only if you want it) ####

dat_fish <- dat %>% filter(!is.na(fish_sc))

m_full_fish <- glmmTMB(
  bites_min ~ CoTS_sc + MAC_sc + TUR_sc + fish_sc + family + (1 | site),
  family = Gamma(link = "log"),
  data = dat_fish
)

summary(m_full_fish) # Higher macroalgal cover is associated with lower feeding rates, after accounting for family and CoTS. 
AIC(m_full, m_full_fish)


#### 5) Diagnostics (run for the main model you will report) ####

sim_full <- simulateResiduals(m_full, n = 500)
plot(sim_full)

testDispersion(sim_full)
testZeroInflation(sim_full)


#### 6) Quick effect-size helper (Gamma log link => multiplicative) ####

get_rr <- function(model, term){
  b <- fixef(model)$cond[term]
  tibble(
    term = term,
    beta = b,
    rate_ratio = exp(b)
  )
}
# this shows 
bind_rows(
  get_rr(m_full, "CoTS_sc"), # A 1 SD increase in site-level CoTS density is associated with an ~10% increase in feeding rate
  get_rr(m_full, "MAC_sc"), # A 1 SD increase in site-level macroalgal cover is associated with an ~10% decrease in feeding rate 
  get_rr(m_full, "TUR_sc") # A 1 SD increase in site-level CoTS density is associated with an ~10% increase in feeding rate
) %>% print()



#### RESULTS ######## 
summary(m_CoTS_to_MAC) 
# Call:
#   lm(formula = MAC_site_mean ~ cots_site_mean_ha, data = benthic_site %>% 
#        left_join(cots_site, by = "site"))
# 
# Residuals:
#   1        2        3        4        5        6 
# 0.12973 -0.44168  1.07237 -0.11586 -0.03638 -0.60818 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)        1.75164    0.49440   3.543    0.024 *
#   cots_site_mean_ha  0.01712    0.02251   0.760    0.489  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6608 on 4 degrees of freedom
# Multiple R-squared:  0.1263,	Adjusted R-squared:  -0.09215 
# F-statistic: 0.5781 on 1 and 4 DF,  p-value: 0.4894

summary(m_CoTS_to_TUR)
# Call:
#   lm(formula = TUR_site_mean ~ cots_site_mean_ha, data = benthic_site %>% 
#        left_join(cots_site, by = "site"))
# 
# Residuals:
#   1       2       3       4       5       6 
# 11.8529 -0.5111 -1.8416 -4.9296 -6.1190  1.5484 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)        28.8672     5.3991   5.347   0.0059 **
#   cots_site_mean_ha   0.3627     0.2458   1.476   0.2141   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 7.216 on 4 degrees of freedom
# Multiple R-squared:  0.3525,	Adjusted R-squared:  0.1906 
# F-statistic: 2.177 on 1 and 4 DF,  p-value: 0.2141


summary(m_total)
# Family: Gamma  ( log )
# Formula:          bites_min ~ CoTS_sc + family + (1 | site)
# Data: dat_nofish
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# 1048.9    1067.3    -518.4    1036.9       152 
# 
# Random effects:
#   
#   Conditional model:
#   Groups Name        Variance Std.Dev.
# site   (Intercept) 0.007123 0.0844  
# Number of obs: 158, groups:  site, 6
# 
# Dispersion estimate for Gamma family (sigma^2): 0.245 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.94276    0.07708  25.204  < 2e-16 ***
#   CoTS_sc       0.02703    0.05252   0.515    0.607    
# familyparrot  1.29879    0.09846  13.191  < 2e-16 ***
#   familyrabbit  0.78097    0.09717   8.037 9.18e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(m_benthic)
# Family: Gamma  ( log )
# Formula:          bites_min ~ MAC_sc + TUR_sc + family + (1 | site)
# Data: dat_nofish
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# 1048.4    1069.8    -517.2    1034.4       151 
# 
# Random effects:
#   
#   Conditional model:
#   Groups Name        Variance Std.Dev.
# site   (Intercept) 0.001525 0.03904 
# Number of obs: 158, groups:  site, 6
# 
# Dispersion estimate for Gamma family (sigma^2): 0.245 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.942935   0.070842  27.426  < 2e-16 ***
#   MAC_sc       -0.081569   0.043077  -1.894   0.0583 .  
# TUR_sc        0.004154   0.044089   0.094   0.9249    
# familyparrot  1.298999   0.098372  13.205  < 2e-16 ***
#   familyrabbit  0.781825   0.097422   8.025 1.01e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(m_full) 
# Family: Gamma  ( log )
# Formula:          bites_min ~ CoTS_sc + MAC_sc + TUR_sc + family + (1 | site)
# Data: dat_nofish
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# 1047.3    1071.8    -515.7    1031.3       150 
# 
# Random effects:
#   
#   Conditional model:
#   Groups Name        Variance  Std.Dev. 
# site   (Intercept) 4.626e-10 2.151e-05
# Number of obs: 158, groups:  site, 6
# 
# Dispersion estimate for Gamma family (sigma^2): 0.242 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.94474    0.06855  28.369  < 2e-16 ***
#   CoTS_sc       0.09262    0.05159   1.795   0.0726 .  
# MAC_sc       -0.10427    0.04201  -2.482   0.0131 *  
#   TUR_sc       -0.04834    0.04909  -0.985   0.3247    
# familyparrot  1.29891    0.09644  13.469  < 2e-16 ***
#   familyrabbit  0.77156    0.09699   7.955 1.79e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(m_full_fish) 
# Family: Gamma  ( log )
# Formula:          bites_min ~ CoTS_sc + MAC_sc + TUR_sc + fish_sc + family + (1 |  
#                                                                                 site)
# Data: dat_fish
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# 454.7     474.9    -218.4     436.7        61 
# 
# Random effects:
#   
#   Conditional model:
#   Groups Name        Variance  Std.Dev. 
# site   (Intercept) 1.653e-10 1.286e-05
# Number of obs: 70, groups:  site, 6
# 
# Dispersion estimate for Gamma family (sigma^2): 0.203 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   1.82189    0.09893  18.416  < 2e-16 ***
#   CoTS_sc       0.17529    0.09380   1.869   0.0617 .  
# MAC_sc       -0.11439    0.05895  -1.940   0.0523 .  
# TUR_sc       -0.09730    0.08831  -1.102   0.2706    
# fish_sc       0.05242    0.06705   0.782   0.4344    
# familyparrot  1.30434    0.14289   9.128  < 2e-16 ***
#   familyrabbit  0.95744    0.13705   6.986 2.83e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# => No evidence that higher CoTS density is associated with increase Macroalgae cover.
# => Turf cover shows a positive but *non significant* association with cots density.
# => Higher macroalgal cover is associated with lower feeding rates, after accounting for family and CoTS.
# => Turf cover does not explain variation in bite rates once macroalgae and family are accounted for.
# => CoTS density was not directly associated with reduced feeding rates; however, the estimated effect of CoTS became more positive once benthic state was accounted for, suggesting indirect or context-dependent relationships rather than direct suppression.
# => Feeding rates are not strongly driven by local fish density in this dataset, likely because behaviour is more substrate-limited than abundance-limited.


write.csv(
  dat_nofish,
  file = "finaldata_scar.csv",
  row.names = FALSE
)

write.csv(
  dat,
  file = "finaldata_withfish_scar.csv",
  row.names = FALSE
)
sink("model_outputs.txt")

summary(m_CoTS_to_MAC)
summary(m_CoTS_to_TUR)
summary(m_total)
summary(m_benthic)
summary(m_full)
summary(m_full_fish)

AIC(m_total, m_benthic, m_full, m_full_fish)

sink()
cat("\n==============================\n")
cat("FINAL MODEL OUTPUTS\n")
cat("==============================\n\n")

cat("\n--- CoTS → Macroalgae ---\n")
print(summary(m_CoTS_to_MAC))

cat("\n--- CoTS → Turf ---\n")
print(summary(m_CoTS_to_TUR))

cat("\n--- Feeding ~ CoTS (total effect) ---\n")
print(summary(m_total))

cat("\n--- Feeding ~ Benthos ---\n")
print(summary(m_benthic))

cat("\n--- Feeding ~ CoTS + Benthos ---\n")
print(summary(m_full))

cat("\n--- Feeding ~ CoTS + Benthos + Fish ---\n")
print(summary(m_full_fish))

cat("\n--- AIC comparison ---\n")
print(AIC(m_total, m_benthic, m_full, m_full_fish))






library(tidyverse)

# trasformazione long
substrate_long <- substrate_site %>%
  pivot_longer(cols = c(HC, MAC, TUR),
               names_to = "substrate",
               values_to = "cover")

# nomi per i facet
substrate_long$substrate <- factor(
  substrate_long$substrate,
  levels = c("HC","MAC","TUR"),
  labels = c("Hard coral (%)",
             "Macroalgae (%)",
             "Turf algae (%)")
)

# grafico
p_sub_cots <- ggplot(substrate_long,
                     aes(x = cots_density,
                         y = cover)) +
  
  geom_smooth(method = "lm",
              se = TRUE,
              color = "black",
              fill = "grey70") +
  
  geom_point(size = 4) +
  
  facet_wrap(~ substrate,
             scales = "free_y") +
  
  labs(
    x = expression("CoTS density (ind ha"^-1*")"),
    y = "Mean benthic cover (%)"
  ) +
  
  theme_classic() +
  
  theme(
    text = element_text(family = "Times New Roman", size = 38),
    strip.background = element_rect(color = "black",
                                    fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

print(p_sub_cots)

ggsave(
  "Fig10_substrate_vs_CoTS.png",
  plot = p_sub_cots,
  width = 20,
  height = 12,
  units = "cm",
  dpi = 300
)

library(dplyr)
library(ggplot2)

#-----------------------------------------------------------
# 1. Calcolo macroalgal cover: media dei 4 transetti per sito
#-----------------------------------------------------------

macro_site <- substrate_transect %>%
  filter(substrate_grp == "MAC") %>%
  group_by(site, transect) %>%
  summarise(
    mac_transect = mean(pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(site) %>%
  summarise(
    MAC = mean(mac_transect, na.rm = TRUE),
    .groups = "drop"
  )

#-----------------------------------------------------------
# 2. Unione con feeding dataset
#-----------------------------------------------------------

feed_macro <- feed_long %>%
  left_join(macro_site, by = "site") %>%
  mutate(
    family = tools::toTitleCase(as.character(family))
  )

feed_macro$family <- factor(
  feed_macro$family,
  levels = c("Butterfly","Parrot","Rabbit")
)

#-----------------------------------------------------------
# 3. Grafico feeding vs macroalgae
#-----------------------------------------------------------

p_feed_macro <- ggplot(feed_macro,
                       aes(x = MAC,
                           y = bites_min)) +
  
  geom_point(alpha = 0.6, size = 2) +
  
  geom_smooth(method = "lm",
              color = "#3b82f6",
              fill = "grey70") +
  
  facet_wrap(~ family) +
  
  labs(
    x = "Macroalgal cover (%)",
    y = "Feeding rate (bites per minute)"
  ) +
  
  theme_classic() +
  
  theme(
    text = element_text(family = "Times New Roman", size = 38),
    strip.background = element_rect(color = "black", fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

print(p_feed_macro)

#-----------------------------------------------------------
# 4. Salvataggio figura
#-----------------------------------------------------------

ggsave(
  "Fig11_feeding_vs_macroalgae.png",
  plot = p_feed_macro,
  width = 20,
  height = 12,
  units = "cm",
  dpi = 300
)

print(macro_site)

substrate_transect %>%
  filter(substrate_grp == "MAC") %>%
  group_by(site, transect) %>%
  summarise(mac_transect = mean(pct, na.rm = TRUE))












feed_turf <- feed_long %>%
  mutate(
    family = tools::toTitleCase(as.character(family)),
    prop_turf = tur / bites_total
  )

feed_turf$family <- factor(
  feed_turf$family,
  levels = c("Butterfly","Parrot","Rabbit")
)

p_turf_cots <- ggplot(feed_turf,
                      aes(x = cots_density,
                          y = prop_turf)) +
  
  geom_point(alpha = 0.6, size = 2) +
  
  geom_smooth(method = "lm",
              color = "#3b82f6",
              fill = "grey70") +
  
  facet_wrap(~ family) +
  
  labs(
    x = expression("CoTS density (ind ha"^-1*")"),
    y = "Proportion of bites on turf"
  ) +
  
  theme_classic() +
  
  theme(
    text = element_text(family = "Times New Roman", size = 38),
    strip.background = element_rect(color = "black", fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

print(p_turf_cots)




library(dplyr)
library(ggplot2)

# densità CoTS per sito
cots_site <- data.frame(
  site = c("TB","SI","AL","GR","TW","RR"),
  cots_density = c(4.0,6.5,7.2,28.9,29.3,29.4)
)

# unione con dataset feeding
feed_turf <- feed_long %>%
  left_join(cots_site, by = "site") %>%
  mutate(
    family = tools::toTitleCase(as.character(family)),
    prop_turf = tur / bites_total
  )

feed_turf$family <- factor(
  feed_turf$family,
  levels = c("Butterfly","Parrot","Rabbit")
)

# grafico
p_turf_cots <- ggplot(feed_turf,
                      aes(x = cots_density,
                          y = prop_turf)) +
  
  geom_point(alpha = 0.6, size = 2) +
  
  geom_smooth(method = "lm",
              color = "#3b82f6",
              fill = "grey70") +
  
  facet_wrap(~ family) +
  
  labs(
    x = expression("CoTS density (ind ha"^-1*")"),
    y = "Proportion of bites on turf"
  ) +
  
  theme_classic() +
  
  theme(
    text = element_text(family = "Times New Roman", size = 38),
    strip.background = element_rect(color = "black", fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

print(p_turf_cots)

# salvataggio
ggsave(
  "Fig12_turf_feeding_vs_cots.png",
  plot = p_turf_cots,
  width = 20,
  height = 12,
  units = "cm",
  dpi = 300
)

