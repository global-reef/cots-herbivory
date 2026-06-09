### 01_CLEAN ####

#### helpers ####

site_code_lookup <- function(site) {
  case_when(
    str_detect(site, regex("^aow\\s*leuk|^ao\\s*leuk", ignore_case = TRUE)) ~ "AL",
    str_detect(site, regex("^shark\\s*island", ignore_case = TRUE)) ~ "SI",
    str_detect(site, regex("^tanote", ignore_case = TRUE)) ~ "TB",
    str_detect(site, regex("^red\\s*rock", ignore_case = TRUE)) ~ "RR",
    str_detect(site, regex("^green\\s*rock", ignore_case = TRUE)) ~ "GR",
    str_detect(site, regex("^twins", ignore_case = TRUE)) ~ "TW",
    TRUE ~ NA_character_
  )
}

clean_family <- function(x) str_to_title(str_trim(x))

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
} # for cots surveys that we did but found no specimens

#### 01. LOAD RAW DATA ####

fish_bites_raw <- read.csv(file.path(data_raw_dir, "2026.06.09_Fish_Bites.csv"), check.names = FALSE)
fish_abund_raw <- read.csv(file.path(data_raw_dir, "2026.06.09_Fish_Abund.csv"), check.names = FALSE)
cots_abund_raw <- read.csv(file.path(data_raw_dir, "2026.06.09_COTS_Abund.csv"), check.names = FALSE)
### 01_CLEAN ####

#### helpers ####

site_code_lookup <- function(site) {
  case_when(
    str_detect(site, regex("^aow\\s*leuk|^ao\\s*leuk", ignore_case = TRUE)) ~ "AL",
    str_detect(site, regex("^shark\\s*island", ignore_case = TRUE)) ~ "SI",
    str_detect(site, regex("^tanote", ignore_case = TRUE)) ~ "TB",
    str_detect(site, regex("^red\\s*rock", ignore_case = TRUE)) ~ "RR",
    str_detect(site, regex("^green\\s*rock", ignore_case = TRUE)) ~ "GR",
    str_detect(site, regex("^twins", ignore_case = TRUE)) ~ "TW",
    TRUE ~ NA_character_
  )
}

clean_family <- function(x) str_to_title(str_trim(x))

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
} # for cots surveys that we did but found no specimens

#### 01. LOAD RAW DATA ####


fish_bites_file <- file.path(data_raw_dir, paste0(analysis_date, "_Fish_Bites.csv"))
fish_abund_file <- file.path(data_raw_dir, paste0(analysis_date, "_Fish_Abund.csv"))
cots_abund_file <- file.path(data_raw_dir, paste0(analysis_date, "_COTS_Abund.csv"))
cpce_raw  <- file.path(data_raw_dir, paste0(analysis_date, "_cpce_long.csv"))

str(fish_bites_raw)
str(fish_abund_raw)
str(cots_abund_raw)
str(cpce_raw)


#### 02. CLEAN FISH BITES ####

bite_cols <- c("macro_algal_bites", "turf_filamentous_algal_bites", "live_coral_bites",
               "dead_coral_with_algae", "co_ts_scars", "sponge_bites", "sediment", "unknown")

fish_bites <- fish_bites_raw %>%
  clean_names() %>%
  mutate(across(c(site, buddy, family, species), str_trim),
         site_code = site_code_lookup(site),
         family = clean_family(family),
         species = str_trim(species),
         date = mdy(date),
         bites_total = rowSums(across(all_of(bite_cols)), na.rm = TRUE),
         minutes_obs = time_s / 60,
         bites_min = bites_total / minutes_obs) %>%
  left_join(fish_bite_species_lookup, by = c("family", "species")) %>%
  filter(!is.na(site_code), !is.na(date), time_s > 0) %>%
  mutate(across(c(site, buddy, family, species, site_code, trophic_group,
                  taxon_family, genus, species_epithet, sci_name, taxon_notes),
                as.factor))

str(fish_bites)

fish_bites %>% count(site_code, family, species, sci_name)

fish_bites %>%
  filter(is.na(sci_name)) %>%
  distinct(family, species)

write.csv(fish_bites, file.path(data_clean_dir, "fish_bites_clean.csv"), row.names = FALSE)


#### 03. CLEAN FISH ABUNDANCE ####

fish_families <- c("parrotfish", "rabbitfish", "butterflyfish")

fish_abund <- fish_abund_raw %>%
  clean_names() %>%
  select(survey_id = surveycode, site, date = date_mm_dd_yyyy, avg_depth_m, time_24hr,
         transect_m, visibility_m, weather = weather_1_5, current, local_boats,
         researcher = researchers, all_of(fish_families)) %>%
  mutate(across(c(survey_id, site, researcher), str_trim),
         site_code = site_code_lookup(site),
         date = mdy(date),
         area_m2 = 1000,
         area_ha = area_m2 / 10000) %>%
  pivot_longer(all_of(fish_families), names_to = "family", values_to = "fish_count") %>%
  mutate(family = clean_family(family),
         fish_density_ha = fish_count / area_ha) %>%
  left_join(functional_taxa, by = "family") %>%
  filter(!is.na(site_code), !is.na(date), !is.na(fish_count)) %>%
  mutate(across(c(survey_id, site, site_code, time_24hr, current, researcher,
                  family, trophic_group, taxon_family, genus, species_epithet, sci_name),
                as.factor)) %>%
  arrange(site_code, date, survey_id, researcher, family)

str(fish_abund)

fish_abund %>% count(site_code, family, sci_name)

fish_abund %>%
  filter(is.na(sci_name)) %>%
  distinct(family)

write.csv(fish_abund, file.path(data_clean_dir, "fish_abund_clean.csv"), row.names = FALSE)

#### 04. CLEAN COTS ABUNDANCE ####

cots_indiv <- cots_abund_raw %>%
  clean_names() %>%
  rename(date = date_mm_dd_yyyy, growth_form = growth_form, removal = removal) %>%
  mutate(across(c(site, survey_id, researcher, active, substrate, growth_form, site_type, removal), str_trim),
         site_code = site_code_lookup(site), date = mdy(date),
         substrate = str_to_upper(substrate), growth_form = str_to_upper(growth_form),
         active = str_to_upper(active), removal = str_to_upper(removal),
         size_cm = readr::parse_number(size_cm),
         depth_m = readr::parse_number(depth_m),
         area_m2 = 2000, area_ha = area_m2 / 10000) %>% # 200m * 5m inclusion each diver * 2 dives
  filter(!is.na(site_code), !is.na(date), !is.na(survey_id))

cots_survey <- cots_indiv %>%
  group_by(site_code, site, site_type, survey_id, date, avg_depth_m, time,
           duration_min, vis_m, researcher, area_m2, area_ha) %>%
  summarise(cots_count = sum(specimen > 0, na.rm = TRUE),
            mean_size_cm = safe_mean(size_cm),
            mean_depth_m = safe_mean(depth_m),
            .groups = "drop") %>%
  mutate(cots_density_ha = cots_count / area_ha) %>%
  arrange(site_code, date, survey_id)

str(cots_indiv)
str(cots_survey)
# check
cots_survey %>%
  group_by(site_code, site_type) %>%
  summarise(
    n_surveys = n(),
    n_zero = sum(cots_count == 0),
    prop_zero = mean(cots_count == 0),
    mean_density_ha = mean(cots_density_ha, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(cots_indiv, file.path(data_clean_dir, "cots_indiv_clean.csv"), row.names = FALSE)
write.csv(cots_survey, file.path(data_clean_dir, "cots_survey_clean.csv"), row.names = FALSE)


#### 05. CLEAN SUBSTRATE ####

substrate <- substrate_raw %>%
  clean_names() %>%
  mutate(across(c(site, transect, substrate, genus, morph, scar), str_trim),
         site_code = site_code_lookup(site),
         substrate = str_to_upper(substrate),
         morph = str_to_upper(morph),
         scar = str_to_upper(scar)) %>%
  filter(!is.na(site_code), !is.na(substrate))

#### 05A. SUBSTRATE TRANSECT COVER ####

substrate_transect <- substrate %>%
  count(site_code, site, transect, substrate, name = "n_points") %>%
  group_by(site_code, site, transect) %>%
  mutate(prop = n_points / sum(n_points),
         pct = prop * 100) %>%
  ungroup()

#### 05B. SUBSTRATE WIDE ####

substrate_wide <- substrate_transect %>%
  select(site_code, site, transect, substrate, pct) %>%
  pivot_wider(names_from = substrate, values_from = pct, values_fill = 0)


str(fish_bites_raw)
str(fish_abund_raw)
str(cots_abund_raw)
# str(substrate_raw)


#### 02. CLEAN FISH BITES ####

bite_cols <- c("macro_algal_bites", "turf_filamentous_algal_bites", "live_coral_bites",
               "dead_coral_with_algae", "co_ts_scars", "sponge_bites", "sediment", "unknown")

fish_bites <- fish_bites_raw %>%
  clean_names() %>%
  mutate(across(c(site, buddy, family, species), str_trim),
         site_code = site_code_lookup(site),
         family = clean_family(family),
         species = str_trim(species),
         date = mdy(date),
         bites_total = rowSums(across(all_of(bite_cols)), na.rm = TRUE),
         minutes_obs = time_s / 60,
         bites_min = bites_total / minutes_obs) %>%
  left_join(fish_bite_species_lookup, by = c("family", "species")) %>%
  filter(!is.na(site_code), !is.na(date), time_s > 0) %>%
  mutate(across(c(site, buddy, family, species, site_code, trophic_group,
                  taxon_family, genus, species_epithet, sci_name, taxon_notes),
                as.factor))

str(fish_bites)

fish_bites %>% count(site_code, family, species, sci_name)

fish_bites %>%
  filter(is.na(sci_name)) %>%
  distinct(family, species)

write.csv(fish_bites, file.path(data_clean_dir, "fish_bites_clean.csv"), row.names = FALSE)


#### 03. CLEAN FISH ABUNDANCE ####

fish_families <- c("parrotfish", "rabbitfish", "butterflyfish")

fish_abund <- fish_abund_raw %>%
  clean_names() %>%
  select(survey_id = surveycode, site, date = date_mm_dd_yyyy, avg_depth_m, time_24hr,
         transect_m, visibility_m, weather = weather_1_5, current, local_boats,
         researcher = researchers, all_of(fish_families)) %>%
  mutate(across(c(survey_id, site, researcher), str_trim),
         site_code = site_code_lookup(site),
         date = mdy(date),
         area_m2 = 1000,
         area_ha = area_m2 / 10000) %>%
  pivot_longer(all_of(fish_families), names_to = "family", values_to = "fish_count") %>%
  mutate(family = clean_family(family),
         fish_density_ha = fish_count / area_ha) %>%
  left_join(functional_taxa, by = "family") %>%
  filter(!is.na(site_code), !is.na(date), !is.na(fish_count)) %>%
  mutate(across(c(survey_id, site, site_code, time_24hr, current, researcher,
                  family, trophic_group, taxon_family, genus, species_epithet, sci_name),
                as.factor)) %>%
  arrange(site_code, date, survey_id, researcher, family)

str(fish_abund)

fish_abund %>% count(site_code, family, sci_name)

fish_abund %>%
  filter(is.na(sci_name)) %>%
  distinct(family)

write.csv(fish_abund, file.path(data_clean_dir, "fish_abund_clean.csv"), row.names = FALSE)

#### 04. CLEAN COTS ABUNDANCE ####

cots_indiv <- cots_abund_raw %>%
  clean_names() %>%
  rename(date = date_mm_dd_yyyy, growth_form = growth_form, removal = removal) %>%
  mutate(across(c(site, survey_id, researcher, active, substrate, growth_form, site_type, removal), str_trim),
         site_code = site_code_lookup(site), date = mdy(date),
         substrate = str_to_upper(substrate), growth_form = str_to_upper(growth_form),
         active = str_to_upper(active), removal = str_to_upper(removal),
         size_cm = readr::parse_number(size_cm),
         depth_m = readr::parse_number(depth_m),
         area_m2 = 2000, area_ha = area_m2 / 10000) %>% # 200m * 5m inclusion each diver * 2 dives
  filter(!is.na(site_code), !is.na(date), !is.na(survey_id))

cots_survey <- cots_indiv %>%
  group_by(site_code, site, site_type, survey_id, date, avg_depth_m, time,
           duration_min, vis_m, researcher, area_m2, area_ha) %>%
  summarise(cots_count = sum(specimen > 0, na.rm = TRUE),
            mean_size_cm = safe_mean(size_cm),
            mean_depth_m = safe_mean(depth_m),
            .groups = "drop") %>%
  mutate(cots_density_ha = cots_count / area_ha) %>%
  arrange(site_code, date, survey_id)

str(cots_indiv)
str(cots_survey)
# check
cots_survey %>%
  group_by(site_code, site_type) %>%
  summarise(
    n_surveys = n(),
    n_zero = sum(cots_count == 0),
    prop_zero = mean(cots_count == 0),
    mean_density_ha = mean(cots_density_ha, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(cots_indiv, file.path(data_clean_dir, "cots_indiv_clean.csv"), row.names = FALSE)
write.csv(cots_survey, file.path(data_clean_dir, "cots_survey_clean.csv"), row.names = FALSE)


#### 05. CLEAN SUBSTRATE ####

substrate <- substrate_raw %>%
  clean_names() %>%
  mutate(across(c(site, transect, substrate, genus, morph, scar), str_trim),
         site_code = site_code_lookup(site),
         substrate = str_to_upper(substrate),
         morph = str_to_upper(morph),
         scar = str_to_upper(scar)) %>%
  filter(!is.na(site_code), !is.na(substrate))

#### 05A. SUBSTRATE TRANSECT COVER ####

substrate_transect <- substrate %>%
  count(site_code, site, transect, substrate, name = "n_points") %>%
  group_by(site_code, site, transect) %>%
  mutate(prop = n_points / sum(n_points),
         pct = prop * 100) %>%
  ungroup()

#### 05B. SUBSTRATE WIDE ####

substrate_wide <- substrate_transect %>%
  select(site_code, site, transect, substrate, pct) %>%
  pivot_wider(names_from = substrate, values_from = pct, values_fill = 0)
