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
raw_file_names <- c(
  fish_bites = "Fish_Bites.csv",
  fish_abund = "Fish_Abund.csv",
  cots_abund = "COTS_Abund.csv",
  cpce = "cpce_long.csv"
)

raw_files <- file.path(data_raw_dir, paste0(analysis_date, "_", raw_file_names)) %>%
  set_names(names(raw_file_names))

raw_data <- raw_files %>%
  map(~ read_csv(.x, show_col_types = FALSE) %>%
        select(-matches("^\\.\\.\\.[0-9]+$"))) %>%
  set_names(names(raw_files))

iwalk(raw_data, ~ {
  cat("\n---", .y, "---\n")
  str(.x)
})

fish_bites_raw <- raw_data$fish_bites
fish_abund_raw <- raw_data$fish_abund
cots_abund_raw <- raw_data$cots_abund
cpce_raw <- raw_data$cpce



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

# Fix likely raw-entry date typo
fish_bites <- fish_bites %>%
  mutate(
    date = case_when(
      site_code == "AL" & date == as.Date("2029-09-29") ~ as.Date("2025-09-29"),
      TRUE ~ date
    )
  )

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
  arrange(site_code, date, survey_id)  %>% 
  mutate(across(c(site_code, site, site_type, survey_id, researcher), as.factor))

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
#### 05. CLEAN SUBSTRATE ####

valid_substrate <- c("HC", "SC", "TUR", "MAC", "SP", "OB", "AB", "AN", "UKN")
valid_gf <- c("MA", "SMA", "CAE", "CB", "TB", "ARB", "DI", "SOL", "FOL", "LM", "ENC")

substrate <- cpce_raw %>%
  clean_names() %>%
  mutate(
    across(c(site, transect, substrate, gf, scar), ~ str_to_upper(str_trim(.x))),
    site = str_to_title(str_to_lower(site)),
    site_code = site_code_lookup(site),
    quadrat = as.integer(quadrat),
    cpce_point = as.integer(cpce_point),
    
    substrate = case_when(
      substrate %in% valid_substrate ~ substrate,
      substrate %in% c("UK", "??", "@@@@", "0", "NUMBER NOT IN IMAGE", "UKN/OB") ~ "UKN",
      substrate %in% c("TUIR", "TURB", "TURF", "TURTUR", "ATUR", "UR") ~ "TUR",
      substrate %in% c("A B", "ABA", "RAB", "SAB", "HAB") ~ "AB",
      substrate %in% c("OB?", "OB (SP)", "O") ~ "OB",
      substrate %in% c("HHC", "HCC", "HCH", "HCMA", "HCENC", "HCTUR") ~ "HC",
      TRUE ~ NA_character_
    ),
    
    gf = case_when(
      gf %in% valid_gf ~ gf,
      gf %in% c("MAS", "M", "MAMA") ~ "MA",
      gf %in% c("SUB", "SM", "MS") ~ "SMA",
      gf == "DIG" ~ "DI",
      gf == "LAM" ~ "LM",
      gf == "TAB" ~ "TB",
      gf %in% c("FOLF", "FO") ~ "FOL",
      gf %in% c("SM/ENC", "ENC/SM") ~ "ENC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(site_code), !is.na(substrate)) %>%
  select(site_code, site, transect, quadrat, cpce_point, substrate, gf, scar) %>%
  mutate(across(c(site_code, site, transect, substrate, gf, scar), as.factor))

str(substrate)

#### 05A. SUBSTRATE TRANSECT COVER ####

substrate_transect <- substrate %>%
  count(site_code, site, transect, substrate, name = "n_points") %>%
  group_by(site_code, site, transect) %>%
  mutate(
    total_points = sum(n_points),
    prop = n_points / total_points,
    pct = prop * 100
  ) %>%
  ungroup() %>%
  mutate(across(c(site_code, site, transect, substrate), as.factor))

str(substrate_transect)
#### 05B. SUBSTRATE WIDE ####

substrate_wide <- substrate_transect %>%
  select(site_code, site, transect, substrate, pct) %>%
  pivot_wider(
    names_from = substrate,
    values_from = pct,
    values_fill = 0
  ) %>%
  mutate(across(c(site_code, site, transect), as.factor))

str(substrate_wide)


write_csv(substrate, file.path(data_clean_dir, "substrate_point_clean.csv"))
write_csv(substrate_transect, file.path(data_clean_dir, "substrate_transect_cover_clean.csv"))
write_csv(substrate_wide, file.path(data_clean_dir, "substrate_transect_cover_wide_clean.csv"))



# we now have cleaned 
str(fish_bites)
str(fish_abund)

str(cots_indiv)
str(cots_survey)

str(substrate)
str(substrate_transect)
str(substrate_wide)




# fish_bites          # 336 feeding observations
# fish_abund          # 327 observations (only parrot, rabbit and butterfly), 8826 fish counted 
# cots_indiv          # 698 indivudal cots records
# cots_survey         # 147 surveys
# substrate           # 50,239 CPCE points
# substrate_transect  # 255 substrate-by-transect quadrats



