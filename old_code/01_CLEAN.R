suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(purrr)
})

raw_dir <- "RAW_data"

# ----------------------
# constants + helpers
# ----------------------
sites <- c("AL","SI","TB","RR","GR","TW")

clean_site <- function(x) x %>% str_trim() %>% str_to_upper()
clean_family <- function(x) x %>% str_trim() %>% str_to_lower()
as_num <- function(x) readr::parse_number(as.character(x))

# area rules
# area helpers
fish_area_m2 <- function(area_m2) area_m2
cots_area_m2 <- function(area_m2) area_m2

# density conversion
to_density_ha <- function(count, area_m2) count / (area_m2 / 10000)

# ----------------------
# helper: read all sheets from a workbook
# assumes each sheet name is a site code
# ----------------------
read_sites_xlsx <- function(path, sites_keep) {
  
  sheet_names <- readxl::excel_sheets(path)
  keep <- base::intersect(sheet_names, sites_keep)
  
  if (length(keep) == 0) {
    stop("No matching site sheets found in: ", basename(path),
         "\nSheets present: ", paste(sheet_names, collapse = ", "))
  }
  
  purrr::map_dfr(
    keep,
    ~ readxl::read_excel(path, sheet = .x) %>%
      dplyr::mutate(site = .x)
  )
}

cots_path <- file.path(raw_dir, "Cots_abundance_RAW.xlsx")
cots_raw <- read_excel(cots_path)
str(cots_raw)

  
  ## COTS #### 
cots <- cots_raw %>%
  rename(
    site        = SITE,
    survey_id   = SURVEY,
    cots_count  = SPECIMEN,
    area_m2     = AREA
  ) %>%
  mutate(
    site             = clean_site(site),
    survey_id        = str_trim(survey_id),
    cots_count       = as.numeric(cots_count),
    area_m2=as.numeric(area_m2), # fix here
    cots_density_ha  = cots_count / (area_m2 / 10000)
  ) %>%
  select(site, survey_id, cots_count, area_m2, cots_density_ha)

## 02 FEEDING ######
feed_path <- file.path(raw_dir, "Feeding_behaviour_RAW.xlsx")
feed_raw <- read_sites_xlsx(feed_path, sites_keep = sites)
str(feed_raw)

# clean 
# helper
as_num <- function(x) readr::parse_number(as.character(x))

feed_long <- feed_raw %>%
  mutate(site = str_trim(site) %>% str_to_upper()) %>%
  mutate(row_id = row_number())

# 1) Parrot block
feed_parrot <- feed_long %>%
  transmute(
    site,
    row_id,
    family = "parrot",
    species = .data$`Species...2`,
    date    = as.Date(.data$`Date...3`),
    time_s  = as_num(.data$`Time...4`),
    mac = as_num(.data$`Mac...5`),
    tur = as_num(.data$`Tur...6`),
    hc  = as_num(.data$`Hc...7`),
    spo = as_num(.data$`Spo...8`),
    sed = as_num(.data$`Sed...9`),
    uk  = as_num(.data$`Uk...10`)
  )

# 2) Rabbit block
feed_rabbit <- feed_long %>%
  transmute(
    site,
    row_id,
    family = "rabbit",
    species = .data$`Species...12`,
    date    = as.Date(.data$`Date...13`),
    time_s  = as_num(.data$`Time...14`),
    mac = as_num(.data$`Mac...15`),
    tur = as_num(.data$`Tur...16`),
    hc  = as_num(.data$`Hc...17`),
    spo = as_num(.data$`Spo...18`),
    sed = as_num(.data$`Sed...19`),
    uk  = as_num(.data$`Uk...20`)
  )

# 3) Butterfly block
feed_butterfly <- feed_long %>%
  transmute(
    site,
    row_id,
    family = "butterfly",
    species = .data$`Species...22`,
    date    = as.Date(.data$`Date...23`),
    time_s  = as_num(.data$`Time...24`),
    mac = as_num(.data$`Mac...25`),
    tur = as_num(.data$`Tur...26`),
    hc  = as_num(.data$`Hc...27`),
    spo = as_num(.data$`Spo...28`),
    sed = as_num(.data$`Sed...29`),
    uk  = as_num(.data$`Uk...30`)
  )

# combine + compute rates
feed_long <- bind_rows(feed_parrot, feed_rabbit, feed_butterfly) %>%
  filter(!is.na(date), !is.na(time_s), time_s > 0) %>%
  mutate(
    bites_total = rowSums(across(c(mac, tur, hc, spo, sed, uk)), na.rm = TRUE),
    minutes_obs = time_s / 60,
    bites_min   = bites_total / minutes_obs
  )

# quick check
feed_long %>% count(site, family)
feed_long %>% summarise(n = n(), min_time = min(time_s), max_time = max(time_s))

### 03. FISH ####

# --- load + normalise headers ---

# rha fish transects
fish_raw <- read.csv("RAW_data/fish_abund_RAW.csv", header = TRUE)
names(fish_raw) <- tolower(names(fish_raw))

# --- main fish transects ---
fish <- fish_raw %>%
  select(!matches("^x(\\.|$)")) %>%
  select(where(~ !all(is.na(.x)))) %>%
  rename(
    survey_id  = surveycode,
    site_name  = site,
    date_chr   = date..mm.dd.yyyy.,
    transect_m = transect..m.,
    researcher = researchers
  ) %>%
  mutate(
    survey_id  = str_trim(survey_id),
    site_name  = str_trim(site_name),
    researcher = str_trim(researcher),
    
    site = case_when(
      str_detect(site_name, regex("^twins", ignore_case = TRUE)) ~ "TW",
      str_detect(site_name, regex("^shark\\s*island", ignore_case = TRUE)) ~ "SI",
      str_detect(site_name, regex("^tanote\\s*bay", ignore_case = TRUE)) ~ "TB",
      str_detect(site_name, regex("^green\\s*rock", ignore_case = TRUE)) ~ "GR",
      str_detect(site_name, regex("^red\\s*rock", ignore_case = TRUE)) ~ "RR",
      str_detect(site_name, regex("^aow\\s*leuk|^ao\\s*leuk", ignore_case = TRUE)) ~ "AL",
      TRUE ~ NA_character_
    ),
    
    date = mdy(date_chr),
    area_m2 = 1000 # updated
  ) %>%
  select(
    survey_id, site, site_name, date, researcher, area_m2,
    parrotfish, rabbitfish, butterflyfish
  ) %>%
  pivot_longer(
    c(parrotfish, rabbitfish, butterflyfish),
    names_to = "family",
    values_to = "fish_count"
  ) %>%
  mutate(
    family = recode(
      family,
      parrotfish = "parrot",
      rabbitfish = "rabbit",
      butterflyfish = "butterfly"
    ),
    fish_count = as.integer(fish_count),
    fish_density_ha = fish_count / (area_m2 / 10000)
  ) %>%
  filter(!is.na(site), !is.na(date), !is.na(fish_count))

# --- append RR + GR timed fish swims ---

timed_raw <- read.csv("RAW_data/TimedFishSurveys_Shallow_MASTER - data.csv", header = TRUE)
names(timed_raw) <- tolower(names(timed_raw))

timed_add <- timed_raw %>%
  mutate(
    site_clean = str_trim(site),
    site = case_when(
      str_detect(site_clean, regex("^red\\s*rock", ignore_case = TRUE)) ~ "RR",
      str_detect(site_clean, regex("^green\\s*rock", ignore_case = TRUE)) ~ "GR",
      TRUE ~ NA_character_
    ),
    date = mdy(date_mm.dd.yy), # fix
    researcher = str_trim(researcher)
  ) %>%
  filter(site %in% c("RR", "GR")) %>%
  transmute(
    survey_id  = paste0(site, "TIMED", format(date, "%Y%m%d"), "_", make.names(researcher)),
    site,
    site_name  = if_else(site == "RR", "Red Rock", "Green Rock"),
    date,
    researcher,
    area_m2    = 5000,  # updated
    parrotfish,
    rabbitfish,
    butterflyfish
  ) %>%
  pivot_longer(
    c(parrotfish, rabbitfish, butterflyfish),
    names_to = "family",
    values_to = "fish_count"
  ) %>%
  mutate(
    family = recode(
      family,
      parrotfish = "parrot",
      rabbitfish = "rabbit",
      butterflyfish = "butterfly"
    ),
    fish_count = as.integer(fish_count),
    fish_density_ha = fish_count / (area_m2 / 10000)
  ) %>%
  filter(!is.na(date), !is.na(fish_count))

fish <- bind_rows(fish, timed_add) %>%
  arrange(site, date, family, survey_id, researcher)

# sanity checks
fish %>% count(site, family)
fish %>% summarise(
  n = n(),
  n_surveys = n_distinct(survey_id),
  n_sites = n_distinct(site)
)
fish %>% group_by(site) %>%
  summarise(area_vals = paste(sort(unique(area_m2)), collapse = ", "))
### 04 substrate #### 
sub_path <- file.path(raw_dir, "Substrate_RAW.xlsx")
substrate_raw <- read_sites_xlsx(sub_path, sites_keep = c("AL","SI","TB","RR","GR","TW"))

str(substrate_raw)

substrate_transect <- substrate_raw %>%
  mutate(
    site = str_trim(site) %>% str_to_upper(),
    transect = str_trim(transect) %>% str_to_lower(),
    substrate = str_trim(substrate) %>% str_to_upper()
  ) %>%
  # If you have multiple hard-coral codes, map them to HC here
  mutate(
    substrate_grp = case_when(
      substrate %in% c("HC", "HCC", "HARD_CORAL") ~ "HC",
      substrate %in% c("MAC", "MA", "MACROALGAE") ~ "MAC",
      substrate %in% c("TUR", "TURF") ~ "TUR",
      TRUE ~ substrate
    )
  ) %>%
  group_by(site, transect) %>%
  mutate(
    n_points = n()
  ) %>%
  count(site, transect, substrate_grp, name = "n_sub") %>%
  group_by(site, transect) %>%
  mutate(
    prop = n_sub / sum(n_sub),
    pct  = 100 * prop
  ) %>%
  ungroup() %>%
  select(site, transect, substrate_grp, n_sub, prop, pct)

# pivot 
substrate_wide <- substrate_transect %>%
  filter(substrate_grp %in% c("MAC", "TUR", "HC")) %>%
  select(site, transect, substrate_grp, pct) %>%
  pivot_wider(
    names_from = substrate_grp,
    values_from = pct,
    values_fill = 0
  )

# sanity checks 
substrate_raw %>%
  count(site, transect, name = "n_points") %>%
  summarise(
    min_points = min(n_points),
    max_points = max(n_points),
    n_transects = n()
  )

### MERGING #### 
str(cots)
str(fish)
str(feed_long)
str(substrate_wide)

# scale everything to Ha 
cots <- cots %>%
  mutate(
    area_ha = area_m2 / 10000,
    cots_per_ha = cots_count / area_ha
  ) %>%
  select(-cots_density_ha) %>% 
  rename(cots_density_ha = cots_per_ha)
fish <- fish %>%
  mutate(
    area_ha = area_m2 / 10000,
    fish_per_ha = fish_count / area_ha
  ) %>%
  select(-fish_density_ha) %>%
  rename(fish_density_ha = fish_per_ha)
feed_long <- feed_long %>%
  mutate(
    feeding_rate = bites_min   # bites per minute, area-invariant
  )



# WHAT WE HAVE NOW 
  # | Dataset       | Native unit                        | Variation we want to preserve.            |
  # | ------------- | ---------------------------------- | ----------------------------------------- |
  # | **Feeding**   | site × date × family × observation | within-site, within-date feeding variance |
  # | **Fish**      | site × date × family × researcher  | replicate fish counts                     |
  # | **CoTS**      | site × survey_id                   | within-site CoTS variability              |
  # | **Substrate** | site × transect                    | fine-scale benthic structure              |
str(timed_raw)
str(cots)

