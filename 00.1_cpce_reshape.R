### cpce data reshaping #### 
library(tidyverse)
library(readxl)
library(janitor)

site_lookup <- tibble(
  sheet_name = c(
    "AOWLEUK",
    "TANOTE BAY",
    "TWINS WALL",
    "GREEN ROCK",
    "RED ROCK",
    "SHARK ISLAND"
  ),
  site_code = c(
    "AL",
    "TB",
    "TW",
    "GR",
    "RR",
    "SI"
  ),
  Site = c(
    "Aow Leuk",
    "Tanote Bay",
    "Twins",
    "Green Rock",
    "Red Rock",
    "Shark Island"
  )
)
clean_cpce_sheet <- function(file, sheet_name) {
  
  site_info <- site_lookup %>% 
    filter(sheet_name == !!sheet_name)
  
  if (nrow(site_info) != 1) {
    stop("Sheet name not found in site_lookup: ", sheet_name)
  }
  
  raw <- read_excel(file, sheet = sheet_name, col_names = FALSE)
  
  point_rows <- raw %>%
    mutate(
      point_check = suppressWarnings(as.integer(...1))
    ) %>%
    filter(point_check %in% 1:25)
  
  if (nrow(point_rows) %% 25 != 0) {
    warning("CPCE point rows are not divisible by 25 in sheet: ", sheet_name)
  }
  
  point_rows <- point_rows %>%
    mutate(
      transect_num = ((row_number() - 1) %/% 25) + 1,
      point_row = ((row_number() - 1) %% 25) + 1,
      Transect = paste0(site_info$site_code, transect_num)
    )
  
  point_rows %>%
    select(Transect, point_row, starts_with("...")) %>%
    pivot_longer(
      cols = starts_with("..."),
      names_to = "col",
      values_to = "value"
    ) %>%
    mutate(
      col_num = as.integer(str_remove(col, "^\\.\\.\\.")),
      Quadrat = ((col_num - 1) %/% 4) + 1,
      field = case_when(
        (col_num - 1) %% 4 == 0 ~ "CPCE Point",
        (col_num - 1) %% 4 == 1 ~ "Substrate",
        (col_num - 1) %% 4 == 2 ~ "GF",
        (col_num - 1) %% 4 == 3 ~ "Scar"
      )
    ) %>%
    filter(Quadrat <= 50, !is.na(field)) %>%
    select(Transect, point_row, Quadrat, field, value) %>%
    pivot_wider(
      id_cols = c(Transect, point_row, Quadrat),
      names_from = field,
      values_from = value
    ) %>%
    mutate(
      Site = site_info$Site,
      Quadrat = as.integer(Quadrat),
      `CPCE Point` = as.integer(`CPCE Point`)
    ) %>%
    select(Site, Transect, Quadrat, `CPCE Point`, Substrate, GF, Scar) %>%
    arrange(Transect, Quadrat, `CPCE Point`)
}

cpce_excel_file <- file.path(data_raw_dir, paste0(analysis_date, "_CPCE.xlsx"))

cpce_long <- excel_sheets(file) %>%
  map_dfr(~ clean_cpce_sheet(file, .x))

write_csv(cpce_long, file.path(data_raw_dir, paste0(analysis_date, "_cpce_long.csv")))
str(cpce_long)

# checks 
cpce_long %>%
  count(Site, Transect)

cpce_long %>%
  count(Site, Transect, Quadrat)

