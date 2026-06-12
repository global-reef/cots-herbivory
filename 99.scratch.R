### 99_SCRATCH_ARCHIVE.R ####
### Parking lot for abandoned or exploratory code.
### Anything needed to reproduce manuscript methods, results, figures, or tables has been moved into numbered scripts.

#### CHECK TOTAL BITE COUNTS ####

fish_bites_mod %>%
  summarise(
    n_obs = n(),
    total_bites = sum(bites_total, na.rm = TRUE),
    mean_bites = mean(bites_total, na.rm = TRUE),
    median_bites = median(bites_total, na.rm = TRUE),
    max_bites = max(bites_total, na.rm = TRUE),
    mean_bites_min = mean(bites_min, na.rm = TRUE),
    median_bites_min = median(bites_min, na.rm = TRUE),
    max_bites_min = max(bites_min, na.rm = TRUE)
  )

fish_bites_mod %>%
  arrange(desc(bites_total)) %>%
  select(site_code, date, family, species, minutes_obs, bites_total, bites_min,
         macro_algal_bites, turf_filamentous_algal_bites, live_coral_bites,
         dead_coral_with_algae, sponge_bites, sediment, unknown) %>%
  slice_head(n = 20)
