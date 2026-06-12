# 🪸 CoTS effects on grazing fish feeding behaviour  🐠

Using long-term crown-of-thorns seastar density, benthic substrate, fish abundance, and fish feeding data to test whether high CoTS densities are associated with ecological change on reefs around Koh Tao, Thailand.  
This repository supports a Global Reef manuscript examining links between **CoTS density**, **benthic condition**, and **grazer fish responses**.

---

## 🐠 Overview  
This project uses reef monitoring data from six sites around Koh Tao, Thailand, spanning a gradient of long-term CoTS density.
The analysis tests whether mean site-level CoTS density is associated with changes in hard coral cover, turf algae, macroalgae, grazer fish abundance, and substrate-specific feeding behaviour.

The main analysis uses continuous mean CoTS density as the predictor.
Additional sensitivity analyses test the local 15 CoTS ha⁻¹ outbreak threshold and the influence of Aow Leuk as a potentially distinct control site.

---

## 🎯 Objectives  
- Quantify site-level CoTS density and benthic substrate composition across six reef sites.  
- Test whether hard coral, turf algae, and macroalgae cover vary across the CoTS density gradient.  
- Test whether abundance of parrotfish, rabbitfish, and butterflyfish changes with CoTS density.  
- Test whether substrate-specific feeding behaviour shifts across the CoTS density gradient.  
- Compare continuous CoTS density models with categorical high vs low CoTS-density sensitivity models.  
- Assess whether results are robust to excluding Aow Leuk.  

---

## 📁 Repository structure   
```
data_raw/                       # raw field datasets
data_clean/                     # cleaned datasets produced by 01_CLEAN.R
docs/                           # notes, manuscript drafts, supporting documents

00_SETUP.R                      # packages, paths, palettes, helper functions
01_CLEAN.R                      # data cleaning and preprocessing
02_EXPLORE.R                    # exploratory checks following Zuur et al. (2010)
03_ANALYSIS_TABLES.R            # builds model-ready analysis tables
04_MODELS.R                     # final manuscript models and sensitivity models
05_RESULTS.R                    # manuscript tables, result values, predictions
06_PLOTS.R                      # publication figures
99_SCRATCH_ARCHIVE.R            # unused exploration / testing / notes

Analysis_YYYY.MM.DD/
  eda/                          # exploratory outputs
  fits/                         # saved model objects
  plots/                        # saved manuscript figures
  stats/                        # model summaries and diagnostics
  summaries/                    # saved summary objects
  tables/                       # model-ready tables and manuscript tables
```

Outputs are generated within each script and saved into the dated `Analysis_YYYY.MM.DD/` folder.

---

## 🧼 Data cleaning  
`01_CLEAN.R` cleans the raw datasets and writes cleaned CSV files to `data_clean/`.

Cleaned outputs include:

```
fish_bites_clean.csv
fish_abund_clean.csv
cots_indiv_clean.csv
cots_survey_clean.csv
substrate_point_clean.csv
substrate_transect_cover_clean.csv
substrate_transect_cover_wide_clean.csv
```

The cleaning script keeps the main datasets separate:

- CoTS individual observations  
- CoTS survey-level densities  
- Fish abundance by survey and family  
- Fish feeding observations  
- CPCe point-level substrate classifications  
- Transect-level substrate cover  

---

## 🔎 Exploratory analysis  
`02_EXPLORE.R` keeps the exploratory workflow visible and reproducible.

Checks include:

- missing values  
- impossible values  
- outliers  
- response distributions  
- zero inflation  
- overdispersion  
- substrate collinearity  
- temporal coverage  
- sampling effort by site  
- join feasibility across CoTS, substrate, fish abundance, and feeding datasets  

This script is exploratory, and based on Zuur et al. (2010)

---

## 📊 Model-ready tables  
`03_ANALYSIS_TABLES.R` builds the final analysis tables used by the models.

Main outputs include:

```
cots_site_context.csv
substrate_site_context.csv
benthos_mod.csv
fish_abund_mod.csv
fish_bites_mod.csv
feeding_pref_mod.csv
```

These tables join site-level CoTS context and substrate context onto the benthic, fish abundance, and feeding datasets.

---

## 🧪 Models  
`04_MODELS.R` contains the final manuscript models and retained sensitivity analyses.

Main model sets:

- Benthic substrate models  
- Fish abundance model  
- Substrate-specific feeding model  
- Secondary total bite-rate models  
- 15 CoTS ha⁻¹ threshold sensitivity models  
- Aow Leuk exclusion sensitivity models  

The main predictor is mean site-level CoTS density.

### Benthic models  
Separate Gaussian mixed-effects models are used for hard coral, macroalgae, and turf algae cover:

```r
substrate_cover ~ cots_mean_ha_sc + (1 | site_code)
```

### Fish abundance model  
Fish counts are modelled directly using a negative binomial mixed-effects model:

```r
fish_count ~ family * cots_mean_ha_sc +
  offset(log(area_ha)) +
  (1 | site_code) +
  (1 | date)
```

### Feeding preference model  
Substrate-specific feeding is modelled using a beta-binomial mixed-effects model:

```r
cbind(bites_sub, bites_other) ~ substrate * family * cots_mean_ha_sc +
  avail_pct_sc +
  (1 | obs_id)
```

Total bite-rate models are retained as secondary feeding intensity analyses, but the main feeding results focus on substrate-specific bite allocation.

---

## 📈 Results and figures  
`05_RESULTS.R` extracts manuscript-ready values, including:

- survey effort summaries  
- Table 1 site context values  
- model fixed effects  
- predicted fish densities  
- predicted feeding probabilities  
- sensitivity analysis summaries  

`06_PLOTS.R` generates the main manuscript figures:

- site-level benthic composition and CoTS density  
- CoTS density vs hard coral, turf algae, and macroalgae  
- fish abundance predictions across the CoTS gradient  
- substrate-specific feeding predictions across the CoTS gradient  

---

## ▶️ Running the workflow  
Run scripts in order:

```r
source("00_SETUP.R")
source("01_CLEAN.R")
source("02_EXPLORE.R")
source("03_ANALYSIS_TABLES.R")
source("04_MODELS.R")
source("05_RESULTS.R")
source("06_PLOTS.R")
```

Raw files should be stored in `data_raw/` and named using the `analysis_date` prefix set in `00_SETUP.R`.

Example:

```
2026.06.09_Fish_Bites.csv
2026.06.09_Fish_Abund.csv
2026.06.09_COTS_Abund.csv
2026.06.09_cpce_long.csv
```

---

## 📦 R packages  
Main packages used:

```r
tidyverse
ggplot2
glmmTMB
broom.mixed
janitor
DHARMa
purrr
```

Optional:

```r
patchwork
```

---

## 🚧 Status  
Active manuscript analysis for Global Reef’s CoTS, benthic condition, and grazer fish feeding project on Koh Tao reefs.

## Notes
-----

- Survey data are collected by Global Reef researchers based in Koh Tao, Thailand.
- Fieldwork and data processing are ongoing; results may be updated as more surveys are completed.

## License
-------

This project is private and not licensed for redistribution. For collaboration inquiries, please contact scarlett@global-reef.com.
---

**Affiliation:** [Global Reef](https://global-reef.com), Koh Tao, Thailand  

