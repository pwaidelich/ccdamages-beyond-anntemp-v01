# clean out the environment
rm(list = ls())

# load packages
library(ggpattern)       # for hatching in maps created with geom_sf()
library(rnaturalearth)   # for geodata to create map charts
library(tidyverse)       # for general data wrangling & plotting
library(ggpubr)          # for grid charts via ggarrange()
library(sf)              # for GIS analysis
library(readxl)          # for importing Excel files
library(rnaturalearthdata) # for low-resolution ADM0-level shapefiles
library(furrr)            # for parallel processing
library(feather)          # for importing .feather files
library(ggcorrplot)      # for correlation matrix plots
library(xtable)          # for exporting tables
library(dtplyr)          # to speed up data wrangling for large tibbles

# packages that must be installed but are not fully loaded (due to conflicting namespace with tidyverse): Hmisc
if(!"Hmisc" %in% rownames(installed.packages())) stop("Please install the Hmisc package before running this script")

# set global options to recalculate data objects - if this is FALSE, RDS objects from previous executions will be loaded
# NOTE: re-calculating data objects can take multiple hours of runtime even on a high-performance cluster
redo_data_objects <- F

# source utils functions
source(file.path("src", "utils", "analyse_impacts.R"))
source(file.path("src", "utils", "project_impacts.R"))
source(file.path("src", "utils", "extract_climate_indices.R"))

# # map the directory where climate projection files are stored (each vector has N_countries x N_modelscenariopairs length)
dir_imp <- "/net/cfc/landclim1/pwaidelich"
# ADM0-level projections using point estimates (= MC sample size of 1)
files_adm0_point <- file.path(dir_imp, "adm0_new", "pointestimates") %>% list.files(full.names = T) %>% str_subset("_mc1_.+\\.feather$")
files_adm0_bc_point <- file.path(dir_imp, "adm0_bc_new", "pointestimates") %>% list.files(full.names = T) %>% str_subset("_mc1_.+\\.feather$")
# main results for all six climate indicators and for status quo approach
files_adm0_mc <- file.path(dir_imp, "adm0_new", "montecarlo") %>% list.files(full.names = T) %>% str_subset("_mc1000_.+\\.feather$")
files_adm0_mc_T_Pt_Pt2zero <- file.path(dir_imp, "adm0_new", "T_Pt_Pt2zero", "montecarlo") %>% list.files(full.names = T) %>% str_subset("_mc1000_.+\\.feather$")
# main results using bias-corrected CMIP6 data - ramp this up to 1,000 Monte Carlo runs as well for final results XX
files_adm0_bc_mc <- file.path(dir_imp, "adm0_bc_new", "montecarlo") %>% list.files(full.names = T) %>% str_subset("_mc1000_.+\\.feather$")
files_adm0_bc_mc_T_Pt_Pt2zero <- file.path(dir_imp, "adm0_bc_new", "T_Pt_Pt2zero", "montecarlo") %>% list.files(full.names = T) %>% str_subset("_mc1000_.+\\.feather$")

# throw an error if these vectors do not have the same length (= different no. of projection files)
if((map_int(list(files_adm0_point, files_adm0_bc_point, files_adm0_mc,
            files_adm0_mc_T_Pt_Pt2zero, files_adm0_bc_mc, files_adm0_bc_mc_T_Pt_Pt2zero), length) %>% 
   unique() %>% length()) > 1) stop("The number of filenames in the different files_... object are not identical")

# create a small-scale polygon of countries from rnaturalearthdata
world <- ne_countries(scale = "medium", returnclass = "sf") %>% select(sovereignt, type, admin, GID_0 = "iso_a3", abbrev, 
                                                                       pop_est, pop_year, lastcensus, gdp_md_est, gdp_year, 
                                                                       economy, income_grp, continent, region_un, subregion, region_wb, geometry) %>%
  # aggregate up the economic region (currently not used)
  mutate(economic_groups = case_when(str_detect(economy, "Developed region") ~ "Developed region",
                                     str_detect(economy, "Emerging region") ~ "Emerging region",
                                     str_detect(economy, "Developing region") ~ "Developing region",
                                     str_detect(economy, "Least developed region") ~ "LDC",
                                     TRUE ~ "not classified - please inspect"),
         income_grp_groups = case_when(str_detect(income_grp, "High income") ~ "High income",
                                       str_detect(income_grp, "middle income") ~ "Middle income",
                                       str_detect(income_grp, "Low income") ~ "Low income"))

# create a vector of all sovereign GID_0 territories
sovereign_countries_gid0 <- world %>% filter(type %in% c("Country", "Sovereign country") & admin == sovereignt & !is.na(GID_0)) %>% pull(GID_0)

# inspect tiny countries omitted in the shapefile pulled above
world_tinycountries <- ne_countries(type = "tiny_countries", returnclass = "sf", scale = "medium") %>% rename(GID_0 = "iso_a3")
world_tinycountries %>% filter(!GID_0 %in% world$GID_0)

# add Tuvalu to the list of sovereign countries and the world shapefile (too small to feature on the map)
sovereign_countries_gid0 <- c(sovereign_countries_gid0, "TUV") %>% unique()
world <- bind_rows(world, world_tinycountries %>% filter(GID_0 == "TUV"))
rm(world_tinycountries)

# load the GADM shapefile at the region1 and region0 (= country) level
sf_gadm0_simple <- readRDS(file.path("data", "sf_gadm0_simple.rds"))

# World Bank data on income-based country classifications
df_worldbank <- read_excel(file.path("data", "input", "230108 World Bank Country Classification.xlsx"),
                           sheet = "List of economies") %>%
  # discard all aggregates (e.g. EU) that have an NA in the region column
  filter(!is.na(Region))

# load the global warming levels and the SSP data on GDP created in 'src/00 Prepare Analysis.R'
df_gwl <- readRDS(file.path("data", "global_warming_levels.rds"))
df_ssp <- readRDS(file.path("data", "df_ssp.rds")) %>% ungroup()

# we use SSP GDP weighting for ADM0-to-global aggregation for all warming levels, which could be distorted if the SSP data does not cover some early years
# therefore, we test if the lowest warming level for which we calculate global results (+1C) is fully covered by data in df_ssp
if(min(df_gwl$year[df_gwl$warming_level == 1]) < min(df_ssp$year)) stop("SSP database does not cover some years for +1C warming level")

# calculate ADM0-level GDP weights
df_ssp <- df_ssp %>%
  group_by(ssp, year) %>%
  mutate(adm0_weight = gdp_annual/sum(gdp_annual)) %>%
  ungroup()
  
# throw an error if the lowest warming level window analysed (+1C) for some model-scenario pair dates further back than the years in df_ssp
if(min(df_gwl$year[df_gwl$warming_level == 1]) < min(df_ssp$year)) stop("For some model-scenario pair, +1C GWL includes years not covered by df_ssp. This will result in wrong global aggregate figures since these years are omitted from the GWL")

# create a data frame hosting all files in files_adm0_mc and some info on the underlying country and model-scenario
create_mcfiles_tibble <- function(filename_vector = NULL, folder_name = "montecarlo") {
  tibble(filename = filename_vector) %>%
    mutate(GID_0 = str_extract(filename, "[A-Z][A-Z][A-Z](?=\\.feather$)"),
           has_sspdata = GID_0 %in% unique(df_ssp$GID_0),
           modelscen_string = str_extract(filename, paste0("(?<=", folder_name, "\\/).+(?=_mc[:digit:])"))) %>%
    left_join(sf_gadm0_simple %>% as_tibble() %>% select(GID_0, NAME_0), by = "GID_0")
}

# create tibbles hosting all the Monte Carlo GDP projection files
df_mcfiles <- create_mcfiles_tibble(files_adm0_mc)
df_mcfiles_bc <- create_mcfiles_tibble(files_adm0_bc_mc)
df_mcfiles_T_Pt_Pt2zero <- create_mcfiles_tibble(files_adm0_mc_T_Pt_Pt2zero)
df_mcfiles_bc_T_Pt_Pt2zero <- create_mcfiles_tibble(files_adm0_bc_mc_T_Pt_Pt2zero)
df_pointfiles <- create_mcfiles_tibble(files_adm0_point, folder_name = "pointestimates") 
df_pointfiles_bc <- create_mcfiles_tibble(files_adm0_bc_point, folder_name = "pointestimates")

# print out all ADM0 regions and only sovereign ADM0 countries that are not in the SSP database
df_mcfiles %>% filter(!has_sspdata) %>% pull(NAME_0) %>% unique() %>% sort()
tibble(Countries = c(paste(df_mcfiles %>% filter(!has_sspdata & GID_0 %in% sovereign_countries_gid0) %>% pull(NAME_0) %>% unique() %>% sort(), collapse = ", "),
                     paste(df_mcfiles %>% filter(!GID_0 %in% sovereign_countries_gid0) %>% pull(NAME_0) %>% unique() %>% sort(), collapse = ", ")),
       reason = c("Not included in the SSP database",
                                 "Listed as a non-sovereign country in Natural Earth shapefiles"),
       implication = c("Included in country-level results but omitted from (GDP-weighted) global aggregate results",
                                 "Omitted from all projections")) %>%
  setNames(names(.) %>% str_replace("reason", "Reason for omission") %>% str_replace("implication", "Implication of omission")) %>%
  xtable(label = "tab:omitted_adm0",
                 caption = "ADM0-level territories omitted from (some) analyses") %>%
  print(file = file.path("..", "sharepoint", "Tables", paste0(Sys.Date(), " omitted_adm0.tex")),
        size="\\small")

## GCM model weights
# calculate the weights for each model-scenario-year-warminglevel obs as # of obs for respective GCM
# based on one example ADM0-level territory
df_gcm_weights <- df_mcfiles$filename[df_mcfiles$GID_0 == "UKR"] %>% map_dfr(read_feather) %>%
  left_join(df_gwl %>% select(-ensemble), by = c("model", "scenario", "year")) %>%
  filter(!is.na(warming_level), monte_carlo_draw == 1) %>%
  group_by(model, warming_level) %>%
  summarise(gcm_weight = 1/n()) %>% ungroup()

# repeat for bias-corrected data and ensure that this yields the same results (= same model-scenario pairs used for raw & bias-corrected data)
df_gcm_weights_bc <- df_mcfiles_bc$filename[df_mcfiles_bc$GID_0 == "AFG"] %>% map_dfr(read_feather) %>%
  left_join(df_gwl %>% select(-ensemble), by = c("model", "scenario", "year")) %>%
  filter(!is.na(warming_level), monte_carlo_draw == 1) %>%
  group_by(model, warming_level) %>%
  summarise(gcm_weight = 1/n()) %>% ungroup()

# if weights differ, throw an error; otherwise, delete df_gcm_weights_bc again
if(identical(df_gcm_weights_bc, df_gcm_weights)) {
  rm(df_gcm_weights_bc)
} else {
  stop("df_gcm_weights and df_gcm_weights_bc are not identical. Please inspect")
}

# plot the # of obs by CMIP6 for each global warming level as a bar chart - colour models for which +0.61 is cut off due to 1979 threshold
df_gcm_weights %>%
  filter(warming_level %in% c(1, 1.5, 2, 3, 4)) %>%
  mutate(model_gwl = paste0(model, " - ", warming_level),
         warming_level = factor(paste0("GWL: +", warming_level, "C"),
                                levels = c("GWL: +0.84C",
                                           "GWL: +1C",
                                           "GWL: +1.5C",
                                           "GWL: +2C",
                                           "GWL: +3C",
                                           "GWL: +4C"))) %>%
  ggplot(aes(reorder(model, -gcm_weight), 1/gcm_weight)) + geom_col() +
  #geom_text(aes(label = 1/gcm_weight, y = -3)) +
  labs(y = "# of model-scenario-years", x = "CMIP6 model") +
  theme_classic() +
  coord_flip() + 
  facet_wrap(~ warming_level, nrow = 1) +
  guides(fill = "none")
ggsave(file.path("..", "sharepoint", "Figures", paste0(Sys.Date(), " no_obs_by_CMIP6_gwl.png")), width = 10, height = 6)

# initiate parallel processing and adjust the limit of how large globals can be passed to workers
plan(multisession, gc = T)
options(future.globals.maxSize= 1500*1024^2)

# set the options for ggpattern to avoid a bug, see here: https://stackoverflow.com/questions/70394227/pattern-in-ggpattern-doesnt-follow-polygon-borders-after-update
options(ggpattern_use_R4.1_features = FALSE)
sf_use_s2(FALSE)


################################################################################
############################## ADM1 IMPACTS ####################################
################################################################################

# NOTE: all charts for ADM1-level results are based only on point estimates since
# they only serve for the SI. If Monte Carlo results were available at ADM1-level,
# this would increase storage requirements for outputs by > 1 order of magnitude

# read in ADM1-level results created with point estimates for damage function parameters & merge in warming levels
df_adm1 <- read_feather(file.path("/net", "cfc", "landclim1", "pwaidelich", "adm1_new", "pointestimates", "df_adm1.feather")) %>% 
  lazy_dt() %>%
  left_join(df_gwl %>% select(-ensemble), by = c("model", "scenario", "year")) %>%
  filter(!is.na(warming_level)) %>%
  as_tibble()

# repeat for bias-corrected data
df_adm1_bc <- read_feather(file.path("/net", "cfc", "landclim1", "pwaidelich", "adm1_bc_new", "pointestimates", "df_adm1_bc.feather")) %>% 
  lazy_dt() %>%
  left_join(df_gwl %>% select(-ensemble), by = c("model", "scenario", "year")) %>%
  filter(!is.na(warming_level)) %>%
  as_tibble()

# load the data from Kotz et al 2022
df_kotz2022 <- read.csv(file.path("..", "sharepoint", "Data", "Kotz et al - data & code", "Data_code_for_Rainfall_economic_production_Kotz_2021", "zenodo", "secondary_data", "PT_master.csv")) %>%
  # create the same ADM1-level identifier
  mutate(GID_1 = paste0(iso, ".", id_1, "_1"))

# create a subset of df_adm1 that only features territories & years included in Kotz et al 2022 data over 1979-2019
# NOTE: to do this, we map the ERA-5 1979-2019 years to the corresponding years in the 0.84C global warming level (= baseline) for each model-scenario
df_adm1_in_kotz <- df_adm1 %>%
  lazy_dt() %>%
  filter(is_baseline) %>%
  # for each model, territory and year, order by SSP and extract the first to ensure we have only one scenario
  arrange(model, GID_1, year, desc(scenario)) %>% group_by(model, GID_1, year) %>% slice_head(n=1) %>% ungroup() %>%
  # make years equal to 1979-2019 (ERA5 baseline) - we can use min(year) to determine the baseline start since we filtered to baseline
  # and data is grouped by model-scenario
  group_by(model, GID_1, scenario) %>%
  mutate(baseline_startyear = min(year),
         year_kotz = year + (1979 - baseline_startyear)) %>% ungroup() %>%
  as_tibble() %>%
  # now merge in the Kotz et al 2022 data
  left_join(df_kotz2022 %>% select(GID_1, year, Tmean, Tstd, Pt, wet_days_1, vwet_days1_am_99p9, Wn) %>% mutate(is_in_kotz2022 = T),
            by = c("GID_1", "year_kotz" = "year"), suffix = c("", "_kotz")) %>%
  filter(is_in_kotz2022) %>%
  select(model, scenario, GID_1, year, year_kotz, baseline_startyear,
         Tmean, Tstd, Pt, wet_days_1, vwet_days1_am_99p9, Wn,
         Tmean_kotz, Tstd_kotz, Pt_kotz, wet_days_1_kotz, vwet_days1_am_99p9_kotz, Wn_kotz)

# repeat for bias-corrected data
df_adm1_bc_in_kotz <- df_adm1_bc %>%
  lazy_dt() %>%
  filter(is_baseline) %>%
  # for each model, territory and year, order by SSP and extract the first to ensure we have only one scenario
  arrange(model, GID_1, year, desc(scenario)) %>% group_by(model, GID_1, year) %>% slice_head(n=1) %>% ungroup() %>%
  # make years equal to 1979-2019 (ERA5 baseline) - we can use min(year) to determine the baseline start since we filtered to baseline
  # and data is grouped by model-scenario
  group_by(model, GID_1, scenario) %>%
  mutate(baseline_startyear = min(year),
         year_kotz = year + (1979 - baseline_startyear)) %>% ungroup() %>%
  as_tibble() %>%
  # now merge in the Kotz et al 2022 data
  left_join(df_kotz2022 %>% select(GID_1, year, Tmean, Tstd, Pt, wet_days_1, vwet_days1_am_99p9, Wn) %>% mutate(is_in_kotz2022 = T),
            by = c("GID_1", "year_kotz" = "year"), suffix = c("", "_kotz")) %>%
  filter(is_in_kotz2022) %>%
  select(model, scenario, GID_1, year, year_kotz, baseline_startyear,
         Tmean, Tstd, Pt, wet_days_1, vwet_days1_am_99p9, Wn,
         Tmean_kotz, Tstd_kotz, Pt_kotz, wet_days_1_kotz, vwet_days1_am_99p9_kotz, Wn_kotz)

# correlation of ADM1-level climate indicators in CMIP6 projections (using only obs included in Kotz et al)
cor_climate_adm1 <- df_adm1_in_kotz %>% select(Tmean, Tstd, Pt, Wn, wet_days_1, vwet_days1_am_99p9) %>%
  setNames(names(.) %>% rewrite_climate_label()) %>%
  cor(use = "complete.obs")

# repeat for the Kotz et al 2022 data (using only the 1979-2019 baseline period)
cor_climate_kotz2022 <- df_kotz2022 %>% filter(year %in% 1979:2019) %>% select(Tmean, Tstd, Pt, Wn, wet_days_1, vwet_days1_am_99p9) %>%
  setNames(names(.) %>% rewrite_climate_label()) %>%
  cor(use = "complete.obs")

# plot the two corrplots next to each for comparisons and save out
ggarrange(
  ggcorrplot(cor_climate_kotz2022, outline.col = "white", type = "lower", hc.order = T,
                         ggtheme = theme_classic, lab = T, pch.col = "grey", pch.cex =  10, digits = 2) +
    labs(subtitle = "ERA5-based sample from Kotz et al., 2022"),
  ggcorrplot(cor_climate_adm1, outline.col = "white", type = "lower", hc.order = T,
                       ggtheme = theme_classic, lab = T, pch.col = "grey", pch.cex =  10, digits = 2) +
    labs(subtitle = "CMIP6 sample for +0.84\u00B0C baseline"),
  nrow = 2, ncol = 1
)
ggsave(file.path("..", "sharepoint", "Figures", "adm1",
                 paste0(Sys.Date(), " ADM1correlation_era5_vs_cmip6.png")),
       width = 10, height = 10)

# correlation between models and Kotz
df_adm1_in_kotz %>% group_by(model) %>%
  summarise(Tmean = cor(Tmean, Tmean_kotz),
            Pt = cor(Pt, Pt_kotz),
            Tstd = cor(Tstd, Tstd_kotz),
            Wn = cor(Wn, Wn_kotz),
            wet_days_1 = cor(wet_days_1, wet_days_1_kotz),
            vwet_days1_am_99p9 = cor(vwet_days1_am_99p9, vwet_days1_am_99p9_kotz)
            ) %>%
 pivot_longer(cols = !starts_with("model"), names_to = "variable", values_to = "cor") %>%
  mutate(variable = rewrite_climate_label(variable)) %>%
  ggplot(aes(model, cor)) + geom_col() +
  labs(y = "Correlation coefficient between\n CMIP6 baseline (+0.84\u00B0C) & ERA-5 1979-2019",
       x = NULL) +
  facet_wrap(~variable) +
  coord_flip()
ggsave(file.path("..", "sharepoint", "Figures", "adm1",
                 paste0(Sys.Date(), " ADM1correlation_cmip6_with_era5.png")),
       width = 7, height = 8)

# baseline summary stats for Kotz et al (using 1979-2019 data only)
df_adm1_baseline_sumstat <- bind_rows(bind_cols(tibble(Statistic = c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."),
                 Model = "Kotz et al., 2022"),
          df_kotz2022 %>% filter(year %in% 1979:2019) %>% 
            select(Tmean, Tstd, Pt, wet_days_1, vwet_days1_am_99p9, Wn) %>%
            map_dfc(~ summary(.x) %>% as.numeric())),
  # add summary stats calculations for all models during the baseline period
  map_dfr(unique(df_adm1_in_kotz$model), .f = function(x) {
    bind_cols(tibble(Statistic = c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."),
                     Model = x),
              map_dfc(df_adm1_in_kotz %>% filter(model == x) %>%
                        select(Tmean, Tstd, Pt, wet_days_1, vwet_days1_am_99p9, Wn),
                      ~ summary(.x) %>% as.numeric()))
    })
)

# repeat for bias-corrected
df_adm1_bc_baseline_sumstat <- bind_rows(bind_cols(tibble(Statistic = c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."),
                                                       Model = "Kotz et al., 2022"),
                                                df_kotz2022 %>% filter(year %in% 1979:2019) %>% 
                                                  select(Tmean, Tstd, Pt, wet_days_1, vwet_days1_am_99p9, Wn) %>%
                                                  map_dfc(~ summary(.x) %>% as.numeric())),
                                      # add summary stats calculations for all models during the baseline period
                                      map_dfr(unique(df_adm1_bc_in_kotz$model), .f = function(x) {
                                        bind_cols(tibble(Statistic = c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."),
                                                         Model = x),
                                                  map_dfc(df_adm1_bc_in_kotz %>% filter(model == x) %>%
                                                            select(Tmean, Tstd, Pt, wet_days_1, vwet_days1_am_99p9, Wn),
                                                          ~ summary(.x) %>% as.numeric()))
                                      })
)

# plot distribution summary stats for each variable over the baseline period as a boxplot
# (whiskers for min & max)
df_adm1_baseline_sumstat %>% pivot_longer(cols = c(Tmean, Tstd, Pt, wet_days_1, vwet_days1_am_99p9, Wn),
                                          names_to = "variable") %>%
  pivot_wider(names_from = "Statistic", values_from = "value") %>%
  mutate(variable = rewrite_climate_label(variable)) %>%
  ggplot(aes(reorder(Model, Model == "Kotz et al., 2022"))) +
  geom_boxplot(aes(ymin = `Min.`,
                   ymax = `Max.`,
                   lower = `1st Qu.`,
                   upper = `3rd Qu.`,
                   middle = Median,
                   fill = Model == "Kotz et al., 2022"),
               stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  labs(x = NULL) +
  coord_flip() +
  guides(fill = "none")
ggsave(file.path("..", "sharepoint", "Figures", "adm1",
                 paste0(Sys.Date(), " ADM1distribution_cmip6_vs_kotz.png")),
       width = 11, height = 8)

# repeat for bias-corrected data
df_adm1_bc_baseline_sumstat %>% pivot_longer(cols = c(Tmean, Tstd, Pt, wet_days_1, vwet_days1_am_99p9, Wn),
                                          names_to = "variable") %>%
  pivot_wider(names_from = "Statistic", values_from = "value") %>%
  mutate(variable = rewrite_climate_label(variable)) %>%
  ggplot(aes(reorder(Model, Model == "Kotz et al., 2022"))) +
  geom_boxplot(aes(ymin = `Min.`,
                   ymax = `Max.`,
                   lower = `1st Qu.`,
                   upper = `3rd Qu.`,
                   middle = Median,
                   fill = Model == "Kotz et al., 2022"),
               stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  labs(x = NULL) +
  coord_flip() +
  guides(fill = "none")
ggsave(file.path("..", "sharepoint", "Figures", "adm1",
                 paste0(Sys.Date(), " ADM1distribution_cmip6_vs_kotz_bc.png")),
       width = 11, height = 8)

# clean up
rm(cor_climate_adm1, cor_climate_kotz2022, df_adm1, df_adm1_bc, df_kotz2022,
   df_adm1_baseline_sumstat, df_adm1_bc_baseline_sumstat, df_adm1_in_kotz, df_adm1_bc_in_kotz)
gc()


################################################################################
############################## GLOBAL IMPACTS ##################################
################################################################################

# load & aggregate the distribution
import_global_distr <- function(file_data = NULL) {
  unique(file_data$modelscen_string) %>%
    future_map_dfr( ~ file_data %>% filter(modelscen_string == .x & has_sspdata) %>% pull(filename) %>% map_dfr(read_feather) %>%
                      # extract the SSP number from the scenario (RCP-SSP) column
                      mutate(ssp = paste0("SSP", str_extract(scenario, "(?<=^ssp)[:digit:]"))) %>%
                      # merge in GDP weights for the respective SSP
                      left_join(df_ssp %>% select(ssp, GID_0, year, adm0_weight), by = c("ssp", "GID_0", "year")) %>%
                      # exclude ADM0-level territories for which there is no GDP data in SSP Database and hence no weight
                      filter(!is.na(adm0_weight)) %>%
                      # group by model-scenario-year as well as MC run and calculate weighted averages (NOTE: data is in %GDP, not in log-scale)
                      group_by(model, scenario, year, monte_carlo_draw) %>%
                      summarise_at(vars(starts_with("imp")), ~ Hmisc::wtd.mean(.x, weights = adm0_weight, normwt = T)),
                    .progress = T, .options = furrr_options(seed = NULL)) %>%
    ungroup()
}

# import the global results or load the RDS files, depending on the 'redo_data_objects' option
if(redo_data_objects) {
  # read in the full distributions from the relevant specifications
  df_global_pointdistr <- import_global_distr(df_pointfiles)
  df_global_bc_pointdistr <- import_global_distr(df_pointfiles_bc)
  df_global_fulldistr <- import_global_distr(df_mcfiles)
  df_global_bc_fulldistr <- import_global_distr(df_mcfiles_bc)
  df_global_fulldistr_T_Pt_Pt2zero <- import_global_distr(df_mcfiles_T_Pt_Pt2zero)
  df_global_bc_fulldistr_T_Pt_Pt2zero <- import_global_distr(df_mcfiles_bc_T_Pt_Pt2zero)
  
  # save out
  saveRDS(df_global_pointdistr, file.path("data", "df_global_pointdistr.rds"))
  saveRDS(df_global_bc_pointdistr, file.path("data", "df_global_bc_pointdistr.rds"))
  saveRDS(df_global_fulldistr, file.path("data", "df_global_fulldistr.rds"))
  saveRDS(df_global_bc_fulldistr, file.path("data", "df_global_bc_fulldistr.rds"))
  saveRDS(df_global_fulldistr_T_Pt_Pt2zero, file.path("data", "df_global_fulldistr_T_Pt_Pt2zero.rds"))
  saveRDS(df_global_bc_fulldistr_T_Pt_Pt2zero, file.path("data", "df_global_bc_fulldistr_T_Pt_Pt2zero.rds"))
} else {
  df_global_pointdistr <- readRDS(file.path("data", "df_global_pointdistr.rds"))
  df_global_bc_pointdistr <- readRDS(file.path("data", "df_global_bc_pointdistr.rds"))
  df_global_fulldistr <- readRDS(file.path("data", "df_global_fulldistr.rds"))
  df_global_bc_fulldistr <- readRDS(file.path("data", "df_global_bc_fulldistr.rds"))
  df_global_fulldistr_T_Pt_Pt2zero <- readRDS(file.path("data", "df_global_fulldistr_T_Pt_Pt2zero.rds"))
  df_global_bc_fulldistr_T_Pt_Pt2zero <- readRDS(file.path("data", "df_global_bc_fulldistr_T_Pt_Pt2zero.rds"))
}

# write a function to aggregate up (weighting all summary stats by # of obs per CMIP6 model for each warming level)
aggregate_global_distr <- function(data = NULL, gcm_weights_used = NULL) {
  data %>%
    lazy_dt() %>%
    # merge in GWL and then the CMIP6 model weights
    left_join(df_gwl %>% filter(warming_level %in% c(1, 1.5, 2, 3, 4)) %>%
                select(-ensemble), by = c("model", "scenario", "year")) %>%
    filter(!is.na(warming_level)) %>%
    left_join(gcm_weights_used, by = c("model", "warming_level")) %>%
    group_by(warming_level) %>%
    # calculate all relevant summary stats for all impact channels using CMIP6 model weights
    summarise_at(vars(starts_with("imp")), .funs = list(
      mean = ~ Hmisc::wtd.mean(.x, weights = gcm_weight, normwt = T),
      var = ~ Hmisc::wtd.var(.x, weights = gcm_weight, normwt = T),
      min = min,
      max = max,
      perc05 = ~ Hmisc::wtd.quantile(.x, probs = 0.05, weights = gcm_weight, normwt = T),
      perc10 = ~ Hmisc::wtd.quantile(.x, probs = 0.1, weights = gcm_weight, normwt = T),
      median = ~ Hmisc::wtd.quantile(.x, probs = 0.5, weights = gcm_weight, normwt = T),
      perc90 = ~ Hmisc::wtd.quantile(.x, probs = 0.9, weights = gcm_weight, normwt = T),
      perc95 = ~ Hmisc::wtd.quantile(.x, probs = 0.95, weights = gcm_weight, normwt = T),
      agreement_meansign = ~ if_else(Hmisc::wtd.mean(.x, weights = gcm_weight, normwt = T) < 0,
                                     Hmisc::wtd.mean(.x < 0, weights = gcm_weight, normwt = T),
                                     Hmisc::wtd.mean(.x >= 0, weights = gcm_weight, normwt = T)),
      signalnoise = ~ Hmisc::wtd.mean(.x, weights = gcm_weight, normwt = T)/(Hmisc::wtd.var(.x, weights = gcm_weight, normwt = T))^0.5
    )) %>% ungroup() %>% as_tibble()
}

# aggregate up
if(redo_data_objects) {
  df_global_fulldistr_aggr <- aggregate_global_distr(df_global_fulldistr, gcm_weights_used = df_gcm_weights)
  df_global_fulldistr_T_Pt_Pt2zero_aggr <- aggregate_global_distr(df_global_fulldistr_T_Pt_Pt2zero, gcm_weights_used = df_gcm_weights)
  df_global_bc_fulldistr_aggr <- aggregate_global_distr(df_global_bc_fulldistr, gcm_weights_used = df_gcm_weights)
  df_global_bc_fulldistr_T_Pt_Pt2zero_aggr <- aggregate_global_distr(df_global_bc_fulldistr_T_Pt_Pt2zero, gcm_weights_used = df_gcm_weights)
  
  saveRDS(df_global_fulldistr_aggr, file.path("data", "df_global_fulldistr_aggr.rds"))
  saveRDS(df_global_fulldistr_T_Pt_Pt2zero_aggr, file.path("data", "df_global_fulldistr_T_Pt_Pt2zero_aggr.rds"))
  saveRDS(df_global_bc_fulldistr_aggr, file.path("data", "df_global_bc_fulldistr_aggr.rds"))
  saveRDS(df_global_bc_fulldistr_T_Pt_Pt2zero_aggr, file.path("data", "df_global_bc_fulldistr_T_Pt_Pt2zero_aggr.rds"))
} else {
  df_global_fulldistr_aggr <- readRDS(file.path("data", "df_global_fulldistr_aggr.rds"))
  df_global_fulldistr_T_Pt_Pt2zero_aggr <- readRDS(file.path("data", "df_global_fulldistr_T_Pt_Pt2zero_aggr.rds"))
  df_global_bc_fulldistr_aggr <- readRDS(file.path("data", "df_global_bc_fulldistr_aggr.rds"))
  df_global_bc_fulldistr_T_Pt_Pt2zero_aggr <- readRDS(file.path("data", "df_global_bc_fulldistr_T_Pt_Pt2zero_aggr.rds"))
}

# print out means as well as upper/lower deciles for all impact channels & global warming levels
df_global_fulldistr_aggr %>%
  mutate(imp_total = paste0(100*round(imp_total_mean, 4), "% (",
                            100*round(imp_total_perc10, 4), " to ",
                            100*round(imp_total_perc90, 4), "%)"),
         imp_temp_diff = paste0(100*round(imp_temp_diff_mean, 4), "% (",
                                100*round(imp_temp_diff_perc10, 4), " to ",
                                100*round(imp_temp_diff_perc90, 4), "%)"),
         imp_Pt_diff = paste0(100*round(imp_Pt_diff_mean, 4), "% (",
                              100*round(imp_Pt_diff_perc10, 4), " to ",
                              100*round(imp_Pt_diff_perc90, 4), "%)"),
         imp_Tstd_diff = paste0(100*round(imp_Tstd_diff_mean, 4), "% (",
                                100*round(imp_Tstd_diff_perc10, 4), " to ",
                                100*round(imp_Tstd_diff_perc90, 4), "%)"),
         imp_Wn_diff = paste0(100*round(imp_Wn_diff_mean, 4), "% (",
                              100*round(imp_Wn_diff_perc10, 4), " to ",
                              100*round(imp_Wn_diff_perc90, 4), "%)"),
         imp_wet_days_1_diff = paste0(100*round(imp_wet_days_1_diff_mean, 4), "% (",
                                      100*round(imp_wet_days_1_diff_perc10, 4), " to ",
                                      100*round(imp_wet_days_1_diff_perc90, 4), "%)"),
         imp_vwet_days1_am_99p9_diff = paste0(100*round(imp_vwet_days1_am_99p9_diff_mean, 4), "% (",
                                              100*round(imp_vwet_days1_am_99p9_diff_perc10, 4), " to ",
                                              100*round(imp_vwet_days1_am_99p9_diff_perc90, 4), "%)"),
         # write the warming level as "+3C" instead of "3"
         warming_level = paste0("+", warming_level, "\u00B0C")) %>%
  # reformat and write out as a .tex table
  select(warming_level, imp_total, ends_with("diff")) %>%
  pivot_longer(-warming_level) %>% pivot_wider(names_from=warming_level, values_from=value) %>%
  mutate(name = name %>% rewrite_impact_label()) %>%
  setNames(names(.) %>% str_replace("^name$", "Impact channel")) %>%
  xtable(label = "tab:si_globalresults",
        caption = "Global mean GDP impacts by warming level and climate indicator (upper and lower decile in parentheses)") %>%
  print(file = file.path("..", "sharepoint", "Tables", paste0(Sys.Date(), " si_globalresults_overview.tex")),
        size="\\tiny")

# repeat for bias-corrected data
# print out means as well as upper/lower deciles for all impact channels & global warming levels
df_global_bc_fulldistr_aggr %>%
  mutate(imp_total = paste0(100*round(imp_total_mean, 4), "% (",
                            100*round(imp_total_perc10, 4), " to ",
                            100*round(imp_total_perc90, 4), "%)"),
         imp_temp_diff = paste0(100*round(imp_temp_diff_mean, 4), "% (",
                                100*round(imp_temp_diff_perc10, 4), " to ",
                                100*round(imp_temp_diff_perc90, 4), "%)"),
         imp_Pt_diff = paste0(100*round(imp_Pt_diff_mean, 4), "% (",
                              100*round(imp_Pt_diff_perc10, 4), " to ",
                              100*round(imp_Pt_diff_perc90, 4), "%)"),
         imp_Tstd_diff = paste0(100*round(imp_Tstd_diff_mean, 4), "% (",
                                100*round(imp_Tstd_diff_perc10, 4), " to ",
                                100*round(imp_Tstd_diff_perc90, 4), "%)"),
         imp_Wn_diff = paste0(100*round(imp_Wn_diff_mean, 4), "% (",
                              100*round(imp_Wn_diff_perc10, 4), " to ",
                              100*round(imp_Wn_diff_perc90, 4), "%)"),
         imp_wet_days_1_diff = paste0(100*round(imp_wet_days_1_diff_mean, 4), "% (",
                                      100*round(imp_wet_days_1_diff_perc10, 4), " to ",
                                      100*round(imp_wet_days_1_diff_perc90, 4), "%)"),
         imp_vwet_days1_am_99p9_diff = paste0(100*round(imp_vwet_days1_am_99p9_diff_mean, 4), "% (",
                                              100*round(imp_vwet_days1_am_99p9_diff_perc10, 4), " to ",
                                              100*round(imp_vwet_days1_am_99p9_diff_perc90, 4), "%)"),
         # write the warming level as "+3C" instead of "3"
         warming_level = paste0("+", warming_level, "\u00B0C")) %>%
  # reformat and write out as a .tex table
  select(warming_level, imp_total, ends_with("diff")) %>%
  pivot_longer(-warming_level) %>% pivot_wider(names_from=warming_level, values_from=value) %>%
  mutate(name = name %>% rewrite_impact_label()) %>%
  setNames(names(.) %>% str_replace("^name$", "Impact channel")) %>%
  xtable(label = "tab:si_globalresults_bc",
         caption = "Global mean GDP impacts by warming level and climate indicator (upper and lower decile in parentheses)") %>%
  print(file = file.path("..", "sharepoint", "Tables", paste0(Sys.Date(), " si_globalresults_overview_bc.tex")),
        size="\\tiny")
  
# function to plot distribution & summary stats of global GDP impacts by impact variable
figure_impdistr_global <- function(data_aggr = NULL, data_distr = NULL, var_selection = c("imp_Tstd_diff",
                                                                                          "imp_Wn_diff",
                                                                                          "imp_wet_days_1_diff",
                                                                                          "imp_vwet_days1_am_99p9_diff"),
                                   warming_levels_shown = c(1.5, 2, 3)) {
  data_aggr %>%
    setNames(names(.) %>% str_replace("agreement_meansign", "agreementmeansign")) %>%
    pivot_longer(cols = starts_with("imp"), names_to = "variable") %>%
    # extract the summary stat label in the 'variable' column, then remove it
    mutate(summary_stat = str_extract(variable, "_[A-Za-z0-9]+$") %>% str_remove("^_"),
           variable = str_remove(variable, "_[A-Za-z0-9]+$")) %>%
    # subset to the climate indicators & warming levels selected
    filter(variable %in% var_selection) %>%
    filter(warming_level %in% warming_levels_shown) %>%
    # convert to wide format and use chart-compatible labels for the climate indicators
    pivot_wider(names_from = "summary_stat", values_from = "value") %>%
    mutate(variable = rewrite_impact_label_linebreaks(variable)) %>%
    ggplot(aes(paste0(warming_level, "\u00B0C"), colour = variable)) +
    # distribution as a violin chart - this requires the full data, not the aggregate summary stats
    geom_violin(data = data_distr %>% left_join(df_gwl %>% select(-ensemble),
                                                by = c("model", "scenario", "year")) %>%
                  pivot_longer(cols = starts_with("imp"), names_to = "variable") %>%
                  filter(variable %in% var_selection,
                         !is.na(warming_level),
                         warming_level %in% warming_levels_shown) %>%
                  mutate(variable = rewrite_impact_label_linebreaks(variable)),
                aes(y = value, colour = variable, fill = variable), position = position_dodge(width = 0.6), alpha = 0.3,
                linewidth = 0.1) +
    # horizontal line to mark zero
    geom_hline(yintercept = 0, colour = "black", linetype = "dashed", alpha = 0.8)  +
    # point for the average impact, errobar indicating upper-to-lower-decile range
    geom_point(aes(y = mean), position = position_dodge(width = 0.6)) +
    geom_errorbar(aes(ymin = perc10, ymax = perc90), position = position_dodge(width = 0.6), width =0.3) +
    guides(fill = NULL, alpha = NULL) + 
    labs(y = "Global GDP impact", x = "Global warming level",colour = NULL, fill = NULL ) +
    scale_y_continuous(labels = scales::percent) +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.background = element_blank(),
          legend.box.background = element_blank())
}

# panel for total impacts, annual temp, annual precip and variability & extremes (joint impacts) - with manual GDP shock comparisons
# NOTE: selected year-to-year GDP contractions are sourced from 'data/input/230503 World Bank GDP_2015USD.csv'
figure1_1 <- figure_impdistr_global(df_global_fulldistr_aggr, df_global_fulldistr, var_selection = c("imp_total", "imp_temp_diff", "imp_varextremes",
                                                                                                     "imp_Pt_diff")) +
  geom_hline(yintercept = -0.263, linetype = "dotted", colour = "grey") +
  annotate("text", x = 0.5, y = -0.253, label = "Syrian GDP contraction in 2012 (civil war)", size = 3, fontface = 2, hjust = 0, colour = "azure4") +
  geom_hline(yintercept = -0.1015, linetype = "dotted", colour = "grey") +
  annotate("text", x = 0.5, y = -0.0915, label = "Greece GDP contraction in 2011 (Euro crisis)", size = 3, fontface = 2, hjust = 0, colour = "azure4") +
  scale_colour_manual(values = c("#E64B35FF", "#4DBBD5FF", "#3C5488FF", "#229954")) +
  scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#3C5488FF", "#229954"))

# same chart for the four variability & extremes indicators disaggregated - again, with manual points of comparison based on World Bank WDI data (for 'World')
figure1_2 <- figure_impdistr_global(df_global_fulldistr_aggr, df_global_fulldistr) +
  geom_hline(yintercept = -0.0134, linetype = "dotted", colour = "grey") +
  annotate("text", x = 0.5, y = -0.0114, label = "Global GDP contraction in 2009 (financial crisis)", size = 3, fontface = 2, hjust = 0, colour = "azure4") +
  geom_hline(yintercept = -0.0311, linetype = "dotted", colour = "grey") +
  annotate("text", x = 0.5, y = -0.0291, label = "Global GDP contraction in 2020 (Covid pandemic)", size = 3, fontface = 2, hjust = 0, colour = "azure4") +
  scale_colour_manual(values = c("#229954", "#229954", "#229954", "#229954"))

# repeat for bias-corrected data w/o horizontal lines for GDP shocks
figure1_1_bc <- figure_impdistr_global(df_global_bc_fulldistr_aggr, df_global_bc_fulldistr, var_selection = c("imp_total", "imp_temp_diff", "imp_varextremes",
                                                                                                     "imp_Pt_diff")) +
  scale_colour_manual(values = c("#E64B35FF", "#4DBBD5FF", "#3C5488FF", "#229954")) +
  scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#3C5488FF", "#229954"))
figure1_2_bc <- figure_impdistr_global(df_global_bc_fulldistr_aggr, df_global_bc_fulldistr) +
  scale_colour_manual(values = c("#229954", "#229954", "#229954", "#229954"))

# compare mean temperature damages for main spec & status quo (= Figure 3, panel 1)
figure_global_kotzvsstatusquo <- function(data = NULL, data_statusquo = NULL, warming_levels_shown = c(1.5, 2, 3)) {
  # bind together the results for the main & the status quo specification
  bind_rows(data %>% select(warming_level, starts_with("imp_total")) %>% mutate(spec = "All climate indicators (= Kotz et al., 2022)"),
            data_statusquo %>% select(warming_level, starts_with("imp_total")) %>% mutate(spec = "Annual temperature only (= status quo)")) %>%
    filter(warming_level %in% warming_levels_shown) %>%
    mutate(spec = factor(spec, levels = c("Annual temperature only (= status quo)", "All climate indicators (= Kotz et al., 2022)"))) %>%
    # create labelled bar chart with total impacts at each warming level & errobar for upper-to-lower-decile range
    ggplot(aes(paste0(warming_level, "\u00B0C"), imp_total_mean, group = spec)) + geom_col(aes(fill = spec), position = position_dodge()) +
    geom_text(aes(y = 0.005, label = paste0(100*round(imp_total_mean, 3), "%")), position = position_dodge(width = 1), size = 3) +
    geom_errorbar(aes(ymin = imp_total_perc10, ymax = imp_total_perc90), position = position_dodge(width = 1), width = 0.5) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
    labs(y = "Global GDP impact", x = "Global warming level", fill = "Impact projections based on...") +
    theme_classic() +
    theme(legend.position = c(0.385, 0.14), legend.background = element_rect(fill='transparent'))
}

# create the chart for raw and bias-corrected data
figure3_1 <- figure_global_kotzvsstatusquo(df_global_fulldistr_aggr, df_global_fulldistr_T_Pt_Pt2zero_aggr)
figure3_1_bc <- figure_global_kotzvsstatusquo(df_global_bc_fulldistr_aggr, df_global_bc_fulldistr_T_Pt_Pt2zero_aggr)


## marginal effects of annual mean temperature for different controls

# write a function to calculate marginal effects of annual mean temp. following Kalkuhl & Wenz 2020 (Footnote 6)
marginal_effect_temperature <- function(T_0 = NULL, T_shock = 1, coef_used = NULL) {
  out <- coef_used["TmeanD"] + coef_used["TmeanLD"] + T_0*(coef_used["TmeanD.Tmean"] + coef_used["TmeanLD.TmeanL"])
  names(out) <- NULL
  return(out)
}
# vectorize the function and write a wrapper for Monte Carlo analysis
marginal_effect_temperature <- Vectorize(marginal_effect_temperature, vectorize.args = "T_0")
marginal_effect_temperature_mc <- function(T_0_range = seq(-10, 30, by = 0.01), T_shock = 1, coef_mc_matrix = NULL) {
  map_dfr(1:ncol(coef_mc_matrix), ~ tibble(T_0 = T_0_range,
                                           damages = marginal_effect_temperature(T_0 = T_0_range, T_shock = T_shock, coef_used = coef_mc_matrix[, .x]),
                                           monte_carlo_draw = .x)
  )
}

# use the function to calculate uncertainty ranges by drawing coefficients from multivariate Gaussian for 1,000 MC draws
if(redo_data_objects) {
  df_margeffect_uncertainty <- bind_rows(marginal_effect_temperature_mc(coef_mc_matrix = get_coefficients(draw_from_vcov = T,
                                                                                                     model_used = "main",
                                                                                                     n=1000)) %>% mutate(spec = "main"),
                                        marginal_effect_temperature_mc(coef_mc_matrix = get_coefficients(draw_from_vcov = T,
                                                                                                         model_used = "T_Pt",
                                                                                                         n=1000)) %>% mutate(spec = "T_Pt"),
                                        marginal_effect_temperature_mc(coef_mc_matrix = get_coefficients(draw_from_vcov = T,
                                                                                                         model_used = "T_Pt_Tstd",
                                                                                                         n=1000)) %>% mutate(spec = "T_Pt_Tstd")
  )
  
  saveRDS(df_margeffect_uncertainty, file = file.path("data", "df_margeffect_uncertainty.rds"))
} else {
  df_margeffect_uncertainty <- readRDS(file.path("data", "df_margeffect_uncertainty.rds"))
}


# get marginal impact of +1C and the uncertainty ranges based on Monte Carlo draws from multivariate Gaussian
figure3_2 <- df_margeffect_uncertainty %>%
  mutate(spec = factor(spec %>% str_replace("^T_Pt$", "Annual precipitation only\n (= status quo)") %>%
                                 str_replace("^T_Pt_Tstd$", "+ Temperature variability") %>%
                                 str_replace("^main$", "+ Other precipitation indicators\n (= Kotz et al. 2022)"),
                       levels = c("Annual precipitation only\n (= status quo)",
                                  "+ Temperature variability",
                                  "+ Other precipitation indicators\n (= Kotz et al. 2022)"
                       )
  )) %>%
  # NOTE: marginal_effect_temperature() gives results in log-scale, so we convert to %GDP
  mutate(damages = log_to_gdp_impacts(damages)) %>%
  ggplot(aes(T_0, damages, colour = spec)) +
  stat_summary(aes(linetype = spec), geom="line", fun=mean) +
  # 95% CI as a geom_ribbon
  stat_summary(data = df_margeffect_uncertainty %>% filter(spec %in% c("main", "T_Pt")),
               aes(fill = spec), geom="ribbon", fun.min=~quantile(.x, 0.025), fun.max = ~quantile(.x, 0.975), alpha = 0.15, colour = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4")) +
  scale_colour_manual(values = c("#00BFC4", "black", "#F8766D")) +
  scale_linetype_manual(values = c("solid", "dotted", "solid")) +
  guides(fill = "none") +
  labs(colour = "Controlling for ...", y = "Marginal GDP impact of annual temperature\n(+1\u00B0C increase)",
       fill = "Controlling for ...", linetype = "Controlling for ...",
       x = "Annual temperature (\u00B0C)") +
  theme_classic() +
  theme(legend.position = c(0.3, 0.2), legend.background=element_blank())


### Variance decomposition following Hsiang et al 2017

# merge warming levels and model weights into the full distribution
df_global_fulldistr_gwl <- df_global_fulldistr %>% left_join(df_gwl %>% select(-ensemble), by = c("model", "scenario", "year")) %>%
  filter(warming_level >= 1) %>%
  pivot_longer(cols = starts_with("imp"), names_to = "variable", values_to = "imp") %>%
  mutate(scenarioyear = paste0(scenario, year)) %>%
  left_join(df_gcm_weights, by = c("model", "warming_level"))

# repeat for bias-corrected data
df_global_bc_fulldistr_gwl <- df_global_bc_fulldistr %>% left_join(df_gwl %>% select(-ensemble), by = c("model", "scenario", "year")) %>%
  filter(warming_level >= 1) %>%
  pivot_longer(cols = starts_with("imp"), names_to = "variable", values_to = "imp") %>%
  mutate(scenarioyear = paste0(scenario, year)) %>%
  left_join(df_gcm_weights, by = c("model", "warming_level"))

# write a function for variance decomposition
decompose_variance_hsiang <- function(data = NULL, is_data_global = T, warming_levels_shown = c(1.5, 2, 3),
                                      dir_charts = file.path("..", "sharepoint", "Figures", "global", "socioeconomic"),
                                      charts_suffix = NULL, medianmodel_type = NULL) {
  # check inputs
  if(!medianmodel_type %in% c("by_indicator_3deg", "totaldamages_3deg", unique(data$model))) stop("medianmodel_type has a value that is not supported")
  
  # if the data is global, we add a GID_0 column, so we can use the same steps for global results & country-level decompositions
  if(is_data_global) data <- data %>% mutate(GID_0 = "Global")
  
  # identify the median model
  df_medianmodel <- data %>%
    # NOTE: we take means per model, so no need for using GCM weights (which are for across-model summary stats)
    group_by(GID_0, warming_level, model, variable) %>% summarise(imp = mean(imp)) %>% ungroup() %>%
    group_by(GID_0, warming_level, variable) %>%
    arrange(GID_0, warming_level, variable, imp) %>%
    mutate(n = n(),
           ind = 1:n(),
           is_median_model = if_else(n %% 2 == 0, ind == n/2, ind == n/2 + 0.5),
           variable_clean = rewrite_impact_label(variable)) %>%
    ungroup()
  
  # chart the median-model for each climate indicator & warming level
  if(!is.null(charts_suffix)) {
    for(gwl_jj in warming_levels_shown) {
      ggarrange(
        plotlist = levels(df_medianmodel$variable_clean) %>% str_subset("except|OLD\\sVERSION", negate = T) %>%
          map(~ df_medianmodel %>% filter(variable_clean == .x, warming_level == gwl_jj) %>%
                mutate(median_label = if_else(is_median_model, "Median model", "Other CMIP6")) %>%
                ggplot(aes(reorder(model, imp), imp)) +
                geom_hline(yintercept = 0, linetype = "dashed", colour = "black", alpha = 0.7) +
                geom_col(aes(fill = median_label)) + facet_wrap(~ variable_clean) +
                scale_y_continuous(labels = scales::percent, n.breaks = 3) +
                labs(x = NULL, y = paste0("Average global GDP impact\nat +", gwl_jj, "\u00B0C in %"), fill = NULL) +
                coord_flip()),
        nrow = 2, ncol = 4, legend = "bottom", common.legend = T
      )
      ggsave(file.path(dir_charts, paste0(Sys.Date(), " HsiangEtAl_medianmodel_", gwl_jj, "_", charts_suffix, ".png")), width = 14, height = 10)
    }
  }
  
  # merge in median model for +3C
  if(medianmodel_type == "by_indicator_3deg") {
    data <- data %>%
      left_join(df_medianmodel %>% filter(warming_level == 3, is_median_model) %>% select(GID_0, model, variable, is_median_model), by = c("GID_0", "model", "variable")) %>%
      mutate(is_median_model = replace_na(is_median_model, F))
  } else if(medianmodel_type == "totaldamages_3deg") {
    data <- data %>%
      left_join(df_medianmodel %>% filter(warming_level == 3, variable == "imp_total", is_median_model) %>% select(GID_0, model, is_median_model), by = c("GID_0", "model")) %>%
      mutate(is_median_model = replace_na(is_median_model, F))
      #mutate(is_median_model = model == df_medianmodel %>% filter(warming_level == 3, variable == "imp_total", is_median_model) %>% pull(model))
  # if medianmodel_type provides a model name directly, we use it to set up the dummy
  } else if(medianmodel_type %in% unique(data$model)) {
    data <- data %>% mutate(is_median_model = model == medianmodel_type)
  } else {
    stop("Provided input for argument medianmodel_type is not supported")
  }
  
  # identify the median scenario-year for each model
  df_medianscenarioyear <- data %>% 
    # again, we take model-specific means, so no need for GCM weights
    group_by(GID_0, model, warming_level, scenarioyear, variable) %>% summarise(imp = mean(imp)) %>% ungroup() %>%
    group_by(GID_0, warming_level, variable, model) %>%
    arrange(GID_0, warming_level, variable, imp) %>%
    mutate(n = n(),
           ind = 1:n(),
           is_median_scenarioyear = if_else(n %% 2 == 0, ind == n/2, ind == n/2 + 0.5),
           variable_clean = rewrite_impact_label(variable)) %>%
    ungroup()
  
  # plot
  if(!is.null(charts_suffix)) {
    for(gwl_jj in warming_levels_shown) {
      ggarrange(
        plotlist = levels(df_medianscenarioyear$variable_clean) %>% str_subset("except|OLD\\sVERSION", negate = T) %>%
          map(~ df_medianscenarioyear %>% filter(variable_clean == .x, warming_level == gwl_jj) %>%
                mutate(median_label = if_else(is_median_scenarioyear, "Median scenario-year", "Other"),
                       scenario = str_extract(scenarioyear, "ssp[:digit:][:digit:][:digit:]"),
                       year = str_remove(scenarioyear, scenario),
                       scenarioyear_clean = paste0(year, " (", scenario, ")")) %>%
                filter(model == "MPI-ESM1-2-HR") %>%
                ggplot(aes(reorder(scenarioyear_clean, imp), imp)) +
                geom_hline(yintercept = 0, linetype = "dashed", colour = "black", alpha = 0.7) +
                geom_col(aes(fill = median_label)) + facet_wrap(~ variable_clean) +
                scale_y_continuous(labels = scales::percent, n.breaks = 3) +
                labs(x = NULL, y = paste0("Average global GDP impact\nat +", gwl_jj, "\u00B0C in %"), fill = NULL) +
                coord_flip()),
        nrow = 2, ncol = 4, legend = "bottom", common.legend = T
      )
      ggsave(file.path(dir_charts, paste0(Sys.Date(), " HsiangEtAl_medianscenarioyear_", gwl_jj, "_", charts_suffix, ".png")), width = 14, height = 14)
    }
  }
  
  # merge in medians for all models & warming levels (since scenario-years are GWL-specific)
  data <- data %>%
    left_join(df_medianscenarioyear %>% filter(is_median_scenarioyear) %>% select(GID_0, warming_level, model, scenarioyear, variable, is_median_scenarioyear),
              by = c("GID_0", "model", "scenarioyear", "variable", "warming_level")) %>%
    mutate(is_median_scenarioyear = replace_na(is_median_scenarioyear, F))

  # identify the median monte carlo draw
  df_medianmc <- data %>%
    # we take means per MC draw across models, so we need to use GCM weights
    group_by(GID_0, warming_level, monte_carlo_draw, variable) %>% summarise(imp = Hmisc::wtd.mean(imp, weights = gcm_weight, normwt = T)) %>% ungroup() %>%
    group_by(GID_0, warming_level, variable) %>%
    arrange(GID_0, warming_level, variable, imp) %>%
    mutate(n = n(),
           ind = 1:n(),
           is_median_mc = if_else(n %% 2 == 0, ind == n/2, ind == n/2 + 0.5),
           variable_clean = rewrite_impact_label(variable)) %>%
    ungroup()
  
  # plot
  if(!is.null(charts_suffix)) {
    for(gwl_jj in warming_levels_shown) {
      ggarrange(
        plotlist = levels(df_medianmc$variable_clean) %>% str_subset("except|OLD\\sVERSION", negate = T) %>%
          map(~ df_medianmc %>% filter(variable_clean == .x, warming_level == gwl_jj) %>%
                mutate(median_label = if_else(is_median_mc, "Median parameter draw", "Other draws")) %>%
                ggplot(aes(reorder(monte_carlo_draw, imp), imp)) +
                geom_hline(yintercept = 0, linetype = "dashed", colour = "black", alpha = 0.7) +
                geom_col(aes(fill = median_label)) + facet_wrap(~ variable_clean) +
                scale_y_continuous(labels = scales::percent, n.breaks = 3) +
                labs(x = NULL, y = paste0("Average global GDP impact\nat +", gwl_jj, "\u00B0C in %"), fill = NULL) +
                theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
                coord_flip()),
        nrow = 2, ncol = 4, legend = "bottom", common.legend = T
      )
      ggsave(file.path(dir_charts, paste0(Sys.Date(), " HsiangEtAl_medianmc_", gwl_jj, "_", charts_suffix, ".png")), width = 14, height = 10)
    }
  }
  
  # merge in median for +3C
  data <- data %>%
    left_join(df_medianmc %>% filter(is_median_mc, warming_level == 3) %>%
                select(GID_0, monte_carlo_draw, variable, is_median_mc), by = c("GID_0", "monte_carlo_draw", "variable")) %>%
    mutate(is_median_mc = replace_na(is_median_mc, F))
  
  # calculate variance over MC draws
  df_vardisaggr_hsiang <- list()
  df_vardisaggr_hsiang$between_mc <- data %>% filter(is_median_model & is_median_scenarioyear) %>%
    group_by(GID_0, warming_level, variable) %>%
    summarise(variance = var(imp),
              variance_source = "monte_carlo_draw",
              n = n()) %>% ungroup()
  # between CMIP6 models
  df_vardisaggr_hsiang$between_model <- data %>% filter(is_median_scenarioyear & is_median_mc) %>%
    group_by(GID_0, warming_level, variable) %>%
    # NOTE: we summarize across models, hence we use GCM weights here
    summarise(variance = Hmisc::wtd.var(imp, weights = gcm_weight, normwt = T),
              variance_source = "model",
              n = n()) %>% ungroup()

  # between scenario-years
  df_vardisaggr_hsiang$between_scenarioyear <- data %>% filter(is_median_model & is_median_mc) %>%
    group_by(GID_0, warming_level, variable) %>%
    summarise(variance = var(imp),
              variance_source = "scenarioyear",
              n = n()) %>% ungroup()
  
  #
  df_vardisaggr_hsiang <- df_vardisaggr_hsiang %>% reduce(bind_rows) %>%
    pivot_wider(values_from = "variance", names_from = "variance_source", -n) %>%
    left_join(data %>% group_by(GID_0, warming_level, variable) %>% summarise(variance_total = var(imp)) %>% ungroup(),
              by = c("GID_0", "warming_level", "variable")) %>%
    mutate(variance_explained = monte_carlo_draw + model + scenarioyear,
           interaction = variance_total - variance_explained)
  
  return(df_vardisaggr_hsiang)
}

# plot variances
figure_vardisaggr_global_hsiang <- function(data = NULL, warming_levels_shown = c(1.5, 2, 3), nrow_used = NA, ncol_used = NA,
                                            show_interaction = TRUE, interaction_as_bar = TRUE,
                                            absolute_value_rescaling = TRUE) {
  data %>%
    filter(warming_level %in% warming_levels_shown) %>%
    # if absolute_value_rescaling is TRUE, we rescale the variance if interaction term is negative and then take absolute values
    {if(absolute_value_rescaling) mutate(.,
                                         variance_total = if_else(interaction < 0, variance_total - 2*interaction, variance_total),
                                         interaction = abs(interaction)) else .} %>%
  {if(interaction_as_bar) pivot_longer(., cols = c("monte_carlo_draw", "model", "scenarioyear", "interaction"),
                                       names_to = "variance_source", values_to = "variance") else pivot_longer(., cols = c("monte_carlo_draw", "model", "scenarioyear"),
                                                      names_to = "variance_source", values_to = "variance")
                                         }%>%
    group_by(GID_0, warming_level, variable) %>%
    {if(interaction_as_bar) mutate(., variance_share = variance/variance_total) else mutate(., variance_share = variance/variance_explained)} %>%
  mutate(variance_total_relative = variance_total/variance_explained) %>%
  {if(!interaction_as_bar) filter(., variance_source != "interaction") else . } %>%
  mutate(variable = rewrite_impact_label_linebreaks(variable),
         variance_source = variance_source %>% str_replace("^model$", "Climate model\nuncertainty") %>%
           str_replace("^year$", "Inter-annual\nvariability") %>%
           str_replace("^scenarioyear$", "Inter-annual\nvariability") %>%
           str_replace("^monte_carlo_draw$", "Damage function\nuncertainty") %>%
           str_replace("^scenario$", "Scenario\nuncertainty") %>%
           str_replace("^interaction$", "Interaction/\nshared variance")) %>%
  mutate(variance_source = factor(variance_source, levels = c("Damage function\nuncertainty", "Inter-annual\nvariability", "Climate model\nuncertainty",
                                                              "Scenario\nuncertainty", "Interaction/\nshared variance"))) %>%
  ungroup() %>%
  ggplot(aes(paste0(warming_level, "\u00B0C"), variance_share)) +
    {if(show_interaction) geom_hline(yintercept = 1, alpha = 0.3, linetype = "dashed")} +
  geom_bar(aes(fill = variance_source), stat = "identity", position = position_stack(vjust = 0.5, reverse = F), alpha = 0.7) +
  geom_label(aes(fill = variance_source, label = paste0(100*(round(variance_share, ifelse(round(variance_share, 2) == 0, 3, 2))), "%")),
             position = position_stack(vjust = 0.5, reverse = F),
             size = 2.5,
             show.legend = F) +
  {if(show_interaction & !interaction_as_bar) geom_point(aes(y = variance_total_relative, shape = "Total variance\nplus interaction"), colour = "black", alpha = 0.5)} +

  labs(x = "Global warming level", y = "Share in variance", colour = NULL, fill = NULL, shape = NULL) +
  scale_fill_manual(values = c("#F8766D", "#619CFF", "#00BA38", "grey")) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~ variable, nrow = 2, ncol = 4, scales = ifelse(show_interaction & !absolute_value_rescaling, "free_y", "fixed")) +
  theme_classic() +
  theme(legend.position = "bottom")
}

# # make the variance share charts as stacked bar charts
figure1_3_bar <- figure_vardisaggr_global_hsiang(decompose_variance_hsiang(df_global_fulldistr_gwl, medianmodel_type = "ACCESS-CM2",
                                                                           charts_suffix = "ACCESS-CM2") %>%
                                                   filter(variable %in% c("imp_total", "imp_temp_diff", "imp_Pt_diff", "imp_varextremes")),
                                                     nrow_used = 1, ncol_used = 4)
figure1_4_bar <- figure_vardisaggr_global_hsiang(decompose_variance_hsiang(df_global_fulldistr_gwl, medianmodel_type = "ACCESS-CM2") %>%
                                                   filter(variable %in% c("imp_Tstd_diff", "imp_Wn_diff", "imp_wet_days_1_diff", "imp_vwet_days1_am_99p9_diff")),
                                                 nrow_used = 1, ncol_used = 4)

figure1_3_bar_bc <- figure_vardisaggr_global_hsiang(decompose_variance_hsiang(df_global_bc_fulldistr_gwl, medianmodel_type = "ACCESS-CM2",
                                                                              charts_suffix = "bc_ACCESS-CM2") %>%
                                                   filter(variable %in% c("imp_total", "imp_temp_diff", "imp_Pt_diff", "imp_varextremes")),
                                                 nrow_used = 1, ncol_used = 4)
figure1_4_bar_bc <- figure_vardisaggr_global_hsiang(decompose_variance_hsiang(df_global_bc_fulldistr_gwl, medianmodel_type = "ACCESS-CM2") %>%
                                                   filter(variable %in% c("imp_Tstd_diff", "imp_Wn_diff", "imp_wet_days_1_diff", "imp_vwet_days1_am_99p9_diff")),
                                                 nrow_used = 1, ncol_used = 4)


# SI versions
figure_vardisaggr_global_hsiang(decompose_variance_hsiang(df_global_fulldistr_gwl, medianmodel_type = "ACCESS-CM2"),
                                absolute_value_rescaling = F,
                                nrow_used = 2, ncol_used = 4)
ggsave(file.path("..", "sharepoint", "Figures", paste0(Sys.Date(), " Figure2_var_withnegativeinteraction_ACCESS-CM2.png")),
       height = 7, width = 10)
figure_vardisaggr_global_hsiang(decompose_variance_hsiang(df_global_fulldistr_gwl, medianmodel_type = "totaldamages_3deg"),
                                absolute_value_rescaling = F,
                                nrow_used = 2, ncol_used = 4)
ggsave(file.path("..", "sharepoint", "Figures", paste0(Sys.Date(), " Figure2_var_withnegativeinteraction_medianmodel_totaldamages_3deg.png")),
       height = 7, width = 10)
figure_vardisaggr_global_hsiang(decompose_variance_hsiang(df_global_fulldistr_gwl, medianmodel_type = "EC-Earth3"),
                                absolute_value_rescaling = F,
                                nrow_used = 2, ncol_used = 4)
ggsave(file.path("..", "sharepoint", "Figures", paste0(Sys.Date(), " Figure2_var_withnegativeinteraction_medianmodel_EC-Earth3.png")),
       height = 7, width = 10)


# create Figure 1 on global results (distribution summary stats & variance decomposition by impact channels)
figure1 <- ggarrange(ggarrange(figure1_1, figure1_3_bar, nrow = 2, ncol = 1,
                    labels = c("a", "c")),
          ggarrange(figure1_2, figure1_4_bar, nrow = 2, ncol = 1,
                    labels = c("b", "d")),
          nrow = 1, ncol = 2)
ggsave(file.path("..", "sharepoint", "Figures",
                 paste0(Sys.Date(), " Figure2.pdf")),
       plot = figure1,
       width = 11, height = 7)

# make a bias-corrected version (using stacked bar charts)
figure1_bc <- ggarrange(ggarrange(figure1_1_bc, figure1_3_bar_bc, nrow = 2, ncol = 1,
                    labels = c("a", "c")),
          ggarrange(figure1_2_bc, figure1_4_bar_bc, nrow = 2, ncol = 1,
                    labels = c("b", "d")),
          nrow = 1, ncol = 2)
ggsave(file.path("..", "sharepoint", "Figures",
                 paste0(Sys.Date(), " Figure2_bc.pdf")),
       plot = figure1_bc,
       width = 11, height = 7)

# clean up the memory
gc()



################################################################################
############################## ADM0 IMPACTS ####################################
################################################################################


######################### IMPORT & AGGREGATE RESULTS  ##########################

# write a function to read in & aggregate each observation by 1/N_gcm (= # of GCM obs for the respective warming level = 20 * the # of SSPS for which model reaches GWL)
import_and_aggregate_adm0_distr <- function(file_data = NULL, include_intramodel_variation = TRUE, gcm_weights_used = NULL) {
  unique(file_data$GID_0) %>%
    future_map_dfr(.f = function(iso) {
      df_distr_adm0 <- file_data$filename[file_data$GID_0 == iso] %>%
        map_dfr(read_feather) %>%
        left_join(df_gwl %>% filter(warming_level >= 1) %>% select(-ensemble), by = c("model", "scenario", "year")) %>%
        filter(!is.na(warming_level)) %>%
        left_join(gcm_weights_used, by = c("model", "warming_level"))
      
      # calculate worst-in-20-years effects
      df_onein20_adm0 <- df_distr_adm0 %>%
        # for each global warming level window (20 years), get the worst impact (no weighting needed as this is within-GCM)
        group_by(model, scenario, monte_carlo_draw, warming_level, GID_0) %>%
        summarise_at(vars(starts_with("imp")), .funs = list(onein20 = ~ min(.x))) %>%
        ungroup() %>%
        # aggregate across the whole distribution
        left_join(gcm_weights_used, by = c("model", "warming_level")) %>%
        group_by(warming_level, GID_0) %>%
        summarise_at(vars(starts_with("imp")),
                     .funs = list(
                       mean = ~ Hmisc::wtd.mean(.x, weights = gcm_weight, normwt = T),
                       var = ~ Hmisc::wtd.var(.x, weights = gcm_weight, normwt = T),
                       agreement_meansign = ~ if_else(Hmisc::wtd.mean(.x, weights = gcm_weight, normwt = T) < 0,
                                                      Hmisc::wtd.mean(.x < 0, weights = gcm_weight, normwt = T),
                                                      Hmisc::wtd.mean(.x >= 0, weights = gcm_weight, normwt = T)
                       ))) %>% ungroup()
      
      df_out <- df_distr_adm0 %>%
        # if include_intramodel_variation = FALSE, we compute weighted average across scenario-years for each model & MC draw (= exclude inter-annual variability)
        { if(!include_intramodel_variation) {
          group_by(., warming_level, model, monte_carlo_draw, GID_0) %>%
            summarise_at(vars(starts_with("imp")), ~ Hmisc::wtd.mean(.x, weights = gcm_weight, normwt = T)) %>%
            ungroup() %>%
            # we add a uniform weight of 1 to all variables because in later part of the pipe, we use weighted mean formulas
            mutate(gcm_weight = 1)
        } else . } %>%
        group_by(warming_level, GID_0) %>%
        summarise_at(vars(starts_with("imp")), .funs = list(
          # calculate certainty equivalents using src/utils/analyse_impacts.R
          # NOTE: this is currently commented out because we would need to implement GCM weights into the CE function, low priority
          #ce_eta1 = ~ get_certainty_equivalent(.x, eta = 1),
          #ce_eta1p35 = ~ get_certainty_equivalent(.x, eta = 1.35),
          mean = ~ Hmisc::wtd.mean(.x, weights = gcm_weight, normwt = T),
          var = ~ Hmisc::wtd.var(.x, weights = gcm_weight, normwt = T),
          min = min,
          max = max,
          perc05 = ~ Hmisc::wtd.quantile(.x, probs = 0.05, weights = gcm_weight, normwt = T),
          perc10 = ~ Hmisc::wtd.quantile(.x, probs = 0.1, weights = gcm_weight, normwt = T),
          median = ~ Hmisc::wtd.quantile(.x, probs = 0.5, weights = gcm_weight, normwt = T),
          perc90 = ~ Hmisc::wtd.quantile(.x, probs = 0.9, weights = gcm_weight, normwt = T),
          perc95 = ~ Hmisc::wtd.quantile(.x, probs = 0.95, weights = gcm_weight, normwt = T),
          agreement_meansign = ~ if_else(Hmisc::wtd.mean(.x, weights = gcm_weight, normwt = T) < 0,
                                         Hmisc::wtd.mean(.x < 0, weights = gcm_weight, normwt = T),
                                         Hmisc::wtd.mean(.x >= 0, weights = gcm_weight, normwt = T)
          ),
          signalnoise = ~ Hmisc::wtd.mean(.x, weights = gcm_weight, normwt = T)/(Hmisc::wtd.var(.x, weights = gcm_weight, normwt = T))^0.5
        )) %>% ungroup() %>%
        # merge in the worst-in-20-years results
        left_join(df_onein20_adm0, by = c("warming_level", "GID_0"))
      
      return(df_out)
    }, .progress = T, .options = furrr_options(seed = NULL))
}

### read in & aggregate distributions for different specifications

if(redo_data_objects) {
  # point estimate results - no bias correction
  df_adm0_pointdistr <- import_and_aggregate_adm0_distr(df_pointfiles, gcm_weights_used = df_gcm_weights)
  df_adm0_bc_pointdistr <- import_and_aggregate_adm0_distr(df_pointfiles_bc, gcm_weights_used = df_gcm_weights)
  df_adm0_fulldistr <- import_and_aggregate_adm0_distr(df_mcfiles, gcm_weights_used = df_gcm_weights) %>%
    mutate(spec = "kotz2022")
  df_adm0_bc_fulldistr <- import_and_aggregate_adm0_distr(df_mcfiles_bc, gcm_weights_used = df_gcm_weights) %>%
    mutate(spec = "kotz2022")
  df_adm0_fulldistr_T_Pt_Pt2zero <- import_and_aggregate_adm0_distr(df_mcfiles_T_Pt_Pt2zero, gcm_weights_used = df_gcm_weights) %>%
    mutate(spec = "T_Pt_Pt2zero")
  df_adm0_bc_fulldistr_T_Pt_Pt2zero <- import_and_aggregate_adm0_distr(df_mcfiles_bc_T_Pt_Pt2zero, gcm_weights_used = df_gcm_weights) %>%
    mutate(spec = "T_Pt_Pt2zero")
  
  # save out
  saveRDS(df_adm0_pointdistr, file = file.path("data", "df_adm0_pointdistr.rds"))
  saveRDS(df_adm0_bc_pointdistr, file = file.path("data", "df_adm0_bc_pointdistr.rds"))
  saveRDS(df_adm0_fulldistr, file = file.path("data", "df_adm0_fulldistr.rds"))
  saveRDS(df_adm0_bc_fulldistr, file = file.path("data", "df_adm0_bc_fulldistr.rds"))
  saveRDS(df_adm0_fulldistr_T_Pt_Pt2zero, file = file.path("data", "df_adm0_fulldistr_T_Pt_Pt2zero.rds"))
  saveRDS(df_adm0_bc_fulldistr_T_Pt_Pt2zero, file = file.path("data", "df_adm0_bc_fulldistr_T_Pt_Pt2zero.rds"))
} else {
  # read in (to save runtime)
  df_adm0_pointdistr <- readRDS(file.path("data", "df_adm0_pointdistr.rds"))
  df_adm0_fulldistr <- readRDS(file.path("data", "df_adm0_fulldistr.rds"))
  df_adm0_bc_fulldistr <- readRDS(file.path("data", "df_adm0_bc_fulldistr.rds"))
  df_adm0_fulldistr_T_Pt_Pt2zero <- readRDS(file.path("data", "df_adm0_fulldistr_T_Pt_Pt2zero.rds"))
  df_adm0_bc_fulldistr_T_Pt_Pt2zero <- readRDS(file.path("data", "df_adm0_bc_fulldistr_T_Pt_Pt2zero.rds"))
}


######################### COUNTRY-LEVEL MEAN IMPACTS & UNCERTAINTY #############

### make hatchmaps - CURRENTLY NOT USED IN THE PAPER

# create a agreement hatchmap using the function defined in src/utils/analyse_impacts.R
# apply the summarystat_agreement_hatchmap() function to all vars & GWLs, then save out charts
crossing(var = names(df_adm0_fulldistr) %>% str_subset("_mean$") %>%
           str_subset("onein20", negate = T) %>%
           str_remove_all("\\_mean$"),
         gwl = unique(df_adm0_fulldistr$warming_level)) %>%
  future_pmap(function(var, gwl) {
    out_agree <- summarystat_agreement_hatchmap(data = df_adm0_fulldistr,
                                                var_selected = var, gwl_selected = gwl,
                                                use_agreement = TRUE,
                                                summarystat_selected = "mean")
    
    ggsave(file.path("..", "sharepoint", "Figures", "adm0", "socioeconomic", "montecarlo", "agreement_maps",
                     paste0(Sys.Date(), " ", var, "_fulldistr_gwl", gwl, ".pdf")),
           plot = out_agree,
           width = 6, height = 4, dpi = 50)
    
    out_signal <- summarystat_agreement_hatchmap(data = df_adm0_fulldistr,
                                                 var_selected = var, gwl_selected = gwl,
                                                 use_agreement = FALSE,
                                                 summarystat_selected = "mean")
    
    ggsave(file.path("..", "sharepoint", "Figures", "adm0", "socioeconomic", "montecarlo", "signalnoise_maps",
                     paste0(Sys.Date(), " ", var, "_fulldistr_gwl", gwl, ".pdf")),
           plot = out_signal,
           width = 6, height = 4, dpi = 50)
    
    return("Map saved out")
  }, .progress = T, .options = furrr_options(seed = NULL))

# plot total impacts as a map
figure_map_impacts <- function(data = NULL, var_selected = "imp_total_mean", warming_level_selected = 3,
                               subtitle_selected = NULL) {
  # rename the variable selected
  names(data)[names(data) == var_selected] <- "var_for_plot"
  
  data %>%
    filter(warming_level == warming_level_selected) %>%
    left_join(world, by = "GID_0") %>%  st_as_sf() %>%
    ggplot() + geom_sf(aes(fill = ifelse(GID_0 %in% sovereign_countries_gid0, var_for_plot, NA)), size = 0.1) +
    scale_fill_gradient2(low = "red", mid = "gray88", high = "blue",
                         n.breaks = 4,
                         #breaks = fillscale_breaks,
                         labels = scales::percent,
                         # set the upper bound of the fill scale legend to zero if all values in the data are negative
                         limits = c(NA, ifelse((data %>% filter(warming_level == warming_level_selected,
                                                                GID_0 %in% sovereign_countries_gid0) %>%
                                                  pull(var_for_plot) %>% max()) < 0, 0, NA))) +
    labs(fill = NULL, subtitle = subtitle_selected) +
    theme_classic() +
    theme(legend.position ="bottom",
          legend.direction = "horizontal",
          axis.line = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank())
} 

# make all map of mean GDP impacts for total impacts & variability and extremes
figure2_1 <- figure_map_impacts(df_adm0_fulldistr,
                                subtitle_selected = "Total GDP impacts (all indicators) at +3\u00B0C")
figure2_1_bc <- figure_map_impacts(df_adm0_bc_fulldistr,
                                   subtitle_selected = "Total GDP impacts (all indicators) at +3\u00B0C")
figure2_2 <- figure_map_impacts(df_adm0_fulldistr, "imp_varextremes_mean",
                                subtitle_selected = "GDP impacts of variability & extremes at +3\u00B0C")
figure2_2_bc <- figure_map_impacts(df_adm0_bc_fulldistr, "imp_varextremes_mean",
                                   subtitle_selected = "GDP impacts of variability & extremes at +3\u00B0C")

### calculate impact of including all climate indicators by comparing to T_Pt_Pt2zero
# first, ensure that order of variables and observations are equal
are_impdata_subtractable <- function(df1 = NULL, df2 = NULL) {
  # we check if the warming_level and GID_0 vectors are identical (row identifiers) and if the column names match (col identifier)¨
  # if both is true, then the values in row i col j in both data frames correspond to each other
  identical(df1$warming_level, df2$warming_level) & identical(df1$GID_0, df2$GID_0) & identical(names(df1), names(df2))
}
if(!are_impdata_subtractable(df_adm0_fulldistr, df_adm0_fulldistr_T_Pt_Pt2zero)) stop("df_adm0_fulldistr and df_adm0_fulldistr_T_Pt_Pt2zero are not subtractable")
if(!are_impdata_subtractable(df_adm0_bc_fulldistr, df_adm0_bc_fulldistr_T_Pt_Pt2zero)) stop("df_adm0_bc_fulldistr and df_adm0_bc_fulldistr_T_Pt_Pt2zero are not subtractable")

# calculate ADM0-level differences in temperature impacts & total GDP impacts between i) full Kotz et al., 2022, ii) status quo, and iii) status quo + temp. variability
df_adm0_fulldistr_diffstatusquo <- as_tibble((df_adm0_fulldistr %>% select(-warming_level, -GID_0, -spec)) - (df_adm0_fulldistr_T_Pt_Pt2zero %>% select(-warming_level, -GID_0, -spec))) %>%
  bind_cols(df_adm0_fulldistr %>% select(warming_level, GID_0)) %>%
  mutate(spec = "kotz_vs_T_Pt_Pt2zero") %>%
  select(warming_level, GID_0, everything())

# repeat for bias corrected data
df_adm0_bc_fulldistr_diffstatusquo <- as_tibble((df_adm0_bc_fulldistr %>% select(-warming_level, -GID_0, -spec)) - (df_adm0_bc_fulldistr_T_Pt_Pt2zero %>% select(-warming_level, -GID_0, -spec))) %>%
  bind_cols(df_adm0_bc_fulldistr %>% select(warming_level, GID_0)) %>%
  mutate(spec = "kotz_vs_T_Pt_Pt2zero") %>%
  select(warming_level, GID_0, everything())

# make a map to compare impacts for all climate indicators and status quo
map_change_vs_statusquo <- function(data = NULL, var_selected = "imp_total_mean", warming_level_selected = 3,
                                    subtitle_selected = NULL) {
  names(data)[names(data) == var_selected] <- "var_for_plot"
  
  df_impactdiff_map <- data %>% filter(warming_level == warming_level_selected) %>%
    left_join(world, by = "GID_0") %>% st_as_sf() %>%
    mutate(var_for_plot = ifelse(GID_0 %in% sovereign_countries_gid0, var_for_plot, NA))
  
  df_impactdiff_map %>%
    ggplot() +
    geom_sf(aes(fill = var_for_plot), size = 0.1) +
    scale_fill_gradient2(midpoint = 0, 
                         #n.breaks = 3,
                         breaks = c(-0.02, 0),
                         labels = scales::percent,
                         low = "red",
                         # expand the fill scale to zero if necessary
                         limits = c(NA, ifelse(max(df_impactdiff_map$var_for_plot, na.rm = T) < 0, 0, NA)),
                         high = "blue",
                         mid = "gray88") +
    labs(fill = paste0("Difference in %-pts"), subtitle = subtitle_selected) +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          axis.line = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank())
}

# make the maps for total impacts
figure3_3 <- map_change_vs_statusquo(df_adm0_fulldistr_diffstatusquo,
                                     subtitle_selected = "   Change in GDP impact by using all climate indicators (instead of status quo) at +3\u00B0C global warming")
figure3_3_bc <- map_change_vs_statusquo(df_adm0_bc_fulldistr_diffstatusquo,
                                     subtitle_selected = "   Change in GDP impact by using all climate indicators (instead of status quo) at +3\u00B0C global warming")

# make the same maps only for how the damages of annual temperature change and save out
figure3_3_si <- map_change_vs_statusquo(df_adm0_fulldistr_diffstatusquo, var_selected = "imp_temp_diff_mean",
                                     subtitle_selected = "   Change in annual temperature impacts by using all climate indicators (instead of status quo)\n    at +3\u00B0C global warming")
figure3_3_si_bc <- map_change_vs_statusquo(df_adm0_bc_fulldistr_diffstatusquo, var_selected = "imp_temp_diff_mean",
                                           subtitle_selected = "   Change in annual temperature impacts by using all climate indicators (instead of status quo)\n    at +3\u00B0C global warming")
ggsave("../sharepoint/Figures/Figure3_3_SI.pdf",
       plot = figure3_3_si,
       width = 8, height = 8, dpi = 300)
ggsave("../sharepoint/Figures/Figure3_3_SI_bc.pdf",
       plot = figure3_3_si_bc,
       width = 8, height = 8, dpi = 300)

# scatter plot using share of distribution agreeing with the mean's sign
figure_agreementshare_adm0 <- function(data = NULL, data_global = NULL, warming_level_selected = 3) {
  data %>%
    # merge in income group info from df_worldbank as well as the global results
    left_join(df_worldbank %>% rename(GID_0 = "Code"), by = "GID_0") %>%
    bind_rows(data_global %>% mutate(GID_0 = "Global",
                                     Region = "Global",
                                     `Income group` = "Global",
                                     spec = "kotz2022")) %>%
    select(warming_level, GID_0, ends_with("_agreement_meansign"), ends_with("_mean"), Region, `Income group`) %>%
    # collapse "Lower middle income" and "Upper middle income" into "Middle income"
    mutate(`Income group` =  case_when(str_detect(`Income group`, "middle income") ~ "Middle income",
                                       TRUE ~ `Income group`)) %>%
    pivot_longer(cols = c(ends_with("agreement_meansign"), ends_with("_mean"))) %>%
    mutate(variable = str_remove(name, "_agreement_meansign$|_mean$"),
           summary_stat = str_remove(name, paste0(variable, "_"))) %>%
    select(-name) %>%
    pivot_wider(names_from = "summary_stat", values_from = "value") %>%
    # select the six individual climate indicators (= excluding total impacts and total variability & extremes)
    filter(variable %in% c("imp_temp_diff", "imp_Pt_diff", "imp_Tstd_diff", "imp_Wn_diff",
                           "imp_wet_days_1_diff", "imp_vwet_days1_am_99p9_diff")) %>% 
    mutate(variable = rewrite_impact_label(variable),
           is_global = GID_0 == "Global") %>%
    filter(warming_level == warming_level_selected, GID_0 %in% c("Global", sovereign_countries_gid0)) %>%
    # Income group in the World Bank file is missing for Venezuela, which is either upper or lower middle income (source: https://publications.iadb.org/en/venezuela-still-upper-middle-income-country-estimating-gni-capita-2015-2021#:~:text=In%20the%202022%20World%20Bank,income%20(GNI)%20of%202013.)
    # therefore, we impute this value
    mutate(`Income group` = if_else(GID_0 == "VEN" & is.na(`Income group`), "Middle income", `Income group`)) %>%
    ggplot(aes(mean, agreement_meansign)) +
    geom_vline(xintercept = c(0), linetype = "dashed") +
    # add IPCC likelihood thresholds
    geom_hline(yintercept = c(0.9, 0.666), linetype = "dotted", colour = "grey") +
    geom_point(aes(colour = `Income group`, shape = `Income group`, size = `Income group`), alpha = 0.6) +
    facet_wrap(~ variable, scales = "free_x", nrow = 3, ncol =  2) +
    scale_x_continuous(labels = scales::percent,
                       breaks = scales::pretty_breaks(3)) +
    scale_colour_manual(values = c("black", scales::hue_pal()(3))) +
    scale_shape_manual(values = c(18, 16, 16, 16)) +
    scale_size_manual(values = c(3, 1.5, 1.5, 1.5)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = "Mean GDP impact", y = "Confidence in sign of mean GDP impact",
         colour = NULL, shape = NULL, size = NULL) +
    guides(shape = guide_legend(nrow = 2), size = guide_legend(nrow = 2), colour =  guide_legend(nrow = 2)) +
    theme_classic() +
    theme(legend.position = "bottom")
}

# make the chart
figure2_3 <- figure_agreementshare_adm0(df_adm0_fulldistr, df_global_fulldistr_aggr, warming_level_selected = 3)
figure2_3_bc <- figure_agreementshare_adm0(df_adm0_bc_fulldistr, df_global_bc_fulldistr_aggr, warming_level_selected = 3)


############ TAIL RISK EXPOSURE CHARTS #########################################

# write a function to subset a data frame to countries with their 5th percentile of damages below a threshold
# NOTE: function also allows to further subset based on i) mean impact exceeding a threshold or ii) signalnoise ratio exceeding a threshold
subset_risk_countries <- function(data = NULL,
                                  impact_selected = NULL,
                                  warming_level_selected = 3,
                                  perc05_threshold = -0.01,
                                  mean_threshold = NULL,
                                  signalnoise_threshold = NULL) {
  
  # check inputs
  if(!paste0(impact_selected, "_perc05") %in% names(data)) stop(paste0("Could not find the perc05 variable for argument impact_selected set to ", impact_selected))
  if(!paste0(impact_selected, "_mean") %in% names(data)) stop(paste0("Could not find the mean variable for argument impact_selected set to ", impact_selected))
  if(!paste0(impact_selected, "_signalnoise") %in% names(data)) stop(paste0("Could not find the signalnoise variable for argument impact_selected set to ", impact_selected))
  
  # subset
  data %>%
    filter(warming_level == warming_level_selected, GID_0 %in% sovereign_countries_gid0) %>%
    {if(!is.null(perc05_threshold)) filter(., !!as.name(paste0(impact_selected, "_perc05")) <= perc05_threshold) else .} %>%
    {if(!is.null(mean_threshold)) filter(., !!as.name(paste0(impact_selected, "_mean")) <= mean_threshold) else .} %>%
    {if(!is.null(signalnoise_threshold)) filter(., !!as.name(paste0(impact_selected, "_signalnoise")) <= signalnoise_threshold) else .} %>%
    mutate(var_used_for_subsetting = impact_selected) %>%
    left_join(world, by = "GID_0") %>%
    select(var_used_for_subsetting, GID_0, admin, warming_level, everything())
}

# calculate the total population in sovereign countries based on the data in the world shapefile
total_pop <- world %>% filter(GID_0 %in% df_adm0_fulldistr$GID_0 & GID_0 %in% sovereign_countries_gid0) %>%
  pull(pop_est) %>% sum()

# write a function to define the exposure data
create_pop_exposure_data <- function(data = NULL, warming_levels_used = c(1, 1.5, 2, 3, 4),
                                     perc05_thresholds_used = round(seq(-0.2, 0, by = 0.01), 2),
                                     spec_used = "kotz") {
  pmap_dfr(.l = crossing(impact_selected = names(data) %>% str_subset("_mean$") %>%
                           # omit the onein20 variable (if it exists) for which we have not calculated percentiles
                           str_subset("onein20", negate = T) %>% str_remove("_mean$"),
                         warming_level_selected = warming_levels_used,
                         perc05_threshold = perc05_thresholds_used),
           .f = ~ subset_risk_countries(data = data,
                                        impact_selected = ..1,
                                        warming_level_selected = ..2,
                                        perc05_threshold = ..3) %>%
             summarise(variable = ..1,
                       warming_level = ..2,
                       perc05_threshold = ..3,
                       # we add up the population which exceeds the thresholds - if this does not hold for ANY ADM0 country, the tibble is empty and yields an NA
                       # in this case, we replace this NA with a value of zero
                       pop_under_threshold = sum(pop_est, na.rm = T)/total_pop %>% replace_na(0)) %>%
             mutate(spec = spec_used)
  )
}

# create the exposure risk data objects for main spec & status quo
df_percrisk_fulldistr <- bind_rows(create_pop_exposure_data(df_adm0_fulldistr, spec_used = "kotz"),
                                   create_pop_exposure_data(df_adm0_fulldistr_T_Pt_Pt2zero, spec_used = "status_quo"))
# repeat for bias-corrected data
df_percrisk_fulldistr_bc <- bind_rows(create_pop_exposure_data(df_adm0_bc_fulldistr, spec_used = "kotz"),
                                      create_pop_exposure_data(df_adm0_bc_fulldistr_T_Pt_Pt2zero, spec_used = "status_quo"))

# write a function to make a line chart of global population shares of exposure to 5\% risks (= using 5th percentile of impacts)
figure_exposure_linechart <- function(data_exposure = NULL, warming_levels_selected = c(1.5, 2, 3)) {
  
  # subset the exposure chart and relabel the specification
  df_linechart_exposure <- data_exposure %>% filter(warming_level %in% warming_levels_selected, variable == "imp_total") %>%
    mutate(spec_label = factor(case_when(spec == "kotz" ~ "Impacts from all climate indicators\n(= Kotz et al., 2022)",
                                         spec == "status_quo" ~ "Impacts from ann. temp. only\n(= status quo)"),
                               levels = c("Impacts from ann. temp. only\n(= status quo)",
                                          "Impacts from all climate indicators\n(= Kotz et al., 2022)")))
  
  # make the linechart
  df_linechart_exposure %>%
    # NOTE: we flip the sign of perc05_threshold because impacts are labelled as 'damages'
    ggplot(aes(-perc05_threshold, pop_under_threshold)) +
    # one line for each warming level, different linetypes by specification
    geom_line(aes(linetype = spec_label, colour = paste0(warming_level, "\u00B0C"))) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(labels = scales::percent) +
    # we add an errorbar for the +3C warming level & 15%GDP threshold, marking the difference between the specs
    annotate("errorbar", x = 0.15, width = 0.005, colour = "black", alpha = 0.9,
             ymin = df_linechart_exposure %>% filter(perc05_threshold == -0.15, warming_level == 3, spec == "status_quo") %>% pull(pop_under_threshold),
             ymax = df_linechart_exposure %>% filter(perc05_threshold == -0.15,  warming_level == 3, spec == "kotz") %>% pull(pop_under_threshold)) +
    # via annotated text label, we explain how to read the errorbar
    annotate("text", x = 0.155, y = 0.44, size = 7.5/.pt,
             label = paste0("Increase for risk of\n15%GDP loss by\nincluding all climate\nindicators: +",
                            100*round(df_linechart_exposure %>% filter(perc05_threshold == -0.15, warming_level == 3) %>%
                                        summarise(diff = pop_under_threshold[spec == "kotz"] - pop_under_threshold[spec == "status_quo"]) %>% pull(diff), 2), "%POP"),
             hjust = 0) +
    labs(linetype = NULL, y = "Global POP share of countries with\n>= 5% chance of damages above threshold",
         x = "Threshold (GDP loss)", colour = "Global warming level") +
    guides(colour = guide_legend(byrow = T),
           linetype = guide_legend(byrow = T)) +
    scale_linetype_manual(values = c("dotted", "solid")) +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          axis.title = element_text(size = 9),
          legend.title = element_text(size = 9),
          legend.spacing.y = unit(-0.2, 'cm'),
          legend.background = element_rect(fill = "transparent", colour = "transparent"))
}

# create the charts
figure5_1 <- figure_exposure_linechart(df_percrisk_fulldistr)
figure5_1_bc <- figure_exposure_linechart(df_percrisk_fulldistr_bc)

# write a function to create the population exposure table chart
figure_exposure_table <- function(data_exposure = NULL, warming_levels_selected = c(1.5, 2, 3)) {
  data_exposure %>%
    filter(variable %in% c("imp_total"),
           warming_level %in% warming_levels_selected,
           perc05_threshold %in% c(0, -0.05, -0.1, -0.15, -0.2)) %>%
    mutate(variable = rewrite_impact_label(variable),
           perc05_threshold_label = paste0(100*round(perc05_threshold, 2), "%"),
           spec_label = factor(case_when(spec == "kotz" ~ "Impacts from all climate indicators\n(= Kotz et al., 2022)",
                                         spec == "status_quo" ~ "Impacts from ann. temp. only\n (= status quo)"),
                               levels = c("Impacts from ann. temp. only\n (= status quo)",
                                          "Impacts from all climate indicators\n(= Kotz et al., 2022)")),
           warming_level_label = paste0(warming_level, "\u00B0C")) %>%
    ggplot(aes(reorder(warming_level_label, warming_level), reorder(perc05_threshold_label, perc05_threshold))) + geom_tile(aes(fill = pop_under_threshold)) +
    facet_wrap(~ spec_label) +
    geom_text(aes(label = paste0(100*round(pop_under_threshold, 2), "%")), size = 3) +
    scale_fill_gradient2(labels = scales::percent, high = "#F8766D", low = "white") +
    guides(fill = "none") +
    labs(x = "Global warming level", y = "Threshold (GDP loss)",
         subtitle = "Selected values of global population exposed to >=5% risk (from Panel a)") +
    theme_classic() +
    theme(plot.subtitle=element_text(size = 9, face = "bold"),
          axis.title = element_text(size = 9))
}

# make the exposure table figures
figure5_2 <- figure_exposure_table(df_percrisk_fulldistr)
figure5_2_bc <- figure_exposure_table(df_percrisk_fulldistr_bc)

# create figure 5 and save out
ggarrange(figure5_1,
          figure5_2,
          heights = c(1.9, 1),
          labels = "auto",
          nrow = 2, ncol = 1)
ggsave(paste0("../sharepoint/Figures/",
              Sys.Date(), " Figure5.pdf"),
       width = 5.2, height = 6)

# repeat for bias-corrected data
ggarrange(figure5_1_bc,
          figure5_2_bc,
          heights = c(1.9, 1),
          labels = "auto",
          nrow = 2, ncol = 1)
ggsave(paste0("../sharepoint/Figures/",
              Sys.Date(), " Figure5_bc.png"),
       width = 5.2, height = 6)

# make an SI version that is more detailed
figure_exposure_table_siversion <- function(data_exposure = NULL, warming_levels_selected = c(1.5, 2, 3)) {
  data_exposure %>%
  filter(warming_level %in% warming_levels_selected) %>%
  filter(!variable %in% c("imp_nontemp", "imp_total", "imp_varextremes", "imp_temp"),
         spec == "kotz",
         perc05_threshold > -0.05) %>%
  mutate(variable = rewrite_impact_label(variable)) %>%
  ggplot(aes(as.character(warming_level), perc05_threshold )) + geom_tile(aes(fill = pop_under_threshold)) + facet_wrap(~ variable) +
  geom_text(aes(label = paste0(100*round(pop_under_threshold, 2), "%")), size = 3) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_gradient2(labels = scales::percent, high = "#F8766D", low = "white") +
  guides(fill = "none") +
  labs(x = "Global warming level (\u00B0C)", y = "Threshold (in %GDP)",
       subtitle = "% of global POP living in countries with >= 5% chance of GDP impacts worse than...") +
  theme_classic()
}

# create SI version
ggsave(paste0("../sharepoint/Figures/",
              Sys.Date(), " population_exposed_perc05_fulldistr.png"),
       plot = figure_exposure_table_siversion(df_percrisk_fulldistr),
       width = 6.35, height = 4)

# clean up
rm(list = ls() %>% str_subset("^df_percrisk"))
rm(list = ls() %>% str_subset("^figure5"))
gc()



###### COMBINE GRID CHARTS TO FINAL FIGURES

# maps of total country-level impacts and variability & extremes, as well as uncertainty levels
figure2 <- ggarrange(ggarrange(ggarrange(figure2_1 + theme(legend.position = "none"), NULL,
                                         as_ggplot(get_legend(figure2_1)), nrow = 3, heights = c(6, -1.2, 1)),
                               ggarrange(figure2_2 + theme(legend.position = "none"), NULL,
                                         as_ggplot(get_legend(figure2_2)), nrow = 3, heights = c(6, -1.2, 1)), nrow = 2, ncol = 1, align = "hv", labels = "auto"),
                     ggarrange(figure2_3, labels = "c"),
                     nrow = 1, ncol = 2, widths = c(1.3, 1))
ggsave(file.path("..", "sharepoint", "Figures",
                 paste0(Sys.Date(), " Figure3.pdf")),
       plot = figure2,
       width = 8.75, height = 5.75)

# bias-corrected version
figure2_bc <- ggarrange(ggarrange(ggarrange(figure2_1_bc + theme(legend.position = "none"), NULL,
                                         as_ggplot(get_legend(figure2_1_bc)), nrow = 3, heights = c(6, -1.2, 1)),
                               ggarrange(figure2_2_bc + theme(legend.position = "none"), NULL,
                                         as_ggplot(get_legend(figure2_2_bc)), nrow = 3, heights = c(6, -1.2, 1)), nrow = 2, ncol = 1, align = "hv", labels = "auto"),
                     ggarrange(figure2_3_bc, labels = "c"),
                     nrow = 1, ncol = 2, widths = c(1.3, 1))
ggsave(file.path("..", "sharepoint", "Figures",
                 paste0(Sys.Date(), " Figure3_bc.pdf")),
       plot = figure2_bc,
       width = 8.75, height = 5.75)

# grid figure on Kotz et al 2022 vs status quo
figure3 <- ggarrange(ggarrange(figure3_1, figure3_2 + labs(y = "Marginal effect of annual temperature \n(+1\u00B0C increase) on GDP"), nrow = 1, ncol = 2, align = "hv", labels = c("a", "b")),
                     ggarrange(figure3_3, labels = "c"),
                     heights = c(1, 1.2),
                     nrow = 2, ncol = 1)
ggsave(file.path("..", "sharepoint", "Figures",
                 paste0(Sys.Date(), " Figure4.pdf")),
       plot = figure3,
       width = 9, height = 7.5)

# bias-corrected version - NOTE: figure3_2 does not depend on CMIP6 data and hence there is no separate version under bias correction
figure3_bc <- ggarrange(ggarrange(figure3_1_bc, figure3_2 + labs(y = "Marginal effect of annual temperature \n(+1\u00B0C increase) on GDP"), nrow = 1, ncol = 2, align = "hv", labels = c("a", "b")),
                     ggarrange(figure3_3_bc, labels = "c"),
                     heights = c(1, 1.2),
                     nrow = 2, ncol = 1)
ggsave(file.path("..", "sharepoint", "Figures",
                 paste0(Sys.Date(), " Figure4_bc.png")),
       plot = figure3_bc,
       width = 9, height = 7.5)


################################################################################
######################### SUMMARY STATS ########################################
################################################################################

# number of models
unique(df_global_fulldistr$model) %>% length()
df_global_fulldistr %>% count(model, scenario) %>% nrow()

# look at summary stats for global GDP impacts
df_global_fulldistr_aggr %>% select(warming_level, starts_with("imp_total")) %>%
  pivot_longer(cols = starts_with("imp")) %>% view()
df_global_fulldistr_aggr %>% select(warming_level, starts_with("imp_temp_diff")) %>%
  pivot_longer(cols = starts_with("imp")) %>% view()
df_global_fulldistr_aggr %>% select(warming_level, starts_with("imp_Pt_diff")) %>%
  pivot_longer(cols = starts_with("imp")) %>% view()
df_global_fulldistr_aggr %>% select(warming_level, starts_with("imp_varextremes")) %>%
  pivot_longer(cols = starts_with("imp")) %>% view()
df_global_fulldistr_aggr %>% select(warming_level, starts_with("imp_vwet")) %>%
  pivot_longer(cols = starts_with("imp")) %>% view()
df_global_fulldistr_aggr %>% select(warming_level, starts_with("imp_Tstd")) %>%
  pivot_longer(cols = starts_with("imp")) %>% view()

# signal-noise ratio of global impacts at different warming levels
df_global_fulldistr_aggr %>% select(warming_level, starts_with("imp_total_signalnoise")) %>%
  pivot_longer(cols = starts_with("imp"))
df_global_fulldistr_T_Pt_Pt2zero_aggr %>% select(warming_level, starts_with("imp_total_signalnoise")) %>%
  pivot_longer(cols = starts_with("imp"))

# status-quo comparison
df_global_fulldistr_T_Pt_Pt2zero_aggr %>% select(warming_level, starts_with("imp_total")) %>%
  pivot_longer(cols = starts_with("imp")) %>% filter(warming_level == 3)
df_global_fulldistr_aggr %>% select(warming_level, starts_with("imp_total")) %>%
  pivot_longer(cols = starts_with("imp")) %>% filter(warming_level == 3)

# temperature impacts
df_global_fulldistr_aggr %>% select(warming_level, starts_with("imp_temp_diff_mean")) %>%
  pivot_longer(cols = starts_with("imp"))

# same code for country-level results
df_adm0_fulldistr %>% select(warming_level, GID_0, starts_with("imp_total_mean")) %>%
  pivot_longer(cols = starts_with("imp")) %>% view()

# find the worst impacted countries (based on MC mean) for each warming level
df_adm0_fulldistr %>% select(warming_level, GID_0, imp_total_mean, imp_temp_diff_mean) %>%
  filter(GID_0 %in% sovereign_countries_gid0) %>%
  group_by(warming_level) %>% slice_min(imp_total_mean, n = 3, with_ties = T)

df_adm0_fulldistr %>% filter(warming_level == 3, imp_wet_days_1_diff_signalnoise > 1, GID_0 %in% sovereign_countries_gid0)

## how much does the marginal effect of temperature rise by including variability & exgtremes?
# to inspect this, we simply compare the marginal effect at T_0 = 0 which equals the first two coefficients (as interactions with T_0 equal zero)
exp((get_coefficients()["TmeanD", 1] + get_coefficients()["TmeanLD", 1])) - exp((get_coefficients(model_used = "T_Pt_Pt2zero")["TmeanD", 1] + get_coefficients(model_used = "T_Pt_Pt2zero")["TmeanLD", 1]))

## who are the countries with the highest extreme precip damages in Figure 3?
df_adm0_fulldistr %>% filter(warming_level == 3) %>% slice_max(n=5, order_by = imp_vwet)

