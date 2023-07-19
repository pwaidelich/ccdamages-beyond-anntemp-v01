# clean out the environment
rm(list = ls())

# load packages
library(sf)
library(tidyverse)       # for general data wrangling & plotting
library(furrr)           # for parallelizing
library(feather)         # for reading in and writing out large files faster

# load helper functions
source(file.path("src", "utils", "extract_climate_indices.R"))
source(file.path("src", "utils", "project_impacts.R"))

# load helper objects
vars <- readRDS(file.path("data", "vars.rds"))
vars_bc <- readRDS(file.path("data", "vars_bc.rds"))

# GDP weights for ADM1-to-ADM0 aggregation
gdp_weights <- read_csv(file.path("data", "adm1-weights.csv"))  %>%
  group_by(GID_0) %>% mutate(adm1_weight = GDP/sum(GDP)) %>% ungroup() %>%
  select(GID_0, GID_1, GDP, adm1_weight)

# draw coefficients
# NOTE: get_coefficients() uses a default seed for Monte Carlo draws from multi-variate Gaussian, so we do not need to fix the seed here)
coefs_point <- get_coefficients()
coefs_mc <- get_coefficients(n = 1000, draw_from_vcov = T)
coefs_mc_T_Pt_Pt2zero <- get_coefficients(model_used = "T_Pt_Pt2zero", draw_from_vcov = T, n = 1000)

# load global warming level data for model-scenario-specific baseline periods based on global warming level (defaulting to 0.84 in src/utils/project_impacts() )
df_gwl <- readRDS(file.path("data", "global_warming_levels.rds"))

# map the Sharepoint directory where climate projection files are stored
dir_climate_data_adm1 <- file.path("..", "sharepoint", "Data", "Climate projections", "cmip6_ng_adm1_new")

# map where GDP projection files (ADM0-level) will be stored
dir_output_data <- file.path("/net", "cfc", "landclim1", "pwaidelich")

# load the tibble objects storing model-scenarios for which we have all 6 climate indicators available
df_modelscen <- readRDS(file.path("data", "df_modelscen.rds"))

# initiate parallel processing
# NOTE: adjust this to your system for optimal trade-off between memory load & computing power
# cl <- parallel::makeCluster(2)
# plan(cluster, workers = cl, gc = T)
#plan(sequential, gc = T)
plan(multisession(workers = 36, gc = TRUE))


################################################################################
############### CALCULATE CEILING FOR BIAS-CORRECTED PRECIP DEVIATION ##########
################################################################################

# extract the maximum monthly precip deviation for raw CMIP6 data (used to cap bias-corrected values)
ceiling_Wn_bc <- map2_dbl(.x = df_modelscen$model, .y = df_modelscen$scenario,
                          .f = ~ extract_adm1_var(varname = "new_mon_dev_regrid",
                                                  vars_object = vars,
                                                  model = .x,
                                                  scenario = .y,
                                                  dir = dir_climate_data_adm1) %>%
                            pull(Wn) %>% max()) %>% max()

# determine all observations affected by the ceiling
adm1regions_capped_by_ceiling <- map2_dfr(.x = df_modelscen$model, .y = df_modelscen$scenario,
                                          .f = ~ extract_adm1_var(varname = "bc_new_mon_dev_regrid",
                                                                  vars_object = vars_bc,
                                                                  model = .x,
                                                                  scenario = .y,
                                                                  dir = dir_climate_data_adm1)) %>%
  filter(Wn > ceiling_Wn_bc) %>%
  # reorder the data
  select(model, scenario, GID_1, year, Wn) %>% arrange(model, scenario, GID_1, year)
saveRDS(adm1regions_capped_by_ceiling, "adm1regions_capped_by_ceiling.rds")


################################################################################
############### PROJECTIONS FOR NON-BIAS-CORRECTED DATA ########################
################################################################################

# run produce_gdp_projections in parallel using only the point estimates for damage functions - this takes ~ 3min for 60 cores
future_map2(.x = df_modelscen$model, .y = df_modelscen$scenario,
            .f = ~ produce_gdp_projections(model = .x,
                                           scenario = .y,
                                           dir_climate_adm1 = dir_climate_data_adm1,
                                           dir_out_mc = file.path(dir_output_data, "adm0_new", "pointestimates"),
                                           vars_object = vars,
                                           coefs_matrix_used = coefs_point,
                                           adm1_weights_used = gdp_weights,
                                           df_gwl_baseline = df_gwl,
                                           use_gwl_baseline = T,
                                           overwrite_previous_files = TRUE,
                                           ceiling_precipdeviation = NULL),
            .options = furrr_options(seed = NULL),
            .progress = TRUE)

# repeat for the Monte Carlo draws of the damage function coefficients
future_map2(.x = df_modelscen$model, .y = df_modelscen$scenario,
            .f = ~ produce_gdp_projections(model = .x,
                                           scenario = .y,
                                           dir_climate_adm1 = dir_climate_data_adm1,
                                           dir_out_mc = file.path(dir_output_data, "adm0_new", "montecarlo"),
                                           vars_object = vars,
                                           coefs_matrix_used = coefs_mc,
                                           adm1_weights_used = gdp_weights,
                                           df_gwl_baseline = df_gwl,
                                           use_gwl_baseline = T,
                                           overwrite_previous_files = F,
                                           ceiling_precipdeviation = NULL),
            .options = furrr_options(seed = NULL),
            .progress = TRUE)

# repeat for the Monte Carlo draws of the damage function coefficients of the "status quo" approach
future_map2(.x = df_modelscen$model, .y = df_modelscen$scenario,
            .f = ~ produce_gdp_projections(model = .x,
                                           scenario = .y,
                                           dir_climate_adm1 = dir_climate_data_adm1,
                                           dir_out_mc = file.path(dir_output_data, "adm0_new", "T_Pt_Pt2zero", "montecarlo"),
                                           vars_object = vars,
                                           coefs_matrix_used = coefs_mc_T_Pt_Pt2zero,
                                           adm1_weights_used = gdp_weights,
                                           df_gwl_baseline = df_gwl,
                                           use_gwl_baseline = T,
                                           overwrite_previous_files = TRUE,
                                           ceiling_precipdeviation = NULL),
            .options = furrr_options(seed = NULL),
            .progress = TRUE)

# create the ADM1-level projections using point estimates
df_adm1 <- future_map2_dfr(.x = df_modelscen$model, .y = df_modelscen$scenario,
            .f = ~  extract_adm1_allvars(model = .x, scenario = .y, dir = dir_climate_data_adm1,
                                                    vars_object = vars) %>%
              project_impacts(coefs = coefs_point[, 1], return_helper_columns = T,
                              df_gwl_baseline = df_gwl,
                              use_gwl_baseline = T),
            .options = furrr_options(seed = NULL),
            .progress = TRUE) %>%
  select(-c(imp_temp_absolutescale_base, Tmean_fd, dimp_temp))
write_feather(df_adm1, file.path(dir_output_data, "adm1_new", "pointestimates", "df_adm1.feather"))


################################################################################
############### PROJECTIONS FOR BIAS-CORRECTED DATA ############################
################################################################################

# run produce_gdp_projections in parallel using only the point estimates for damage functions - this takes ~ 1min for 60 cores
future(future_map2(.x = df_modelscen$model, .y = df_modelscen$scenario,
            .f = ~ produce_gdp_projections(model = .x,
                                           scenario = .y,
                                           dir_climate_adm1 = dir_climate_data_adm1,
                                           dir_out_mc = file.path(dir_output_data, "adm0_bc_new", "pointestimates"),
                                           vars_object = vars_bc,
                                           coefs_matrix_used = coefs_point,
                                           adm1_weights_used = gdp_weights,
                                           df_gwl_baseline = df_gwl,
                                           use_gwl_baseline = T,
                                           overwrite_previous_files = TRUE,
                                           ceiling_precipdeviation = ceiling_Wn_bc),
            .options = furrr_options(seed = NULL)))

# repeat for the Monte Carlo draws of the damage function coefficients
future_map2(.x = df_modelscen$model, .y = df_modelscen$scenario,
            .f = ~ produce_gdp_projections(model = .x,
                                           scenario = .y,
                                           dir_climate_adm1 = dir_climate_data_adm1,
                                           dir_out_mc = file.path(dir_output_data, "adm0_bc_new", "montecarlo"),
                                           vars_object = vars_bc,
                                           coefs_matrix_used = coefs_mc,
                                           adm1_weights_used = gdp_weights,
                                           df_gwl_baseline = df_gwl,
                                           use_gwl_baseline = T,
                                           overwrite_previous_files = T,
                                           ceiling_precipdeviation = ceiling_Wn_bc),
            .options = furrr_options(seed = NULL),
            .progress = TRUE)

# repeat for the Monte Carlo draws of the damage function coefficients of the "status quo" approach
future_map2(.x = df_modelscen$model, .y = df_modelscen$scenario,
            .f = ~ produce_gdp_projections(model = .x,
                                           scenario = .y,
                                           dir_climate_adm1 = dir_climate_data_adm1,
                                           dir_out_mc = file.path(dir_output_data, "adm0_bc_new", "T_Pt_Pt2zero", "montecarlo"),
                                           vars_object = vars_bc,
                                           coefs_matrix_used = coefs_mc_T_Pt_Pt2zero,
                                           adm1_weights_used = gdp_weights,
                                           df_gwl_baseline = df_gwl,
                                           use_gwl_baseline = T,
                                           overwrite_previous_files = F,
                                           ceiling_precipdeviation = ceiling_Wn_bc),
            .options = furrr_options(seed = NULL),
            .progress = TRUE)

# create the ADM1-level projections using point estimates
df_adm1_bc <- future_map2_dfr(.x = df_modelscen$model, .y = df_modelscen$scenario,
                           .f = ~  extract_adm1_allvars(model = .x, scenario = .y, dir = dir_climate_data_adm1,
                                                        vars_object = vars_bc) %>%
                             mutate(Wn = if_else(Wn > ceiling_Wn_bc, ceiling_Wn_bc, Wn)) %>%
                             project_impacts(coefs = coefs_point[, 1], return_helper_columns = T,
                                             df_gwl_baseline = df_gwl,
                                             use_gwl_baseline = T),
                           .options = furrr_options(seed = NULL),
                           .progress = TRUE) %>%
  select(-c(imp_temp_absolutescale_base, Tmean_fd, dimp_temp))
write_feather(df_adm1_bc, file.path(dir_output_data, "adm1_bc_new", "pointestimates", "df_adm1_bc.feather"))
