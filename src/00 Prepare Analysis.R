# clean the environment
rm(list = ls())

# load packages
library(tidyverse)  # for general data wrangling
library(readxl)     # for importing Excel workbook data
library(xtable)     # for exporting Latex tables

# source functions to rewrite variable names (used for chart labels below)
source(file.path("src", "utils", "analyse_impacts.R"))

# NOTE: for this script the following scripts have to be executed already:
# - aggregate-025.ipynb (aggregates 0.25-degree annual CMIP6 climate indicators to ADM1 level)
# - rainfall_deviations.ipynb (calculates the standardized monthly rainfall deviation measure at the ADM1-level)
# - prepare_weights_and_shapefiles.R (calculates ADM1 region weights and prepares shapefiles used in the analysis)

# map the directory for the ADM1-level climate data
dir_adm1 <- file.path("..", "sharepoint", "Data", "Climate projections", "cmip6_ng_adm1_new")


################################################################################
########## STORE VARIABLE INFO & AVAILABILITY INTO TIBBLE OBJECTS ##############
################################################################################

# define all variables & the short name used for them in the netCDF projection files (tas = temperature, pr = precipitation)
# also define the names in the Kotz et al 2022 code
vars <- tibble(names = c("ann_mean_regrid", "dtd_var_regrid", "ann_totd_regrid",
                         "wet_days_regrid", "btstrp_pr999_regrid",
                         "new_mon_dev_regrid"),
               short = c(rep("tas", 2), rep("pr", 4)),
               names_kotz = c("Tmean", "Tstd", "Pt",
                              "wet_days_1", "vwet_days1_am_99p9", "Wn"),
               # add empty columns where models for which the variable is available are stored (separately for each RCP-SSP pair)
               models_ssp119 = list(NA),
               models_ssp126 = list(NA),
               models_ssp370 = list(NA))

# repeat for the bias corrected variables (everything identical except for the bc_ prefix in variable names)
vars_bc <- vars %>% mutate(names = paste0("bc_", names))

# check for which model-scenario we have which climate variable and store that info in the vars/vars_bc objects
check_available_models <- function(vars_object = NULL) {
  for(jj in 1:nrow(vars_object)) {
    
    # create a RegEx to extract the model name
    pattern_temp <- paste0("(?<=", vars_object$short[jj], "_).+(?=_ssp)")
    
    # identify the available files for this variable on the FTP
    filenames_temp <- file.path(dir_adm1, vars_object$names[jj]) %>% list.files()
    
    # extract all matches and store them, nested into a list, into the models columns in the vars_object data frame
    vars_object$models_ssp119[jj] <- filenames_temp %>% str_subset("ssp119") %>% str_extract(pattern_temp) %>% list()
    
    vars_object$models_ssp126[jj] <- filenames_temp %>% str_subset("ssp126") %>% str_extract(pattern_temp) %>% list()
    
    vars_object$models_ssp370[jj] <- filenames_temp %>% str_subset("ssp370") %>% str_extract(pattern_temp) %>% list()
    
    # clean up
    rm(filenames_temp, pattern_temp)
  }
  return(vars_object)
}
vars <- check_available_models(vars)
vars_bc <- check_available_models(vars_bc)

# read in the CMIP6 model data from Hausfather et al 2022
df_hausfather <- read_excel(file.path("data", "input", "HausfatherEtAl2022 - SI Data - CMIP6 Models.xlsx"),
                 sheet = "ECSTCR", range = "A3:H61")

# save out the vars data frame
saveRDS(vars, file.path("data", "vars.rds"))
saveRDS(vars_bc, file.path("data", "vars_bc.rds"))

# chart for each model-SSP pair which indicators are (not) available and denote which CMIP6 models are missing altogether
df_for_modelchart <- vars %>% 
  pivot_longer(cols = starts_with("models_ssp"), names_to = "scenario", values_to = "models") %>%
  mutate(scenario = str_extract(scenario, "ssp[:digit:]+$")) %>%
  unnest(models) %>%
  arrange(models, scenario) %>%
  group_by(models, scenario) %>%
  summarise(Tmean = "Tmean" %in% names_kotz,
            Tstd = "Tstd" %in% names_kotz,
            Pt = "Pt" %in% names_kotz,
            wet_days_1 = "wet_days_1" %in% names_kotz,
            vwet_days1_am_99p9 = "vwet_days1_am_99p9" %in% names_kotz,
            Wn = "Wn" %in% names_kotz,
            ) %>% ungroup() %>%
  pivot_longer(cols = c(Tmean, Tstd, Pt, wet_days_1, vwet_days1_am_99p9, Wn),
               names_to = "var", values_to = "is_provided") %>%
  mutate(modelscen = paste0(models, " - ", scenario),
         var = rewrite_climate_label(var))

# create the chart (and annotate info about exclusion of CAMS-CSM1-0 due to inconcistent realizations by indicator)  
df_for_modelchart %>%
ggplot(aes(modelscen, var)) + geom_tile(aes(fill = is_provided), colour = "white", size = 0.8) +
  coord_flip() +
  labs(subtitle = paste0("CAMS-CSM1-0 excluded due to inconsistent realizations \n(r1 for ann & monthly values, r2 for daily)",
                         "\nCMIP6 models not covered:\n ", paste0(df_hausfather$`Model Name`[!df_hausfather$`Model Name` %in% df_for_modelchart$models], collapse = ", ") %>%
                           str_replace("E3SM-1-0", "\nE3SM-1-0") %>%
                           str_replace("GISS-E2", "\nGISS_E2")),
       fill = "Climate indicator available?",
       x = "CMIP6 model & SSP") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(colour = ifelse(df_for_modelchart %>% group_by(modelscen) %>% summarise(is_all_provided = mean(is_provided) == 1) %>%
                                                                            pull(is_all_provided), 'black', 'darkred')),
        plot.subtitle=element_text(size=8, face = "italic"))
ggsave(file.path("..", "sharepoint", "Figures", "cmip6", paste0(Sys.Date(), " data_availability_by_CMIP6-scenario.png")), width = 6, height = 12)


################################################################################
########## PREPARE GLOBAL WARMING LEVELS DATA ##################################
################################################################################

###  create a tibble storing all model-scenario-ensemble combination we have available somewhere
# first, detect all netcdf files
# NOTE: when executing the full script, this step can sometimes file and lead to errors further below
#       if the object is empty after (failing to) execute the whole script, run the script line by line
nc_files_adm1 <- dir_adm1 %>% list.files(recursive = TRUE) %>% str_subset("\\.nc$") %>%
  unique()

# then extract the indicative part of the filename and the relevant information using regex
df_models <- tibble(string = map_chr(nc_files_adm1, ~ str_extract(.x, "(?<=_pr_|_tas_).+(?=_g025\\.nc)"))) %>%
  mutate(model = string %>% str_extract(".+(?=_ssp)"),
         scenario = string %>% str_extract("(?<=_)ssp[:digit:][:digit:][:digit:](?=_r)"),
         ensemble = string %>% str_extract("(?<=ssp[:digit:][:digit:][:digit:]_).+$")
  ) %>%
  distinct() %>%
  # drop the raw string
  select(-string)

# confirm that for each model, we have only one model-scenario-ensemble combination
df_models %>% group_by(model, scenario) %>% mutate(n = n()) %>% filter(n > 1)
# -> actually, for CAMS-CSM1-0 different climate indicators come from different ensemble id runs
# -> therefore, we remove the model below
df_models <- df_models %>% filter(!(model == "CAMS-CSM1-0" & ensemble != "r1i1p1f1"))
df_models %>%  group_by(model, scenario) %>% mutate(n = n()) %>% filter(n > 1)

# identify only those model-scenario-ensemble combinations for which we have all the required data
create_df_modelscen <- function(vars_object = NULL) {
  # dynamically create a tibble with all model-scenario pairs for which we have ALL climate indicators using the vars object
  # step 1: extract all scenario names from the vars object column names
  out <- tibble(scenario = names(vars_object) %>% str_subset("ssp[:digit:][:digit:][:digit:]") %>%
                  str_extract("ssp[:digit:][:digit:][:digit:]") %>%
                  unique())
  
  # store all available models for the respective scenario in a new column as a list, then unnest into LONG format
  out$model <- NA
  for(jj in 1:nrow(out)) {
    out$model[jj] <- vars_object %>% 
      pull(paste0("models_", out$scenario[jj])) %>% reduce(intersect) %>% list()
  }
  out <- unnest(out, cols = "model")
  
  # remove the CAMS-CSM1-0 model for which we have data only for inconsistent ensemble ids (r1 v r2) and merge
  # in the ensemble information
  out <- out %>% filter(model != "CAMS-CSM1-0") %>% left_join(df_models, by = c("model", "scenario"))
  
  return(out)
}
df_modelscen <- create_df_modelscen(vars)
df_modelscen_bc <- create_df_modelscen(vars_bc)

# remove the df_models object to avoid confusions
rm(df_models)

# ensure that we get the same model-scenario pairings for raw and bias-corrected data
if(!identical(df_modelscen, df_modelscen_bc)) stop("Available model-scenario pairs are not identical for raw and bias-corrected CMIP6 data. Please inspect")

# read in the global warming levels and transform to long format
df_gwl_nested <- read_excel(file.path("data", "input", "cmip6_wl_for_tas_annmean_inclGWL4.xlsx"), sheet = "final") %>%
  mutate(model = str_extract(modelscenarioensemble, "^.+(?=_ssp[:digit:][:digit:][:digit:])"),
         scenario = str_extract(modelscenarioensemble, "(?<=_)ssp[:digit:][:digit:][:digit:](?=_r[:digit:]i[:digit:])"),
         ensemble = str_extract(modelscenarioensemble, "(?<=_ssp[:digit:][:digit:][:digit:]_)r[:digit:]+i[:digit:]+p[:digit:]+f[:digit:]+")) %>%
  select(-modelscenarioensemble) %>%
  pivot_longer(cols = c(starts_with("0"), starts_with("1"), starts_with("2"), starts_with("3"), starts_with("4")),
               names_to = "warming_level", values_to = "year") %>%
  group_by(model, scenario, ensemble, warming_level) %>% summarise(beg = year[beg_or_end == "beg"],
                                                                   end = year[beg_or_end == "end"]) %>% ungroup() %>%
  filter(!(is.na(beg) & is.na(end))) %>%
  mutate(gwl_years = NA, gwl_year_no = NA)

# for each model-scenario & warming level, create a nested vector with all the years in the respective GWL window
for(jj in 1:nrow(df_gwl_nested)) {
  df_gwl_nested$gwl_years[jj] <- list(df_gwl_nested$beg[jj]:df_gwl_nested$end[jj])
  df_gwl_nested$gwl_year_no[jj] <- length(df_gwl_nested$beg[jj]:df_gwl_nested$end[jj])
}

# create an unnested version
df_gwl <- df_gwl_nested %>% select(model, scenario, ensemble, warming_level, year = "gwl_years") %>%
  unnest(year) %>% mutate(warming_level = as.numeric(warming_level))

# check that for with available climate data we do not have GWL info (= NA values after a left_join)
if((df_modelscen %>% left_join(df_gwl %>% select(-year) %>% distinct(), by = c("model", "scenario", "ensemble")) %>%
  filter(is.na(warming_level)) %>% nrow()) > 0) stop("For some model-scenarios in df_modelscen, there is no GWL data available")
if((df_modelscen_bc %>% left_join(df_gwl %>% select(-year) %>% distinct(), by = c("model", "scenario", "ensemble")) %>%
    filter(is.na(warming_level)) %>% nrow()) > 0) stop("For some model-scenarios in df_modelscen_bc, there is no GWL data available")

# check that we have exactly one ensemble run per model-scenario pair
if((bind_rows(df_modelscen, df_modelscen_bc) %>% distinct() %>%
   count(model, scenario) %>% pull(n) %>% max()) > 1) stop("For at least one model in df_modelscen, we have multiple ensembles. Please inspect")

# merge df_gwl into the tibble with the available models
# NOTE: this has the advantage that we can simply merge data on model-scenario pairs and do not have to carry around the ensemble id in all our results
df_gwl <- df_modelscen %>% distinct() %>% left_join(df_gwl, by = c("model", "scenario", "ensemble"))

# save this out as an RDS file, together with the df_modelscen and df_modelscen_bc objects
saveRDS(df_gwl, file.path("data", "global_warming_levels.rds"))
saveRDS(df_modelscen, file.path("data", "df_modelscen.rds"))
saveRDS(df_modelscen_bc, file.path("data", "df_modelscen_bc.rds"))

# visualize the warming levels
unique(df_gwl_nested$warming_level[!is.na(df_gwl_nested$warming_level)]) %>%
  map(.f = function(x){
    df_gwl_nested %>% filter(warming_level == x) %>%
        mutate(modelscen = paste0(model, " - ", scenario)) %>%
        left_join(df_for_modelchart %>% group_by(modelscen) %>% summarise(all_provided = mean(is_provided) == 1) %>% filter(all_provided),
                  by = "modelscen") %>%
        filter(all_provided) %>%
        ggplot(aes(model)) +
        geom_errorbar(aes(ymin = beg, ymax = end, colour = scenario)) +
      {if(x == 0.84) geom_hline(yintercept = c(1979, 2019), linetype = "dashed")} +
        labs(y = NULL, x = NULL, subtitle = paste0("GWL: +", x, "C"), colour = NULL) +
        coord_flip() +
      theme(legend.position = "bottom")
      
        ggsave(file.path("..", "sharepoint", "Figures", "cmip6", paste0(Sys.Date(), " gwlperiod_by_gcm_gwl", x, ".png")),
               width = 4, height = 10)
  }
)

# write out the models/scenarios/ensembles used in a Latex table
df_gwl %>% left_join(df_gwl_nested %>% mutate(window = paste0(beg, "-", end)) %>% select(model, scenario, ensemble, warming_level, window) %>%
                       # ensure that GWL is a numeric
                       mutate(warming_level = as.numeric(warming_level)),
                     by = c("model", "scenario", "ensemble", "warming_level")) %>%
  select(-year) %>% distinct() %>%
  pivot_wider(names_from = "warming_level", values_from = "window") %>%
  mutate_at(vars("1", "1.5", "2", "3", "4"), ~ replace_na(.x, "")) %>%
  # add degC and + symbols to the warming levels in the column names
  setNames(names(.) %>% str_replace("(?<=[:digit:])$", "\u00B0C") %>%
             str_replace("^(?=[:digit:])", "\\+")) %>%
  # rename and rearrange columns
  select(GCM = "model", Scenario = "scenario", Ensemble = "ensemble", everything()) %>% arrange(GCM, Scenario) %>%
  # export to a .tex file
  xtable(label = "tab:gcm_ssp_list",
                 digits = 0,
                 caption = "Models, scenarios and ensembles used for climatic projections") %>%
  print(file = file.path("..", "sharepoint", "Tables", paste0(Sys.Date(), " GCM_SSP_used_list.tex")),
        size="\\tiny")


################################################################################
###################### PREPARE SSP DATA ########################################
################################################################################

# read in country-level SSP GDP
df_ssp_raw <- read_excel(file.path("data", "input", "220803 SSP-Database GDP PPP.xlsx"),
                     sheet = "data",
                     range = "A1:X921") %>%
  pivot_longer(cols = starts_with("2"), names_to = "year", values_to = "gdp") %>%
  mutate(is_historical = FALSE)

# read in the historic data starting in row 922
df_ssp_historic <- read_excel(file.path("data", "input", "220803 SSP-Database GDP PPP.xlsx"),
                       sheet = "data",
                       range = "A922:Q1103") %>%
              pivot_longer(cols = c(starts_with("2"), starts_with("1")), names_to = "year", values_to = "gdp") %>%
              rename(Scenario = "Scenario (History)") %>%
              mutate(is_historical = TRUE)

# combine the two
df_ssp <- bind_rows(df_ssp_raw, map_dfr(paste0("SSP", 1:5), ~ df_ssp_historic %>% mutate(Scenario = .x))) %>%
  # make year a numeric and harmonize column names
          rename(GID_0 = "Region", ssp = "Scenario") %>%
          mutate(year = as.integer(year)) %>%
          # filter out 2010 scenario values since we have actuals to avoid two values for the same year
          filter(!(year == 2010 & !is_historical)) %>%
          # create a log version of GDP for log-linear interpolation
          mutate(gdp_ln = log(gdp))

# check that we have exactly one value per country, year and scenario
if((df_ssp %>% group_by(ssp, GID_0, year) %>% summarise(n = n()) %>% filter(n != 1) %>% nrow()) > 0) stop("There are year-country-SSP pairings with zero or > 1 values") 

# create the data set to be filled
# step 1: create a data frame with all possible combinations of SSPs, ADM0-level countries and years
df_ssp_out <- crossing(ssp = unique(df_ssp$ssp), GID_0 = unique(df_ssp$GID_0),
                       year = min(df_ssp$year):max(df_ssp$year)) %>%
  # merge in the SSP data
  left_join(df_ssp %>% select(ssp, GID_0, year, gdp, is_historical, gdp_ln),
            by = c("ssp", "GID_0", "year")) %>%
  # for each country and SSP, order by year and perform linear interpolation on the log-GDP via na.approx()
  group_by(ssp, GID_0) %>% arrange(ssp, GID_0, year) %>%
  mutate(gdp_ln_interpolated = zoo::na.approx(gdp_ln, maxgap = 4)) %>%
  # convert back to normal scale, then discard rows with missing GDP information after interpolation
  mutate(gdp_annual = exp(gdp_ln_interpolated)) %>%
  filter(!is.na(gdp_annual))

# check summary stats
summary(df_ssp_out)

# save this out as an RDS file
saveRDS(df_ssp_out %>% select(ssp, GID_0, year, gdp_annual) %>% ungroup(), file.path("data", "df_ssp.rds"))

