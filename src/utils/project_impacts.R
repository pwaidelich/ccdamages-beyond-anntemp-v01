# write a function to load the coefficients for the dose-response function
get_coefficients <- function(draw_from_vcov = FALSE, cluster_national = TRUE, seed = 2022,
                             n = 1, model_used = "main") {
  
  if(!is.logical(cluster_national)) stop("Argument cluster_national must be TRUE or FALSE")
  if(n != 1 & !draw_from_vcov) stop("Argument n can only be unequal to 1 if draw_from_vcov is TRUE")
  if(!model_used %in% c("main", "T_Pt", "T_Pt_Tstd", "T_Pt_Tstd_wetdays", "T_Pt_Tstd_wetdays_Wn", "T_Pt_Pt2zero", "T_Pt_Tstd_Pt2zero" )) {
    stop("model_used must be either main, T_Pt, T_Pt_Tstd, T_Pt_Tstd_wetdays, T_Pt_Tstd_wetdays_Wn, T_Pt_Pt2zero or T_Pt_Tstd_Pt2zero")
  }
  
  # load the regression coefficients - if specification is not main, we load the list of coefficient vectors and extract the specified element by name (= model_used)
  if(model_used == "main") {
    coefs <- readRDS(file.path("data", "kotzetal2022_coefs_main.rds"))
  } else {
    coefs <- readRDS(file.path("data", "coef_list.rds"))[[model_used]]
  }
  
  # ensure that interactions are shown with .
  names(coefs) <- gsub(":", ".", names(coefs))
  
  # return the coefficients if no vcov-drawing is required
  if(!draw_from_vcov) return(coefs %>% as.matrix())
  
  # otherwise load the correct, draw the required sample and return it
  if(draw_from_vcov) {
    if(cluster_national) {
      if(model_used == "main") {
        vcov <- readRDS(file.path("data", "kotzetal2022_cov_clustnational_main.rds"))
      } else {
        vcov <- readRDS(file.path("data", "cov_iso_list.rds"))[[model_used]]
      }
    } else {
      if(model_used == "main") {
        vcov <- readRDS(file.path("data", "kotzetal2022_cov_clustregional_main.rds"))
      } else {
        vcov <- readRDS(file.path("data", "cov_id_list.rds"))[[model_used]]
      }
    }
    
    # round the matrix to avoid asymmetry due to precision, then draw
    vcov <- round(vcov, 18)
    set.seed(seed)
    out <- MASS::mvrnorm(n = n, mu = coefs %>% setNames(NULL), Sigma = vcov) %>% t()
    
    # NOTE: due to imprecision, mvrnorm() can return extremely small non-zero values for variables with coefficient & covariance set to zero
    # to address this, we manually overwrite these values with a clean zero
    vars_excluded <- names(coefs)[coefs == 0]
    out[vars_excluded, ] <- 0
    
    return(out)
  }
}

# function to calculate impacts in a given year relative to the baseline period
project_impacts <- function(data = NULL, coefs = NULL, df_gwl_baseline = NULL, baseline_start_manual = 1979, baseline_end_manual = 2019,
                            return_helper_columns = FALSE, use_gwl_baseline = TRUE, gwl_baseline_level = 0.84) {
  ## check inputs
  # names of the elements in the coefs object
  if(mean(c("Tstd", "TmeanD", "TmeanLD", "Pt", "Pt2", "Wn", "Wn_2", "wet_days_1", "wet_days_1_2","vwet_days1_am_99p9",
           "TmeanD.Tmean",  "TmeanLD.TmeanL", "Tmean.vwet_days1_am_99p9") %in% names(coefs)) != 1) stop("Names of elements in coefs do not match the required variables in Kotz et al 2022 regression")
  # column names in data
  if(mean(c("model", "scenario", "GID_1", "year", "Tmean", "Tstd", "Pt", "wet_days_1", "vwet_days1_am_99p9", "Wn") %in% names(data)) != 1) stop("Argument data must have columns named model, scenario, GID_1, year, Tmean, Tstd, Pt, wet_days_1, vwet_days1_am_99p9, Wn. Please inspect")
  # column names in df_gwl_baseline
  if(mean(c("model", "scenario", "warming_level", "year") %in% names(df_gwl_baseline)) != 1) stop("Argument df_gwl_baseline must feature columns called model, scenario, warming_level, year")
  # argument types
  if(!"numeric" %in% class(coefs)) stop("Argument coefs must be of type numeric")
  if(!is.logical(return_helper_columns)) stop("Argument return_helper_columns must of type logical")
  # does data contain data for exactly one model-scenario pair?
  if((data %>% select(model, scenario) %>% distinct() %>% nrow()) != 1) stop("Argument data features more than one model-scenario pair. project_impacts() only works for data for a single model-scenario pair")
  # is the GWL baseline value a single numeric
  if(!gwl_baseline_level %in% df_gwl_baseline$warming_level) stop("Argument gwl_baseline_level must be a numeric that features in the warming_level column of df_gwl_baseline")
  
  # depending on user input, we extract baseline period definition either from user input or from the baseline global warming level window in df_gwl_baseline  
  if(!use_gwl_baseline) {
    
    baseline_start_used <- baseline_start_manual
    baseline_end_used <- baseline_end_manual
    
  } else {

    baseline_gwl_window <- data %>% select(model, scenario) %>% distinct() %>%
                            left_join(df_gwl_baseline %>% filter(warming_level == gwl_baseline_level), by = c("model", "scenario")) %>%
                            summarise(start = min(year), end = max(year))
    
    baseline_start_used <- baseline_gwl_window$start
    baseline_end_used <- baseline_gwl_window$end
    
  }
  
  if(!baseline_start_used %in% data$year | !baseline_end_used %in% data$year) stop("Argument data does not feature the specified baseline start/end years in column year")
  
  # identify baseline years via a dummy column
  data <- data %>% mutate(is_baseline = year %in% baseline_start_used:baseline_end_used)
  
  # calculate the impacts in the baseline period
  df_base <- data %>% filter(is_baseline) %>%
    # calculate the mean annual temperature across the baseline period (used to calculate extreme precip damages to avoid implicit adaptation through warming temperatures)
    group_by(GID_1) %>% mutate(Tmean_base = mean(Tmean)) %>% ungroup() %>%
    mutate(imp_Pt = Pt*coefs["Pt"] + (Pt)^2 * coefs["Pt2"],
           imp_Wn = Wn*coefs["Wn"] + (Wn)^2 * coefs["Wn_2"],
           imp_Tstd = Tstd * coefs["Tstd"],
           imp_wet_days_1 = wet_days_1 * coefs["wet_days_1"] + (wet_days_1)^2 * coefs["wet_days_1_2"],
           imp_vwet_days1_am_99p9 = vwet_days1_am_99p9 * coefs["vwet_days1_am_99p9"] + vwet_days1_am_99p9*Tmean_base*coefs["Tmean.vwet_days1_am_99p9"]
    )
  
  # take average (NOTE: we convert to %GDP scale for this to avoid adding implicit risk premia by averaging the log impacts)
  df_base_mean <- df_base %>% group_by(GID_1) %>% summarise_at(vars(starts_with("imp")), ~ log(mean(exp(.x)))) %>%
    # add the Tmean_base (which is already averaged, see above)
    left_join(df_base %>% select(GID_1, Tmean_base) %>% distinct(), by = "GID_1")
  
  # calculate impacts for all years
  df_out <- data %>%
    # discard all years prior to the baseline period as we do not use them anyways to save runtime
    filter(year >= baseline_start_used) %>%
    # select variables of interest
    select(model, scenario, GID_1, year, is_baseline,
                                   Tmean, Tstd, Pt, wet_days_1, Wn, vwet_days1_am_99p9) %>%
    # merge in the baseline average temperature (required to project out extreme precipitation damages)
    left_join(df_base_mean %>% select(Tmean_base, GID_1), by = "GID_1") %>%
    # calculate current impacts (in log-scale)
    mutate(imp_Pt = Pt*coefs["Pt"] + (Pt)^2 * coefs["Pt2"],
           imp_Wn = Wn*coefs["Wn"] + (Wn)^2 * coefs["Wn_2"],
           imp_Tstd = Tstd * coefs["Tstd"],
           imp_wet_days_1 = wet_days_1 * coefs["wet_days_1"] + (wet_days_1)^2 * coefs["wet_days_1_2"],
           imp_vwet_days1_am_99p9 = vwet_days1_am_99p9 * coefs["vwet_days1_am_99p9"]  + vwet_days1_am_99p9 * Tmean_base * coefs["Tmean.vwet_days1_am_99p9"]
    ) %>%
    # now merge in the average impacts during the baseline period (in log-scale) and denote them via a "_base" suffix
    left_join(df_base_mean %>% select(GID_1, starts_with("imp")), by = "GID_1", suffix = c("", "_base")) %>%
    # calculate differences (still in log-scale)
    mutate(imp_Tstd_diff = imp_Tstd - imp_Tstd_base,
           imp_Pt_diff = imp_Pt - imp_Pt_base,
           imp_Wn_diff = imp_Wn - imp_Wn_base,
           imp_wet_days_1_diff = imp_wet_days_1 - imp_wet_days_1_base,
           imp_vwet_days1_am_99p9_diff = imp_vwet_days1_am_99p9 - imp_vwet_days1_am_99p9_base
    ) %>%
    # temperature damages
    group_by(GID_1) %>% arrange(year, .by_group=T) %>%
    mutate(Tmean_fd = c(NA, diff(Tmean))) %>%
    mutate(dimp_temp = coefs["TmeanD"]*Tmean_fd + coefs["TmeanLD"]*lag(Tmean_fd) +
             coefs["TmeanD.Tmean"]*Tmean_fd*Tmean + coefs["TmeanLD.TmeanL"]*lag(Tmean_fd)*lag(Tmean)) %>%
    # impact is the cumulative sum for all years since the baseline start year
    mutate(imp_temp = cumsum(ifelse(is.na(dimp_temp), 0, dimp_temp))) %>%
    # impact over baseline is imp_temp minus the average over the baseline period (in absolute scale, not log scale)
    mutate(imp_temp_absolutescale_base = if_else(is_baseline, exp(imp_temp), NA_real_),
           imp_temp_diff = imp_temp - log(mean(imp_temp_absolutescale_base, na.rm = T))) %>%
    # create total impacts and varextremes impact
    mutate(imp_total = imp_temp_diff + imp_Pt_diff + imp_Tstd_diff + imp_Wn_diff + imp_wet_days_1_diff + imp_vwet_days1_am_99p9_diff,
           imp_varextremes = imp_Tstd_diff + imp_Wn_diff + imp_wet_days_1_diff + imp_vwet_days1_am_99p9_diff) %>%
    # convert all impacts also from log-scale to %GDP, so we can average and aggregate them up
    mutate_at(vars(starts_with("imp")), ~ exp(.x) - 1) %>%
    ungroup()
  
  # write out results (incl. helper columns if this is specified by the user)
  if(return_helper_columns) {
    return(df_out)
  } else {
    return(df_out %>% select(model, scenario, GID_1, year,
                             imp_total, imp_varextremes, imp_temp_diff, imp_Tstd_diff, imp_Pt_diff, imp_Wn_diff,
                             imp_wet_days_1_diff, imp_vwet_days1_am_99p9_diff))
  }
}

# wrapper to calculate impacts across Monte Carlo draws for the damage function
project_impacts_mc <- function(data = NULL, coefs_mc_matrix = NULL, df_gwl_baseline = NULL, baseline_start_manual = 1979, baseline_end_manual = 2019,
                               use_gwl_baseline = TRUE, gwl_baseline_level = 0.84) {
  # NOTE: get_coefficients() returns a matrix with one column per parameter draw. Therefore, here we loop through columns using map()
  
  if(!"matrix" %in% class(coefs_mc_matrix)) stop("coefs_mc_matrix must be of type matrix")
  
  map_dfr(1:ncol(coefs_mc_matrix), ~ project_impacts(data = data,
                                                     coefs = coefs_mc_matrix[, .x],
                                                     df_gwl_baseline = df_gwl_baseline,
                                                     baseline_start_manual = baseline_start_manual,
                                                     baseline_end_manual = baseline_end_manual,
                                                     use_gwl_baseline = use_gwl_baseline,
                                                     gwl_baseline_level = gwl_baseline_level) %>%
            mutate(monte_carlo_draw = .x)) %>%
    select(model, scenario, monte_carlo_draw, GID_1, year, everything())
}

# function to aggregate ADM1-level results to ADM0 using user-defined weights
aggregate_adm1_to_adm0 <- function(data = NULL, adm1_weights_used = NULL, is_mc_output = NULL) {
  if(mean(unique(data$GID_1) %in% adm1_weights_used$GID_1) < 1) stop("Not all ADM1 regions in data are covered by the ADM1 weights data provided")
  if(mean(c("GID_0", "GID_1", "adm1_weight") %in% names(adm1_weights_used)) != 1) stop("adm1_weights_used must contain columns GID_0, GID_1 and adm1_weight")
  if(is.null(is_mc_output)) stop("is_mc_output must be set to TRUE or FALSE depending on whether argument data is created via Monte Carlo analysis or not")
  
  df_out <- data %>% left_join(adm1_weights_used %>% select(GID_1, GID_0, adm1_weight), by = "GID_1") %>%
            dtplyr::lazy_dt() %>%
           {if(is_mc_output) group_by(., model, scenario, monte_carlo_draw, year, GID_0) else group_by(., model, scenario, year, GID_0)} %>%
            mutate(test_dummy = 1) %>%
            summarise_at(vars(starts_with("imp"), "test_dummy"), ~ sum(.x * adm1_weight)) %>%
            ungroup() %>%
            as_tibble()
  
  # the aggregated test_dummy reflects to what weights sum up - throw an error if this is not one (rounding here to avoid precision issues)
  if(mean(round(df_out$test_dummy, 5) ==1) != 1) stop(paste0("For some ADM0-level regions, weights sum up to the following values:\n",
                                                              paste0(df_out %>% filter(round(test_dummy, 10) != 1) %>% pull(test_dummy) %>% unique(), collapse = ", "),
                                                              "\nThis applies to model ", unique(df_out$model), " scenario ", unique(df_out$scenario)))
  
  # return the data frame without the test dummy variable
  return(df_out %>% select(-test_dummy))
}

# write a function to produce GDP impact projections by country and save them out directly and in separate files
produce_gdp_projections <- function(model = NULL, scenario = NULL,
                                    dir_climate_adm1 = NULL,
                                    vars_object = NULL,
                                    dir_out_mc = NULL, 
                                    overwrite_previous_files = TRUE,
                                    coefs_matrix_used = NULL,
                                    adm1_weights_used = NULL,
                                    df_gwl_baseline = NULL,
                                    baseline_start_manual = 1979,
                                    baseline_end_manual = 2019,
                                    use_gwl_baseline = TRUE,
                                    gwl_baseline_level = 0.84,
                                    ceiling_precipdeviation = NULL) {
  # load the climate variables
  df_adm1 <- extract_adm1_allvars(model = model,
                                  scenario = scenario,
                                  dir = dir_climate_adm1,
                                  vars_object = vars_object) %>%
            mutate(GID_0 = str_extract(GID_1, "^[A-Z][A-Z][A-Z]"))
  
  # apply the ceiling to monthly precip deviation ('Wn') if specified 
  if(!is.null(ceiling_precipdeviation)) df_adm1 <- df_adm1 %>% mutate(Wn = if_else(Wn > ceiling_precipdeviation, ceiling_precipdeviation, Wn))
  
  # read in all GID_0 codes in the data
  iso_strings <- unique(df_adm1$GID_0)
    
  # apply this function to the different countries using map()
  map(iso_strings, .f = function(iso) {
    file_name_iso <- file.path(dir_out_mc, paste0(model, "_", scenario, "_mc", ncol(coefs_matrix_used), "_", iso, ".feather"))
    
    if(!file.exists(file_name_iso) | overwrite_previous_files) {
      write_feather(df_adm1 %>% filter(GID_0 == iso) %>%
                        project_impacts_mc(coefs_mc_matrix = coefs_matrix_used,
                                           df_gwl_baseline = df_gwl_baseline,
                                           baseline_start_manual = baseline_start_manual,
                                           baseline_end_manual = baseline_end_manual,
                                           use_gwl_baseline = use_gwl_baseline,
                                           gwl_baseline_level = gwl_baseline_level) %>%
                        aggregate_adm1_to_adm0(adm1_weights_used = adm1_weights_used, is_mc_output = TRUE),
                    path = file_name_iso)
  }})
  
  return(TRUE)
}
