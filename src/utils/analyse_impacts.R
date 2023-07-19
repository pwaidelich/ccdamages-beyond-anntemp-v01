library(ggpattern)       # for hatching in maps created with geom_sf()
library(tidyverse)
library(ggpubr)
library(sf)


# write to helper functions for converting between log-scale and %GDP scale
log_to_gdp_impacts <- function(x) exp(x) - 1
gdp_to_log_impacts <- function(x) log(x+1)

# function to write out the climate variable names from Kotz et al 2022 used in our code
rewrite_climate_label <- function(x) {
  factor(x %>% str_replace("^Pt$", "Annual precip.") %>%
           str_replace("^Tstd$", "Day-to-day temp. var.") %>%
           str_replace("^wet_days_1$", "# of wet days") %>%
           str_replace("^Wn", "Monthly precip. deviation") %>%
           str_replace("^Tmean$", "Annual temperature") %>%
           str_replace("^vwet_days1_am_99p9$", "Extreme precip."),
         levels = c("Annual temperature", 
                    "Day-to-day temp. var.",
                    "Annual precip.",
                    "Monthly precip. deviation",
                    "# of wet days",
                    "Extreme precip.")
  )
}

# function to write out the impact labels used in our code (deriving from Kotz et al variable names)
rewrite_impact_label <- function(x) {
  factor(x %>% str_replace("^imp_Pt_diff$", "Annual precip.") %>%
           str_replace("^imp_Tstd_diff$", "Day-to-day temp var.") %>%
           str_replace("^imp_wet_days_1_diff$", "# of wet days") %>%
           str_replace("^imp_Wn_diff$", "Monthly precip. deviation") %>%
           str_replace("^imp_temp$", "Annual temperature (OLD VERSION XX)") %>%
           str_replace("^imp_temp_diff$", "Annual temperature") %>%
           str_replace("^imp_vwet_days1_am_99p9_diff$", "Extreme precip.") %>%
           str_replace("^imp_total$", "All indicators") %>%
           str_replace("^imp_nontemp$", "All indicators except ann. temp.") %>%
           str_replace("^imp_nontemp_diff$", "All indicators except ann. temp.") %>%
           str_replace("^imp_varextremes$", "Variability & extremes"),
         levels = c("All indicators", "All indicators except ann. temp.",
                    "Annual temperature",
                    "Annual temperature (OLD VERSION XX)",
                    "Annual precip.",
                    "Variability & extremes",
                    "Day-to-day temp var.", 
                    "Monthly precip. deviation",
                    "# of wet days", "Extreme precip.")
  )
}

# create another version with line breaks in the impact labels
rewrite_impact_label_linebreaks <- function(x) {
  factor(x %>% str_replace("^imp_Pt_diff$", "Annual\nprecip.") %>%
           str_replace("^imp_Tstd_diff$", "Day-to-day\ntemp. var.") %>%
           str_replace("^imp_wet_days_1_diff$", "# of wet\ndays") %>%
           str_replace("^imp_Wn_diff$", "Monthly precip.\ndeviation") %>%
           str_replace("^imp_temp_diff$", "Annual\ntemperature") %>%
           str_replace("^imp_temp$", "Annual temperature (OLD VERSION XX)") %>%
           str_replace("^imp_vwet_days1_am_99p9_diff$", "Extreme\nprecip.") %>%
           str_replace("^imp_total$", "All indicators") %>%
           str_replace("^imp_nontemp$", "All indicators except ann. temp.") %>%
           str_replace("^imp_nontemp_diff$", "All indicators except ann. temp.") %>%
           str_replace("^imp_varextremes$", "Variability &\nextremes"),
         levels = c("All indicators",
                    "All indicators except ann. temp.",
                    "Annual\ntemperature",
                    "Annual temperature (OLD VERSION XX)",
                    "Annual\nprecip.",
                    "Variability &\nextremes",
                    "Day-to-day\ntemp. var.", 
                    "Monthly precip.\ndeviation",
                    "# of wet\ndays", "Extreme\nprecip.")
  )
}


# function to write out the labels of summary stats used in the code
rewrite_summarystat_label <- function(x) {
  factor(x %>% str_replace("mean", "Mean") %>%
           str_replace("^rp_eta1$", "Risk premium (eta=1)") %>%
           str_replace("^ce_eta1$", "Certainty equivalent\n(eta=1)") %>%
           str_replace("^rp_eta1p35$", "Risk premium (eta=1.35)") %>%
           str_replace("^ce_eta1p35$", "Certainty equivalent\n(eta=1.35)") %>%
           str_replace("^var$", "Variance") %>%
           str_replace("^sd$", "St.Dv.") %>%
           str_replace("^var_mc$", "Socio-economic variance") %>%
           str_replace("^var_climate$", "Climatic variance") %>%
           str_replace("^var_interannual", "Inter-annual clim variance") %>%
           str_replace("^var_otherclimate", "Other clim variance"),
         levels = c("Risk premium (eta=1.35)",
                    "Risk premium \n(eta=1.35)",
                    "Risk premium (eta=1)", "Mean",
                    "Certainty equivalent\n(eta=1)",
                    "Certainty equivalent\n(eta=1.35)",
                    "Variance", "St.Dv.",
                    "Socio-economic variance", "Climatic variance",
                    "Inter-annual clim variance", "Other clim variance")
  )
}

# write a function for mapping mean impacts & agreement shares by ADM0 country
summarystat_agreement_hatchmap <- function(data = NULL,
                                    gwl_selected = 3,
                                    var_selected = "imp_nontemp",
                                    summarystat_selected = "mean",
                                    use_agreement = TRUE,
                                    low_colour = "red",
                                    high_colour = "blue",
                                    mid_colour = "gray88",
                                    fillscale_limits = c(NA, NA),
                                    fillscale_breaks = waiver(),
                                    pattern_fill = "white",
                                    pattern_fill_legend_key = "gray88",
                                    iso_subset = NULL,
                                    area_level = "ADM0",
                                    iso3_sovereign_countries = NULL) {
  
  if(area_level == "ADM1") stop("Currently not implemented (must first deal with ggpatterns inability with missing factor levels)")
  
  if(!area_level %in% c("ADM0", "ADM1")) stop("Argument area_level must be either ADM0 or ADM1")
  if(!summarystat_selected %in% c("mean", "var", "min", "max", paste0("perc", c("05", "10", "25", "75", "90", "95")))) stop("summarystat_selected does not match an accepted abbreviation for summary statistics of the distribution")
  
  var_selected_string <- paste0(var_selected, "_", summarystat_selected)
  agreement_string <- paste0(var_selected, "_agreement_meansign")
  signalnoise_string <- paste0(var_selected, "_signalnoise")
  
  if(mean(c(var_selected_string, agreement_string, signalnoise_string) %in% names(data)) != 1) stop("Either var_selected_string, signalnoise_string or agreement_string does not match a column name in argument data")
  
  data %>%
    filter(warming_level == gwl_selected) %>%
    mutate(var_for_plot = !!as.name(var_selected_string),
           agreement_for_plot = !!as.name(agreement_string),
           signalnoise_for_plot = !!as.name(signalnoise_string)) %>%
    {if(!is.null(iso_subset)) filter(., GID_0 %in% iso_subset) else .} %>%
    {if(area_level == "ADM0") left_join(., world, by = "GID_0") %>% st_as_sf() else .} %>%
    {if(area_level == "ADM1") left_join(., sf_gadm1_augmented, by = "GID_1") %>% st_as_sf() else .} %>%
    {if(use_agreement) mutate(., uncertainty = case_when(agreement_for_plot >= 0.9 ~ ">=90% (Very likely)",
                                                         agreement_for_plot >= 0.66 ~ "66-90% (Likely)",
                                                         TRUE ~ "<66%")) else mutate(.,
                                                                                     uncertainty = case_when(abs(signalnoise_for_plot) >= 1.65 ~ ">=90% (Very likely)",
                                                                                                             abs(signalnoise_for_plot) >= 0.98 ~ "66-90% (Likely)",
                                                                                                                TRUE ~ "<66%"))} %>% 
    {if(area_level == "ADM0" & !is.null(iso3_sovereign_countries)) mutate(., is_sovereign = GID_0 %in% iso3_sovereign_countries,
                                                                          var_for_plot = if_else(is_sovereign, var_for_plot, NA_real_),
                                                                          uncertainty = if_else(is_sovereign, uncertainty, ">=90% (Very likely)")) else .} %>%
    
    # hack-ish solution to the problem that ggpattern cannot deal with missing factor levels: add all potential values for uncertainty in a non-visible multi-polygon
    bind_rows(sf_gadm0_simple %>% filter(GID_0 == "VAT") %>% mutate(uncertainty = ">=90% (Very likely)")) %>%
    bind_rows(sf_gadm0_simple %>% filter(GID_0 == "VAT") %>% mutate(uncertainty = "66-90% (Likely)")) %>%
    bind_rows(sf_gadm0_simple %>% filter(GID_0 == "VAT") %>% mutate(uncertainty = "<66%")) %>%
    # recode uncertainty as a factor
    mutate(uncertainty = factor(uncertainty, levels = c(">=90% (Very likely)", "66-90% (Likely)", "<66%"))) %>%
    ggplot() + geom_sf_pattern(aes(fill = var_for_plot, pattern = uncertainty,
                                   pattern_angle = uncertainty,
                                   pattern_density = uncertainty,
                                   pattern_spacing = uncertainty), size = 0.1,
                               pattern_fill = pattern_fill,
                               pattern_colour = NA) +
    scale_fill_gradient2(low = low_colour, mid = mid_colour, high = high_colour,
                         n.breaks = 3,
                         breaks = fillscale_breaks,
                         labels = scales::percent,
                         limits = fillscale_limits) +
    scale_pattern_angle_manual(name = ifelse(use_agreement, "Agreement\n on mean's sign", "Confidence \n(signal-noise)"),
                               labels = c(">=90% (Very likely)", "66-90% (Likely)", "<66%"),
                               values = c(0, 45, 180)) +
    scale_pattern_spacing_manual(name = ifelse(use_agreement, "Agreement\n on mean's sign", "Confidence \n(signal-noise)"),
                                 labels = c(">=90% (Very likely)", "66-90% (Likely)", "<66%"),
                                 values = c(0, 0.025, 0.01)) +
    scale_pattern_density_manual(name = ifelse(use_agreement, "Agreement\n on mean's sign", "Confidence \n(signal-noise)"),
                                 labels = c(">=90% (Very likely)", "66-90% (Likely)", "<66%"),
                                 values = c(0, 0.2, 0.35)) +
    scale_pattern_manual(name = ifelse(use_agreement, "Agreement\n on mean's sign", "Confidence \n(signal-noise)"),
                         labels = c(">=90% (Very likely)", "66-90% (Likely)", "<66%"),
                         values = c("none", "stripe", "stripe"),
                         guide = guide_legend(override.aes = list(fill = pattern_fill_legend_key))) +
    guides(pattern = guide_legend(nrow = 3, byrow = TRUE)) +
    labs(fill = paste0("GDP impact (", summarystat_selected, ")"),
         title = paste0("GWL: +", gwl_selected, "C - ", var_selected %>% rewrite_impact_label())
    ) +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          axis.line = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank())
}

# write a function for mapping mean impacts & agreement shares by ADM0 country
new_hatchmap <- function(data = NULL,
                                           gwl_selected = 3,
                                           var_selected = "imp_nontemp",
                                           summarystat_selected = "mean",
                                           low_colour = "red",
                                           high_colour = "blue",
                                           mid_colour = "gray88",
                                           fillscale_limits = c(NA, NA),
                                           fillscale_breaks = waiver(),
                                           pattern_fill = "white",
                                           pattern_fill_legend_key = "gray88",
                                           iso_subset = NULL,
                                           area_level = "ADM0") {
  
  if(!area_level %in% c("ADM0", "ADM1")) stop("Argument area_level must be either ADM0 or ADM1")
  if(!summarystat_selected %in% c("mean", "var", "min", "max", paste0("perc", c("05", "10", "25", "75", "90", "95")))) stop("summarystat_selected does not match an accepted abbreviation for summary statistics of the distribution")
  
  var_selected_string <- paste0(var_selected, "_", summarystat_selected)
  mean_string <- paste0(var_selected, "_mean")
  variance_string <- paste0(var_selected, "_var")
  
  if(mean(c(var_selected_string, mean_string, variance_string) %in% names(data)) != 1) stop("Either var_selected_string, mean_string or variance_string does not match a column name in argument data")
  
  data %>%
    filter(warming_level == gwl_selected) %>%
    mutate(var_for_plot = !!as.name(var_selected_string),
           mean_for_signalnoise = !!as.name(mean_string),
           var_for_signalnoise = !!as.name(variance_string)) %>%
    {if(!is.null(iso_subset)) filter(., GID_0 %in% iso_subset) else .} %>%
    mutate(signalnoise = mean_for_signalnoise/(var_for_signalnoise)^0.5) %>%
    mutate(uncertainty = factor(case_when(abs(signalnoise) >= 1.65 ~ ">=90% (Very likely)",
                                          abs(signalnoise) >= 0.98 ~ "66-90% (Likely)",
                                          TRUE ~ "<66%"),
                                levels = c(">=90% (Very likely)", "66-90% (Likely)",  "<66%"))) %>%
    {if(area_level == "ADM0") left_join(., world, by = "GID_0") %>% st_as_sf() else .} %>%
    {if(area_level == "ADM1") left_join(., sf_gadm1_augmented, by = "GID_1") %>% st_as_sf() else .} %>%
    ggplot() + geom_sf_pattern(aes(fill = var_for_plot, pattern = uncertainty,
                                   pattern_angle = uncertainty,
                                   pattern_density = uncertainty,
                                   pattern_spacing = uncertainty), size = 0.1,
                               pattern_fill = pattern_fill,
                               pattern_colour = NA) +
    scale_fill_gradient2(low = low_colour, mid = mid_colour, high = high_colour,
                         n.breaks = 3,
                         breaks = fillscale_breaks,
                         labels = scales::percent,
                         limits = fillscale_limits) +
    scale_pattern_angle_manual(name = "Confidence signal-noise ratio",
                               labels = c(">=90% (Very likely)", "66-90% (Likely)",  "<66%"),
                               values = c(0, 45, 180)) +
    scale_pattern_spacing_manual(name = "Confidence signal-noise ratio",
                                 labels = c(">=90% (Very likely)", "66-90% (Likely)",  "<66%"),
                                 values = c(0, 0.025, 0.01)) +
    scale_pattern_density_manual(name = "Confidence signal-noise ratio",
                                 labels = c(">=90% (Very likely)", "66-90% (Likely)",  "<66%"),
                                 values = c(0, 0.2, 0.35)) +
    scale_pattern_manual(name = "Confidence signal-noise ratio",
                         labels = c(">=90% (Very likely)", "66-90% (Likely)",  "<66%"),
                         values = c("none", "stripe", "stripe"),
                         guide = guide_legend(override.aes = list(fill = pattern_fill_legend_key))) +
    guides(pattern = guide_legend(nrow = 3, byrow = TRUE)) +
    labs(fill = paste0("GDP impact (", summarystat_selected, ")"),
         title = paste0("GWL: +", gwl_selected, "C - ", var_selected %>% rewrite_impact_label()),
         #subtitle = var_label_subtitle
    ) +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          axis.line = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank())
}


# write a function for mapping mean impacts & agreement shares by ADM0 country
mean_agreement_hatchmap <- function(data = NULL,
                                    gwl_selected = 3,
                                    var_selected = "imp_nontemp",
                                    low_colour = "red",
                                    high_colour = "blue",
                                    mid_colour = "gray88",
                                    pattern_fill = "white",
                                    pattern_fill_legend_key = "gray88",
                                    iso_subset = NULL,
                                    area_level = "ADM0") {
  
  if(!area_level %in% c("ADM0", "ADM1")) stop("Argument area_level must be either ADM0 or ADM1")
  
  plot_out <- data %>%
    filter(warming_level == gwl_selected) %>%
    filter(var == var_selected) %>%
    {if(!is.null(iso_subset)) filter(., GID_0 %in% iso_subset) else .} %>%
    mutate(uncertainty = factor(case_when(agreement_share >= 0.9 ~ ">=90% (Very likely)",
                                          agreement_share >= 0.66 ~ "66-90% (Likely)",
                                          TRUE ~ "<66%"),
                                levels = c(">=90% (Very likely)", "66-90% (Likely)",  "<66%"))) %>%
    {if(area_level == "ADM0") left_join(., world, by = "GID_0") %>% st_as_sf() else .} %>%
    {if(area_level == "ADM1") left_join(., sf_gadm1_augmented, by = "GID_1") %>% st_as_sf() else .} %>%
    ggplot() + geom_sf_pattern(aes(fill = mean, pattern = uncertainty,
                                   pattern_angle = uncertainty,
                                   pattern_density = uncertainty,
                                   pattern_spacing = uncertainty), size = 0.1,
                               pattern_fill = pattern_fill,
                               pattern_colour = NA) +
    scale_fill_gradient2(low = low_colour, mid = mid_colour, high = high_colour,
                         n.breaks = 3,
                         breaks = waiver(),
                         labels = scales::percent) +
    scale_pattern_angle_manual(name = "Agreement\n on sign",
                               labels = c(">=90% (Very likely)", "66-90% (Likely)",  "<66%"),
                               values = c(0, 45, 180)) +
    scale_pattern_spacing_manual(name = "Agreement\n on sign",
                                 labels = c(">=90% (Very likely)", "66-90% (Likely)",  "<66%"),
                                 values = c(0, 0.025, 0.01)) +
    scale_pattern_density_manual(name = "Agreement\n on sign",
                                 labels = c(">=90% (Very likely)", "66-90% (Likely)",  "<66%"),
                                 values = c(0, 0.2, 0.35)) +
    scale_pattern_manual(name = "Agreement\n on sign",
                         labels = c(">=90% (Very likely)", "66-90% (Likely)",  "<66%"),
                         values = c("none", "stripe", "stripe"),
                         guide = guide_legend(override.aes = list(fill = pattern_fill_legend_key))) +
    guides(pattern = guide_legend(nrow = 3, byrow = TRUE)) +
    labs(fill = "GDP impact",
         title = paste0("GWL: +", gwl_selected, "C - ", var_selected %>% rewrite_impact_label()),
         #subtitle = var_label_subtitle
    ) +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          axis.line = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank())
  
  plot_out
  # ggsave(file.path("../sharepoint/Figures/adm0/socioeconomic/montecarlo",
  #                  paste0(Sys.Date(), " ", var_selected, "_map_gwl", gwl_selected, ".pdf")),
  #        width = 6, height = 4, dpi = 50)
  # 
  # print("map saved out successfully")
}

########################### CHARTS #############################################

# write a function to plot projection based on averaging CMIP6 together with model agreement share
cmip6_sumstat_map <- function(data = NULL, var_selected = NULL, gwl_selected = NULL, var_label = NULL,
                              var_label_subtitle = "Projection based on CMIP6 average",
                              agreement_label_fill = "% of obs agreeing \nwith sign",
                              agreement_label_subtitle = "CMIP6 inter-model agreement",
                              midpoint_agreement = 0.66,
                              low_colour = "red", high_colour = "blue", level = "GID_1") {
  if(!var_selected %in% data$var) stop("var_selected is not featured in the data")
  
  # create the relevant subset (variable & global warming level of interest) & add polygon shapefile info
  if(level == "GID_1") {
    df_map <- data %>% filter(warming_level == gwl_selected, var == var_selected) %>%
      right_join(sf_gadm1_augmented, by = c("GID_0", "GID_1")) %>%
      st_as_sf()
  } else if(level == "GID_0") {
    df_map <- data %>% filter(warming_level == gwl_selected, var == var_selected) %>%
      right_join(sf_gadm0_simple, by = c("GID_0")) %>%
      st_as_sf()
  } else {
    stop("Argument level must be either GID_1 or GID_0")
  }
  
  # make a 2x1 grid chart of two maps: predicted mean of the selected variable & model agreement, both by ADM1 region
  ggarrange(  
    df_map %>%
      # filter again to throw out territories for which we have polygons but no impact projection
      filter(warming_level == gwl_selected, var == var_selected) %>%
      ggplot() + geom_sf(aes(fill = var_mean), size = 0.1) +
      geom_sf(data = world, fill = NA, colour = "black", size = 0.25) +
      scale_fill_gradient2(low = low_colour, mid = "white", high = high_colour,
                           n.breaks = 3) +
      labs(fill = var_label,
           title = paste0("Global warming level: +", gwl_selected, "C"),
           subtitle = var_label_subtitle) +
      #guides(x = guide_axis(n.dodge=2)) +
      #guides(x = guide_axis(angle = 45)) + #, fill = guide_colorbar(direction = "horizontal")) +
      theme_classic() +
      theme(legend.position = "bottom"),
    
    df_map %>%
      filter(warming_level == gwl_selected, var == var_selected) %>%
      ggplot() + geom_sf(aes(fill = agreement_share), size = 0.1) +
      geom_sf(data = world, fill = NA, colour = "black", size = 0.25) +
      scale_fill_gradientn(colours = c("darkred", "white", "darkgreen"),
                           values = scales::rescale(c(0, midpoint_agreement, 1)),
                           limits = c(0, 1)) +
      labs(fill = agreement_label_fill,
           subtitle = agreement_label_subtitle) +
      theme_classic() +
      theme(legend.position = "bottom"),
    nrow = 2, ncol = 1, align = "hv"
  )
}

utility_isoelastic <- function(x = NULL, eta = 1) {
  if(!is.numeric(x) | !is.numeric(eta)) stop("Arguments x and eta must be numeric")
  if(length(eta) != 1) stop("Argument eta must have length one")
  
  if(eta == 1) {
    return(log(x))
  } else {
    return((x^(1-eta) - 1)/(1 - eta))
  }
}

