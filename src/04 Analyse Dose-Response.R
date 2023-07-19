# clean out the environment
rm(list = ls())

# load packages
library(PupillometryR)   # for raincloud charts
library(tidyverse)       # for general data wrangling & plotting
library(ggpubr)          # for grid charts via ggarrange()
library(furrr)           # for parallel processing
library(feather)         # for importing .feather files

# source utils functions
source(file.path("src", "utils", "project_impacts.R"))
source(file.path("src", "utils", "analyse_impacts.R"))

# set the alpha parameter for dose-response charts
alpha_doseresponse <- 0.1
theme_set(theme_classic())


################################################################################
############################## CONCEPTUAL METHODOLOGY CHART ####################
################################################################################

# read in the point estimate coefficients from Kotz et al 2022 & 1000 MC draws (clustered at country level)
coefs <- get_coefficients()[,1]
coefs_mc <- get_coefficients(draw_from_vcov = T, n = 1000) %>% t() %>% as_tibble()

# load shapefiles for ADM1-level regions
sf_gadm1_augmented <- st_read(file.path("..", "sharepoint", "Data", "GADM", "gadm36_levels_simple_augmented", "adm1_augmented.shp"))

# read in the ADM1-level data that uses point estimates for damage function coefficients
df_adm1 <- read_feather(file.path("/net", "cfc", "landclim1", "pwaidelich", "adm1_new", "pointestimates", "df_adm1.feather"))

# select the example region of interest (New York state)
gid1_example <- sf_gadm1_augmented %>% filter(str_detect(NAME_1, "New York")) %>% pull(GID_1)

# set the example CMIP6 model
example_model <- "ACCESS-CM2"

# extract the data for the example region under SSP3-7.0 & mark the example model
data_chart <- df_adm1 %>% filter(GID_1 == gid1_example, scenario == "ssp370") %>%
  mutate(main_model_factor = factor(model, levels = c(example_model, unique(df_adm1$model) %>% str_subset(example_model, negate= T))),
         is_main_model = model == example_model)

# determine the warming-level window
df_gwl <- readRDS(file.path("data", "global_warming_levels.rds"))

# extracct start and end year for the baseline period (+0.84C) of the example model (under SSP3-7.0) as well as for +3C
warming_level0p84_example_beg <- df_gwl %>% filter(scenario == "ssp370", model == example_model, warming_level == 0.84) %>% summarise(beg = min(year)) %>% pull(beg)
warming_level0p84_example_end <- df_gwl %>% filter(scenario == "ssp370", model == example_model, warming_level == 0.84) %>% summarise(end = max(year)) %>% pull(end)
warming_level_baseline <- warming_level0p84_example_beg:warming_level0p84_example_end
warming_level3_example_beg <- df_gwl %>% filter(scenario == "ssp370", model == example_model, warming_level == 3) %>% summarise(beg = min(year)) %>% pull(beg)
warming_level3_example_end <- df_gwl %>% filter(scenario == "ssp370", model == example_model, warming_level == 3) %>% summarise(end = max(year))%>% pull(end)

# plot the extreme precipitation variable for the example region and the example model since start of baseline window
exampleadm1_cmip6_timeseries <- data_chart %>% filter(is_main_model, year >= warming_level0p84_example_beg) %>%
  ggplot(aes(year, vwet_days1_am_99p9)) +
  # add verticle lines to mark the GWL windows
  geom_vline(aes(xintercept = warming_level0p84_example_beg), linetype = "dashed", colour = "#7CAE00") +
  geom_vline(aes(xintercept = warming_level0p84_example_end), linetype = "dashed", colour = "#7CAE00") +
  geom_vline(aes(xintercept = warming_level3_example_beg),
             linetype = "dashed", colour = "#619CFF") +
  geom_vline(aes(xintercept = warming_level3_example_end),
             linetype = "dashed", colour = "#619CFF") +
  geom_line(data = data_chart %>% filter(!is_main_model, year >= warming_level0p84_example_beg), aes(colour = model), alpha = 0.3) +
  geom_line(colour = "black") +
  # add coloured points for all data points in the GWL windows
  geom_point(aes(y = ifelse(is_baseline, vwet_days1_am_99p9, NA)), colour = "#7CAE00") +
  geom_point(aes(y = ifelse(year %in% warming_level3_example_beg:warming_level3_example_end, vwet_days1_am_99p9, NA)), colour = "#619CFF") +
  # add text labels explaining the example model as well as the vertical lines
  annotate("text", x = mean(warming_level_baseline), y = 200, label = "Baseline period\n(+0.84\u00B0C)", colour = "#7CAE00", fontface = 2, size = 3) +
  annotate("text", x = warming_level3_example_end + 1, y = -10, label = paste0("Example model (", example_model, ")"), hjust = 0, fontface = 2, size = 2) +
  annotate("text", x = warming_level3_example_end + 1, y = -20, label = "Other CMIP6 models", hjust = 0, fontface = 2, colour = "darkgrey", size = 2) +
  annotate("text", x = df_gwl %>% filter(scenario == "ssp370", model == example_model, warming_level == 3) %>% summarise(beg = min(year), end = max(year)) %>%
             mutate(mid = 0.5*(beg + end)) %>% pull(mid), y = 200, label = "+3\u00B0C global\n warming", colour = "#619CFF", fontface = 2, size = 3) +
  labs(y = "Extreme precipitation in New York\n state under SSP3-7.0 (mm)") +
  scale_colour_manual(values = c(rep("lightgrey", length(unique(df_adm1$model)) - 1))) +
  theme(legend.position = "none", axis.title.x = element_blank())

# inspect
exampleadm1_cmip6_timeseries

# extract the impact projections for the example region & model under SSP3-7.0
df_adm1_example <- df_adm1 %>% filter(GID_1 == gid1_example, model == example_model, scenario == "ssp370") %>%
    # identify all years in the baseline period and in +3C window in a new 'period' column, set to NA for other years
    mutate(period = case_when(year %in% warming_level_baseline ~ "Baseline period",
                              year %in% warming_level3_example_beg:warming_level3_example_end ~ "+3\u00B0C global warming",
                              TRUE ~ NA_character_))

# extract the average temperature and the average extreme precipitation impacts in the baseline period
# NOTE: we need baseline mean temp to calculate the 95% CI for this particular region as the extreme precip dose-response function interacts with it
Tmean_base_example <- df_adm1_example %>% filter(year %in% warming_level_baseline) %>% pull(Tmean) %>% mean()
# we use baseline mean impacts & mean extreme precip to illustrate how GDP losses vis-a-vis baseline are calculated
imp_basemean_example <- df_adm1_example %>% filter(year %in% warming_level_baseline) %>% summarise(imp_vwet_days1_am_99p9 = mean(imp_vwet_days1_am_99p9),
                                                                                                   vwet_days1_am_99p9 = mean(vwet_days1_am_99p9))

# select an example year for which the chart illustrates how damages (impact - baseline mean of impact) is calculated
highlighted_year_example <- 2067
df_adm1_example %>% filter(year == highlighted_year_example)

# calculate the 95% CI for a sufficiently large sample size
# NOTE: to do so, calculate the marginal effect (which is constant as the function is linear) for each draw and then take the respective quantiles
ci_bounds_example <- get_coefficients(draw_from_vcov = T, n = 10^5) %>% t() %>% as_tibble() %>% mutate(coef_extreme_precip = vwet_days1_am_99p9 + Tmean_base_example * Tmean.vwet_days1_am_99p9) %>%
  summarize(perc0p025 = quantile(coef_extreme_precip, 0.025),
            perc0p975 = quantile(coef_extreme_precip, 0.975))

# plot the dose-response function together with the projected data points
exampleadm1_impacts <- expand.grid(Tmean = Tmean_base_example,
             # for value range of extreme precip, we go from zero to the maximum in the relevant warming levels and add a minor value for visual reasons only
            vwet_days1_am_99p9 = seq(0, 3 + max(df_adm1_example %>% filter(year %in% warming_level0p84_example_beg:warming_level3_example_end) %>% pull(vwet_days1_am_99p9)))) %>%
  as_tibble() %>%
  # NOTE: raw impact calculations using coefficients & climate indicators are in log-scale, so we need to convert to %GDP
  mutate(y = log_to_gdp_impacts(coefs["vwet_days1_am_99p9"]*vwet_days1_am_99p9 + coefs["Tmean.vwet_days1_am_99p9"]*vwet_days1_am_99p9*Tmean),
         y_0p025 = log_to_gdp_impacts(vwet_days1_am_99p9 * ci_bounds_example$perc0p025),
         y_0p975 = log_to_gdp_impacts(vwet_days1_am_99p9 * ci_bounds_example$perc0p975)) %>%
  ggplot(aes(vwet_days1_am_99p9, y)) +
  # calculate the dose-response functions for different Monte Carlo draws - COMMENTED OUT SINCE WE USE THE 95% CI INSTEAD
  # plyr::alply(as.matrix(coefs_mc), 1, function(coefs) {
  #   stat_function(fun= ~ log_to_gdp_impacts(coefs["vwet_days1_am_99p9"]*.x + coefs["Tmean.vwet_days1_am_99p9"]*.x*Tmean_base_example), colour="grey", alpha = alpha_doseresponse)
  # }) +
  # plot the point estimates as a black line and the 95% CI around it
  geom_line(colour = "black") +
  geom_ribbon(aes(ymin = y_0p025, ymax = y_0p975), fill = "lightgrey", alpha = 0.3) +
  # draw a dashed horizontal line from the baseline average point to the example year
  annotate("segment", x = imp_basemean_example$vwet_days1_am_99p9,
           xend = df_adm1_example %>% filter(year == highlighted_year_example) %>% pull(vwet_days1_am_99p9),
           y = imp_basemean_example$imp_vwet_days1_am_99p9, yend = imp_basemean_example$imp_vwet_days1_am_99p9,
           linetype = "dashed", alpha = 0.2) +
  # draw a vertical errorbar illustrating the difference between the example year impact and the baseline average
  annotate("errorbar", x = df_adm1_example %>% filter(year == highlighted_year_example) %>% pull(vwet_days1_am_99p9),
           ymin = df_adm1_example %>% filter(year == highlighted_year_example) %>% pull(imp_vwet_days1_am_99p9),
           ymax = imp_basemean_example$imp_vwet_days1_am_99p9,
           linetype = "solid",
           colour = "#f8766d",
           width = 3) +
  geom_point(data = df_adm1_example %>% filter(!is.na(period)), aes(vwet_days1_am_99p9, imp_vwet_days1_am_99p9,
                                                                    colour = period), alpha = 0.8) +
  # add an extra point with different shape and larger size for the baseline average
  geom_point(data = imp_basemean_example %>% mutate(period = "Baseline period",
                                                    summary_stat = "Baseline average"), aes(vwet_days1_am_99p9, imp_vwet_days1_am_99p9,
                                                  colour = period, shape = summary_stat), colour = "#7CAE00", size = 6) +
  # # text labels to explain what line corresponds to the point estimates and the MC draws - COMMENTED OUT AS WE USE THE 95% CI
  # annotate("text", x = 50, y = -0.055, label = "Monte Carlo draws", hjust = 0, fontface = 2, colour = "darkgrey") +
  # annotate("text", x = 50, y = -0.05, label = "Central damage function", hjust = 0, fontface = 2, colour = "black") +
  # text labels for an example year to illustrate how damages are derived using the dose-response function & the baseline
  annotate("text", x = 90, y = mean(c(df_adm1_example %>% filter(year == highlighted_year_example) %>% pull(imp_vwet_days1_am_99p9),
                                    imp_basemean_example$imp_vwet_days1_am_99p9)),
           label = paste0("Projected damages\nfor an example\nyear (", highlighted_year_example, "): ",
                          100*round(df_adm1_example %>% filter(year == highlighted_year_example) %>%
                            pull(imp_vwet_days1_am_99p9_diff), 3), "%"), hjust = 0, fontface = 2, colour = "#be1509", size = 3) +
  guides(alpha = "none") +
  scale_linetype_manual(values = c("solid", "dashed", "dashed"), guide = "none") +
  scale_shape_manual(values = c(18)) +
  labs(x = "Extreme precipitation in New York state under SSP3-7.0 (mm)", y = "Change in GDP",
       colour = NULL, shape = NULL) +
  scale_colour_manual(values = c("#7CAE00", "#619CFF"), limits = rev) +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = c(0.25, 0.25),
        legend.text=element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.background = element_blank(),
        legend.margin=unit(0, "cm"))
exampleadm1_impacts

# import the ADM0-level Monte Carlo impact projections created in "src/01 Create Projections.R"
df_adm0_mc_us <- file.path("/net", "cfc", "landclim1", "pwaidelich", "adm0_new", "montecarlo") %>% list.files(full.names = T) %>%
  # extract the right Monte Carlo sample size & example country (US for NY state)
  str_subset("_mc1000_.+\\.feather$") %>% str_subset("USA\\.feather$") %>%
  map_dfr(read_feather) %>%
  left_join(df_gwl %>% select(-ensemble), by = c("model", "scenario", "year")) %>%
  filter(warming_level == 3) %>%
  # for each model-scenario combination, merge in the start and end of the +3C global warming level window
  left_join(df_gwl %>% group_by(model, scenario, ensemble, warming_level) %>% summarise(beg = min(year), end = max(year)) %>% ungroup() %>%
              select(-ensemble),
            by = c("model", "scenario", "warming_level")) %>%
  mutate(model_label = paste0(model, "\n(", beg, "-", end, ")"))

# ensure that we have exactly 1000 (MC draws) x 20 (years) observations for each model
df_adm0_mc_us %>% count(model_label) %>% count(n)

df_adm0_mc_us_boxplot <- df_adm0_mc_us %>% group_by(model_label) %>%
  summarise(median = median(imp_vwet_days1_am_99p9_diff),
            lower_quartile = quantile(imp_vwet_days1_am_99p9_diff, 0.25),
            upper_quartile = quantile(imp_vwet_days1_am_99p9_diff, 0.75),
            lower_decile = quantile(imp_vwet_days1_am_99p9_diff, 0.10),
            upper_decile = quantile(imp_vwet_days1_am_99p9_diff, 0.90))

# plot the ADM0-level impacts for the example region's country (= US)
exampleadm0_impact_distribution <- df_adm0_mc_us %>% group_by(model_label) %>% ggplot() +
  aes(x = model_label, fill = model_label) +
  geom_point(aes(y = imp_vwet_days1_am_99p9_diff),
             color = "#F8766D",
             position = position_jitter(w = .15),
             size = 0.5,
             alpha = 0.05) +
  geom_boxplot(data = df_adm0_mc_us_boxplot %>% group_by(model_label),
               stat = "identity",
               aes(lower = lower_quartile, upper = upper_quartile,
                   middle = median, ymin = lower_decile, ymax = upper_decile),
               width = .25,
               fill = "#F8766D",
               outlier.shape = NA,
               alpha = 0.5) +
  geom_flat_violin(aes(y = imp_vwet_days1_am_99p9_diff),
                   position = position_nudge(x = .2),
                   fill = "#F8766D",
                   alpha = 0.7,
                   adjust = 0.5)  +
  # clip the chart due to some outliers for ACCESS-CM2 & KACE-1-0 for the sake of readability
  coord_flip(ylim = c(-0.025, NA)) +
  guides(fill = "none", color = "none") +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "CMIP6 model (and respective +3\u00B0C window)",
       y = "Projected damages of extreme precipitation\n in the US for years in +3\u00B0C window") +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text.y=element_text(face = c(rep("plain", length(unique(df_adm0_mc_us$model_label)) - 1),
                                          "bold"), size = 7.5))
# inspect the result
exampleadm0_impact_distribution

# combine the three figures for the grid chart
figure_conceptual <- ggarrange(
  ggarrange(exampleadm1_cmip6_timeseries, exampleadm1_impacts, nrow = 2, ncol = 1, 
          labels = "auto", align = "hv"),
  exampleadm0_impact_distribution,
  widths = c(1.2, 1),
  nrow = 1, ncol = 2, labels = c(NA, "c")
)

# print out values exceeding the 2.5% threshold used for the chart
df_adm0_mc_us %>% filter(imp_vwet_days1_am_99p9_diff < -0.025) %>% count(model)
df_adm0_mc_us %>% group_by(model) %>% slice_min(n=1, imp_vwet_days1_am_99p9_diff) %>%
  select(model, scenario, monte_carlo_draw, year, GID_0, imp_vwet_days1_am_99p9_diff) %>%
  arrange(imp_vwet_days1_am_99p9_diff)

# save out
ggsave(file.path("../sharepoint/Figures",
                 paste0(Sys.Date(), " Figure1.png")),
       plot = annotate_figure(figure_conceptual, top = text_grob("Illustration for one example climate indicator (out of six), CMIP6 model and region (NY state)",
                                                                 color = "black", face = "bold", size = 12, hjust = 0.565)),
       width = 10.1, height = 6.5)


################################################################################
############################## DOSE-RESPONSE CHARTS ############################
################################################################################

# read out the 0.1th and 99.9th percentiles of the ADM1-level distribution over the baseline period
# NOTE: we use these to determine value ranges for the x-axes in dose-response charts
df_perc <- df_adm1 %>% filter(is_baseline) %>% select(Tmean, Pt, Tstd, wet_days_1, Wn, vwet_days1_am_99p9) %>%
  summarize_all(.f = list(perc0p001 = ~ quantile(.x, 0.001),
                          perc0p999 = ~ quantile(.x, 0.999)))

# overwrite some of these percentiles with natural bounds (e.g., lower bound of zero for variability)
df_bounds <- tibble(Tmean_lower = df_perc$Tmean_perc0p001,
                    Tmean_upper = df_perc$Tmean_perc0p999,
                    Pt_lower = min(0, df_perc$Pt_perc0p001),
                    Pt_upper = df_perc$Pt_perc0p999,
                    Tstd_lower = 0,
                    Tstd_upper = df_perc$Tstd_perc0p999,
                    Wn_lower = df_perc$Wn_perc0p001,
                    Wn_upper = df_perc$Wn_perc0p999,
                    wet_days_1_lower = 0,
                    wet_days_1_upper = 366,
                    vwet_days1_am_99p9_lower = min(0, df_perc$vwet_days1_am_99p9_perc0p001),
                    vwet_days1_am_99p9_upper = df_perc$vwet_days1_am_99p9_perc0p999)

# create a list that hosts the different dose-response function charts for each climate indicator
dose <- list()

# dose-response for annual mean temperature
# NOTE: for mean temperature, impacts are actually trajectory-dependent - here we use the integral of the marginal effect as a simplification for illustration purposes (see SI)
dose[["Tmean"]] <- tibble(Tmean = seq(df_bounds$Tmean_lower, df_bounds$Tmean_upper, by = 0.1)) %>%
  ggplot(aes(Tmean)) +
  # draw one line for each coefficient draw in coefs_mc
  plyr::alply(as.matrix(coefs_mc), 1, function(x_coef) {
    stat_function(fun= ~ log_to_gdp_impacts((x_coef["TmeanD"] + x_coef["TmeanLD"])*.x + 0.5*(x_coef["TmeanD.Tmean"] + x_coef["TmeanLD.TmeanL"])*.x^2), alpha = alpha_doseresponse*0.1)
  }) +
  geom_function(fun = ~ log_to_gdp_impacts((coefs["TmeanD"] + coefs["TmeanLD"])*.x + 0.5*(coefs["TmeanD.Tmean"] + coefs["TmeanLD.TmeanL"])*.x^2)) +
  guides(alpha = "none") +
  labs(x = "Annual temperature (\u00B0C)", y = "Change in GDP") +
  scale_y_continuous(labels = scales::percent)

# dose-response for annual precip
dose[["Pt"]] <- tibble(Pt = seq(df_bounds$Pt_lower, df_bounds$Pt_upper, by = 100)) %>%
  ggplot(aes(Pt)) +
  plyr::alply(as.matrix(coefs_mc), 1, function(coefs) {
    stat_function(fun= ~ log_to_gdp_impacts(coefs["Pt"]*.x + coefs["Pt2"]*.x^2), colour="grey", alpha = alpha_doseresponse)
  }) +
  geom_function(fun = ~ log_to_gdp_impacts(coefs["Pt"]*.x + coefs["Pt2"]*.x^2)) +
  guides(alpha = "none") +
  labs(x = "Annual precip (mm)", y = "Change in GDP") +
  scale_y_continuous(labels = scales::percent)

# dose-response for monthly precip deviation
dose[["Wn"]] <- tibble(Wn = seq(df_bounds$Wn_lower, df_bounds$Wn_upper, by = 0.1)) %>%
  ggplot(aes(Wn)) +
  plyr::alply(as.matrix(coefs_mc), 1, function(coefs) {
    stat_function(fun= ~ log_to_gdp_impacts(coefs["Wn"]*.x + coefs["Wn_2"]*.x^2), colour="grey", alpha = alpha_doseresponse)
  }) +
  geom_function(fun = ~ log_to_gdp_impacts(coefs["Wn"]*.x + coefs["Wn_2"]*.x^2)) +
  guides(alpha = "none") +
  labs(x = "Monthly precip variation (unitless)", y = "Change in GDP") +
  scale_y_continuous(labels = scales::percent)

# dose-response for # of wet days
dose[["wet_days_1"]] <- tibble(wet_days_1 = seq(df_bounds$wet_days_1_lower, df_bounds$wet_days_1_upper, by = 1)) %>%
  ggplot(aes(wet_days_1)) +
  plyr::alply(as.matrix(coefs_mc), 1, function(coefs) {
    stat_function(fun= ~ log_to_gdp_impacts(coefs["wet_days_1"]*.x + coefs["wet_days_1_2"]*.x^2), colour="grey", alpha = alpha_doseresponse)
  }) +
  geom_function(fun = ~ log_to_gdp_impacts(coefs["wet_days_1"]*.x + coefs["wet_days_1_2"]*.x^2)) +
  guides(alpha = "none") +
  labs(x = "No. of wet days (#)", y = "Change in GDP") +
  scale_y_continuous(labels = scales::percent)

# dose-response for day-to-day temp variability
dose[["Tstd"]] <- tibble(Tstd = seq(df_bounds$Tstd_lower, df_bounds$Tstd_upper, by = 0.05)) %>%
  ggplot(aes(Tstd)) +
  plyr::alply(as.matrix(coefs_mc), 1, function(coefs) {
    stat_function(fun= ~ log_to_gdp_impacts(coefs["Tstd"]*.x), colour="grey", alpha = alpha_doseresponse)
  }) +
  geom_function(fun = ~ log_to_gdp_impacts(coefs["Tstd"]*.x)) +
  guides(alpha = "none") +
  labs(x = "Day-to-day temp. variability (\u00B0C)", y = "Change in GDP") +
  scale_y_continuous(labels = scales::percent)

# dose-response for extreme precip (for three example temperatures: 5C, 15C, and 25C)
dose[["vwet_days1_am_99p9"]] <- expand.grid(Tmean = c(5, 15, 25),
                                            vwet_days1_am_99p9 = seq(df_bounds$vwet_days1_am_99p9_lower, df_bounds$vwet_days1_am_99p9_upper, by = 1)) %>%
  as_tibble() %>%
  mutate(y = log_to_gdp_impacts(coefs["vwet_days1_am_99p9"]*vwet_days1_am_99p9 +
           coefs["Tmean.vwet_days1_am_99p9"]*vwet_days1_am_99p9*Tmean)) %>%
  ggplot(aes(vwet_days1_am_99p9, y)) +
  # carry out the Monte Carlo draw functions only for an annual mean temperature of 15 (example value)
  plyr::alply(as.matrix(coefs_mc), 1, function(coefs) {
    stat_function(fun= ~ log_to_gdp_impacts(coefs["vwet_days1_am_99p9"]*.x + coefs["Tmean.vwet_days1_am_99p9"]*.x*15), colour="grey", alpha = alpha_doseresponse)
  }) +
  geom_line(aes(colour = as.character(Tmean), linetype = as.character(Tmean))) +
  guides(alpha = "none") +
  scale_colour_manual(values = c("black", "#00BA38", "#F8766D")) +
  scale_linetype_manual(values = c("solid", "dashed", "dashed")) +
  labs(x = "Extreme precipitation (mm)", y = "Change in GDP", colour = "Annual temperature (\u00B0C)", linetype = "Annual temperature (\u00B0C)") +
  theme(legend.position = c(0.27, 0.25),
        legend.text=element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.background = element_rect(fill = "transparent", colour = "transparent")) +
  scale_y_continuous(labels = scales::percent)

# write a function to create marginal histograms
make_marginal_distribution <- function(data = df_adm1 %>% filter(is_baseline), variable = NULL) {
  data %>%
    rename(y = variable) %>%
    # use the respective bounds for each variable to ensure that values ranges are identical to the dose-response charts
    filter(y <= df_bounds %>% pull(paste0(variable, "_upper")) %>% setNames(NULL),
           y >= df_bounds %>% pull(paste0(variable, "_lower")) %>% setNames(NULL)) %>%
    ggplot(aes(y)) +
    geom_histogram(position = "identity", alpha = 0.6) +
    labs(fill = NULL, x = variable %>% rewrite_climate_label() %>% str_replace("precip\\.$", "precip. (mm)") %>%
           str_replace("temperature$", "temperature (\u00B0C)") %>% str_replace("var\\.$", "var. (\u00B0C)")) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())
}

# write a function to combine the dose-response function and the marginal histogram
combine_doseresponse_marg <- function(variable = NULL) {
  ggarrange(make_marginal_distribution(variable = variable) + theme(axis.title = element_blank()),
            NULL,
            dose[[variable]],
            nrow = 3, heights = c(1.7, -0.4, 5), align = "hv"
            )
}

# apply the wrapper function to all 6 cliamte indices
figure_dose_marg <- ggarrange(plotlist = map(names(dose), combine_doseresponse_marg), nrow = 2, ncol = 3)

# save out the figure
ggsave(file.path("..", "sharepoint", "Figures", "regressionmodel",
                 paste0(Sys.Date(), " doseresponse_baselinedistribution.png")),
       plot = figure_dose_marg,
       width = 10, height = 7
       )

