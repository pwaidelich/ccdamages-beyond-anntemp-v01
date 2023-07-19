# NOTE: this script is forked from the 'main_results.R' in the Zenodo repo for the replication of Kotz et al, 2022
# several lines of code have been taken directly from this script to ensure that we fully replicate their regression
# results that inform our damage projections

# clean the environment
rm(list = ls())

# load required packages
library(tidyverse)
library(plm)
library(clubSandwich)
library(stargazer)
library(gsubfn)

# source helper functions
source("src/utils/analyse_impacts.R")
source("src/utils/project_impacts.R")

#function to help stargazer package format in scientific format nicely, adapted from: https://stackoverflow.com/a/56924401
replace_numbers = function(x, low=0.01, high=1e3, digits = 3, scipen=-7, ...) {
  x = gsub(mark,'.',x)
  x.num = as.numeric(x)
  if (x.num<0.01){
    ifelse(
      (x.num >= low) & (x.num < high),
      round(x.num, digits = 2),
      prettyNum(x.num, digits=2, scientific = scipen, ...)
    )
  }
  else {
    ifelse(
      (x.num >= low) & (x.num < high),
      round(x.num, digits = digits),
      prettyNum(x.num, digits=digits, scientific = scipen, ...)
    )
  }
}

#load data from Kotz et al 2022
table <- read.csv(file.path("..", "sharepoint", "Data", "Kotz et al - data & code", "Data_code_for_Rainfall_economic_production_Kotz_2021", "zenodo", "secondary_data", "PT_master.csv"))

# generate quadratic terms for annual total, monthly deviations, and the number of wet days
table$Pt2 <- table$Pt^2
table$Wn_2 <- table$Wn^2
table$wet_days_1_2 <- table$wet_days_1^2

# main econometric specification from Kotz et al., 2022
mod <- plm(dlgdp_pc_usd ~ Tstd + TmeanD + TmeanD:Tmean + TmeanLD + TmeanLD:TmeanL + Pt + Pt2 + Wn + Wn_2 + wet_days_1 + wet_days_1_2 + vwet_days1_am_99p9 + vwet_days1_am_99p9:Tmean,data=table,index=c('ID',"year"),model='within',effect='twoways')
coefs <- coef(mod)
cov_id <- vcovCR(mod,cluster=table$ID,type='CR0')
cov_iso <- vcovCR(mod,cluster=table$iso,type='CR0')
se_iso <- sqrt(diag(cov_iso))

# clean up variable names and save out point estimates for coefficients
names(coefs) <- str_replace(names(coefs), ":", ".")
saveRDS(coefs, file = file.path("data", "kotzetal2022_coefs_main.rds"))

# clean up the clubSandwich objects into regular matrices
clean_up_clubsandwich <- function(clubsand_obj = NULL) {
  # set all attributes to NULL
  for(attr_name in c("type", "cluster", "bread", "v_scale",
                     "est_mats", "adjustments", "target", "inverse_var",
                     "ignore_FE")) attr(clubsand_obj, attr_name) <- NULL
  
  # set class to matrix and clean up rownames
  class(clubsand_obj) <- "matrix"
  rownames(clubsand_obj) <- str_replace(rownames(clubsand_obj), ":", ".")
  colnames(clubsand_obj) <- rownames(clubsand_obj)
  
  return(clubsand_obj)
}
# clean up the covariance matrices
cov_id <- clean_up_clubsandwich(cov_id)
cov_iso <- clean_up_clubsandwich(cov_iso)

# inspect the standard errors and compare manually once against Kotz et al., 2022 main regression in the paper
sqrt(diag(cov_id))
sqrt(diag(cov_iso))
# -> checks out

# ensure that the order is consistent between the objects
if(!identical(rownames(cov_id), names(coefs))) stop("Rownames in cov_id do not match the names in the coefs vector")
if(!identical(rownames(cov_iso), names(coefs))) stop("Rownames in cov_iso do not match the names in the coefs vector")

# export as RDS files
saveRDS(cov_id, file = file.path("data", "kotzetal2022_cov_clustregional_main.rds"))
saveRDS(cov_iso, file = file.path("data", "kotzetal2022_cov_clustnational_main.rds"))


### PART II: estimate permutations of the full Kotz et al model and save them out

# write a function to add zero elements to the vcov matrix
add_empty_elements <- function(coefs_input = NULL, cov_input = NULL, coefs_format_model = coefs, cov_format_model = cov_id) {
  
  if(!is.null(coefs_input) & !is.null(cov_input)) stop("Please only provide either coefs_input or cov_input")
  
  if(!is.null(coefs_input)) {
    # create a vector with same length and element names as coefs
    out <- rep(NA, length(coefs_format_model))
    names(out) <- names(coefs_format_model)
    for(name_jj in names(out)) out[name_jj] <- ifelse(name_jj %in% names(coefs_input), coefs_input[name_jj], 0)
    return(out)
  }
  
  if(!is.null(cov_input)) {
    # create a matrix with same # of rows/cols and rownames/column names as cov_format_model
    out <- matrix(NA, nrow = nrow(cov_format_model), ncol = ncol(cov_format_model))
    rownames(out) <- rownames(cov_format_model)
    colnames(out) <- colnames(cov_format_model)
    
    # overwrite ":" with "." in the row/column names of cov_input
    rownames(cov_input) <- rownames(cov_input) %>% str_replace_all("\\:", "\\.")
    colnames(cov_input) <- colnames(cov_input) %>% str_replace_all("\\:", "\\.")
    
    # if rowname and colname of out feature in cov_input, transcribe this value - otherwise set it to zero
    for(rowname_jj in rownames(out)) {
      for(colname_jj in colnames(out)) {
        out[rownames(out) == rowname_jj, colnames(out) == colname_jj] <- ifelse(rowname_jj %in% rownames(cov_input) &
                                                                                colname_jj %in% colnames(cov_input),
                                                                                cov_input[rownames(cov_input) == rowname_jj,
                                                                                          colnames(cov_input) == colname_jj],
                                                                                0)
      }
    }
    
    return(out)
  }
}

# create a list storing the models, their coefficients, and their variance-covariance matrices
formula_list <- list()
mod_list <- list()
coef_list <- list()
cov_id_list <- list()
cov_iso_list <- list()

# create all the formulas of interest by iteratively adding regressors and name them accordingly
tempprecip_element <- "dlgdp_pc_usd ~ TmeanD + TmeanD:Tmean + TmeanLD + TmeanLD:TmeanL + Pt + Pt2"
tempprecip_burke <- "dlgdp_pc_usd ~ Tmean + Tmean2 + Pt + Pt2"
formula_list <- list("T_Pt" = as.formula(tempprecip_element),
                     "T_Pt_Tstd" = as.formula(paste0(tempprecip_element, " + Tstd")),
                     "T_Pt_Tstd_99p9" = as.formula(paste0(tempprecip_element, " + Tstd + vwet_days1_am_99p9 + vwet_days1_am_99p9:Tmean")),
                     "T_Pt_Tstd_Wn" = as.formula(paste0(tempprecip_element, " + Tstd + Wn + Wn_2")),
                     "T_Pt_Tstd_wetdays" = as.formula(paste0(tempprecip_element, " + Tstd + wet_days_1 + wet_days_1_2")),
                     "T_Pt_Tstd_wetdays_Wn" = as.formula(paste0(tempprecip_element, " + Tstd + wet_days_1 + wet_days_1_2 + Wn + Wn_2")),
                     "T_Pt_Tstd_wetdays_Wn_99p9" = as.formula(paste0(tempprecip_element, " + Tstd + wet_days_1 + wet_days_1_2 + Wn + Wn_2 + vwet_days1_am_99p9 + vwet_days1_am_99p9:Tmean")),
                     "T_Pt_Burke" = as.formula(tempprecip_burke))

# estimate models by looping through the formula list
for(specname in names(formula_list)) {
  # NOTE: we add a squared term of annual mean temperature to the data for the Burke-like specification
  mod_list[[specname]] <- plm(formula_list[[specname]],data=table %>% mutate(Tmean2 = Tmean^2), index=c('ID',"year"),model='within',effect='twoways')
  
  # store coefficients and vcov matrices (for national- and regional clustering) into the list objects and add empty elements (= climate indicators that do not feature in the model)
  coef_list[[specname]] <- add_empty_elements(coefs_input = coef(mod_list[[specname]]) %>% setNames(names(.) %>% str_replace(":", ".")))
  cov_id_list[[specname]] <- add_empty_elements(cov_input = vcovCR(mod_list[[specname]],cluster=table$ID,type='CR0')) %>% clean_up_clubsandwich()
  cov_iso_list[[specname]] <- add_empty_elements(cov_input = vcovCR(mod_list[[specname]],cluster=table$iso,type='CR0')) %>% clean_up_clubsandwich()
}

# ensure that the sample size is the same for all regression models
map_dbl(mod_list, .f = ~ length(.x$residuals))
length(mod$residuals)
# -> minor sample size deviation for the Burke-like specification because we do not use first-diffs, so we gain sample size
# -> no adjustment required because the Burke-like specification is not used in the paper as of now

# ensure that the full model is identical to the Kotz model estimated above
if(!all.equal(coef(mod) %>% sort(), coef(mod_list[["T_Pt_Tstd_wetdays_Wn_99p9"]]) %>% sort())) stop("Coefficients for the main specification and the last permutation do no match")
if(!all.equal(cov_iso, cov_iso_list$T_Pt_Tstd_wetdays_Wn_99p9)) stop("Country-clustered variance-covariance matrices for the main specification and the last permutation do no match")

# add modifications of some specs where we set the coef/covariance for annual precip to zero (as this variable is often not projected out but only used as control)
for(specname_jj in c("T_Pt", "T_Pt_Tstd", "T_Pt_Burke")) {
  if(!specname_jj %in% names(coef_list)) stop("specname_jj is not an element of coef_list. Please inspect")
  
  # extract the coefficient vector, set all values for Pt and Pt2 to zero and add to coef_list with the suffix "_Pt2zero"
  coef_jj <- coef_list[[specname_jj]]
  coef_jj[names(coef_jj) %in% c("Pt", "Pt2")] <- 0
  coef_list[[paste0(specname_jj, "_Pt2zero")]] <- coef_jj
  rm(coef_jj)
  
  # repeat for the region-clustered variance-covariance matrix
  cov_id_jj <- cov_id_list[[specname_jj]]
  # set all columns to zero
  cov_id_jj[, colnames(cov_id_jj) %in% c("Pt", "Pt2")] <- 0
  # do the same for the rows
  cov_id_jj[rownames(cov_id_jj) %in% c("Pt", "Pt2"),]  <- 0
  cov_id_list[[paste0(specname_jj, "_Pt2zero")]] <- cov_id_jj
  rm(cov_id_jj)
  
  # repeat for the country-clustered variance-covariance matrix
  cov_iso_jj <- cov_iso_list[[specname_jj]]
  cov_iso_jj[, colnames(cov_iso_jj) %in% c("Pt", "Pt2")] <- 0
  cov_iso_jj[rownames(cov_iso_jj) %in% c("Pt", "Pt2"),]  <- 0
  cov_iso_list[[paste0(specname_jj, "_Pt2zero")]] <- cov_iso_jj
  rm(cov_iso_jj)
}

# save out the coefficient and the national-clustered vcov 
saveRDS(coef_list, file = file.path("data", "coef_list.rds"))
saveRDS(cov_id_list, file = file.path("data", "cov_id_list.rds"))
saveRDS(cov_iso_list, file = file.path("data", "cov_iso_list.rds"))

# export the respective regression tables
mod_for_table <- list(m1 = mod_list$T_Pt,
                      m2 = mod_list$T_Pt_Tstd,
                      m3 = mod_list$T_Pt_Tstd_wetdays,
                      m4 = mod_list$T_Pt_Tstd_wetdays_Wn,
                      m5 = mod_list$T_Pt_Tstd_wetdays_Wn_99p9)
se_for_table <- list(m1 = sqrt(diag(cov_iso_list$T_Pt)) %>% setNames(names(.) %>% str_replace("\\.", ":")),
                     m2 = sqrt(diag(cov_iso_list$T_Pt_Tstd)) %>% setNames(names(.) %>% str_replace("\\.", ":")),
                     m3 = sqrt(diag(cov_iso_list$T_Pt_Tstd_wetdays)) %>% setNames(names(.) %>% str_replace("\\.", ":")),
                     m4 = sqrt(diag(cov_iso_list$T_Pt_Tstd_wetdays_Wn)) %>% setNames(names(.) %>% str_replace("\\.", ":")),
                     m5 = sqrt(diag(cov_iso_list$T_Pt_Tstd_wetdays_Wn_99p9)) %>% setNames(names(.) %>% str_replace("\\.", ":")))

# print out the table in a clean format despite high & heterogeneous number of digits
mark  = '::::'
star = stargazer(mod_for_table,se = se_for_table,dep.var.labels=c("Subnational income growth"),column.labels=NULL,omit.stat="F",align=TRUE,digits=10,decimal.mark  = mark, digit.separator = '')
cat(gsubfn(paste0("([0-9.\\-]+", mark, "[0-9.\\-]+)"), ~replace_numbers(x), star), sep='\n')

# clean up
rm(mod_for_table, se_for_table, star)


################################################################################
##################### ADDITIONAL REGRESSIONS ###################################
################################################################################

## a) use real GDPpc instead of nominal following Callahan & Mankin 2022

# load the GDP_deflator variable for the US and implement following Callahan & Mankin 2022
df_gdpdeflator <- read_csv(file.path("..", "sharepoint", "Data", "CallahanMankin_ExtremeHeatEconomics_2022", "Data", "GDP_Deflator", "API_NY.GDP.DEFL.ZS_DS2_en_excel_v2_4353530.csv"))
df_deflator_us <- df_gdpdeflator %>% filter(`Country Code` == "USA") %>% select(-`Country Name`, -starts_with("Indicator")) %>%
                        gather(key = "year", value = "gdp_deflator", -c("Country Code")) %>% rename(GID_0 = "Country Code")

# extract 2010 value and calculate the factor by which we need to multiply to inflate/deflate to 2010USD
us_gdpdeflator_2010 <- df_deflator_us %>% filter(year == 2010) %>% pull(gdp_deflator)
df_deflator_us <- df_deflator_us %>% mutate(conversion_to_usd2010 = us_gdpdeflator_2010 / gdp_deflator,
                                            year = as.integer(year))

# clean up
rm(df_gdpdeflator, us_gdpdeflator_2010)

# quickly check that we can reproduce the dependent variable of Kotz et al 2022 (dlgdp_pc_usd) by taking ADM1-specific
# first-difference of log-transformed GDPpc
# (serves to rule out that data is not structured such that we cannot use lag() functions)
table %>% group_by(ID) %>% mutate(id_checking = 1:n(),
                                  lag_checking = if_else(id_checking == 1, NA_real_, dplyr::lag(lgdp_pc_usd)),
                                  dlgdp_pc_usd_checking = lgdp_pc_usd - lag_checking,
                                  diff_checking = dlgdp_pc_usd_checking - dlgdp_pc_usd) %>%
  select(id_checking, lgdp_pc_usd, lag_checking, dlgdp_pc_usd, dlgdp_pc_usd_checking, diff_checking) %>%
  summary()
# -> difference is basically zero (although not precisely due to precision issues)

# merge into Kotz et al 2022 data
table <- table %>% left_join(df_deflator_us %>% select(year, conversion_to_usd2010), by = "year") %>%
          # convert nominal GDPpc in current USD to 2010USD using the GDP deflator ratio for the US
          # NOTE: the Kotz et al 2022 data features some zero values for a Turkish region in gdp_pc_usd (with no implications for their analysis)
          #  to avoid numerical issues arising from taking log(0), we set these values to NA instead
          mutate(gdp_pc_2010 = if_else(gdp_pc_usd != 0, gdp_pc_usd*conversion_to_usd2010, NA_real_),
                 # create dependent variable by calculating logs
                 lgdp_pc_2010 = log(gdp_pc_2010)) %>%
          # take first-differences
          group_by(ID) %>% mutate(dlgdp_pc_2010 = lgdp_pc_2010 - dplyr::lag(lgdp_pc_2010, n = 1)) %>% ungroup()


# compare the two measures
table %>% select(dlgdp_pc_2010, dlgdp_pc_usd) %>% summary()
# -> very similar range but median and mean are considerably lower once we adjust for inflation

# plot the difference over time
table %>% filter(year >= 1979) %>% ggplot(aes(year, dlgdp_pc_2010 - dlgdp_pc_usd)) + geom_point() +
  labs(y = "Log-scale effect of inflation-adjustment on ADM1-level GDPpc growth\n (converted to 2010USD via US GDP Deflator)",
       x = NULL) + theme_classic()
ggsave(file.path("..", "sharepoint", "Figures", "regressionmodel",
                 paste0(Sys.Date(), " inflation_adjusted_gdpcgrowth.png")),
       width = 5, height = 5)
# -> basically a constant for each year, so should be addressed via year FEs

# repeat the main specification regression using inflation-adjuste GDPpc growth (in log-scale) as dep. var.
mod_infladjusted <- plm(dlgdp_pc_2010 ~ Tstd + TmeanD + TmeanD:Tmean + TmeanLD + TmeanLD:TmeanL + Pt + Pt2 + Wn + Wn_2 + wet_days_1 + wet_days_1_2 + vwet_days1_am_99p9 + vwet_days1_am_99p9:Tmean,data=table,index=c('ID',"year"),model='within',effect='twoways')
cov_iso_infladjusted <- vcovCR(mod_infladjusted,cluster=table$iso,type='CR0')
se_iso_infladjusted <- sqrt(diag(cov_iso_infladjusted))

# ensure that order of variables is identical across coef()-output and the 'se_...' objects
if(!identical(names(coef(mod)), names(se_iso))) stop("order of climate indicators in coef(mod) and se_iso is not compatible")
if(!identical(names(coef(mod_infladjusted)), names(se_iso_infladjusted))) stop("order of climate indicators in coef(mod_infladjusted) and se_iso_infladjusted is not compatible")

# combine results and plot them incl. 95% CI
bind_rows(tibble(term = names(coef(mod)),
                 estimate = coef(mod),
                 std.error = se_iso,
                 spec = "Kotz et al. 2022\n (main specification)"),
          tibble(term = names(coef(mod_infladjusted)),
                 estimate = coef(mod_infladjusted),
                 std.error = se_iso_infladjusted,
                 spec = "Inflation-adjustment\n (via US GDP Deflator)")
) %>%
  mutate(ci_low = estimate - 1.96*std.error,
         ci_high = estimate + 1.96*std.error) %>%
  ggplot(aes(term, estimate)) + geom_point(aes(colour = spec, group = spec), position = position_dodge2(width = 0.9)) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high, colour = spec), position = position_dodge2()) +
  labs(colour = NULL, y = "Estimated coefficient and 95% CI", x = "Independent variable") +
  scale_colour_manual(values = c("#619CFF", "#F8766D"), guide = guide_legend(reverse = TRUE)) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = c(0.3, 0.2))
ggsave(file.path("..", "sharepoint", "Figures", "regressionmodel",
                  paste0(Sys.Date(), " inflation_adjusted_regression.png")),
        width = 5, height = 5)

# write out the inflation-adjusted regression table
star_infladjusted <- stargazer(mod, mod_infladjusted, se = list(se_iso, se_iso_infladjusted),dep.var.labels=c("Subnational income growth"),column.labels=NULL,omit.stat="F",align=TRUE,digits=10,decimal.mark  = mark, digit.separator = '')
cat(gsubfn(paste0("([0-9.\\-]+", mark, "[0-9.\\-]+)"), ~replace_numbers(x), star_infladjusted), sep='\n')

