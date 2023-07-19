# clean the environment
rm(list = ls())

# load packages
library(tidyverse)
library(sf)         # for geodata operations
library(furrr)      # for running map() in parallel
library(rmapshaper) # to simplify polygon shapefiles using ms_simplify()


################################################################################
###################### CREATE SIMPLIFIED SHAPEFILES ############################
################################################################################

# To speed up plotting ADM0-level maps, we create a simplified shapefile version

# load country-level shapefile
sf_gadm0 <- st_read('../sharepoint/Data/Kotz et al - data & code/Data_code_for_Rainfall_economic_production_Kotz_2021/zenodo/shapefiles/gadm36_levels_gpkg/gadm36_levels.gpkg',
                    layer = "level0")

# simplify the ADM0-level shapefile to accelerate the plotting of maps
sf_gadm0_simple <- future_map_dfr(sf_gadm0 %>% group_split(GID_0), ~  ms_simplify(.x, keep = 0.01, keep_shapes = TRUE),
                                  .progress = T)
# inspect the result visually
sf_gadm0_simple %>% ggplot() + geom_sf()
object.size(sf_gadm0_simple)/10^6

# export as RDS file
saveRDS(sf_gadm0_simple, file.path("data", "sf_gadm0_simple.rds"))


################################################################################
########## CREATE AUGMENTED ADM1-LEVEL DATA OBJECTS ############################
################################################################################

# The ADM1-level (geo-)data from GADM does not feature countries that have no ADM1-level subnational regions.
# Among countries in the SSP database, this is the case for the Maldives (MDV), Aruba (ABW) and Kiribati (KIR).
# Therefore, we create an augmented version of the ADM1-level data that features these countries with a single
# (imputed) ADM1-level region.

# load ADM1 shapefile
sf_gadm1 <- st_read('../sharepoint/Data/Kotz et al - data & code/Data_code_for_Rainfall_economic_production_Kotz_2021/zenodo/shapefiles/gadm36_levels_gpkg/gadm36_levels.gpkg',
                    layer = "level1")

# confirm that the three countries are not in it
if(sum(str_count(sf_gadm1$GID_0, "MDV|ABW|KIR")) != 0) stop("There are pattern matches for MDV, ABW or KIR in the ADM1-level shapefile. Please inspect")

# create additional information
sf_gadm1_augmented <- bind_rows(sf_gadm0 %>% filter(GID_0 %in% c("MDV", "ABW", "KIR")) %>% mutate(GID_1 = paste0(GID_0, ".1_1"),
                                                                                                  NAME_1 = paste0(NAME_0, " (imputed)")),
                                sf_gadm1)

# write out as a shapefile
st_write(sf_gadm1_augmented, file.path("..", "sharepoint", "Data", "GADM", "gadm36_levels_simple_augmented", "adm1_augmented.shp"), delete_layer = TRUE)



################################################################################
############################# WEIGHTS CALCULATION ##############################
################################################################################

## a) GDP weights for ADM1-level regions to aggregate to ADM0

# read in the GDP weights produced in a previous script
gdp_weights <- read_csv(file.path("data", "adm1-weights.csv")) 

# add the three missing countries w/o ADM1-level subnational regions w/ a weight of 1
gdp_weights <- bind_rows(gdp_weights,
                         # add the countries with an imputed GDP of 1 (value is irrelevant as the region with have a weight of 100% anyways)
                         sf_gadm1_augmented %>% filter(GID_0 %in% c("ABW", "KIR", "MDV")) %>% as_tibble() %>% select(GID_0, GID_1, NAME_0, NAME_1) %>% mutate(GDP = 1))

# write out the GDP weights
write_csv(gdp_weights, file.path("data", "adm1-weights.csv"))


## b) area weights for ADM1-level regions to aggregate to ADM0

# calculate share in ADM0 area by ADM1 region
area_weights <- sf_gadm1_augmented %>% mutate(area = st_area(geom)) %>%
  group_by(GID_0) %>% mutate(adm1_weight = as.numeric(area/sum(area))) %>% ungroup() %>%
  as_tibble() %>% select(GID_0, NAME_0, GID_1, NAME_1, adm1_weight)

# export area weights as CSV file
write_csv(area_weights, file.path("data", "adm1-weights_areaweighted.csv"))

