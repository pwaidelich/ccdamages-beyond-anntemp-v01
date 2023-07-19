# Climate damage projections beyond annual temperature
Replication materials for Waidelich, Batibeniz, Rising, Kikstra, and Seneviratne (2023). If you have questions about the code or find any errors/bugs in it, please reach out to Paul Waidelich at paul.waidelich@gess.ethz.ch (corresponding author).

NOTE: THIS CLEAN REPOSITORY SERVES FOR THE REVIEW PROCESS ONLY TO ALLOW FOR TRANSPARENCY REGARDING THE SCRIPTS & ANALYSES UNDERLYING THE SUBMITTED MANUSCRIPT. AN EXTENDED REPOSITORY THAT ALSO COVERS LARGE INPUT & OUTPUT DATA FILES (STORED ON THE `sharepoint` FOLDER AND THE EXTERNAL DRIVE MENTIONED IN THE NEXT SECTION OF THIS README) WILL BE PUBLISHED UPON ACCEPTANCE.

## Organization of the overall project repository
Please note that both the underlying climate data and the subnational GDP projections produced require a substantial amount of storage. Therefore, we use a threefold folder structure:
1. The `ccdamages-beyond-anntemp-v01` repo folder, which hosts i) all code, ii) small input files (e.g., country-level GDP projections from the SSP Database), and iii) intermediate data objects (e.g., RDS files). The latter, however, are not pushed to the repo due to their size
2. A `sharepoint` folder on the same hierarchy level as the 'ccdamages-beyond-anntemp-v01' repo, which hosts i) larger input files (e.g., geospatial data on administrative boundaries from GADM) and ii) final outputs, such as figures and tables
3. External drives that host the GDP projections produced (around 2.6TB for 1,000 parameter draws for dose-response function parameters). This drive must be hard-mapped by the user in script `src/01 Create Projections.R` and should provide sufficient storage to host the resulting projection results

### Organization of the `ccdamages-beyond-anntemp-v01` repository
* **data**: intermediate files created and saved out during the analysis and read in again, primarily to avoid expensive re-runs.
* **data/input**: all small input data
* **src**: all code used throughout the analysis
* **src/utils**: helper functions to load the climate data, produce GDP projections, and create figures/tables. Functions in this folder are sourced in the header of the main scripts where required

Note that the repo includes an `ccdamages-beyond-anntemp-v01.Rproj` file, which will automatically map the working directory to the `ccdamages-beyond-anntemp-v01` repo folder. If you do not use the Rproj file, you must set the working directory manually for all R scripts.

### Organization of the `sharepoint` directory
* **Data**: hosts all larger input files, namely i) administrative boundary shapefiles under `Data/GADM`, ii) replication data from Kotz et al., 2022 under `Data/Kotz et al - data & code`, and iii) the processed climate data under `Data/Climate projections`, which are created by the Jupyter Notebooks in the `ccdamages-beyond-anntemp-v01` repo
* **Figures**: contains all figures, which are created by the R scripts in the `ccdamages-beyond-anntemp-v01` repo
* **Tables**: contains all tables, which are created by the R scripts in the `ccdamages-beyond-anntemp-v01` repo

## System requirements
All scripts were executed on a high-performance cluster with 500GB RAM and 48 cores. Note that due to the large filesizes (particularly resulting from the Monte Carlo runs of the GDP projections), the scripts in their original version are likely to max out the memory of a conventional work computer. Depending on the computational resources available to you, you might need to tweak the code to reduce memory requirements, e.g., via chunking parts of the analyses. External storage should be at least 2.75TB to host all GDP projection files for 1,000 Monte Carlo draws for damage function parameters.

## Scripts
To reproduce the results in the paper, run the following scripts in the specified order:
1. `src/create_gdp_weights_adm1.R`: uses gridded GDP data to calculate ADM1-level weights for aggregating data from ADM1 regions to the ADM0 level using GDP weighting. Runtime is about half a day on a high-performance cluster.
2. `src/prepare_weights_and_shapefiles.R`: adds countries that have no ADM1-level subnational region in the GADM database but do feature in the SSP Database (Maldives, Aruba, Kiribati) to the ADM1-level administrative boundary shapefile and other auxiliary files (GDP & area weights) to ensure that they do feature in our projections.
3. `src/aggregate-0.25.ipynb`: aggregates the downscaled CMIP6 data (0.25deg) to ADM1-level subnational regions in the GADM database using the `xagg` package, GADM shapefiles, and area-weighting. Aggregated values for all six climate indicators are saved as netCDF files under the `sharepoint` directory
4. `src/rainfall_deviations.ipynb`: calculates the monthly precipitation deviation, which Kotz et al., 2022 calculate only at the ADM1 level and not at the grid cell level. Results are saved as netCDF files under the `sharepoint` directory
5. `src/00 Estimate Regressions.R`: estimates the main specification from Kotz et al., 2022 (using their replication code) as well as the permutations of their regression model we use throughout our paper. Saves out regression models, coefficients, and variance-covariance matrices as R objects (RDS files) under `ccdamages-beyond-anntemp-v01/data`
6. `src/00 Prepare Analysis.R`: carries out multiple tasks to prepare the GDP projection calculation, such as formatting the global warming level data accordingly or interpolating the SSP data on GDP. Resulting objects are saved out as RDS files under `ccdamages-beyond-anntemp-v01/data`
7. `src/01 Create Projections.R`: creates GDP projections at the ADM1 level for each CMIP6 model-scenario pair for all years between the start of the +0.84C warming level period (our baseline) and 2100 for 1,000 Monte Carlo draws of damage function parameters, aggregates them up to the ADM0 level, and saves them out as `.feather` files on the external drive. Runtime is multiple days even on a high-performance cluster
8. `src/03 Analyse Projections.R`: creates all figures, tables, and numeric results in the paper as well as the Supplementary Information except for our conceptual Figure 1. The script also includes the variance decomposition and the aggregation of ADM0-level results to the global level using current GDP (as per the SSP Database) as country weights. All figures and visualizations are saved out under the `sharepoint` directory
9. `src/04 Analyse Dose-Response.R`: creates the conceptual Figure 1 and saves it out under the `sharepoint` directory

## Data files from external sources (timestamp = download date)
* `data/input/220803 SSP-Database GDP PPP.xlsx`: taken from the SSP Database
* `data/input/230108 World Bank Country Classification.xlsx`: taken from the World Bank Data Help Desk on the World Bank Country and Lending Groups
* `data/input/230503 World Bank GDP_2015USD.csv`: taken from the World Bank's World Development Indicator database
* `data/input/HausfatherEtAl2022 - SI Data - CMIP6 Modelx.xlsx`: taken from Hausfather, Marvel, Schmidt, Nielsen-Gammon, and Zelinka (2022, Nat Clim Change)
* `data/input/cmip6_wl_for_tas_annmean_inclGWL4.xlsx`: global warming level windows calculated following Batibeniz, Hauser, and Seneviratne (2023, Earth Syst Dyn). Filtering this file to the +0.84C warming level window (our baseline period) and reformatting results in the file `data/input/cmip6_wl_0p84_baselines_41years.csv`

## sessionInfo()
```
R version 4.2.2 (2022-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: openSUSE Leap 15.4

Matrix products: default
BLAS/LAPACK: /usr/local/OpenBLAS-0.3.21/lib/libopenblas_haswellp-r0.3.21.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ggcorrplot_0.1.4        ncdf4.helpers_0.3-6     metR_0.13.0             ncdf4_1.19             
 [5] feather_0.3.5           rmapshaper_0.4.6        furrr_0.3.1             future_1.29.0          
 [9] rnaturalearthdata_0.1.0 ggridges_0.5.4          readxl_1.4.1            sf_1.0-9               
[13] ggpubr_0.5.0            forcats_0.5.2           stringr_1.4.1           dplyr_1.0.10           
[17] purrr_0.3.5             readr_2.1.3             tidyr_1.2.1             tibble_3.1.8           
[21] ggplot2_3.4.0           tidyverse_1.3.2         rnaturalearth_0.1.0     ggExtra_0.10.0         
[25] ggpattern_1.0.1        

loaded via a namespace (and not attached):
  [1] googledrive_2.0.0   colorspace_2.0-3    deldir_1.0-6        ggsignif_0.6.4      ellipsis_0.3.2     
  [6] class_7.3-20        htmlTable_2.4.1     base64enc_0.1-3     fs_1.5.2            rstudioapi_0.14    
 [11] httpcode_0.3.0      proxy_0.4-27        listenv_0.8.0       fansi_1.0.3         lubridate_1.9.0    
 [16] xml2_1.3.3          splines_4.2.2       codetools_0.2-18    memuse_4.2-3        cachem_1.0.6       
 [21] knitr_1.41          geojsonlint_0.4.0   Formula_1.2-4       jsonlite_1.8.3      broom_1.0.1        
 [26] cluster_2.1.4       dbplyr_2.2.1        png_0.1-8           shiny_1.7.3         compiler_4.2.2     
 [31] httr_1.4.4          backports_1.4.1     Matrix_1.5-3        assertthat_0.2.1    fastmap_1.1.0      
 [36] gargle_1.2.1        cli_3.4.1           later_1.3.0         htmltools_0.5.3     tools_4.2.2        
 [41] gtable_0.3.1        glue_1.6.2          V8_4.2.2            Rcpp_1.0.9          carData_3.0-5      
 [46] cellranger_1.1.0    vctrs_0.5.1         crul_1.3            xfun_0.35           globals_0.16.2     
 [51] rvest_1.0.3         timechange_0.1.1    mime_0.12           miniUI_0.1.1.1      lifecycle_1.0.3    
 [56] rstatix_0.7.2       googlesheets4_1.0.1 scales_1.2.1        hms_1.1.2           promises_1.2.0.1   
 [61] parallel_4.2.2      RColorBrewer_1.1-3  curl_4.3.3          gridExtra_2.3       memoise_2.0.1      
 [66] dtplyr_1.2.2        rpart_4.1.19        latticeExtra_0.6-30 stringi_1.7.8       jsonvalidate_1.3.2 
 [71] e1071_1.7-12        checkmate_2.1.0     rlang_1.0.6         pkgconfig_2.0.3     lattice_0.20-45    
 [76] htmlwidgets_1.5.4   cowplot_1.1.1       tidyselect_1.2.0    parallelly_1.32.1   plyr_1.8.8         
 [81] magrittr_2.0.3      R6_2.5.1            generics_0.1.3      Hmisc_4.7-2         DBI_1.1.3          
 [86] foreign_0.8-83      pillar_1.8.1        haven_2.5.1         withr_2.5.0         units_0.8-0        
 [91] nnet_7.3-18         survival_3.4-0      abind_1.4-5         sp_1.5-1            modelr_0.1.10      
 [96] crayon_1.5.2        car_3.1-1           interp_1.1-3        KernSmooth_2.23-20  utf8_1.2.2         
[101] tzdb_0.3.0          jpeg_0.1-10         grid_4.2.2          data.table_1.14.6   reprex_2.0.2       
[106] digest_0.6.30       classInt_0.4-8      xtable_1.8-4        httpuv_1.6.6        munsell_0.5.0
```





