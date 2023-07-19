library(dplyr)
library(raster)
library(PBSmapping)

shp <- importShapefile(file.path("..", "sharepoint", "Data", "GADM", "gadm36_levels_simple", "adm1.shp"))
polydata <- attr(shp, 'PolyData')
rr <- raster(file.path("..", "sharepoint", "Data", "Socioeconomic projections", "ISIMIP2b", "gdp_2005soc_0p5deg_annual_2006-2099.nc4"), band=5)

df <- as.data.frame(rr, xy=T)
names(df) <- c('X', 'Y', 'GDP')

## Apply ceiling
topcode <- max(df$GDP[df$GDP < 1e20])
df$GDP[df$GDP >= 1e20] <- topcode
df$EID <- 1:nrow(df)

found <- findPolys(as.EventData(df), shp, maxRows=nrow(df))

## Ensure that every PID is represented
missingpids <- unique(polydata$PID)[!(unique(polydata$PID) %in% df2$PID)]
centroids <- calcCentroid(subset(shp, PID %in% missingpids), rollup=3)

deg2rad <- function(deg) return(deg*pi/180)

gcd.slc <- function(long0, lat0, longs, lats) {
  long0 <- deg2rad(long0)
  lat0 <- deg2rad(lat0)
  longs <- deg2rad(longs)
  lats <- deg2rad(lats)
  
  R <- 6371 # Earth mean radius [km]
  d <- acos(sin(lat0)*sin(lats) + cos(lat0)*cos(lats) * cos(longs-long0)) * R
  return(d) # Distance in km
}

## Assign closest cell to each: add to found so can split
for (ii in 1:nrow(centroids)) {
  closest <- which.min(gcd.slc(centroids$X[ii], centroids$Y[ii], df$X, df$Y))
  found <- rbind(found, data.frame(EID=closest, PID=centroids$PID[ii], SID=centroids$SID[ii], Bdry=NA))
}

found2 <- found %>% left_join(df, by='EID')

## Divide up GDP
found3 <- found2 %>% group_by(EID) %>% summarize(PID=PID, GDP=GDP / length(PID))
found3$GDP[found3$GDP == 0] <- 1

polydata2 <- polydata %>% left_join(found3 %>% group_by(PID) %>% summarize(GDP=sum(GDP)))

write.csv(polydata2, file.path("data", "adm1-weights.csv"), row.names=F)