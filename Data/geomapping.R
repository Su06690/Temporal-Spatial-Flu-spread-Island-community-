library(tidyverse)
library(ggmap)
library(sp)
library(sf)
library(rnaturalearth)
library(maps)
library(mapdata)

map(col="grey80", border = "grey40", fill = TRUE, xlim = c(10, 36), ylim = c(50, 68), mar = rep(0.1, 4))
box()
points(32.9871639,	129.0810861,col=2,pch=18)

# Load the sp package
library(sp)


geocoded <- data.frame(stringsAsFactors = FALSE)

school <- read_excel("data/institution_population_positive_cases_new.xlsx",sheet = "merged")
sources <- distinct(school,address )
coord <- select(school, lon, lat)
school<- as.data.frame(school)



library(readxl)
school <- read_excel("data/school_location.xlsx",sheet = "merged")
locations_df <- read_excel("data/locations.xlsx",sheet = "merged")
locations <- as_tibble(locations_df)
library(sf)
library(mapview)
locations_sf <- st_as_sf(locations, coords = c("lon", "lat"), crs = 4326)
mapview(locations_sf)



location <- read_excel("data/school_location.xlsx",sheet = "merged")
location <- as.data.frame(location)
coord <- select(location, lon, lat)
points_sp <- SpatialPoints(coords = coord, proj4string = CRS("+proj=longlat +datum=WGS84"))
points_spdf <- SpatialPointsDataFrame(coord = coord,
                                      data = location,  
                                      proj4string = CRS("+proj=longlat +datum=WGS84"))
palette(alpha(c("darkorchid", "darkorange"), 0.7))
par(mar = c(1, 1, 3, 1))
# Plot points
plot(points_spdf,
     pch = 20,
     col = points_spdf$type,
     cex = sqrt(points_spdf$n)/2 + 0.25)




locations_df <- as_tibble(location)
locations_sf <- st_as_sf(locations_df, coords = c("lon", "lat"), crs = 4326)





library(mapview)
mapview(locations_sf)




locations_sf <- st_as_sf(school$address, coords = c("logitude", "latitude"), crs = 4326)
geocode_OSM("shinkamigto, Japan")
