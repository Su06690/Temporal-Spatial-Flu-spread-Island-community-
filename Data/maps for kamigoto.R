library(NipponMap)
library(jpndistrict)
library(divagis)
library(cshapes)


library(maps)
library(sf)
library(ggplot2)
library(rgdal)
library(rgeos)

#Reading the shapefiile
#sourcefile for download; https://www.diva-gis.org/datadown    
sf <- st_read(dsn="/Users/suhan/Downloads/JPN_adm", layer="JPN_adm2")
shape <- readOGR(dsn="/Users/suhan/Downloads/JPN_adm", layer="JPN_adm2")

#To view the attributes
head(shape@data)
summary(sf)

#Plotting the shapefile
plot(shape)
plot(sf)

#Plotting the districts only
plot(sf["NAME_2"], axes = TRUE, main = "Districts")





JapanPrefMap()
if (requireNamespace("RColorBrewer", quietly = TRUE)) { 
  cols <- rev(RColorBrewer::brewer.pal(8,"Set2"))
}else{
  cols <- sample(colours(), 47)
}

JapanPrefMap(col = cols, border = gray(.8), axes = TRUE)

if (requireNamespace("foreign", quietly = TRUE)) {
  dat <- foreign::read.dbf(system.file("shapes/jpn.dbf", package="NipponMap")) 
  op <- par(bg = "skyblue")
  p <- JapanPrefMap(col = "ivory")
  col <- c("olivedrab4", "olivedrab1")
  pop <- dat$population / 1e+7
  symbols(p, circles = sqrt(pop / (2 * pi)), inches = FALSE,
          fg = col[1], bg = col[2], add = TRUE)
  idx <- c(1e+6, 5e+6, 1e+7)
  pos <- legend("bottomright", legend = format(idx, scientific = 10, big.mark = ","),
                title = "Population (2010)", bg = "white", x.intersp = 2, y.intersp = 1.5) symbols(pos$text$x - 1, pos$text$y, circles = sqrt(idx / 1e+7 / (2 * pi)),
                                                                                                   inches = FALSE, fg = col[1], bg = col[2], add = TRUE)
  par(op)
}



nagasaki<- jpn_pref(pref_code = 42, district = TRUE)
kamigoto<- jpn_cities(42411)





