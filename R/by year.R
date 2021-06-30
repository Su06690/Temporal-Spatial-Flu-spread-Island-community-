prop_region_yr1<- finalData[, lapply(.SD, sum), by = .(district, season), .SDcols = "prop_recent_travel1"]
prop_region_yr1<- reshape(prop_region_yr1, idvar = "district", timevar = "season", direction = "wide")
prop_region_yr1<- prop_region_yr1 %>% rename(
  "2009/10" =prop_recent_travel1.2010,
  "2010/11" = prop_recent_travel1.2011,
  "2011/12" = prop_recent_travel1.2012,
  "2012/13" =prop_recent_travel1.2013, 
  " 2013/14" =prop_recent_travel1.2014
)
prop_region_yr1[is.na(prop_region_yr1)] <- 0
names(prop_region_yr1)[1] <- "S_NAME"
prop_region_yr1<- melt(prop_region_yr1, id.vars=c("S_NAME"))

prop_region_yr2<- finalData[, lapply(.SD, sum), by = .(district, season), .SDcols = "prop_recent_travel2"]
prop_region_yr2<- reshape(prop_region_yr2, idvar = "district", timevar = "season", direction = "wide")
prop_region_yr2<- prop_region_yr2 %>% rename(
  "2009/10" =prop_recent_travel2.2010,
  "2010/11" = prop_recent_travel2.2011,
  "2011/12" = prop_recent_travel2.2012,
  "2012/13" =prop_recent_travel2.2013, 
  " 2013/14" =prop_recent_travel2.2014)
prop_region_yr2[is.na(prop_region_yr2)] <- 0
names(prop_region_yr2)[1] <- "S_NAME"
prop_region_yr2<- melt(prop_region_yr2, id.vars=c("S_NAME"))

prop_region_yr3<- finalData[, lapply(.SD, sum), by = .(district, season), .SDcols = "prop_recent_travel3"]
prop_region_yr3<- reshape(prop_region_yr3, idvar = "district", timevar = "season", direction = "wide")
prop_region_yr3<- prop_region_yr3 %>% rename(
  "2009/10" =prop_recent_travel3.2010,
  "2010/11" = prop_recent_travel3.2011,
  "2011/12" = prop_recent_travel3.2012,
  "2012/13" =prop_recent_travel3.2013, 
  " 2013/14" =prop_recent_travel3.2014)
prop_region_yr3[is.na(prop_region_yr3)] <- 0
names(prop_region_yr3)[1] <- "S_NAME"
prop_region_yr3<- melt(prop_region_yr3, id.vars=c("S_NAME"))

library(tidyverse)
#prop_reg_yr1<- prop_region_yr1 %>% remove_rownames %>% column_to_rownames(var="district")
#prop_reg_yr2<- prop_region_yr2 %>% remove_rownames %>% column_to_rownames(var="district")
#prop_reg_yr3<- prop_region_yr3 %>% remove_rownames %>% column_to_rownames(var="district")
model1_map<- merge(map1, prop_region_yr1)
model2_map<- merge(map2, prop_region_yr2)
model3_map<- merge(map3, prop_region_yr3)

maps_combined <- rbind(model1_map, model2_map, model3_map)

library(wesanderson)


ggplot(maps_combined) +  geom_sf(aes(fill = value)) + facet_grid(model~variable)  +     
  scale_fill_gradient(low="blue", high="red")+
  #coord_sf(xlim = c(-83.8, -82.2), ylim = c(40.2, 41.3)) +
  theme_classic(base_size = 8)+  theme(axis.text.x = element_blank(),
                                       axis.text.y = element_blank(),
                                       axis.ticks = element_blank())
#rect = element_blank())


















p1<- ggplot(model1_map) +  geom_sf(aes(fill = value)) + facet_grid(~variable)  +     
  scale_fill_gradient(low="blue", high="red")+
  #coord_sf(xlim = c(-83.8, -82.2), ylim = c(40.2, 41.3)) +
  theme_classic(base_size = 9)+  theme(axis.text.x = element_blank(),
                                       axis.text.y = element_blank(),
                                       axis.ticks = element_blank())
                                       #rect = element_blank())

p1



p<- ggplot(model2_map) +  geom_sf(aes(fill = prop_recent_travel2)) + facet_grid(~season)  +     
  scale_fill_gradient2(na.value = "lightgrey", name = "Nb imports",
                       low = "white", mid = "lightblue", high = "darkblue") + 
  #coord_sf(xlim = c(-83.8, -82.2), ylim = c(40.2, 41.3)) +
  theme_classic(base_size = 9)
p<- p+ theme(
  plot.title = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank())
p


