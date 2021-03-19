library(data.table)
finalData<- as.data.table(finalData)
finalData[, prop_recent_travel1 := summary(out1, burnin = 5000)$tree$import]
# Add the proportion of iterations in model 2  & 3where each case is an import
finalData[, prop_recent_travel2 := summary(out2, burnin = 5000)$tree$import]
finalData[, prop_recent_travel3 := summary(out3, burnin = 5000)$tree$import]

## ----n_import_reg_per_reg, warning = FALSE-----------------------------------------------------------------------
# Number of imports per region in model 1
data_2010_11<- finalData %>% filter (season=="2010/11")
data_2011_12<- finalData %>% filter (season=="2011/12")
data_2012_13<- finalData %>% filter (season=="2012/13")
data_2013_14<- finalData %>% filter (season=="2013/14")


prop_reg_year1 <- data_2010_11[, .(prop_per_reg = sum(prop_recent_travel1)), 
                       by = district]$prop_per_reg
prop_reg_year2 <- data_2011_12[, .(prop_per_reg = sum(prop_recent_travel1)), 
                               by = district]$prop_per_reg
prop_reg_year3 <- data_2012_13[, .(prop_per_reg = sum(prop_recent_travel1)), 
                               by = district]$prop_per_reg
prop_reg_year4 <- data_2013_14[, .(prop_per_reg = sum(prop_recent_travel1)), 
                               by = district]$prop_per_reg

names(prop_reg_year1) <- unique(data_2010_11$district)
names(prop_reg_year2) <- unique(data_2011_12$district)
names(prop_reg_year3) <- unique(data_2012_13$district)
names(prop_reg_year4) <- unique(data_2013_14$district)

library(ggplot2)
library(sf)
mapyr1 <- st_read(dsn="/Users/suhan/Project/kamigoto/Kamigoto shape file", layer="kamigoto")
mapyr1$X_CODE <- as.numeric(mapyr1$X_CODE)
mapyr1$Y_CODE <- as.numeric(mapyr1$Y_CODE)
mapyr2 <- mapyr1
mapyr3 <- mapyr1
mapyr4 <- mapyr1
mapyr1$season <- "2010/11"
mapyr2$season <- "2011/12"
mapyr3$season <- "2012/13"
mapyr4$season <- "2013/14"

mapyr1$prop_reg_year <- prop_reg_year1[as.character(mapyr1$S_NAME)]
mapyr2$prop_reg_year <- prop_reg_year2[as.character(mapyr2$S_NAME)]
mapyr3$prop_reg_year <- prop_reg_year3[as.character(mapyr3$S_NAME)]
mapyr4$prop_reg_year <- prop_reg_year4[as.character(mapyr4$S_NAME)]

#map_years <- rbind(mapyr1, mapyr2, mapyr3, mapyr4)


# Plot: number of imports per region, two panels
yr1<- ggplot(mapyr1) +  geom_sf(aes(fill = prop_reg_year)) +  
  scale_fill_gradient2(na.value = "lightgrey", breaks = seq(0, 40, 5), name = "Nb imports",
                       low = "white", mid = "lightblue", high = "darkblue") + 
  theme_classic(base_size = 9)
yr1


yr2<- ggplot(mapyr2) +  geom_sf(aes(fill = prop_reg_year)) +  
  scale_fill_gradient2(na.value = "lightgrey", breaks = seq(0, 40, 5),name = "Nb imports",
                       low = "white", mid = "lightblue", high = "darkblue") + 
  theme_classic(base_size = 9)
yr2

yr3<- ggplot(mapyr3) +  geom_sf(aes(fill = prop_reg_year)) +  
  scale_fill_gradient2(na.value = "lightgrey",breaks = seq(0, 40, 5), name = "Nb imports",
                       low = "white", mid = "lightblue", high = "darkblue") + 
  theme_classic(base_size = 9)
yr3

yr4<- ggplot(mapyr4) +  geom_sf(aes(fill = prop_reg_year)) +  
  scale_fill_gradient2(na.value = "lightgrey",breaks = seq(0, 40, 5), name = "Nb imports",
                       low = "white", mid = "lightblue", high = "darkblue") + 
  theme_classic(base_size = 9)
yr4


library("ggpubr")
library(patchwork)

figure1 <- yr1+yr2+yr3+yr4
figure1


####

#model2

prop_reg_year1_model2 <- data_2010_11[, .(prop_per_reg = sum(prop_recent_travel2)), 
                               by = district]$prop_per_reg
prop_reg_year2_model2 <- data_2011_12[, .(prop_per_reg = sum(prop_recent_travel2)), 
                               by = district]$prop_per_reg
prop_reg_year3_model2 <- data_2012_13[, .(prop_per_reg = sum(prop_recent_travel2)), 
                               by = district]$prop_per_reg
prop_reg_year4_model2 <- data_2013_14[, .(prop_per_reg = sum(prop_recent_travel2)), 
                               by = district]$prop_per_reg

names(prop_reg_year1_model2) <- unique(data_2010_11$district)
names(prop_reg_year2_model2) <- unique(data_2011_12$district)
names(prop_reg_year3_model2) <- unique(data_2012_13$district)
names(prop_reg_year4_model2) <- unique(data_2013_14$district)

library(ggplot2)
library(sf)
mapyr1_model2 <- st_read(dsn="/Users/suhan/Project/kamigoto/Kamigoto shape file", layer="kamigoto")
mapyr1_model2$X_CODE <- as.numeric(mapyr1_model2$X_CODE)
mapyr1_model2$Y_CODE <- as.numeric(mapyr1_model2$Y_CODE)
mapyr2_model2<- mapyr1_model2
mapyr3_model2 <- mapyr1_model2
mapyr4_model2 <- mapyr1_model2
mapyr1_model2$season <- "2010/11"
mapyr2_model2$season <- "2011/12"
mapyr3_model2$season <- "2012/13"
mapyr4_model2$season <- "2013/14"

mapyr1_model2$prop_reg_year <- prop_reg_year1_model2 [as.character(mapyr1_model2$S_NAME)]
mapyr2_model2$prop_reg_year <- prop_reg_year2_model2 [as.character(mapyr2_model2$S_NAME)]
mapyr3_model2$prop_reg_year <- prop_reg_year3_model2 [as.character(mapyr3_model2$S_NAME)]
mapyr4_model2$prop_reg_year <- prop_reg_year4_model2 [as.character(mapyr4_model2$S_NAME)]

#map_years <- rbind(mapyr1, mapyr2, mapyr3, mapyr4)


# Plot: number of imports per region, two panels
yr1_model2<- ggplot(mapyr1_model2) +  geom_sf(aes(fill = prop_reg_year)) +  
  scale_fill_gradient2(na.value = "lightgrey", breaks = seq(0, 40, 5), name = "Nb imports",
                       low = "white", mid = "lightblue", high = "darkblue") + 
  theme_classic(base_size = 9)
yr1_model2


yr2_model2<- ggplot(mapyr2_model2) +  geom_sf(aes(fill = prop_reg_year)) +  
  scale_fill_gradient2(na.value = "lightgrey", breaks = seq(0, 40, 5),name = "Nb imports",
                       low = "white", mid = "lightblue", high = "darkblue") + 
  theme_classic(base_size = 9)
yr2_model2

yr3_model2<- ggplot(mapyr3_model2) +  geom_sf(aes(fill = prop_reg_year)) +  
  scale_fill_gradient2(na.value = "lightgrey",breaks = seq(0, 40, 5), name = "Nb imports",
                       low = "white", mid = "lightblue", high = "darkblue") + 
  theme_classic(base_size = 9)
yr3_model2

yr4_model2<- ggplot(mapyr4_model2) +  geom_sf(aes(fill = prop_reg_year)) +  
  scale_fill_gradient2(na.value = "lightgrey",breaks = seq(0, 40, 5), name = "Nb imports",
                       low = "white", mid = "lightblue", high = "darkblue") + 
  theme_classic(base_size = 9)
yr4_model2


library("ggpubr")
library(patchwork)

figure2 <- yr1_model2+yr2_model2+yr3_model2+yr4_model2
figure2


##model3

prop_reg_year1_model3 <- data_2010_11[, .(prop_per_reg = sum(prop_recent_travel3)), 
                                      by = district]$prop_per_reg
prop_reg_year2_model3 <- data_2011_12[, .(prop_per_reg = sum(prop_recent_travel3)), 
                                      by = district]$prop_per_reg
prop_reg_year3_model3 <- data_2012_13[, .(prop_per_reg = sum(prop_recent_travel3)), 
                                      by = district]$prop_per_reg
prop_reg_year4_model3 <- data_2013_14[, .(prop_per_reg = sum(prop_recent_travel3)), 
                                      by = district]$prop_per_reg

names(prop_reg_year1_model3) <- unique(data_2010_11$district)
names(prop_reg_year2_model3) <- unique(data_2011_12$district)
names(prop_reg_year3_model3) <- unique(data_2012_13$district)
names(prop_reg_year4_model3) <- unique(data_2013_14$district)

library(ggplot2)
library(sf)
mapyr1_model3 <- st_read(dsn="/Users/suhan/Project/kamigoto/Kamigoto shape file", layer="kamigoto")
mapyr1_model3$X_CODE <- as.numeric(mapyr1_model3$X_CODE)
mapyr1_model3$Y_CODE <- as.numeric(mapyr1_model3$Y_CODE)
mapyr2_model3<- mapyr1_model3
mapyr3_model3 <- mapyr1_model3
mapyr4_model3 <- mapyr1_model3
mapyr1_model3$season <- "2010/11"
mapyr2_model3$season <- "2011/12"
mapyr3_model3$season <- "2012/13"
mapyr4_model3$season <- "2013/14"

mapyr1_model3$prop_reg_year <- prop_reg_year1_model3[as.character(mapyr1_model3$S_NAME)]
mapyr2_model3$prop_reg_year <- prop_reg_year2_model3[as.character(mapyr2_model3$S_NAME)]
mapyr3_model3$prop_reg_year <- prop_reg_year3_model3[as.character(mapyr3_model3$S_NAME)]
mapyr4_model3$prop_reg_year <- prop_reg_year4_model3[as.character(mapyr4_model3$S_NAME)]

#map_years <- rbind(mapyr1, mapyr2, mapyr3, mapyr4)


# Plot: number of imports per region, two panels
yr1_model3<- ggplot(mapyr1_model3) +  geom_sf(aes(fill = prop_reg_year)) +  
  scale_fill_gradient2(na.value = "lightgrey", breaks = seq(0, 40, 5), name = "Nb imports",
                       low = "white", mid = "lightblue", high = "darkblue") + 
  theme_classic(base_size = 9)
yr1_model3


yr2_model3<- ggplot(mapyr2_model3) +  geom_sf(aes(fill = prop_reg_year)) +  
  scale_fill_gradient2(na.value = "lightgrey", breaks = seq(0, 40, 5),name = "Nb imports",
                       low = "white", mid = "lightblue", high = "darkblue") + 
  theme_classic(base_size = 9)
yr2_model3

yr3_model3<- ggplot(mapyr3_model3) +  geom_sf(aes(fill = prop_reg_year)) +  
  scale_fill_gradient2(na.value = "lightgrey",breaks = seq(0, 40, 5), name = "Nb imports",
                       low = "white", mid = "lightblue", high = "darkblue") + 
  theme_classic(base_size = 9)
yr3_model3

yr4_model3<- ggplot(mapyr4_model3) +  geom_sf(aes(fill = prop_reg_year)) +  
  scale_fill_gradient2(na.value = "lightgrey",breaks = seq(0, 40, 5), name = "Nb imports",
                       low = "white", mid = "lightblue", high = "darkblue") + 
  theme_classic(base_size = 9)
yr4_model3


library("ggpubr")
library(patchwork)

figure3 <- yr1_model3+yr2_model3+yr3_model3+yr4_model3
figure3


