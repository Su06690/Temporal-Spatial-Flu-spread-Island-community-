rm(list = ls())

library(o2geosocial)

load ("data/kamigoto_combined_data.RData")

## ----distrib-----------------------------------------------------------------------------------------------------
# Distribution of the latent/incubation period  average 1.63 ± 0.06 days and standard deviation (sd) 0.26 ± 0.08 days A.cori etal (https://doi.org/10.1016/j.epidem.2012.06.001)
# 2 days (1-4 days range)
f_dens <- dgamma(x = 1:100, scale = 0.04147239, shape = 39.30325)
#f_dens<- f_dens/sum(f_dens)

#?distribution of infectious period average 0.99 ± 0.25 days and sd 0.96 ± 0.15 days , A.cori etal (https://doi.org/10.1016/j.epidem.2012.06.001)


# Distribution of the generation time 10.1371/journal.pone.0075339  mean 3.2 (2.4-3.9)
#SD is calculated  from the formula from cochrone http://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm 
w_dens <- dnorm(x = 1:100, mean = 3.2, sd = 0.4)
# a mean serial interval of 3.6 days (95% confidence interval = 2.9–4.3 days), with standard deviation 1.6 days

## ----age, message = FALSE, warning = FALSE-----------------------------------------------------------------------
# either jpoly2 created or new social contacft study of nishiura etal
library(readxl)
#currently using weekday socail contact
social_contact<- read_excel("data/contact_weekday.xlsx", 
                            col_types = c("numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric", 
                                          "numeric", "numeric", "numeric"))

social_contact <-t(social_contact)
social_contact <-data.table::as.data.table(social_contact)
# Compute the proportion of connection to each age group
a_dens <- t(t(social_contact)/colSums(social_contact))


## ----distance_pop, warning = FALSE-------------------------------------------------------------------------------
# Extract all regions in the territory
dt_regions<- read.csv("data/geolocations.csv")
school_location<- read.csv("data/school_location.csv")
school_location<- school_location %>%  rename(S_NAME = missing_name, district_id=area.number.of.the.school,district_id=area.number.of.the.school)

vaccination<- read.csv("data/vaccination.csv")
vaccination <- vaccination %>% filter(Year!= "2009/2010")


# Create the matrices of coordinates for each region (one "from"; one "to")
mat_dist_from <- matrix(c(rep(dt_regions$lon, nrow(dt_regions)),
                          rep(dt_regions$lat, nrow(dt_regions))), ncol = 2)
mat_dist_to <- matrix(c(rep(dt_regions$lon, each = nrow(dt_regions)), 
                        rep(dt_regions$lat, each = nrow(dt_regions))),
                      ncol = 2)

# Compute all the distances between the two matrices
all_dist <- geosphere::distGeo(mat_dist_from, mat_dist_to)
# Compile into a distance matrix
dist_mat <- matrix(all_dist/1000, nrow = nrow(dt_regions))

# Extract the population vector
pop_vect <- dt_regions$population

# Rename the matrix columns and rows, and the population vector
names(pop_vect) <- rownames(dist_mat) <- colnames(dist_mat) <- 
  dt_regions$S_NAME


##age_group accroding with social contact data
labs <- c(paste(seq(0, 65, by = 5),seq(0 + 5 - 1, 70 - 1, by = 5),
                sep = "-"), paste(70, "+", sep = ""))
labs
df_all_positive$age_gp<- cut(df_all_positive$age, breaks = c(seq(0, 70, by = 5), Inf), labels = labs, right = FALSE)

df_all_positive$age_gp<- as.numeric(df_all_positive$age_gp)
## ----first_model, message = FALSE, warning = FALSE---------------------------------------------------------------
# Set movement, likelihood and prior lists to default
moves <- custom_moves()
likelihoods <- custom_likelihoods()
priors <- custom_priors()

##***added the folder district 
###removing the variables with NA for the district

#correct district in data_clean_df_all_positive first
finalData<-subset(df_all_positive,!(is.na(df_all_positive ["district"])))
finalData<- finalData[order(as.Date(finalData$visit_date, format="%Y/%m/%d")),]

n_missing_case<- sum(is.na(finalData$onset_date))

dist_onset_visit<- table(finalData$visit_date - finalData$onset_date)/sum(!is.na(finalData$onset_date))

library(tidyverse)
inferred_duration<- sample(x=names(dist_onset_visit)%>% as.numeric,
                           size= n_missing_case,
                           replace = T,
                           prob=dist_onset_visit)

finalData[is.na(finalData$onset_date),
          "onset_date"]<- finalData[is.na(finalData$onset_date), "visit_date"]-
  inferred_duration

finalData<- finalData[order(as.Date(finalData$onset_date, format="%Y/%m/%d")),]

source(Rutils)
#finalData[, season := lubridate::year(visit_date)]
year_start <- year(min(finalData$visit_date))
year_end <- year(max(finalData$visit_date))
years <- year_start:year_end
mid_year <- as.Date(paste(years,"07-01",sep="-"))
season <- paste(years[-length(years)],str_sub(years,start=3,end=4)[-1],sep='/')
finalData<- mutate(finalData,season=cut(visit_date,breaks=mid_year,labels=season))

# Data and config, model 
data1 <- outbreaker_data(dates = finalData$onset_date, #try with visit dates here since the Onset dates were not recorded for first
                         age_group = finalData$age_gp , #Age group
                         region = finalData$district, #Location
                         #genotype = dt_cases$Genotype, #Genotype
                         w_dens = w_dens, #Serial interval
                         f_dens = f_dens, #Latent period
                         a_dens = a_dens,  #Age stratified contact matrix
                         population = pop_vect, #Population 
                         distance = dist_mat #Distance matrix
)


config1 <- create_config(data = data1, 
                         find_import = TRUE,
                         n_iter = 20000, #Iteration number: main run
                         n_iter_import = 10000, #Iteration number: short run
                         sample_every = 50, 
                         burnin = 5000, #burnin period: first run
                         outlier_relative = T, #Absolute(F) / relative threshold 
                         outlier_threshold =0.95, #Value of the threshold
                         prior_a = c(.2, 5),
                         delta=8, #tmeporal threshold for the pre-clsutering , short 
                         verbatim =TRUE,
                         outlier_plot = T
)

# Run model 1:no inference of the importation status of cases,
out1 <- outbreaker(data = data1, config = config1, moves = moves, 
                   priors = priors, likelihoods = likelihoods)

# Set data and config for model 2
data2 <- outbreaker_data(dates = finalData$onset_date, 
                         age_group = finalData$age_gp,
                         region = finalData$district,
                         #genotype = dt_cases$Genotype, 
                         w_dens = w_dens, 
                         f_dens = f_dens, 
                         a_dens = a_dens,
                         population = pop_vect, 
                         distance = dist_mat,
                         import = finalData$recent_travel=="TRUE"  #Import status of the cases
)

config2 <- create_config(data = data2, 
                         find_import = FALSE, # No inference of import status
                         n_iter = 20000, 
                         sample_every = 50, # 1 in 50 iterations is kept
                         prior_a = c(.2, 5),
                         burnin = 5000,  
                         delta=8, #tmeporal threshold for the pre-clsutering , short 
                         verbatim =TRUE,
                         outlier_plot = T
)

# Run model 2
out2 <- outbreaker(data = data2, config = config2, moves = moves, 
                   priors = priors, likelihoods = likelihoods)




# Set data and config for model 3
data3 <- outbreaker_data(dates = finalData$onset_date, 
                         age_group = finalData$age_gp,
                         region = finalData$district,
                         #genotype = dt_cases$Genotype, 
                         w_dens = w_dens, 
                         f_dens = f_dens, 
                         a_dens = a_dens,
                         population = pop_vect, 
                         distance = dist_mat,
                         import = finalData$recent_travel=="TRUE"  #Import status of the cases
)

config1 <- create_config(data = data1, 
                         find_import = TRUE,
                         n_iter = 20000, #Iteration number: main run
                         n_iter_import = 10000, #Iteration number: short run
                         sample_every = 50, 
                         burnin = 5000, #burnin period: first run
                         outlier_relative = T, #Absolute(F) / relative threshold 
                         outlier_threshold =0.95, #Value of the threshold
                         prior_a = c(.2, 5),
                         delta=8, #tmeporal threshold for the pre-clsutering , short 
                         verbatim =TRUE,
                         outlier_plot = T
)



config3 <- create_config(data = data3, 
                         find_import = TRUE, #inference of import status
                         n_iter = 20000, #Iteration number: main run
                         n_iter_import = 10000, #Iteration number: short run
                         burnin = 5000, #burnin period: first run
                         outlier_relative = F, #Absolute / relative threshold 
                         outlier_threshold = 0.05, #Value of the threshold
                         prior_a = c(.2, 5),
                         delta=8, #tmeporal threshold for the pre-clsutering , short 
                         verbatim =TRUE,
                         #sd_b==0.01,
                         outlier_plot = T
)

# Run model 3
out3 <- outbreaker(data = data3, config = config3, moves = moves, 
                   priors = priors, likelihoods = likelihoods)


## ----params, warning = FALSE-------------------------------------------------------------------------------------
# Summary parameters a and b, removing the burnin-period
#Model 1
print(summary(out1, burnin = 5000)$a)
print(summary(out1, burnin = 5000)$b)
# Model 2
print(summary(out2, burnin = 5000)$a)
print(summary(out2, burnin = 5000)$b)
# Model 2¥3
print(summary(out3, burnin = 5000)$a)
print(summary(out3, burnin = 5000)$b)

## ----prepare_plot, warning = FALSE-------------------------------------------------------------------------------
# We create groups of cluster size: initialise the breaks for each group

##to make the smaller clusters size 
group_cluster <- c(1,2,5,10,50,100,200) - 1


# Grouped cluster size distribution in each run
clust_infer1 <- summary(out1, group_cluster = group_cluster, 
                        burnin = 5000)$cluster
clust_infer2 <- summary(out2, group_cluster = group_cluster, 
                        burnin = 5000)$cluster
clust_infer3 <- summary(out3, group_cluster = group_cluster, 
                        burnin = 5000)$cluster
# Merge inferred and reference cluster size distributions into one matrix
clust_size_matrix <- rbind(clust_infer1["Median",] , clust_infer2["Median",],clust_infer3["Median",])

clust_size_matrix_n <- rbind(clust_infer2["Median",],clust_infer3["Median",])

## ----plot_first, warning = FALSE, fig.width = 6.5, fig.height = 4, fig.cap= "Figure 3: Comparison of inferred cluster size distribution with the reference data"----
# Histogram of the inferred and reference cluster size distributions 
b<- barplot(clust_size_matrix, names.arg = colnames(clust_infer1), las=1,
            ylab = "Number of clusters", xlab = "Cluster size", main = "", 
            beside = T, ylim = c(0, max(c(clust_infer2, clust_infer3))))
# Add the 50% CI
arrows(b[1,], clust_infer1["1st Qu.",], b[1,], clust_infer1["3rd Qu.",], 
       angle = 90, code = 3, length = 0.1)
arrows(b[2,], clust_infer2["1st Qu.",], b[2,], clust_infer2["3rd Qu.",], 
       angle = 90, code = 3, length = 0.1)
arrows(b[3,], clust_infer3["1st Qu.",], b[3,], clust_infer3["3rd Qu.",], 
       angle = 90, code = 3, length = 0.1)
# Add legend
legend("topright", fill = grey.colors(3), bty = "n",
       legend = c("inferred import", "Epi import", "Epi import and inferred import"))


c<- barplot(clust_size_matrix_n, names.arg = colnames(clust_infer2), las=1,
             ylab = "Number of clusters", xlab = "Cluster size", main = "", 
             beside = T, ylim = c(0, max(c(clust_infer2, clust_infer3))))
# Add the 50% CI
arrows(b[1,], clust_infer2["1st Qu.",], b[1,], clust_infer2["3rd Qu.",], 
       angle = 90, code = 3, length = 0.1)
arrows(b[2,], clust_infer3["1st Qu.",], b[2,], clust_infer3["3rd Qu.",], 
       angle = 90, code = 3, length = 0.1)
# Add legend
legend("topright", fill = grey.colors(2), bty = "n",
       legend = c("Epi import", "Epi import and inferred import"))

## ----import_sf_file, message = FALSE, warning = FALSE, fig.width = 6.5, fig.height = 2.5-------------------------
library(ggplot2)
library(sf)
# Read the shapefile and create one map for each model

map1 <- st_read(dsn="/Users/suhan/Project/kamigoto/Kamigoto shape file", layer="kamigoto")
map1$X_CODE <- as.numeric(map1$X_CODE)
map1$Y_CODE <- as.numeric(map1$Y_CODE)
map1 <- map1 %>% group_by(S_NAME) %>% summarise()

map2 <- map1
map3 <- map1
map1$model <- "Model 1"
map2$model <- "Model 2"
map3$model <- "Model 3"

## ----n_import_reg_per_case, warning = FALSE----------------------------------------------------------------------
# Add the proportion of iterations in model 1 where each case is an import
library(data.table)
finalData<- as.data.table(finalData)
finalData[, prop_recent_travel1 := summary(out1, burnin = 5000)$tree$import]
# Add the proportion of iterations in model 2  &3where each case is an import
finalData[, prop_recent_travel2 := summary(out2, burnin = 5000)$tree$import]
finalData[, prop_recent_travel3 := summary(out3, burnin = 5000)$tree$import]
## ----n_import_reg_per_reg, warning = FALSE-----------------------------------------------------------------------
# Number of imports per region in model 1
prop_region_yr1<- finalData[, lapply(.SD, sum), by = .(district, season), .SDcols = "prop_recent_travel1"]
prop_region_yr1<- reshape(prop_region_yr1, idvar = "district", timevar = "season", direction = "wide")
prop_region_yr1<- prop_region_yr1 %>% rename(
  "2010/11" = "prop_recent_travel1.2010/11",
  "2011/12" = "prop_recent_travel1.2011/12",
  "2012/13" ="prop_recent_travel1.2012/13", 
  "2013/14" ="prop_recent_travel1.2013/14")

#prop_region_yr1[is.na(prop_region_yr1)] <- 0
names(prop_region_yr1)[1] <- "S_NAME"
prop_region_yr1<- melt(prop_region_yr1, id.vars=c("S_NAME"))

prop_region_yr2<- finalData[, lapply(.SD, sum), by = .(district, season), .SDcols = "prop_recent_travel2"]
prop_region_yr2<- reshape(prop_region_yr2, idvar = "district", timevar = "season", direction = "wide")
prop_region_yr2<- prop_region_yr2 %>% rename(
  "2010/11" = "prop_recent_travel2.2010/11",
  "2011/12" = "prop_recent_travel2.2011/12",
  "2012/13" ="prop_recent_travel2.2012/13", 
  "2013/14" ="prop_recent_travel2.2013/14")
#prop_region_yr2[is.na(prop_region_yr2)] <- 0
names(prop_region_yr2)[1] <- "S_NAME"
prop_region_yr2<- melt(prop_region_yr2, id.vars=c("S_NAME"))

prop_region_yr3<- finalData[, lapply(.SD, sum), by = .(district, season), .SDcols = "prop_recent_travel3"]
prop_region_yr3<- reshape(prop_region_yr3, idvar = "district", timevar = "season", direction = "wide")
prop_region_yr3<- prop_region_yr3 %>% rename(
  "2010/11" = "prop_recent_travel3.2010/11",
  "2011/12" = "prop_recent_travel3.2011/12",
  "2012/13" ="prop_recent_travel3.2012/13", 
  "2013/14" ="prop_recent_travel3.2013/14")
#prop_region_yr3[is.na(prop_region_yr3)] <- 0
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
maps_combined_n <- rbind(model2_map, model3_map)
library(tmap)
import<- tm_shape(maps_combined) +
  tm_polygons(
    col = "value",
    breaks = c(0,1, 5, 10, 20, 30, 40),
    title = "Number of import per region",
    pal = c("#E1F5C4", "#EDE574", "#F9D423", "#FC913A", "#FF4E50","#E16A86"),
    labels = c("0-1","1-5", "5-10", "10-20", "20-30", "30-40"),
    legend.is.portrait = TRUE) +
  tm_layout(
    frame = TRUE,
    legend.outside = TRUE,
    legend.outside.position = "bottom")+
  tm_facets(by=c("model", "variable"), showNA = FALSE)

tmap_save(import, filename = "import.png")


import_n<- tm_shape(maps_combined_n) +
  tm_polygons(
    col = "value",
    breaks = c(0,1, 5, 10, 20, 30, 40),
    title = "Number of import per region",
    pal = c("#E1F5C4", "#EDE574", "#F9D423", "#FC913A", "#FF4E50","#E16A86"),
    labels = c("0-1","1-5", "5-10", "10-20", "20-30", "30-40"),
    legend.is.portrait = TRUE) +
  tm_layout(
    frame = TRUE,
    legend.outside = TRUE,
    legend.outside.position = "bottom")+
  tm_facets(by=c("model", "variable"), showNA = FALSE)

tmap_save(import_n, filename = "import_new.png")


## ----n_sec_region, warning = FALSE-------------------------------------------------------------------------------
n_sec_per_reg <- function(finalData, out, burnin){
  ## Number of secondary cases per case
  n_sec <- apply(out[out$step > burnin, grep("alpha", colnames(out))], 1, 
                 function(X){
                   X <- factor(X, 1:length(X))
                   return(table(X))})
  #aggregatge by year
  # Vector of each season
 all_season <- unique(finalData$season)
 # Initialise output (list or dataframe). If dataframe, add one column for the season
 tot_n_sec_tot<- list()
 # For each season, compute the number of secondary cases per region
 for(i in seq_along(all_season)){
   season_i <- all_season[i]
   n_sec_i <- n_sec[finalData$season == season_i,]
   regionnames <- dt_regions$S_NAME # vector of all the regions
   
   ## Aggregate by region
   tot_n_sec_reg_i <- aggregate(n_sec_i, list(c(finalData[season == season_i,
                                                          district]#, regionnames
   )), sum)
   ## Divide by the number of cases in each region
   tot_n_sec_reg_i <- cbind(tot_n_sec_reg_i[, 1], 
                            tot_n_sec_reg_i[, -1] / table(finalData[season == season_i,
                                                                    district]))
   # check whether any region have not reported any case this year
   missing_regions <- regionnames[!is.element(regionnames, tot_n_sec_reg_i[,1])]
   # If so, then add these missing rows to the data frame tot_n_sec_reg_i
   if(length(missing_regions) > 0){
     # Matrix of missing data
     matrix_missing <- cbind.data.frame(missing_regions, 
                                        matrix(NA, nrow = length(missing_regions),
                                               ncol = ncol(n_sec_i)))
     colnames(matrix_missing) <- colnames(tot_n_sec_reg_i)
     # Merge the data frame to the missing data
     tot_n_sec_reg_i <- rbind.data.frame(tot_n_sec_reg_i, matrix_missing)
   }
   
   # Merge tot_n_sec_reg_i with tot_n_sec_tot
   tot_n_sec_tot [[i]] <- tot_n_sec_reg_i 
 }
 
 return(tot_n_sec_tot)
}

## Generate the number of secondary cases per case in each region
n_sec_tot1 <- n_sec_per_reg(finalData = finalData, out = out1, burnin = 5000)
n_sec_tot1<- data.table::rbindlist(n_sec_tot1, idcol = TRUE)
n_sec_tot1$.id[n_sec_tot1$.id==1] <- "2010/2011"
n_sec_tot1$.id[n_sec_tot1$.id==2] <- "2011/2012"
n_sec_tot1$.id[n_sec_tot1$.id==3] <- "2012/2013"
n_sec_tot1$.id[n_sec_tot1$.id==4] <- "2013/2014"


names(n_sec_tot1)[names(n_sec_tot1) == ".id"] <- "season"
names(n_sec_tot1)[names(n_sec_tot1) =="tot_n_sec_reg_i[, 1]"] <-"S_NAME"


n_sec_tot2 <- n_sec_per_reg(finalData = finalData, out = out2, burnin = 5000)
n_sec_tot2<- data.table::rbindlist(n_sec_tot2, idcol = TRUE)
n_sec_tot2$.id[n_sec_tot2$.id==1] <- "2010/2011"
n_sec_tot2$.id[n_sec_tot2$.id==2] <- "2011/2012"
n_sec_tot2$.id[n_sec_tot2$.id==3] <- "2012/2013"
n_sec_tot2$.id[n_sec_tot2$.id==4] <- "2013/2014"


names(n_sec_tot2)[names(n_sec_tot2) == ".id"] <- "season"
names(n_sec_tot2)[names(n_sec_tot2) =="tot_n_sec_reg_i[, 1]"] <- "S_NAME"

n_sec_tot3 <- n_sec_per_reg(finalData = finalData, out = out3, burnin = 5000)
n_sec_tot3<- data.table::rbindlist(n_sec_tot3, idcol = TRUE)
n_sec_tot3$.id[n_sec_tot3$.id==1] <- "2010/2011"
n_sec_tot3$.id[n_sec_tot3$.id==2] <- "2011/2012"
n_sec_tot3$.id[n_sec_tot3$.id==3] <- "2012/2013"
n_sec_tot3$.id[n_sec_tot3$.id==4] <- "2013/2014"


names(n_sec_tot3)[names(n_sec_tot3) == ".id"] <- "season"
names(n_sec_tot3)[names(n_sec_tot3) =="tot_n_sec_reg_i[, 1]"] <- "S_NAME"


## Compute the median in each model
library(tidyverse)
n_sec_tot1$n_sec <- apply(n_sec_tot1[,c(-1,-2)], 1, median)
n_sec_tot_model1<- n_sec_tot1 %>% select(season, S_NAME,n_sec)


n_sec_tot2$n_sec <- apply(n_sec_tot2[,c(-1,-2)], 1, median)
n_sec_tot_model2<- n_sec_tot2 %>% select(season, S_NAME,n_sec)

n_sec_tot3$n_sec <- apply(n_sec_tot3[,c(-1,-2)], 1, median)
n_sec_tot_model3<- n_sec_tot3 %>% select(season, S_NAME,n_sec)


## Add to the matrices describing the maps
sec_model1_map<- merge(map1, n_sec_tot_model1)
sec_model2_map<- merge(map2, n_sec_tot_model2)
sec_model3_map<- merge(map3, n_sec_tot_model3)

## ----create_maps2, warning = FALSE, fig.width = 6.5, fig.height = 2.5, fig.cap = "Figure 6: Median number of secondary transmission per case in each census tract"----
# Merge maps

maps_n_sec <- rbind(sec_model1_map, sec_model2_map, sec_model3_map)
maps_n_sec_n <- rbind(sec_model2_map, sec_model3_map)

library(tmap)
secondary<- tm_shape(maps_n_sec) +
  tm_polygons(
    col = "n_sec",
    breaks = c(0, 1, 3, 5),
    title = "Distribution of the number of secondary cases",
    pal = c("#E1F5C4", "#EDE574", "#F9D423", "#FC913A"),
    labels = c("0-1","1-3", "3-5"),
    legend.is.portrait = TRUE) +
  tm_layout(
    frame = FALSE,
    legend.outside = TRUE,
    legend.outside.position = "bottom")+
  tm_facets(by=c("model", "season"), showNA = FALSE)

tmap_save(secondary, filename = "secondary.png")

secondary_n<- tm_shape(maps_n_sec_n) +
  tm_polygons(
    col = "n_sec",
    breaks = c(0, 1, 3, 5),
    title = "Distribution of the number of secondary cases",
    pal = c("#E1F5C4", "#EDE574", "#F9D423", "#FC913A"),
    labels = c("0-1","1-3", "3-5"),
    legend.is.portrait = TRUE) +
  tm_layout(
    frame = FALSE,
    legend.outside = TRUE,
    legend.outside.position = "bottom")+
  tm_facets(by=c("model", "season"), showNA = FALSE)

tmap_save(secondary, filename = "secondary_new.png")


###################-----------------##################

n_sec_out3 <- apply(out3[out3$step >5000, grep("alpha", colnames(out3))], 1, 
                    function(X){
                      X <- factor(X, 1:length(X))
                      return(table(X))})

colnames(n_sec_out3) <- paste("X", colnames(n_sec_out3), sep = "")


n_sec_out2 <- apply(out2[out2$step >5000, grep("alpha", colnames(out2))], 1, 
                    function(X){
                      X <- factor(X, 1:length(X))
                      return(table(X))})

n_sec_out1 <- apply(out2[out1$step >5000, grep("alpha", colnames(out1))], 1, 
                    function(X){
                      X <- factor(X, 1:length(X))
                      return(table(X))})


#n_sec_out3[,1] number of cases per case at iteration 1
require(foreign)
require(ggplot2)
library(MASS)   #negative binomial regression


pop<- read.csv("data/population.csv")  #population by district by year
hh<-read.csv("data/household_size.csv")  #households by district
mean_hh<- read.csv("data/mean_hh.csv")  ## mean household per district calculated from 2010 census

#df<- merge(pop, hh,  all.x=TRUE)
df<-merge(pop,mean_hh, all.x = TRUE)
df_new<- merge(finalData, df,  all.x=TRUE)
#df_new<- cbind(df,n_sec_out3)
df_new$Trans_102[df_new$Ite_102 == 1 | df_new$Ite_102 == 2 | df_new$Ite_102 ==3 ] = "Moderate"
df_new$Trans_102[df_new$Ite_102    >3] = "High"
df_new$Trans_102[df_new$Ite_102  == 0] = "No"

##new age_group (<5, 5-14, 15-64,>=65)
df_new[df_new$age <= 4, "n_age_group"] <- "<5"
df_new[df_new$age > 4 & df_new$age <= 14, "n_age_group"] <- "5-14"
df_new[df_new$age > 14 & df_new$age <= 64, "n_age_group"] <- "15-64"
df_new[df_new$age > 64, "n_age_group"] <- ">64"
df_new$n_age_group<-factor(df_new$n_age_group,level = c("15-64","<5","5-14", ">64"))  
## set the 15-64 years as factor level 1 so that i  can use this ag group as ref, i have to check whether it makes sensese

  

##vaccination status
df_new$vac [df_new$n_dose_vaccine_reported  >= 1] = 1       
df_new$vac [df_new$n_dose_vaccine_reported  == 0] = 0 
sum(is.na(df_new$vac))   ##71 missing values

df_new$vac<- factor(df_new$vac, levels=c(0,1), labels=c("No", "Yes"))
summary(df_new$vac)


#normalization data
#removing outliers


cor(df_new) 


# run n regressions
library(MASS)
age_gp_lms <- lapply(1:ncol(n_sec_out3), function(x) glm.nb(n_sec_out3[,x] ~ df_new$n_age_group))
vac_lms <- lapply(1:ncol(n_sec_out3), function(x) glm.nb(n_sec_out3[,x] ~ df_new$vac,na.action=na.omit))


lms <- lapply(1:ncol(n_sec_out3), function(x) glm.nb(n_sec_out3[,x] ~ df_new$n_age_group+ 
                                                              df_new$vac+ 
                                                       df_new$district_hh+
                                                       df_new$mean_hh+
                                                       df_new$season+
                                                       df_new$district_pop,
                                                     na.action=na.omit))
# extract just coefficients
coef_lms<- sapply(lms, coef)
coef_lms


# if you need more info, get full summary call. now you can get whatever, like:
summaries <- lapply(lms, summary)
# ...coefficents with p values:
p_values<- lapply(summaries, function(x) x$coefficients[, c(1,4)])
# ...or r-squared values
sapply(summaries, function(x) c(r_sq = x$r.squared, 
                                adj_r_sq = x$adj.r.squared))




########################################
##Region Index
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="若松郷" ] = 1    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="桐古里郷"] = 2    
n_sec_tot1$region_index [n_sec_tot1$S_NAME =="宿ノ浦郷" ] = 3    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="荒川郷"] = 4    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="西神ノ浦郷"] = 5    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="日島郷"] = 6    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="有福郷" ] = 7    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="漁生浦郷" ] = 8    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="間伏郷" ] = 9    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  == "榊ノ浦郷"  ] = 10    
n_sec_tot1$region_index [n_sec_tot1$S_NAME =="相河郷" ] = 11    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="青方郷" ] = 12    
n_sec_tot1$region_index [n_sec_tot1$S_NAME =="網上郷" ] = 13    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="飯ノ瀬戸郷"] = 14    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="今里郷" ] = 15    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="続浜ノ浦郷" ] = 16    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="奈摩郷" ] = 17    
n_sec_tot1$region_index [n_sec_tot1$S_NAME =="船崎郷"] = 18    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="三日ノ浦郷" ] = 19    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="道土井郷"] = 20    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="浦桑郷" ] = 21    
n_sec_tot1$region_index [n_sec_tot1$S_NAME =="榎津郷"] = 22    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="丸尾郷" ] = 23    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="似首郷"] = 24    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  == "小串郷" ] = 25    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="曽根郷" ] = 26    
n_sec_tot1$region_index [n_sec_tot1$S_NAME == "立串郷"] = 27    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="津和崎郷"] = 28    
n_sec_tot1$region_index [n_sec_tot1$S_NAME =="有川郷"  ] = 29   
n_sec_tot1$region_index [n_sec_tot1$S_NAME =="小河原郷"  ] = 30   
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="赤尾郷"] = 31   
n_sec_tot1$region_index [n_sec_tot1$S_NAME =="友住郷"] = 32    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="江ノ浜郷" ] = 33    
n_sec_tot1$region_index [n_sec_tot1$S_NAME =="太田郷"] = 34    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="七目郷" ] = 35    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="鯛ノ浦阿瀬津郷" ] = 36    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="東神ノ浦郷"  ] = 37    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="奈良尾郷"  ] = 38    
n_sec_tot1$region_index [n_sec_tot1$S_NAME  =="岩瀬浦郷"  ] = 39  


n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="若松郷" ] = 1    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="桐古里郷"] = 2    
n_sec_tot2$region_index [n_sec_tot2$S_NAME =="宿ノ浦郷" ] = 3    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="荒川郷"] = 4    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="西神ノ浦郷"] = 5    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="日島郷"] = 6    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="有福郷" ] = 7    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="漁生浦郷" ] = 8    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="間伏郷" ] = 9    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  == "榊ノ浦郷"  ] = 10    
n_sec_tot2$region_index [n_sec_tot2$S_NAME =="相河郷" ] = 11    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="青方郷" ] = 12    
n_sec_tot2$region_index [n_sec_tot2$S_NAME =="網上郷" ] = 13    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="飯ノ瀬戸郷"] = 14    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="今里郷" ] = 15    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="続浜ノ浦郷" ] = 16    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="奈摩郷" ] = 17    
n_sec_tot2$region_index [n_sec_tot2$S_NAME =="船崎郷"] = 18    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="三日ノ浦郷" ] = 19    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="道土井郷"] = 20    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="浦桑郷" ] = 21    
n_sec_tot2$region_index [n_sec_tot2$S_NAME =="榎津郷"] = 22    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="丸尾郷" ] = 23    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="似首郷"] = 24    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  == "小串郷" ] = 25    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="曽根郷" ] = 26    
n_sec_tot2$region_index [n_sec_tot2$S_NAME == "立串郷"] = 27    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="津和崎郷"] = 28    
n_sec_tot2$region_index [n_sec_tot2$S_NAME =="有川郷"  ] = 29   
n_sec_tot2$region_index [n_sec_tot2$S_NAME =="小河原郷"  ] = 30   
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="赤尾郷"] = 31   
n_sec_tot2$region_index [n_sec_tot2$S_NAME =="友住郷"] = 32    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="江ノ浜郷" ] = 33    
n_sec_tot2$region_index [n_sec_tot2$S_NAME =="太田郷"] = 34    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="七目郷" ] = 35    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="鯛ノ浦阿瀬津郷" ] = 36    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="東神ノ浦郷"  ] = 37    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="奈良尾郷"  ] = 38    
n_sec_tot2$region_index [n_sec_tot2$S_NAME  =="岩瀬浦郷"  ] = 39    


n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="若松郷" ] = 1    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="桐古里郷"] = 2    
n_sec_tot3$region_index [n_sec_tot3$S_NAME =="宿ノ浦郷" ] = 3    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="荒川郷"] = 4    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="西神ノ浦郷"] = 5    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="日島郷"] = 6    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="有福郷" ] = 7    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="漁生浦郷" ] = 8    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="間伏郷" ] = 9    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  == "榊ノ浦郷"  ] = 10    
n_sec_tot3$region_index [n_sec_tot3$S_NAME =="相河郷" ] = 11    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="青方郷" ] = 12    
n_sec_tot3$region_index [n_sec_tot3$S_NAME =="網上郷" ] = 13    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="飯ノ瀬戸郷"] = 14    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="今里郷" ] = 15    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="続浜ノ浦郷" ] = 16    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="奈摩郷" ] = 17    
n_sec_tot3$region_index [n_sec_tot3$S_NAME =="船崎郷"] = 18    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="三日ノ浦郷" ] = 19    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="道土井郷"] = 20    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="浦桑郷" ] = 21    
n_sec_tot3$region_index [n_sec_tot3$S_NAME =="榎津郷"] = 22    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="丸尾郷" ] = 23    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="似首郷"] = 24    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  == "小串郷" ] = 25    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="曽根郷" ] = 26    
n_sec_tot3$region_index [n_sec_tot3$S_NAME == "立串郷"] = 27    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="津和崎郷"] = 28    
n_sec_tot3$region_index [n_sec_tot3$S_NAME =="有川郷"  ] = 29   
n_sec_tot3$region_index [n_sec_tot3$S_NAME =="小河原郷"  ] = 30   
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="赤尾郷"] = 31   
n_sec_tot3$region_index [n_sec_tot3$S_NAME =="友住郷"] = 32    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="江ノ浜郷" ] = 33    
n_sec_tot3$region_index [n_sec_tot3$S_NAME =="太田郷"] = 34    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="七目郷" ] = 35    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="鯛ノ浦阿瀬津郷" ] = 36    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="東神ノ浦郷"  ] = 37    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="奈良尾郷"  ] = 38    
n_sec_tot3$region_index [n_sec_tot3$S_NAME  =="岩瀬浦郷"  ] = 39    

#library(bayestestR)
ci_hdi <- ci(n_sec_tot1$n_sec, method = "HDI")
ci_eti <- ci(posterior, method = "ETI")



    
#plot with credible interval
n_sec_tot1$low_sec <- apply(n_sec_tot1[,c(-1,-2)], 1, function(X) 
  return(quantile(x = X, probs = 0.025, na.rm =TRUE)))
  
n_sec_tot1$up_sec <- apply(n_sec_tot1[,c(-1,-2)], 1, function(X) 
  return(quantile(x = X, probs = 0.975, na.rm =TRUE)))


#plot with credible interval
n_sec_tot1$low_sec <- apply(n_sec_tot1[,c(-1,-2)], 1, function(X) 
  return(quantile(x = X, probs = 0.025, na.rm =TRUE)))

n_sec_tot1$up_sec <- apply(n_sec_tot1[,c(-1,-2)], 1, function(X) 
  return(quantile(x = X, probs = 0.975, na.rm =TRUE)))

ggplot(n_sec_tot1) +
  geom_bar( aes(x=S_NAME, y=n_sec), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=S_NAME, ymin=low_sec, ymax=up_sec), width=0.4, colour="orange", alpha=0.9, size=1.3)+
  facet_wrap(~season,  ncol=1)


#plot with credible interval
n_sec_tot1$low_sec <- apply(n_sec_tot1[,c(-1,-2)], 1, function(X) 
  return(quantile(x = X, probs = 0.025, na.rm =TRUE)))

n_sec_tot1$up_sec <- apply(n_sec_tot1[,c(-1,-2)], 1, function(X) 
  return(quantile(x = X, probs = 0.975, na.rm =TRUE)))

ggplot(n_sec_tot1, aes(x=region_index, y=n_sec, group=season, color=season)) +
  geom_line() +
  geom_errorbar(aes(ymin=low_sec, ymax=up_sec))

#plot with credible interval
n_sec_tot2$low_sec <- apply(n_sec_tot2[,c(-1,-2)], 1, function(X) 
  return(quantile(x = X, probs = 0.025, na.rm =TRUE)))

n_sec_tot2$up_sec <- apply(n_sec_tot2[,c(-1,-2)], 1, function(X) 
  return(quantile(x = X, probs = 0.975, na.rm =TRUE)))


ggplot(n_sec_tot2, aes(x=region_index, y=n_sec, group=season, color=season)) +
  geom_line() +
  geom_errorbar(aes(ymin=low_sec, ymax=up_sec))



#plot with credible interval
new_sec3<- na.omit(n_sec_tot3) 
new_sec3$low_sec <- apply(new_sec3[,c(-1,-2)], 1, function(X) 
  return(quantile(x = X, probs = 0.025, na.rm =TRUE)))

new_sec3$up_sec <- apply(new_sec3[,c(-1,-2)], 1, function(X) 
  return(quantile(x = X, probs = 0.975, na.rm =TRUE)))

##this plot seems ok
ggplot(new_sec3, aes(x=region_index, y=n_sec, group=season, color=season)) +
  geom_point() +
  geom_errorbar(aes(ymin=low_sec, ymax=up_sec))


ggplot(new_sec3, aes(x=region_index, y=n_sec)) +
  geom_point() +
  geom_errorbar(aes(ymin=low_sec, ymax=up_sec))+
  facet_wrap(~season,  ncol=1)

  
############################3
maps_n_sec <- maps_n_sec[maps_n_sec$INTPTLON > lim_lon[1] &
                           maps_n_sec$INTPTLON < lim_lon[2] &
                           maps_n_sec$INTPTLAT > lim_lat[1] & 
                           maps_n_sec$INTPTLAT < lim_lat[2],]

# Plot the geographical distribution of the number of secondary cases
ggplot(maps_n_sec) +  geom_sf(aes(fill = n_sec)) + facet_grid(~model)  +     
  scale_fill_gradient2(na.value = "lightgrey", mid = "lightblue",
                       low = "white", midpoint = 1, high = "darkblue",
                       breaks = seq(0, 5, 0.5),name = "Sec cases") +
  #coord_sf(xlim = c(-83.8, -82.2), ylim = c(40.2, 41.3)) +
  theme_classic(base_size = 9) 


