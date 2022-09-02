library(readxl)
library(tidyverse)
library(data.table)
library(o2geosocial)
library(tmap)
library(sf)
library(foreign)
library(reshape2)
library(ggplot2)
library(showtext)
library(surveillance)
library(ggpubr)
library(incidence)
library(MASS)

dt_regions <- readRDS("data/dt_regions.rds")
finalData  <- readRDS("data/220626_finaldata.rds")
finalData<- as.data.table(finalData)
## population data of 8 years
dis_pop<- readRDS("data/dis_pop.rds")



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


###weekly total contact (5*weekdays, 2*weekends)
weekday_social_contact<- as.matrix(read_excel("data/contact_weekday.xlsx", 
                                              col_types = c("numeric", "numeric", "numeric", 
                                                            "numeric", "numeric", "numeric", 
                                                            "numeric", "numeric", "numeric", 
                                                            "numeric", "numeric", "numeric", 
                                                            "numeric", "numeric", "numeric")))

weekend_social_contact<-as.matrix(read_excel("data/contact_weekend.xlsx", 
                                             col_types = c("numeric", "numeric", "numeric", 
                                                           "numeric", "numeric", "numeric", 
                                                           "numeric", "numeric", "numeric", 
                                                           "numeric", "numeric", "numeric", 
                                                           "numeric", "numeric", "numeric")))


social_contact<- ((5*weekday_social_contact)+(2*weekend_social_contact))/7

social_contact <-t(social_contact)
social_contact <-data.table::as.data.table(social_contact)
# Compute the proportion of connection to each age group
a_dens <- t(t(social_contact)/colSums(social_contact))

# Extract the population vector
pop_vect <- dis_pop$ave_pop

## ----distance_pop, warning = FALSE-------------------------------------------------------------------------------
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


###Stouffer rank 
##dont sue the absolute distance s in previous method, 
#number of inhabitatns in between the region (the distance doest matter)
#only one parameter
# For every column of the distance matrix, use the cumulative sum of the 
# population vector ordered by the distance. Remove the values where 
# the distance between the regions is above gamma

# Rename the matrix columns and rows, and the population vector
names(pop_vect) <- rownames(dist_mat) <- colnames(dist_mat) <- 
  dt_regions$S_NAME

###Qestion on config
dist_mat_stouffer <- apply(dist_mat, 2, function(X){
  pop_X <- cumsum(pop_vect[order(X)])
  omega_X <- pop_X[names(X)]
  # omega_X is set to -1 if the distance between two regions is above gamma
  omega_X[X > config1$gamma] <- -1
  return(omega_X)
})
# The new value of gamma is equal to the maximum of dist_mat_stouffer + 1
gamma <- max(dist_mat_stouffer) + 1

# The values previously set to -1 are now set to the new value of gamma
dist_mat_stouffer[dist_mat_stouffer == -1] <- max(dist_mat_stouffer) * 2

# Set movement, likelihood and prior lists to default
source("code/RCPP_stouffer.R")
moves3 <- custom_moves(a =cpp_stouffer)

# Set movement, likelihood and prior lists to default
moves <- custom_moves()
likelihoods <- custom_likelihoods()
priors <- custom_priors()

f_null <- function(param) {
  return(0.0)
}

priors3 <- custom_priors(b = f_null)

data_stouffer <- outbreaker_data(dates = finalData$onset_date, 
                                 age_group = finalData$age_gp,
                                 region = finalData$district,
                                 genotype = finalData$rdt, 
                                 w_dens = w_dens, 
                                 f_dens = f_dens, 
                                 a_dens = a_dens,
                                 population = pop_vect, 
                                 distance = dist_mat_stouffer, #Distance matrix
                                 import = finalData$recent_travel=="TRUE" 
)

config_stouffer <- create_config(data = data_stouffer, 
                                 find_import = TRUE, #inference of import status
                                 gamma = gamma,
                                 prior_a = c(0, 3),
                                 move_b = FALSE, # b is not estimated
                                 init_b = 0,
                                 spatial_method = "power-law",
                                 n_iter = 40000, #Iteration number: main run
                                 n_iter_import = 20000, #Iteration number: short run
                                 burnin = 5000, #burnin period: first run
                                 outlier_relative = F, #Absolute / relative threshold
                                 outlier_threshold = 0.05, #Value of the threshold
                                 verbatim =TRUE,
                                 #sd_b==0.01,
                                 outlier_plot = T
)

set.seed(1)
# Run the model using the Stouffer's rank method
out3_stouffer <- outbreaker(data = data_stouffer, config = config_stouffer , moves = moves3, 
                            priors = priors3, likelihoods = likelihoods)


# Histogram of cluster size distributions 
##to make the smaller clusters size 
group_cluster <- c(1,2,5,10,50,100,200,500) - 1
clust_stouffer <- summary(out3_stouffer, group_cluster = group_cluster, 
                          burnin = 5000)$cluster
#clust_stouffer <- as.data.frame(clust_stouffer )
#clust_stouffer <- clust_stouffer %>% mutate(across(where(is.numeric), ~ round(., 0)))
clus<- barplot(clust_stouffer["Median",], las=1,ylim=c(0,600),
               ylab = "Number of clusters", xlab = "Cluster size", main = "", 
               beside = T)
# Add the 50% CI
arrows(clus, clust_stouffer["1st Qu.",], clus,clust_stouffer["3rd Qu.",], angle = 90, code = 3, length = 0.1)


##heatmap
summ_out3_stouffer <- summary(out3_stouffer, burnin = 1000)
median_a_st <- summ_out3_stouffer$a["Median"]
#median_b_st <- summ_out3_stouffer$b["Median"]   ##this doesnt matter for stouffer rank coz of only one parameter

connectivity_matrix_st <- o2geosocial:::cpp_log_like_s(population = pop_vect, distance = dist_mat_stouffer, 
                                                       a = median_a_st, b = - median_a_st,
                                                       spatial = "power-law")[[1]]

colnames(connectivity_matrix_st) <- rownames(connectivity_matrix_st) <- dt_regions$S_NAME


data_melt_st <-reshape2::melt(connectivity_matrix_st) ## convert into a long format
data_melt_st$value_char <- cut(data_melt_st$value, breaks = c(0,0.01, 0.1, 0.5, 1), 
                               labels = c("<0.01", "0.01-0.1", "0.1-0.5", ">0.5")) ## group in categories

showtext_auto()
gg_heatmap_stouffer <- ggplot(data_melt_st, aes(Var1, Var2,fill = value_char)) + 
  geom_tile(colour="white", size=0.25) + 
  scale_fill_manual(values = c("lightgrey", "orange", "red", "black")) +
  theme(axis.title = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),
        axis.ticks = element_blank())+ guides(fill=guide_legend(title="Probability of transmission"))
gg_heatmap_stouffer

#import status
finalData<- as.data.table(finalData)
finalData[, prop_recent_travel1 := summary(out3_stouffer, burnin = 5000)$tree$import]

# Number of imports per region
prop_region_stouf<- finalData[, lapply(.SD, sum), by = .(district, season), .SDcols = "prop_recent_travel1"]
prop_region_stouf<- reshape(prop_region_stouf, idvar = "district", timevar = "season", direction = "wide")
prop_region_stouf<- prop_region_stouf %>% rename(
  "2010/11" = "prop_recent_travel1.2010/11",
  "2011/12" = "prop_recent_travel1.2011/12",
  "2012/13" ="prop_recent_travel1.2012/13", 
  "2013/14" ="prop_recent_travel1.2013/14",
  "2014/15" = "prop_recent_travel1.2014/15",
  "2015/16" = "prop_recent_travel1.2015/16",
  "2016/17" = "prop_recent_travel1.2016/17",
  "2017/18" ="prop_recent_travel1.2017/18")
prop_region_stouf[is.na(prop_region_stouf)] <- 0
names(prop_region_stouf)[1] <- "S_NAME"
prop_region_stouf<- reshape2::melt(prop_region_stouf, id.vars=c("S_NAME"))


# Read the shapefile and create one map for each model
#map <- st_read(dsn="data/Kamigoto shape file", layer="kamigoto")
map1 <- st_read(dsn="data/Kamigoto shape file", layer="kamigoto")
map1$X_CODE <- as.numeric(map1$X_CODE)
map1$Y_CODE <- as.numeric(map1$Y_CODE)
map1 <- map1 %>% group_by(S_NAME) %>% summarise()

stouf_map<- merge(map1, prop_region_stouf)
max_import <- round(max(prop_region_stouf$value ,na.rm = TRUE), 0)
if(max_import  < max(prop_region_stouf$value,na.rm = TRUE)) 
  max_import  <- max_import  + 1

import_stouf<- tm_shape(stouf_map) +
  tm_polygons(
    col = "value",
    breaks = c(0,1, 5, 10, 20, 30, max_import),
    title = "Number of import per region",
    pal = c("#E1F5C4", "#EDE574", "#F9D423", "#FC913A", "#FF4E50"),     ###"#E16A86"
    labels = c("0-1","1-5", "5-10", "10-20", "20-30",paste0("30-",max_import )),
    legend.is.portrait = TRUE) +
  tm_layout(
    frame = TRUE,
    legend.outside = TRUE, 
    #legend.text.size = 10, legend.title.size = 8, panel.label.size = 1.5,
    legend.outside.position = "bottom")+
  tm_facets(by="variable", showNA = FALSE)

tmap_save(import_stouf, filename = "Figures/import_stouf.png")
import_stouf


##nubmer of secdondary cases
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
n_sec_tot1_stouffer <-n_sec_per_reg(finalData = finalData, out = out3_stouffer, burnin = 5000)
n_sec_tot1_stouffer<- data.table::rbindlist(n_sec_tot1_stouffer, idcol = TRUE)
n_sec_tot1_stouffer$.id[n_sec_tot1_stouffer$.id==1] <- "2010/2011"
n_sec_tot1_stouffer$.id[n_sec_tot1_stouffer$.id==2] <- "2011/2012"
n_sec_tot1_stouffer$.id[n_sec_tot1_stouffer$.id==3] <- "2012/2013"
n_sec_tot1_stouffer$.id[n_sec_tot1_stouffer$.id==4] <- "2013/2014"
n_sec_tot1_stouffer$.id[n_sec_tot1_stouffer$.id==5] <- "2014/2015"
n_sec_tot1_stouffer$.id[n_sec_tot1_stouffer$.id==6] <- "2015/2016"
n_sec_tot1_stouffer$.id[n_sec_tot1_stouffer$.id==7] <- "2016/2017"
n_sec_tot1_stouffer$.id[n_sec_tot1_stouffer$.id==8] <- "2017/2018"

names(n_sec_tot1_stouffer)[names(n_sec_tot1_stouffer) == ".id"] <- "season"
names(n_sec_tot1_stouffer)[names(n_sec_tot1_stouffer) =="tot_n_sec_reg_i[, 1]"] <-"S_NAME"

## Compute the median
n_sec_tot1_stouffer$n_sec <- apply(n_sec_tot1_stouffer[,c(-1,-2)], 1, median)
# n_sec_tot1_stouffer$n_sec <- apply(n_sec_tot1_stouffer[,c(-1,-2)], 1, function(X){
#   return(quantile(X, probs = 0.025, na.rm = T))
# })
n_sec_tot1_stouffer_model1<- n_sec_tot1_stouffer %>% select(season, S_NAME,n_sec)


## Add to the matrices describing the maps
n_sec_tot1_stouffer_model1_map<- merge(map1, n_sec_tot1_stouffer_model1)

##chekcing the maximum values
max_value <- round(max(n_sec_tot1_stouffer$n_sec,na.rm = TRUE), 1)
if(max_value < max(n_sec_tot1_stouffer$n_sec,na.rm = TRUE)) 
  max_value <- max_value + 0.1


##making map with median value
secondary_stouffer<- tm_shape(n_sec_tot1_stouffer_model1_map) +
  tm_polygons(
    col = "n_sec",
    breaks = c(0, 0.8, 1.2, max_value),
    title = "Distribution of the number of secondary cases",
    pal = c("#E1F5C4", "#EDE574", "#F9D423", "#FC913A"),
    labels = c("0-0.8","0.8-1.2", paste0("1.2-", max_value)),
    legend.is.portrait = TRUE) +
  tm_layout(
    frame = FALSE,
    legend.outside = TRUE,
    #legend.text.size = 10, legend.title.size = 8, panel.label.size = 1.5,
    legend.outside.position = "bottom")+
  tm_facets(by="season", showNA = FALSE)

tmap_save(secondary_stouffer, filename = "Figures/secondary_stouffer.png")
secondary_stouffer

##displaying the Reff with 95%CI
n_sec <- function(finalData, out, burnin){
  ## Number of secondary cases per case
  n_sec <- apply(out[out$step > burnin, grep("alpha", colnames(out))], 1, 
                 function(X){
                   X <- factor(X, 1:length(X))
                   return(table(X))})
}
sec_distri<- n_sec(finalData = finalData, out = out3_stouffer, burnin = 5000)

median_sec_distri <-as.data.frame(apply(sec_distri, 1, median))
low_sec_distri <- as.data.frame(apply(sec_distri, 1, function(X){
  return(quantile(X, probs = 0.025, na.rm = T))
}))
high_sec_distri <- as.data.frame(apply(sec_distri, 1, function(X){
  return(quantile(X, probs = 0.975, na.rm = T))
}))

sec_case<- cbind(median_sec_distri, low_sec_distri, high_sec_distri)
names(sec_case)[1] <- "sec_case_med"
names(sec_case)[2] <- "sec_case_lower"
names(sec_case)[3] <- "sec_case_upper"

dist<- n_sec_tot1_stouffer %>% select(season,S_NAME,n_sec, low_n_sec,high_n_sec)
##n_secondary cases with 9% CI  all district by season
year_fig<- ggplot(dist) +
  geom_bar( aes(x=S_NAME, y=n_sec), stat="identity", fill="skyblue", alpha=0.4) +
  facet_wrap (~season, ncol=2) + theme_bw() 
year_fig + geom_errorbar( aes(x=S_NAME, ymin=low_n_sec, ymax=high_n_sec), width=0.2, colour="orange", alpha=0.4, size=1.0)+
  facet_wrap (~season, ncol=2) + theme_bw() +
  theme(text = element_text(size=8), 
        axis.text.x = element_text(angle=90, hjust=1), 
        axis.title.x=element_blank())

##getting the quantile of the no of secondary cases in quantiles

n_sec_out3_stouf <- apply(out3_stouffer[out3_stouffer$step >5000, grep("alpha", colnames(out3_stouffer))], 1, 
                          function(X){
                            X <- factor(X, 1:length(X))
                            return(table(X))})


test<- apply(n_sec_out3_stouf , 1, function(X) quantile(X,c(0.025,0.25,0.5,0.75,0.975),))
dat.test<- t(test)
dat.test<- as.data.frame(dat.test)
#creating the unique ID
dat.test$id <- seq_len(nrow(dat.test))
dat.test <- dat.test[order(dat.test$id),]

finalData$id <- seq_len(nrow(finalData))
finalData <- finalData[order(finalData$id),]

#neighbourhood vac coverage
map <- map1
if (class(map)[1] != "SpatialPolygonsDataFrame"){
  map <- as(map,"Spatial")    
}

neighbourhood <- nbOrder(poly2adjmat(map), maxlag = 7)
colnames(neighbourhood) <- rownames(neighbourhood) <- map$S_NAME
neighbourhood[1:5, 1:5]
pop$cov_nei1 <-  pop$cov_nei2 <- pop$cov_nei3 <- -1
# For each region, compute the number of vacc / total population in neighbours of degree 1
# For each region, compute the number of vacc / total population in neighbours of degree 2
# For each region, compute the number of vacc / total population in neighbours of degree 3
for(i in seq_len(nrow(neighbourhood))){
  i_reg <- rownames(neighbourhood)[i]
  neighbours_i <- colnames(neighbourhood)[which(neighbourhood[i,] == 1)]
  neighbours_i2 <- colnames(neighbourhood)[which(neighbourhood[i,] == 2)]
  neighbours_i3 <- colnames(neighbourhood)[which(neighbourhood[i,] == 3)]
  for(j in unique(pop$season)){
    neigh_pop <- sum(pop[pop$season == j & is.element(pop$district, neighbours_i),][, "district_pop"])
    neigh_vac <- sum(pop[pop$season == j & is.element(pop$district, neighbours_i),][, "vac"])
    pop[pop$season == j & pop$district == i_reg, ]$cov_nei1 <- neigh_vac/neigh_pop
    neigh_pop2 <- sum(pop[pop$season == j & is.element(pop$district, neighbours_i2),][, "district_pop"])
    neigh_vac2 <- sum(pop[pop$season == j & is.element(pop$district, neighbours_i2),][, "vac"])
    pop[pop$season == j & pop$district == i_reg, ]$cov_nei2 <- neigh_vac2/neigh_pop2
    neigh_pop3 <- sum(pop[pop$season == j & is.element(pop$district, neighbours_i3),][, "district_pop"])
    neigh_vac3 <- sum(pop[pop$season == j & is.element(pop$district, neighbours_i3),][, "vac"])
    pop[pop$season == j & pop$district == i_reg, ]$cov_nei3 <- neigh_vac3/neigh_pop3
  }
}


df1 <- merge (finalData, pop, by = c( "district", "season"))
df1  <- df1 [order(df1 $id),]
mer_df<-merge(df1,dat.test, by="id")

##new age_group (wr)
mer_df[mer_df$age <= 3, "n_age_group"] <- "<3"
mer_df[mer_df$age > 3 & mer_df$age <= 6, "n_age_group"] <- "4-6"
mer_df[mer_df$age > 6 & mer_df$age <= 12, "n_age_group"] <- "7-12"
mer_df[mer_df$age > 12 & mer_df$age <= 18, "n_age_group"] <- "13-18"
mer_df[mer_df$age > 18 & mer_df$age <= 64, "n_age_group"] <- "19-64"
mer_df[mer_df$age > 64 & mer_df$age <= 74, "n_age_group"] <- "65-74"
mer_df[mer_df$age > 74, "n_age_group"] <- ">74"
mer_df$n_age_group<-factor(mer_df$n_age_group,level = c("19-64","<3","4-6","7-12","13-18","65-74", ">74"))  
## set the 15-64 years as factor level 1 so that i  can use this ag group as ref, i have to check whether it makes sensese

##vaccination status
mer_df<- mer_df %>%
  mutate(vac_his = case_when(n_dose_vaccine_reported  >= 1 ~ 1,
                             n_dose_vaccine_reported  == 0 ~ 0,
                             is.na(n_dose_vaccine_reported) ~ -1))
sum(is.na(mer_df$vac_his))   ##0 missing values
mer_df$vac_his <- factor(mer_df$vac_his, levels=c( 0, 1, -1), labels=c("No","Yes", "Unknown"))
summary(mer_df$vac_his)

# No     Yes    Unknown 
#1428    5296    1639 

#regional vaccination
mer_df$vac_cov<- mer_df$vac/mer_df$district_pop *100
mer_df$cov_cat[mer_df$vac_cov <50] = "<50"
mer_df$cov_cat[mer_df$vac_cov >=50 & mer_df$vac_cov <60] = "50-60"
mer_df$cov_cat[mer_df$vac_cov >=60 & mer_df$vac_cov <65 ] = "60-65"
mer_df$cov_cat[mer_df$vac_cov >=65] = ">65"
mer_df$cov_cat<-factor(mer_df$cov_cat,level = c("50-60", "<50","60-65", ">65"))  
summary(mer_df$cov_cat)

##categorizing the district population
mer_df$pop_cat[mer_df$district_pop <500] = "<500"
mer_df$pop_cat[mer_df$district_pop >=500 & mer_df$district_pop <2000] = "500-2000"
mer_df$pop_cat[mer_df$district_pop >=2000] = ">2000"
mer_df$pop_cat<-factor(mer_df$pop_cat,level = c("<500","500-2000", ">2000")) 
# mer_df$pop_cat<-factor(mer_df$pop_cat,level = c("500-2000", "<500",">2000")) 
# mer_df$pop_cat<-factor(mer_df$pop_cat,level = c("<500","500-2000", ">2000")) 
summary(mer_df$pop_cat)

#for seasonality
t <- as.numeric(mer_df$onset_date - as.Date("2010-01-01"))
mer_df$cost <- cos(2*pi*t/365)
mer_df$sint <- sin(2*pi*t/365)

#regression analysis
#univariate
age_gp_lms <- lapply(1:ncol(n_sec_out3_stouf), function(x) glm.nb(n_sec_out3_stouf[,x] ~ mer_df$n_age_group ))
summary(glm.nb(n_sec_out3_stouf[,500] ~ mer_df$n_age_group,na.action=na.omit))

vac_lms <- lapply(1:ncol(n_sec_out3_stouf), function(x) glm.nb(n_sec_out3_stouf[,x] ~ mer_df$vac_his,na.action=na.omit))
summary(glm.nb(n_sec_out3_stouf[,500] ~ mer_df$vac_his,na.action=na.omit))

hh_lms <- lapply(1:ncol(n_sec_out3_stouf), function(x) glm.nb(n_sec_out3_stouf[,x] ~ log(mer_df$mean_hh),na.action=na.omit))
summary(glm.nb(n_sec_out3_stouf[,500] ~ log(mer_df$mean_hh),na.action=na.omit))


pop_den_lms<- lapply(1:ncol(n_sec_out3_stouf), function(x) glm.nb(n_sec_out3_stouf[,x] ~ mer_df$pop_cat,na.action=na.omit))
summary(glm.nb(n_sec_out3_stouf[,500] ~ mer_df$pop_cat,na.action=na.omit))

seasonality<- lapply(1:ncol(n_sec_out3_stouf), function(x) glm.nb(n_sec_out3_stouf[,x] ~ mer_df$season,na.action=na.omit))
summary(glm.nb(n_sec_out3_stouf[,500] ~ mer_df$season,na.action=na.omit))

reg_cov<- lapply(1:ncol(n_sec_out3_stouf), function(x) glm.nb(n_sec_out3_stouf[,x] ~ mer_df$cov_cat,na.action=na.omit))
summary(glm.nb(n_sec_out3_stouf[,500] ~ mer_df$cov_cat, na.action=na.omit))

neig_cov<- lapply(1:ncol(n_sec_out3_stouf), function(x) glm.nb(n_sec_out3_stouf[,x] ~ mer_df$cov_nei1,na.action=na.omit))
summary(glm.nb(n_sec_out3_stouf[,500] ~ mer_df$cov_nei1, na.action=na.omit))

#regression
m1<-lapply(1:ncol(n_sec_out3_stouf), function(x) glm.nb(n_sec_out3_stouf[,x] ~ mer_df$n_age_group ++ #age group
                                                          mer_df$vac_his  + #individual vaccine history
                                                          log(mer_df$mean_hh)+  #mean hh
                                                          log(mer_df$district_pop)  +  #pop per district
                                                          mer_df$season +    #season
                                                          log(1 - mer_df$cov_nei1) +  #neighbourhood vaccination coverage
                                                          log(1 - mer_df$vac_cov/100)+     #regional vaccination coverage
                                                          mer_df$rdt+    #type of influenza
                                                          mer_df$cost+   #seasonality
                                                          mer_df$sint,   #seasonality
                                                        na.action=na.omit))
#taking the coefficients
coef_m1<-exp(t(apply(sapply(m1, coef), 1, function(X) quantile(X,c(0.025,0.25,0.5,0.75,0.975),))))
library(mitml)
##https://cran.r-project.org/web/packages/mitml/mitml.pdf code from this link
est <- testEstimates(m1)
##taking the confidence interval
conf<- exp(confint(est))



#calculating for IRR and CI
results_exp<- exp(results)
results_exp_t<- t(results_exp)
#get the confidence intervals for the coefficients by profiling the likelihood function.
#est <- cbind(Estimate = coef(m1), confint(m1))
#getting IRR and their CIs
#exp(est)

#getting p-values
p_values<- sapply(m1, function(x) summary(x)$coefficients[,4])
p_values_med<-as.data.frame(apply(p_values,1,median)) 


#plotting the proportion
p <- ggplot(data = dfr_perc, 
            aes(x = n_age_group, y = perc, fill = sec_trans)) + 
  geom_bar(stat = 'identity', position = 'dodge', alpha = 2/3) +  
  #scale_y_continuous(labels = perc) +    
  #scale_fill_few("medium", drop = FALSE) +  
  labs(x = "Age Group", y = "Number of transmission")  + 
  theme(panel.background= element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="bottom",legend.title = element_blank())
p


# Number of imports per year
prop_stouf<- finalData[, lapply(.SD, sum), by = .(season), .SDcols = "prop_recent_travel1"]
prop_stouf<- prop_stouf %>% mutate(across(where(is.numeric), ~ round(., 0)))
names(prop_stouf)[2] <- "Stouffer_import"
Epi_import<- c( 57, 82, 85 , 99, 183, 108, 153, 192) 
prop_yr_stouf<- as.data.frame(cbind(prop_stouf, Epi_import)) 

import <- melt(prop_yr_stouf,id.vars = "season")

p1<- ggplot(import,aes(x = season,y = value)) + 
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge") + 
  scale_y_log10()+
  labs(x = "Season", y = "Number of imported cases") +
  scale_fill_discrete(name = " ")  + 
  theme(panel.background= element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="bottom",legend.title = element_blank())
p1


figure <- ggarrange(p, p1,
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2)

figure 

##epicurve
ili.type <- incidence(mer_df$onset_date, interval = "week", groups = mer_df$rdt)
ili_plot<- plot(ili.type, stack = TRUE, border = "white", n_breaks = 20)+
  labs (x =" ", y = "weekly incidence")+
  theme(panel.background= element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="bottom",legend.title = element_blank(),text=element_text(size=18))
ili_plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



#demographic

mer_df$n_age_group<-factor(mer_df$n_age_group,level = c("<3","4-6","7-12","13-18","19-64","65-74", ">74"))  
age<- ggplot(mer_df, aes(x = n_age_group)) +  
  geom_bar(aes(y = (..count..)/sum(..count..)), fill="steelblue")+
  xlab("Age group") + 
  ylab("Proportion")+ 
  theme_bw()
age


sex<-ggplot(mer_df, aes(x =sex, na.rm = TRUE)) +  
  geom_bar(aes(y = (..count..)/sum(..count..)), fill="steelblue",na.rm = TRUE)+
  xlab("Gender") + 
  ylab("Proportion")+
  theme_bw()
sex

vac<- ggplot(mer_df, aes(x =vac_his, na.rm = TRUE)) +  
  geom_bar(aes(y = (..count..)/sum(..count..)), fill="steelblue",na.rm = TRUE)+
  xlab("Vaccination status") + 
  ylab("Proportion")+
  theme_bw()
vac


flu_type<-ggplot(mer_df, aes(x =mer_df$rdt, na.rm = TRUE)) +  
  geom_bar(aes(y = (..count..)/sum(..count..)), fill="steelblue",na.rm = TRUE)+
  xlab("Vaccination status") + 
  ylab("Proportion")+
  theme_bw()
flu_type

a<- ggarrange(age, sex,vac,flu_type,
              labels = c("A", "B", "C", "D"),
              ncol = 2, nrow = 2)
a

demo<- mer_df%>%select(n_age_group,sex,vac_his,rdt)
