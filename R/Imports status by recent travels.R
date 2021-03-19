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
                         n_iter = 20000, #Iteration number: main run
                         n_iter_import = 10000, #Iteration number: short run
                         burnin = 5000, #burnin period: first run
                         outlier_relative = F, #Absolute / relative threshold 
                         outlier_threshold =0.05, #Value of the threshold
                         prior_a = c(.2, 5),
                         delta=8, #tmeporal threshold for the pre-clsutering , short 
                         verbatim =TRUE
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
                         verbatim =TRUE
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
                         import = finalData$recent_travel #Import status of the cases
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
                         verbatim =TRUE
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
group_cluster <- c(1,2,5,10,50,80,100) - 1


# Grouped cluster size distribution in each run
clust_infer1 <- summary(out1, group_cluster = group_cluster, 
                        burnin = 5000)$cluster
clust_infer2 <- summary(out2, group_cluster = group_cluster, 
                        burnin = 5000)$cluster
clust_infer3 <- summary(out3, group_cluster = group_cluster, 
                        burnin = 5000)$cluster
# Merge inferred and reference cluster size distributions into one matrix
clust_size_matrix <- rbind(clust_infer1["Median",], clust_infer2["Median",],clust_infer3["Median",])


## ----plot_first, warning = FALSE, fig.width = 6.5, fig.height = 4, fig.cap= "Figure 3: Comparison of inferred cluster size distribution with the reference data"----
# Histogram of the inferred and reference cluster size distributions 
b <- barplot(clust_size_matrix, names.arg = colnames(clust_infer1), las=1,
             ylab = "Number of clusters", xlab = "Cluster size", main = "", 
             beside = T, ylim = c(0, max(c(clust_infer1, clust_infer2, clust_infer3))))
# Add the 50% CI
arrows(b[1,], clust_infer1["1st Qu.",], b[1,], clust_infer1["3rd Qu.",], 
       angle = 90, code = 3, length = 0.1)
arrows(b[2,], clust_infer2["1st Qu.",], b[2,], clust_infer2["3rd Qu.",], 
       angle = 90, code = 3, length = 0.1)
arrows(b[3,], clust_infer3["1st Qu.",], b[3,], clust_infer3["3rd Qu.",], 
       angle = 90, code = 3, length = 0.1)
# Add legend
legend("topright", fill = grey.colors(3), bty = "n",
       legend = c("Inferred import status",  "Epi import", "Epi import and inferred import"))

## ----import_sf_file, message = FALSE, warning = FALSE, fig.width = 6.5, fig.height = 2.5-------------------------
library(ggplot2)
library(sf)
# Read the shapefile and create one map for each model

map1 <- st_read(dsn="/Users/suhan/Project/kamigoto/Kamigoto shape file", layer="kamigoto")
map1$X_CODE <- as.numeric(map1$X_CODE)
map1$Y_CODE <- as.numeric(map1$Y_CODE)
map2 <- map1
map3 <- map1
map1$model <- "Model 1"
map2$model <- "Model 2"
map3$model <- "Model 3"

## ----n_import_reg_per_case, warning = FALSE----------------------------------------------------------------------
# Add the proportion of iterations in model 1 where each case is an import
finalData<- as.data.table(finalData)
finalData[, prop_recent_travel1 := summary(out1, burnin = 5000)$tree$import]
# Add the proportion of iterations in model 2  &3where each case is an import
finalData[, prop_recent_travel2 := summary(out2, burnin = 5000)$tree$import]
finalData[, prop_recent_travel3 := summary(out3, burnin = 5000)$tree$import]
## ----n_import_reg_per_reg, warning = FALSE-----------------------------------------------------------------------
# Number of imports per region in model 1
prop_reg1 <- finalData[, .(prop_per_reg = sum(prop_recent_travel1)), 
                      by = district]$prop_per_reg
# Number of imports per region in model 2
prop_reg2 <- finalData[, .(prop_per_reg = sum(prop_recent_travel2)), 
                      by = district]$prop_per_reg
# Number of imports per region in model 3
prop_reg3 <- finalData[, .(prop_per_reg = sum(prop_recent_travel3)), 
                       by = district]$prop_per_reg
names(prop_reg1) <- names(prop_reg2) <- names(prop_reg3) <- unique(finalData$district)

# Add the number of imports in each region to the maps
map1$prop_reg <- prop_reg1[as.character(map1$S_NAME)]
map2$prop_reg <- prop_reg2[as.character(map2$S_NAME)]
map3$prop_reg <- prop_reg3[as.character(map3$S_NAME)]

## ----create_map1, warning = FALSE, fig.width = 6.5, fig.height = 2.5, fig.cap = "Figure 5: Average number of imported cases per census tract, regions where no case was reported are shown in grey."----
# Merge maps
maps <- rbind(map1, map2, map3)

#maps <- maps[maps$X_CODE > lim_lon[1] & maps$X_CODE < lim_lon[2] & maps$X_CODE < lim_lon[3]& 
#               maps$Y_CODE > lim_lat[1] & maps$Y_CODE< lim_lat[2] & maps$Y_CODE< lim_lat[3],]


# Plot: number of imports per region, two panels
ggplot(maps) +  geom_sf(aes(fill = prop_reg)) + facet_grid(~model)  +     
  scale_fill_gradient2(na.value = "lightgrey", name = "Nb imports",
                       low = "white", mid = "lightblue", high = "darkblue") + 
  #coord_sf(xlim = c(-83.8, -82.2), ylim = c(40.2, 41.3)) +
  theme_classic(base_size = 9)

## ----n_sec_region, warning = FALSE-------------------------------------------------------------------------------
#' Title: Compute the number of secondary cases per case in each region
#'
#' @param dt_cases: reference dataset
#' @param out: Matrix output of outbreaker()
#' @param burnin: Numeric, length of the burnin phase
#'
#' @return A numeric matrix: the first column is the census tract ID, the
#' other columns show the number of secondary cases per case. Each row 
#' corresponds to a different iteration.
n_sec_per_reg <- function(finalData, out, burnin){
  ## Number of secondary cases per case
  n_sec <- apply(out[out$step > burnin, grep("alpha", colnames(out))], 1, 
                 function(X){
                   X <- factor(X, 1:length(X))
                   return(table(X))})
  ## Aggregate by region
  tot_n_sec_reg <- aggregate(n_sec, list(finalData$district), sum)
  ## Divide by the number of cases in each region
  tot_n_sec_reg <- cbind(tot_n_sec_reg[, 1], 
                         tot_n_sec_reg[, -1] / table(finalData$district))
  return(tot_n_sec_reg)
}


## Generate the number of secondary cases per case in each region
n_sec_tot1 <- n_sec_per_reg(finalData = finalData, out = out1, burnin = 5000)
n_sec_tot2 <- n_sec_per_reg(finalData = finalData, out = out2, burnin = 5000)
n_sec_tot3 <- n_sec_per_reg(finalData = finalData, out = out3, burnin = 5000)
## Compute the median in each model
n_sec1 <- apply(n_sec_tot1[,-1], 1, median)
n_sec2 <- apply(n_sec_tot2[,-1], 1, median)
n_sec3 <- apply(n_sec_tot3[,-1], 1, median)
names(n_sec1) <- names(n_sec2) <- names(n_sec3) <- unique(finalData$district)
## Add to the matrices describing the maps
map1$n_sec <- as.numeric(n_sec1[as.character(map1$S_NAME)])
map2$n_sec <- as.numeric(n_sec2[as.character(map2$S_NAME)])
map3$n_sec <- as.numeric(n_sec3[as.character(map3$S_NAME)])


## ----create_maps2, warning = FALSE, fig.width = 6.5, fig.height = 2.5, fig.cap = "Figure 6: Median number of secondary transmission per case in each census tract"----
# Merge maps
maps_n_sec <- rbind(map1, map2, map3)
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


