rm(list = ls())

library(o2geosocial)
library(data.table)

load ("data/kamigoto_combined_data.RData")

## ----distrib-----------------------------------------------------------------------------------------------------
# Distribution of the latent/incubation period  average 1.63 ± 0.06 days and standard deviation (sd) 0.26 ± 0.08 days A.cori etal (https://doi.org/10.1016/j.epidem.2012.06.001)
# 2 days (1-4 days range)
f_dens <- dgamma(x = 1:100, scale = 0.04147239, shape = 39.30325)

#?distribution of infectious period average 0.99 ± 0.25 days and sd 0.96 ± 0.15 days , A.cori etal (https://doi.org/10.1016/j.epidem.2012.06.001)
# Distribution of the generation time 10.1371/journal.pone.0075339  mean 3.2 (2.4-3.9)

#SD is calculated  from the formula from cochrone http://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm 
w_dens <- dnorm(x = 1:100, mean = 3.2, sd = 4.0)


## ----age, message = FALSE, warning = FALSE-----------------------------------------------------------------------
# either jpoly2 created or new social contacft study of nishiura etal
library(readxl)
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


# Data and config, model 
data1 <- outbreaker_data(dates = finalData$visit_date , #try with visit dates here since the Onset dates were not recorded for first
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
                         outlier_relative = T, #Absolute / relative threshold 
                         outlier_threshold = 0.9 #Value of the threshold
)

# Run model 1:no inference of the importation status of cases,
out1 <- outbreaker(data = data1, config = config1, moves = moves, 
                   priors = priors, likelihoods = likelihoods)

# Set data and config for model 2
data2 <- outbreaker_data(dates = finalData$visit_date, 
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
                         burnin = 5000)
# Run model 2
out2 <- outbreaker(data = data2, config = config2, moves = moves, 
                   priors = priors, likelihoods = likelihoods)

## ----params, warning = FALSE-------------------------------------------------------------------------------------
# Summary parameters a and b, removing the burnin-period
#Model 1
print(summary(out1, burnin = 5000)$a)
print(summary(out1, burnin = 5000)$b)
# Model 2
print(summary(out2, burnin = 5000)$a)
print(summary(out2, burnin = 5000)$b)


## ----prepare_plot, warning = FALSE-------------------------------------------------------------------------------
# We create groups of cluster size: initialise the breaks for each group
group_cluster <- c(1,2, 5,50, 100, 1000) - 1

# Grouped cluster size distribution in each run
clust_infer1 <- summary(out1, group_cluster = group_cluster, 
                        burnin = 5000)$cluster
clust_infer2 <- summary(out2, group_cluster = group_cluster, 
                        burnin = 5000)$cluster
# Merge inferred and reference cluster size distributions into one matrix
clust_size_matrix <- rbind(clust_infer1["Median",], clust_infer2["Median",])


## ----plot_first, warning = FALSE, fig.width = 6.5, fig.height = 4, fig.cap= "Figure 3: Comparison of inferred cluster size distribution with the reference data"----
# Histogram of the inferred and reference cluster size distributions 
b <- barplot(clust_size_matrix, names.arg = colnames(clust_infer1), las=1,
             ylab = "Number of clusters", xlab = "Cluster size", main = "", 
             beside = T, ylim = c(0, max(c(clust_infer1, clust_infer2))))
# Add the 50% CI
arrows(b[1,], clust_infer1["1st Qu.",], b[1,], clust_infer1["3rd Qu.",], 
       angle = 90, code = 3, length = 0.1)
arrows(b[2,], clust_infer2["1st Qu.",], b[2,], clust_infer2["3rd Qu.",], 
       angle = 90, code = 3, length = 0.1)
# Add legend
legend("topright", fill = grey.colors(2), bty = "n",
       legend = c("Inferred import status",  "Simulated dataset"))


## ----index_infer, warning = FALSE--------------------------------------------------------------------------------
#' Title: Compute the proportion of iterations in the outbreaker() output 
#` where the inferred index matches the actual index in dt_cases
#'
#' @param finalData: reference dataset
#' @param out: Matrix output of outbreaker()
#' @param burnin: Numeric, length of the burnin phase
#'
#' @return Numeric vector showing the proportion of iterations pointing to
#' the correct index case
index_infer <- function(finalData, out, burnin){
  ## Generate the data frame listing every infector:
  # Select rows above burnin, and columns describing who infected whom
  out_index <- out[out$step > burnin, grep("alpha", colnames(out))]
  # ID of each infector
  ID_index <- matrix(finalData[unlist(out_index), ID], ncol = nrow(dt_cases))
  # Match inferred (ID_index) and actual infector (column infector_ID)
  match_infer_data <- t(ID_index) == dt_cases$infector_ID
  # If a case is rightly inferred as an ancestor, set match to TRUE
  match_infer_data[is.na(t(ID_index)) & is.na(dt_cases$infector_ID)] <- TRUE
  prop_correct <- rowSums(match_infer_data, na.rm = T)/ncol(match_infer_data)
  
  return(prop_correct)
}
# Same as index_infer, except it returns the proportion of inferred indexes
# who are in the same reference cluster as the case
index_clust <- function(dt_cases, out, burnin){
  ## Generate the data frame listing every infector:
  # Select rows above burnin, and columns describing who infected whom
  out_index <- out[out$step > burnin, grep("alpha", colnames(out))]
  # cluster of each infector
  clust_index <- matrix(dt_cases[unlist(out_index), cluster], 
                        ncol = nrow(dt_cases))
  # Match inferred (cluster_index) and actual cluster (column cluster)
  match_infer_data <- t(clust_index) == dt_cases$cluster
  # Exclude ancestors
  match_infer_data <- match_infer_data[!is.na(dt_cases$infector_ID),]
  
  
  prop_correct <- rowSums(match_infer_data, na.rm = T)/ncol(match_infer_data)
  
  return(prop_correct)
}
# Run index_infer for each model
index_infer1 <- index_infer(dt_cases = dt_cases, out = out1, burnin = 5000)
index_infer2 <- index_infer(dt_cases = dt_cases, out = out2, burnin = 5000)
# Run index_clust for each model
index_clust1 <- index_clust(dt_cases = dt_cases, out = out1, burnin = 5000)
index_clust2 <- index_clust(dt_cases = dt_cases, out = out2, burnin = 5000)


## ----plot index_infer, warning = FALSE, fig.width = 7, fig.height = 3, fig.cap = "Figure 4: Panel A: Proportion of iterations with the correct index for each case; Panel B: Proportion of iterations where the index is from the correct cluster"----
# Plot the sorted proportion in each model
par(bty = "n", mfrow = c(1, 2), mar = c(5,4,2,0), oma = c(0, 0, 0, 0))
# Panel A: Perfect match
plot(sort(index_infer1), type = "l", ylab = "Proportion of iterations", xlab = "Case", 
     main =  "A", las=1, col = grey.colors(3)[1], lwd = 3)
lines(sort(index_infer2), col = grey.colors(3)[2], lwd = 3)

# Panel B: Close match
plot(sort(index_clust1), type = "l", xlab = "Case", ylab = "", 
     main =  "B", las=1, col = grey.colors(3)[1], lwd = 3)
lines(sort(index_clust2), col = grey.colors(3)[2], lwd = 3)
legend("bottomright", col = grey.colors(3)[1:2], lwd = 3, bty = "n",
       legend = c("Inferred import status","Known import status"))


## ----import_sf_file, message = FALSE, warning = FALSE, fig.width = 6.5, fig.height = 2.5-------------------------
library(ggplot2)
library(sf)
# Read the shapefile and create one map for each model

map1 <- st_read(dsn="/Users/suhan/Project/kamigoto/Kamigoto shape file", layer="kamigoto")
map1$X_CODE <- as.numeric(map1$X_CODE)
map1$Y_CODE <- as.numeric(map1$Y_CODE)
map2 <- map1
map1$model <- "Model 1"
map2$model <- "Model 2"


## ----n_import_reg_per_case, warning = FALSE----------------------------------------------------------------------
# Add the proportion of iterations in model 1 where each case is an import
finalData<- as.data.table(finalData)
finalData[, prop_recent_travel1 := summary(out1, burnin = 5000)$tree$import]
# Add the proportion of iterations in model 2 where each case is an import
finalData[, prop_recent_travel2 := summary(out2, burnin = 5000)$tree$import]


## ----n_import_reg_per_reg, warning = FALSE-----------------------------------------------------------------------
# Number of imports per region in model 1
prop_reg1 <- finalData[, .(prop_per_reg = sum(prop_recent_travel1)), 
                      by = district]$prop_per_reg
# Number of imports per region in model 2
prop_reg2 <- finalData[, .(prop_per_reg = sum(prop_recent_travel2)), 
                      by = district]$prop_per_reg
names(prop_reg1) <- names(prop_reg2) <- unique(finalData$district)

# Add the number of imports in each region to the maps
map1$prop_reg <- prop_reg1[as.character(map1$S_NAME)]
map2$prop_reg <- prop_reg2[as.character(map2$S_NAME)]


## ----create_map1, warning = FALSE, fig.width = 6.5, fig.height = 2.5, fig.cap = "Figure 5: Average number of imported cases per census tract, regions where no case was reported are shown in grey."----
# Merge maps
maps <- rbind(map1, map2)

# Crop map to area of interest
lim_lon <- c(-84, -82)
lim_lat <- c(40, 41.5)
maps <- maps[maps$X_CODE > lim_lon[1] & maps$X_CODE < lim_lon[2] & 
               maps$Y_CODE > lim_lat[1] & maps$Y_CODE< lim_lat[2],]

# Plot: number of imports per region, two panels
ggplot(maps) +  geom_sf(aes(fill = prop_reg)) + facet_grid(~model)  +     
  scale_fill_gradient2(na.value = "lightgrey", midpoint = 0.8, 
                       breaks = c(0, 0.5, 1, 1.5), name = "Nb imports",
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
## Compute the median in each model
n_sec1 <- apply(n_sec_tot1[,-1], 1, median)
n_sec2 <- apply(n_sec_tot2[,-1], 1, median)
names(n_sec1) <- names(n_sec2) <- unique(finalData$district)
## Add to the matrices describing the maps
map1$n_sec <- as.numeric(n_sec1[as.character(map1$S_NAME)])
map2$n_sec <- as.numeric(n_sec2[as.character(map2$S_NAME)])


## ----create_maps2, warning = FALSE, fig.width = 6.5, fig.height = 2.5, fig.cap = "Figure 6: Median number of secondary transmission per case in each census tract"----
# Merge maps
maps_n_sec <- rbind(map1, map2)
# Crop map to area of interest
#lim_lon <- c(-84, -82)
#lim_lat <- c(40, 41.5)
#maps_n_sec <- maps_n_sec[maps_n_sec$X_CODE > lim_lon[1] &
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

## ----create_stouff_matrix, message = FALSE, warning = FALSE------------------------------------------------------
# For every column of the distance matrix, use the cumulative sum of the 
# population vector ordered by the distance. Remove the values where 
# the distance between the regions is above gamma
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


## // [[Rcpp::depends(o2geosocial)]]

## #include <Rcpp.h>

## #include <Rmath.h>

## #include <o2geosocial.h>

## // This function is used to estimate new values of the spatial parameter.

## // It is based on the structure as cpp_move_a in o2geosocial,

## // [[Rcpp::export()]]

## Rcpp::List cpp_stouffer(Rcpp::List param, Rcpp::List data, Rcpp::List config,

##                         Rcpp::RObject custom_ll, Rcpp::RObject custom_prior){

##   // Import parameters

##   Rcpp::List new_param = clone(param);

##   double gamma = config["gamma"];

##   int max_kappa = config["max_kappa"];

##   Rcpp::List new_log_s_dens = new_param["log_s_dens"];

##   Rcpp::NumericMatrix dist = data["distance"], probs = new_log_s_dens[0];

##   Rcpp::NumericMatrix ances = data["can_be_ances_reg"];

##   Rcpp::NumericVector pop = data["population"], limits = config["prior_a"];

##   // Size of the probability matrix

##   int nb_cases = pow(probs.size(), 0.5);

##   // Draw new value of a

##   Rcpp::NumericVector new_a = new_param["a"];

##   double sd_a = static_cast<double>(config["sd_a"]);

##   double old_logpost = 0.0, new_logpost = 0.0, p_accept = 0.0;

##   // proposal (normal distribution with SD: config$sd_a)

##   new_a[0] += R::rnorm(0.0, sd_a); // new proposed value

##   if (new_a[0] < limits[0] || new_a[0] > limits[1]) return param;

##   // Generate new probability matrix

##   new_param["log_s_dens"] =

##     o2geosocial::cpp_log_like(pop, dist, ances, new_a[0], new_a[0],

##                               max_kappa, gamma, "power-law", nb_cases);

##   // Compare old and new likelihood values

##   old_logpost = o2geosocial::cpp_ll_space(data, config, param,

##                                           R_NilValue, custom_ll);

##   new_logpost = o2geosocial::cpp_ll_space(data, config, new_param,

##                                           R_NilValue, custom_ll);

##   // Add prior values

##   old_logpost += o2geosocial::cpp_prior_a(param, config, custom_prior);

##   new_logpost += o2geosocial::cpp_prior_a(new_param, config, custom_prior);

##   // Accept or reject proposal

##   p_accept = exp(new_logpost - old_logpost);

##   if (p_accept < unif_rand()) return param;

##   return new_param;

## }


## ----moves_priors_data_conf, message = FALSE, warning = FALSE----------------------------------------------------
# Edit the lists of movements and priors
moves3 <- custom_moves(a = cpp_stouffer)
# Define null function
f_null <- function(param) {
  return(0.0)
}
priors3 <- custom_priors(b = f_null)


# Set data and config lists
data3 <- outbreaker_data(dates = dt_cases$Date, #Onset dates
                         age_group = dt_cases$age_group, #Age group
                         region = dt_cases$Cens_tract, #Location
                         genotype = dt_cases$Genotype, #Genotype
                         w_dens = w_dens, #Serial interval
                         f_dens = f_dens, #Latent period
                         a_dens = a_dens, #Age stratified contact matrix
                         population = pop_vect, #Population 
                         distance = dist_mat_stouffer #Distance matrix
)
config3 <- create_config(data = data3, 
                         gamma = gamma,
                         init_b = 0, move_b = FALSE, # b is not estimated
                         n_iter = 20000, #Iteration number: main run
                         n_iter_import = 10000, #Iteration number: short run
                         burnin = 5000, #burnin period: first run
                         outlier_relative = T, #Absolute / relative threshold
                         outlier_threshold = 0.9 #Value of the threshold
)
# Run the model using the Stouffer's rank method
out_stouffer <- outbreaker(data = data3, config = config3, moves = moves3, 
                           priors = priors3, likelihoods = likelihoods)


## ----stouffer_cluster, message = FALSE, warning = FALSE, fig.width = 6.5, fig.height = 4, fig.cap = "Figure 7: Comparison of inferred cluster size distribution with the reference data"----
# Grouped cluster size distribution in the Stouffer's rank model
clust_infer_stouf <- summary(out_stouffer, burnin = 5000, 
                             group_cluster = group_cluster)$cluster
# Merge inferred and reference cluster size distributions
clust_size_matrix <- rbind(clust_infer_stouf["Median",], h$counts) 
# Plot the two distributions
b <- barplot(clust_size_matrix, names.arg = colnames(clust_infer_stouf), 
             beside = T, ylab = "Number of clusters", xlab = "Cluster size", 
             main = "", las = 1)
# Add CIs
arrows(b[1,], clust_infer_stouf["1st Qu.",], b[1,], 
       clust_infer_stouf["3rd Qu.",], angle = 90, code = 3, length = 0.1)
legend("topright", fill = grey.colors(2), bty = "n",
       legend = c("Inferred import status, Stouffer's rank method", 
                  "Simulated dataset"))


## ----match_stouffer, message = FALSE, warning = FALSE, fig.width = 7, fig.height = 3, fig.cap = "Figure 8: Panel A: Proportion of iterations with the correct index for each case; Panel B: Proportion of iterations where the index is from the correct cluster"----
# Generate the proportion of perfect and close match for each case in out3
index_infer_stouf <- index_infer(dt_cases = dt_cases, out = out_stouffer, 
                                 burnin = 5000)
index_clust_stouf <- index_clust(dt_cases = dt_cases, out = out_stouffer, 
                                 burnin = 5000)
# Plot the sorted proportion in each model
par(bty = "n", mfrow = c(1, 2), mar = c(5,4,2,0), oma = c(0, 0, 0, 0))
# Panel A: Perfect match
plot(sort(index_infer_stouf), main = "A", col = grey.colors(2)[1], lwd = 3,
     xlab = "Case", ylab = "Proportion of iterations", type = "l", las=1)
# Panel B: Close match
plot(sort(index_clust_stouf), type = "l", ylab = "", xlab = "Case", 
     main =  "B", las=1, col = grey.colors(2)[1], lwd = 3)

