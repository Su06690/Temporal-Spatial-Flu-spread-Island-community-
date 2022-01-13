rm(list = ls())

library(o2geosocial)
library(MASS)
library(tidyverse)
library(data.table)
library(sf)
library(ggplot2)
library(foreign)
library(tmap)
library(readxl)


## ----distrib-----------------------------------------------------------------------------------------------------###
# Distribution of the latent/incubation period  average 1.63 ± 0.06 days and standard deviation (sd) 0.26 ± 0.08 days A.cori etal (https://doi.org/10.1016/j.epidem.2012.06.001)
# 2 days (1-4 days range)
f_dens <- dgamma(x = 1:100, scale = 0.04147239, shape = 39.30325)

# Distribution of the generation time 10.1371/journal.pone.0075339  mean 3.2 (2.4-3.9)
#SD is calculated  from the formula from cochrone http://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm 
w_dens <- dnorm(x = 1:100, mean = 3.2, sd = 0.4)

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


## ----first_model, message = FALSE, warning = FALSE---------------------------------------------------------------
# Set movement, likelihood and prior lists to default
moves <- custom_moves()
likelihoods <- custom_likelihoods()
priors <- custom_priors()

load ("data/finalData.RData")

data3 <- outbreaker_data(dates = finalData$onset_date, 
                         age_group = finalData$age_gp,
                         region = finalData$district,
                         genotype = finalData$rdt, 
                         w_dens = w_dens, 
                         f_dens = f_dens, 
                         a_dens = a_dens,
                         population = pop_vect, 
                         distance = dist_mat,
                         import = finalData$recent_travel=="TRUE"  #Import status of the cases
)


config3 <- create_config(data = data3, 
                         find_import = TRUE, #inference of import status
                         n_iter = 40000, #Iteration number: main run
                         n_iter_import = 20000, #Iteration number: short run
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

