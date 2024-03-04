# Run using R version 4.3.2 "Eye Holes"

# 02 - Construct Spatial Overlap Networks
# Author: __MY NAME__
# Associated Publication:
# __TITLE__  
  # __AUTHORS__

# This code sources the functions in 02-1_functions_construct_spatial_overlap_networks
# to estimate parameters describing bank vole space use, construct spatial overlap networks representing
# populations of bank voles at a given study site, and calculate network metrics from these networks


# load packages
library(here) #VERSION ___
library(tidyverse) #VERSION ___
library(igraph) #VERSION ___
library(lubridate) #VERSION ___
library(janitor) #VERSION ___

#clear environment
rm(list = ls())


#######----------- GENERATE SPACE USE DISTRIBUTION PARAMETERS ----------###############
#######---------------- CREATE SPATIAL OVERLAP NETWORKS ----------------###############
#######------------------- CALCULATE NETWORK METRICS -------------------###############

#call the functions
source(here("02-1_functions_construct_overlap_networks.R"))



########### RUN FOR 2021 DATA ################

#load the 2021 capture data from RDS
fulltrap21 <- readRDS(here("fulltrap21_03.04.24.rds"))

# #alternatively, load capture data from csv file and format data columns
# fulltrap21 <- read.csv(here("fulltrap21_03.04.24.csv")) %>%
#   mutate(year = as.numeric(year),
#          month = factor(month, levels=c("june", "july", "aug", "sept", "oct")),
#          season = factor(season, levels=c("summer", "fall")),
#          occ.sess = as.character(occ.sess),
#          occasion = as.numeric(occasion),
#          session = as.numeric(session),
#          site = as.character(site),
#          trt = factor(trt, levels=c("unfed_control", "unfed_deworm", "fed_control", "fed_deworm")),
#          food_trt = factor(food_trt, levels=c("unfed", "fed")),
#          helm_trt = factor(helm_trt, levels=c("control", "deworm")),
#          tag = as.character(tag),
#          firstcap = factor(firstcap),
#          trap = as.character(trap),
#          x = as.numeric(x),
#          y = as.numeric(y),
#          sex = as.factor(sex),
#          season_breeder = factor(season_breeder, levels=c("breeder", "nonbreeder")))


#generate space use distribution parameters
generate_params(data = fulltrap21,
                params_file = "spaceuse_parameters21.rds")
#view output file:
# params21 <- readRDS(here("spaceuse_parameters21.rds"))


#create spatial overlap networks from capture data and parameters
create_overlap_networks(data = fulltrap21,
                        centroids_file = "monthly_centroids21.rds",
                        params_file = "spaceuse_parameters21.rds",
                        networks_file = "overlap_networks21.rds")
#view output files:
# centroids21 <- readRDS(here("monthly_centroids21.rds"))
# overlapnets21 <- readRDS(here("overlap_networks21.rds"))


#calculate network metrics (weighted degree)
calculate_network_metrics(data=fulltrap21,
                          networks_file = "overlap_networks21.rds",
                          netmets_file = "network_metrics21.rds")
#view output file:
# netmets21 <- readRDS(here("network_metrics21.rds"))



########### REPEAT FOR 2022 DATA ################


#load the 2022 capture data from RDS
fulltrap22 <- readRDS(here("fulltrap22_03.04.24.rds"))

# #alternatively, load capture data from csv file and format data columns
# fulltrap22 <- read.csv(here("fulltrap22_03.04.24.csv")) %>%
#   mutate(year = as.numeric(year),
#          month = factor(month, levels=c("june", "july", "aug", "sept", "oct")),
#          season = factor(season, levels=c("summer", "fall")),
#          occ.sess = as.character(occ.sess),
#          occasion = as.numeric(occasion),
#          session = as.numeric(session),
#          site = as.character(site),
#          trt = factor(trt, levels=c("unfed_control", "unfed_deworm", "fed_control", "fed_deworm")),
#          food_trt = factor(food_trt, levels=c("unfed", "fed")),
#          helm_trt = factor(helm_trt, levels=c("control", "deworm")),
#          tag = as.character(tag),
#          firstcap = factor(firstcap),
#          trap = as.character(trap),
#          x = as.numeric(x),
#          y = as.numeric(y),
#          sex = as.factor(sex),
#          season_breeder = factor(season_breeder, levels=c("breeder", "nonbreeder")))


#generate space use distribution parameters
generate_params(data = fulltrap22,
                params_file = "spaceuse_parameters22.rds")
#view output file:
# params22 <- readRDS(here("spaceuse_parameters22.rds"))


#create spatial overlap networks from capture data and parameters
create_overlap_networks(data = fulltrap22,
                        centroids_file = "monthly_centroids22.rds",
                        params_file = "spaceuse_parameters22.rds",
                        networks_file = "overlap_networks22.rds")
#view output files:
# centroids22 <- readRDS(here("monthly_centroids22.rds"))
# overlapnets22 <- readRDS(here("overlap_networks22.rds"))


#calculate network metrics (weighted degree)
calculate_network_metrics(data=fulltrap22,
                          networks_file = "overlap_networks22.rds",
                          netmets_file = "network_metrics22.rds")
#view output file:
netmets22 <- readRDS(here("network_metrics22.rds"))
