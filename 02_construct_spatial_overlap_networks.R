### 03 - Construct Spatial Overlap Networks
### AUTHOR
### 23 March 2024
### this code accompanies the manuscript: "Ecological factors alter how spatial overlap predicts viral 
  # infection dynamics in wild rodent populations"
### Run using R version 4.3.2 (2023-10-31) -- "Eye Holes"

### PURPOSE: 
# This code sources the functions in 02-1_functions_construct_spatial_overlap_networks
# to estimate parameters describing bank vole space use, construct spatial overlap networks representing
# populations of bank voles at a given study site, and calculate network metrics from these networks

###------------------------------------------------------------------------------


# load packages
library(here) #v 1.0.0
library(tidyverse) #v 2.0.0
library(igraph) #v1.6.0
library(lubridate) #v 1.9.3
library(janitor) #v 2.2.0
library(ggridges) #v 0.5.5
library(cowplot) #v 1.1.2

#clear environment
rm(list = ls())


###-----------------------------------------------------------


#######----------- GENERATE SPACE USE DISTRIBUTION PARAMETERS ----------###############
#######---------------- CREATE SPATIAL OVERLAP NETWORKS ----------------###############
#######------------------- CALCULATE NETWORK METRICS -------------------###############

#call the functions
source(here("02-1_functions_construct_overlap_networks.R"))


########### RUN FOR 2021 DATA ################

#load the 2021 capture data from RDS
fulltrap21 <- readRDS(here("fulltrap21_03.04.24.rds"))


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


####-----------------------------------------------------------------------------



##################### DEGREE DISTRIBUTION BY TRT / YEAR #################################

## INCLUDES ALL animals from networks 
  #(i.e., 170 more than netmets_puuv because the voles that are lacking PUUV data are in the networks and are here)

### LOAD INDIVIDUAL VOLE METADATA by year 
fulltrap21 <- readRDS(file="fulltrap21_03.04.24.rds")
fulltrap22 <- readRDS(file="fulltrap22_03.04.24.rds")

fulltrap21.22 <- rbind(fulltrap21, fulltrap22) %>% 
  group_by(year, month, tag) %>% slice(1) %>% #one entry per vole per month
  select(year, site, trt, month, tag)

### LOAD NETWORK METRICS data
netmets21 <- readRDS(file="network_metrics21.rds") %>% mutate(year=as.numeric(2021))
netmets22 <- readRDS(file="network_metrics22.rds") %>% mutate(year=as.numeric(2022))

netmets21.22 <- rbind(netmets21, netmets22) %>% select(!focal_id) #drop duplicate PIT tag # column

#combine METADATA and NETWORK METRICS
netmets_full <- left_join(netmets21.22, fulltrap21.22, by=c("year", "site", "month", "tag")) %>%
  mutate(site=as.factor(site),
         year=as.factor(year),
         month = as.factor(month),
         month = factor(month, levels=c("june", "july", "aug", "sept", "oct")))


### SUMMARISE VALUES for main text ###

# #mean wtdeg by month, trt, in 2021 and 2022
# netmets_full %>% group_by(year, month, trt) %>%
#   summarise(mean=mean(wt.deg),
#             sd=sd(wt.deg))

#mean wtdeg across ALL months, by treatment in 2021, 2022
netmets_full %>% group_by(year, trt) %>%
  summarise(mean=mean(wt.deg),
            sd=sd(wt.deg))

#### THESE VALUES REPORTED in manuscript
#mean wtdeg across ALL trt, ALL months in 2021, 2022
netmets_full %>% group_by(year) %>%
  summarise(mean=mean(wt.deg),
            sd=sd(wt.deg))


### VISUALIZE weighted degree distributions, generate Figure for supplement ###

#labeller for facets
trt_labs <- as_labeller(c(unfed_control="Unfed Control", 
                          unfed_deworm="Unfed Deworm", 
                          fed_control="Fed Control", 
                          fed_deworm="Fed Deworm"))

#increase axis ticks: https://stackoverflow.com/questions/11335836/increase-number-of-axis-ticks
#add mean for each facet: https://stackoverflow.com/questions/44196384/how-to-produce-different-geom-vline-in-different-facets-in-r

#create dataframe
meandata <- netmets_full %>% group_by(year, trt) %>%
  summarize(mean_x = mean(wt.deg))

#visualize weighted degree summarised by trt, year - vertical bars for mean
netmets_full %>%
  ggplot(aes(x=wt.deg)) +
  geom_histogram(stat="bin") +
  geom_vline(data=meandata, aes(xintercept=mean_x, color=trt), linewidth=1, 
             show.legend = FALSE) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  facet_grid(trt~year, labeller=labeller(trt=trt_labs)) +
  xlab("Weighted Degree") + ylab("Count")

#visualize wtdeg by month per trt, year - SUPPLEMENT FIGURE S1
png(filename = "FIG_wt.deg_by_month-trt.png", width=10, height=8, units="in", res=300)
netmets_full %>%
  mutate(month.rev = factor(month, levels=c("oct", "sept", "aug", "july", "june"))) %>%
  ggplot(aes(x=wt.deg, y=month.rev, fill=trt, color=trt)) +
  scale_fill_manual(values = c("#B2DF8A", "#33A02C", "#CAB2D6", "#6A3D9A")) +
  scale_color_manual(values = c("#B2DF8A", "#33A02C", "#CAB2D6", "#6A3D9A")) +
  geom_density_ridges(stat = "binline", alpha=0.5,
                      scale=0.9, rel_min_height = 0.001) +
  geom_vline(data=meandata, aes(xintercept=mean_x, color=trt), linewidth=1, 
             show.legend = FALSE) +
  facet_grid(trt~year, labeller=labeller(trt=trt_labs)) +
  labs(x="Weighted Degree") +
  scale_y_discrete(labels=c("june" = "June", "july" = "July",
                            "aug" = "Aug", "sept" = "Sept", "oct" = "Oct")) +
  theme_half_open() + panel_border() + background_grid() +
  theme(legend.position = "none",
        axis.title.x = element_text(size=18),
        axis.title.y = element_blank(),
        axis.text = element_text(size=14),
        strip.text.x = element_text(size=16),
        strip.text.y = element_text(size=14))
dev.off()

############################################################################
