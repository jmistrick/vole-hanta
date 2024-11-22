### 4 - Plot Space Use Kernels
### AUTHOR
### 23 March 2024
### this code accompanies the manuscript: "Ecological factors alter how spatial overlap predicts viral 
  # infection dynamics in wild rodent populations"
### Run using R version 4.3.2 (2023-10-31) -- "Eye Holes"

### PURPOSE: 
# THIS CODE constructs and plots visualizations of all space use kernel estimates for bank voles captured at each site 
# in a given month. Then combines plots for months June-October at a given site in a given year as a composite figure

###------------------------------------------------------------------------------

# load packages
library(here) #v 1.0.1
# library(propagate) ### #v1.0-6 THIS WILL MASK 'SELECT' in tidyverse, load it first
library(tidyverse) #v 2.0.0
library(janitor) #v 2.2.0
library(lubridate) #v 1.9.3
library(igraph) #v 1.6.0
library(gridExtra) #v 2.3
library(ggforce) #v 0.4.1 #for geom_circle in ggplot 
library(cowplot) #v 1.1.2


#clear environment
rm(list = ls())

###------------------------------------------------------------------------------

######## ABOUT SPACE USE KERNELS ###########
#Space use is estimated using a negative sigmoidal curve following Wanelik & Farine 2022: DOI 10.1007/s00265-022-03222-5
#where the declining probability (P) of an individual being detected at 
#a distance (d) from the centroid of its space use kernel is given by:

# P(d) = 1 / (1 + e^(-a-bd))

#where (a) describes the overall size of the home range
#(b) describes the steepness of the edge of the home range
#and (d) is the logarithmic distance from the centroid

#To visualize space use, I'm using the distance from the centroid (i.e. the radius) where the probability of detection (p) is 1%
#and creating a circle with this radius to represent the space use kernel of the vole

# radius = (log((1/p)-1) + a) / (-b)

###------------------------------------------------------------------------------



####----------- LOAD DATA -----------------

#a and b parameters of negative sigmoidal curve to estimate space use
#params files generated in file: "02_construct_spatial_overlap_networks.R"
params21 <- readRDS(here("spaceuse_parameters21.rds")) %>% mutate(year=2021) 
  # %>% select(!c(aSE, bSE)) #remove SE if it was calculated in 02-1_functions_construct_overlap_networks
params22 <- readRDS(here("spaceuse_parameters22.rds")) %>% mutate(year=2022) 
  # %>% select(!c(aSE, bSE)) #remove SE if it was calculated in 02-1_functions_construct_overlap_networks

#MONTHLY centroid locations for each vole
#centroids files generated in file: "02_construct_spatial_overlap_networks.R"
centroids21 <- readRDS(here("monthly_centroids21.rds")) %>% rename(tag = Tag_ID)
centroids22 <- readRDS(here("monthly_centroids22.rds")) %>% rename(tag = Tag_ID)

#load fulltrap data (trapping metadata for each capture) - pull PIT tag, month, site
ft21 <- readRDS(here("fulltrap21_03.04.24.rds"))
trapdata21 <- ft21 %>% select(c(year, season, trt, site, month, tag, sex, season_breeder)) %>%
  group_by(tag, month) %>% slice(1)

ft22 <- readRDS(here("fulltrap22_03.04.24.rds"))
trapdata22 <- ft22 %>% select(c(year, season, trt, site, month, tag, sex, season_breeder)) %>%
  group_by(tag, month) %>% slice(1)


####--------------------------------------


### SUMMARIZE mean space use kernel size (core and peripheral) by functional group, trt, year ####
### using the a and b params for each fxnl group and plotting space use size in area (meters^2) ####

#core space use - 50%
spaceusedata_core <- rbind(params21, params22) %>%
  mutate(rad_trap = (log((1/0.5)-1) + a) / (-b)) %>% #calculate the radius ("trap units") where probability of detection is 50%
  mutate(rad_m = rad_trap*10) %>% #multiply by 10 since a,b parameters are calculated in 'trap units' and traps are 10m apart
  mutate(area = 2*pi*(rad_m^2)) %>%
  separate_wider_delim(stsb, delim="_", names=c("season", "foodtrt", "helmtrt", "sex", "breeder")) %>%
  unite(trt, foodtrt, helmtrt) %>%
  unite(fxnl, sex, breeder) %>%
  group_by(year, season, trt, fxnl) %>%
  summarize(mean = mean(area)) %>%
  mutate(season = factor(season, levels=c("summer", "fall"))) %>%
  mutate(trt = factor(trt, levels=c("unfed_control", "unfed_deworm",
                                    "fed_control", "fed_deworm")))

#peripheral space use - 95%
spaceusedata_periph <- rbind(params21, params22) %>%
  mutate(rad_trap = (log((1/0.05)-1) + a) / (-b)) %>% #calculate the radius ("trap units") where probability of detection is 95%
  mutate(rad_m = rad_trap*10) %>% #multiply by 10 since a,b parameters are calculated in 'trap units' and traps are 10m apart
  mutate(area = 2*pi*(rad_m^2)) %>%
  separate_wider_delim(stsb, delim="_", names=c("season", "foodtrt", "helmtrt", "sex", "breeder")) %>%
  unite(trt, foodtrt, helmtrt) %>%
  unite(fxnl, sex, breeder) %>%
  group_by(year, season, trt, fxnl) %>%
  summarize(mean = mean(area)) %>%
  mutate(season = factor(season, levels=c("summer", "fall"))) %>%
  mutate(trt = factor(trt, levels=c("unfed_control", "unfed_deworm",
                                    "fed_control", "fed_deworm")))

#assign labellers for figure
trt.labs <- as_labeller(c("unfed_control" = "Unfed-Control",
                          "unfed_deworm" = "Unfed-Deworm",
                          "fed_control" = "Fed-Control",
                          "fed_deworm" = "Fed-Deworm"))
year.labs <- as_labeller(c("2021" = "2021",
                          "2022" = "2022"))

#plot mean core space use area by fxnl group, treatment, and year and save figure as png
png(filename = here("Figure_2_spaceuse_by_fxnl_core.png"), height=6, width = 12, units = "in", res=600)
ggplot(aes(x=season, y=mean, color=trt, shape=fxnl, group=fxnl), data=spaceusedata_core) +
  # geom_jitter(size=5, width=0.15) +
  geom_point(size=5, stroke=1.5) +
  # geom_errorbar(aes(ymin=mean-rad_sd_sq, ymax=mean+rad_sd_sq), width=.5,
  #               position=position_dodge(0.05)) + #add this in if you want error bars (see R script "std dev of HR radii")
  geom_line(size=0.7) +
  scale_x_discrete(labels=c("summer"="Summer", "fall"="Autumn")) +
  scale_shape_manual(values=c(19, 1, 17, 2),
                     labels=c("Reproductive Female", "Non-Reproductive Female", "Reproductive Male", "Non-Reproductive Male")) +
  scale_color_manual(values = c("#B2DF8A", "#33A02C", "#CAB2D6", "#6A3D9A"),
                     name = "Treatment",
                     labels = c("Unfed-Control", "Unfed-Deworm", "Fed-Control",  "Fed-Deworm")) +
  facet_grid(year ~ trt, labeller = labeller(trt=trt.labs, year=year.labs)) +
  labs(y = paste("Mean Space Use", "(m\u00B2)"), x=NULL, shape="Functional Group:") +
  guides(color="none") +
  theme_bw() +
  theme(strip.text = element_text(size=16),
        axis.title.y = element_text(size=20, margin=margin(0,15,0,0)),
        axis.text = element_text(size=14, color="#808080"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=15),
        legend.position = "bottom")
dev.off()

#plot mean peripheral space use area by fxnl group, treatment, and year and save figure as png
png(filename = here("Figure_2_spaceuse_by_fxnl_periph.png"), height=6, width = 12, units = "in", res=600)
ggplot(aes(x=season, y=mean, color=trt, shape=fxnl, group=fxnl), data=spaceusedata_periph) +
  # geom_jitter(size=5, width=0.15) +
  geom_point(size=5, stroke=1.5) +
  # geom_errorbar(aes(ymin=mean-rad_sd_sq, ymax=mean+rad_sd_sq), width=.5,
  #               position=position_dodge(0.05)) + #add this in if you want error bars (see R script "std dev of HR radii")
  geom_line(size=0.7) +
  scale_x_discrete(labels=c("summer"="Summer", "fall"="Autumn")) +
  scale_shape_manual(values=c(19, 1, 17, 2),
                     labels=c("Reproductive Female", "Non-Reproductive Female", "Reproductive Male", "Non-Reproductive Male")) +
  scale_color_manual(values = c("#B2DF8A", "#33A02C", "#CAB2D6", "#6A3D9A"),
                     name = "Treatment",
                     labels = c("Unfed-Control", "Unfed-Deworm", "Fed-Control",  "Fed-Deworm")) +
  facet_grid(year ~ trt, labeller = labeller(trt=trt.labs, year=year.labs)) +
  labs(y = paste("Mean Space Use", "(m\u00B2)"), x=NULL, shape="Functional Group:") +
  guides(color="none") +
  theme_bw() +
  theme(strip.text = element_text(size=16),
        axis.title.y = element_text(size=20, margin=margin(0,15,0,0)),
        axis.text = element_text(size=14, color="#808080"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=15),
        legend.position = "bottom")
dev.off()


####------------------------------ END --------------------------------




######## PLOT SPACE USE KERNELS FOR ALL VOLES PER SITE in a given MONTH #############
############## VISUALIZE MONTHS JUNE-OCTOBER PER SITE, PER YEAR #####################


####----------- CREATE 'CIRCLES' df for PLOTTING SPACE USE KERNELS -----------------

#### INSTEAD - (Nov 21, 2024) USE THE 'homerange2#' files from '07_percent_overlap' file - with mixed size HRs
#NOV 21, 2024 - using the 'homerange21' file instead of circles21 just to try it
circles21 <- homerange21

circles22 <- homerange22

# #for 2021 data - 95% peripheral HR
# circles21 <- left_join(centroids21, trapdata21, by=c("tag", "month", "site")) %>%
#   unite(stsb, season, trt, sex, season_breeder) %>% left_join(params21, by="stsb") %>%
#   mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #make sure month is a factor
#   separate_wider_delim(stsb, delim="_", names=c("season", "food_trt", "helm_trt", "sex", "season_breeder")) %>%
#   unite(trt, food_trt, helm_trt) %>%
#   unite(fxnl_grp, sex, season_breeder, remove=FALSE) %>%
#   mutate(rad_0.05 = (log((1/0.05)-1) + a) / (-b))
# 
# #for 2022 data - 95% peripheral HR
# circles22 <- left_join(centroids22, trapdata22, by=c("tag", "month", "site")) %>%
#   unite(stsb, season, trt, sex, season_breeder) %>% left_join(params22, by="stsb") %>%
#   mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #make sure month is a factor
#   separate_wider_delim(stsb, delim="_", names=c("season", "food_trt", "helm_trt", "sex", "season_breeder")) %>%
#   unite(trt, food_trt, helm_trt) %>%
#   unite(fxnl_grp, sex, season_breeder, remove=FALSE) %>%
#   mutate(rad_0.05 = (log((1/0.05)-1) + a) / (-b))


##---------------- CREATE A NESTED LIST of CIRCLES (nested by site, month) --------------------

#save a vector of the site names (alphabetical order)
site_names <- unique(circles21$site) %>% sort()

#split() makes a list consisting of individual data.frames based on a condition ('site' in this case)
circles21_site_list <- split(circles21, circles21$site)

#create new list to hold nested site, month capture
circles21_list <- list()

for(i in 1:length(circles21_site_list)){
  # print(i)
  temp_list <- split(circles21_site_list[[i]], circles21_site_list[[i]]$month)
  circles21_list[[i]] <- temp_list
}

#name 1e list element as site (2e list elements are months May-Oct)
names(circles21_list) <- site_names

### OUTPUT: CIRCLES##_LIST is a nested list of length 12
#1e level is all the sites (12)
#2e level is all the months per site (5) excluding May

####--------------------repeat for 2022 data----------------------

#save a vector of the site names (alphabetical order)
site_names <- unique(circles22$site) %>% sort()

#split() makes a list consisting of individual data.frames based on a condition ('site' in this case)
circles22_site_list <- split(circles22, circles22$site)

#create new list to hold nested site, month capture
circles22_list <- list()

for(i in 1:length(circles22_site_list)){
  # print(i)
  temp_list <- split(circles22_site_list[[i]], circles22_site_list[[i]]$month)
  circles22_list[[i]] <- temp_list
}

#name 1e list element as site (2e list elements are months May-Oct)
names(circles22_list) <- site_names

### OUTPUT: CIRCLES##_LIST is a nested list of length 12
#1e level is all the sites (12)
#2e level is all the months (June-October) per site (5) excluding May


####----------------------------END----------------------------------



####------------ PLOT SPACE USE KERNELS of all VOLES/MONTH/SITE/YEAR IN A LOOP ---------------------------

#define colors for each fxnl group
#because number of fxnl groups per plot varies, want to be sure that each group is always same color
#https://stackoverflow.com/questions/17180115/manually-setting-group-colors-for-ggplot2
fxnl.colors <- c(F_breeder="#c9184a", F_nonbreeder="#ffa9b9", M_breeder="#023e8a", M_nonbreeder="#a3d5ff")

####------------------ for 2021 ---------------------------

for(i in 1:length(circles21_list)) {

  png(filename = paste("spaceuse_kernel_", "HR_rad_", names(circles21_list)[[i]], "_2021", ".png", sep = ""),
      width=18 , height=5, units="in", res=600)

  p <- list()

  for(j in 1:length(circles21_list[[i]])){

    data <- circles21_list[[i]][[j]]

    #plot
    p[[j]] <- data %>%
      ggplot() +
      geom_point(aes(x=x, y=y, color=fxnl_grp), show.legend=FALSE) +
      xlim(-3,15) + ylim(-3,15) +
      geom_circle( aes(x0=x, y0=y, r=HR_rad, fill=fxnl_grp), alpha=0.5) +
      scale_color_manual(values=fxnl.colors) +
      scale_fill_manual(values=fxnl.colors) +
      geom_rect(aes(xmin = 0, xmax = 11, ymin = 0, ymax = 11),
                fill=NA, alpha = 0.4, color = "#444444", linetype=2) +
      theme_void() +
      theme(plot.title = element_text(size= 16, hjust = 0.5),
            legend.position = "none",
            axis.ticks = element_blank(),
            axis.text  = element_blank(),
            axis.title = element_blank()) +
      labs(title=paste(names(circles21_list[[i]])[j])) +
      coord_fixed()

  }

  do.call(grid.arrange, c(p, ncol=5, top=paste( c(names(circles21_list)[[i]]), "2021") )) #plot all the plots in list p

  dev.off()

}

####------------------repeat for 2022---------------------------


for(i in 1:length(circles22_list)) {

  png(filename = paste("spaceuse_kernel_", "HR_rad_", names(circles22_list)[[i]], "_2022", ".png", sep = ""),
      width=18 , height=5, units="in", res=600)

  p <- list()

  for(j in 1:length(circles22_list[[i]])){

    data <- circles22_list[[i]][[j]]

    #plot
    p[[j]] <- data %>%
      ggplot() +
      geom_point(aes(x=x, y=y, color=fxnl_grp), show.legend=FALSE) +
      xlim(-4,14) + ylim(-4,14) +
      geom_circle( aes(x0=x, y0=y, r=HR_rad, fill=fxnl_grp), alpha=0.5) +
      scale_color_manual(values=fxnl.colors) +
      scale_fill_manual(values=fxnl.colors) +
      geom_rect(aes(xmin = 0, xmax = 11, ymin = 0, ymax = 11),
                fill=NA, alpha = 0.4, color = "#444444", linetype=2) +
      theme_void() +
      theme(plot.title = element_text(size= 16, hjust = 0.5),
            legend.position = "none",
            axis.ticks = element_blank(),
            axis.text  = element_blank(),
            axis.title = element_blank()) +
      labs(title=paste(names(circles22_list[[i]])[j])) +
      coord_fixed()

  }

  do.call(grid.arrange, c(p, ncol=5, top=paste( c(names(circles22_list)[[i]]), "2022") )) #plot all the plots in list p

  dev.off()

}

#####--------------------------------END--------------------------------------------




# #### July 16 2024 - MATT M-S wrote some code to help generate standard deviations for the space use estimates
# ## It works, but the errors are so small on most of the fxnl group/trt/season combos that the error
# ## bars don't show up in the figure OR they're NA for some of the nonbreeders OR there are 3-4 
# ## where they're really big and tbh it just makes the whole figure uglier. So for now, I'm skipping this
# ## unless a reviewer really wants it
# 
# ## this code would be best suited to be run right before the "SUMMARIZE mean space use kernel size..." figure
# 
# # library(propagate) #this library fucks with dplyr::separate() so make sure to load it before tidyverse
# 
# #load the spaceuse parameters data and keep the SE estimates for a and b parameters
# paramsFULL21 <- readRDS(here("spaceuse_parameters21.rds")) %>% mutate(year=2021) #with SE
# paramsFULL22 <- readRDS(here("spaceuse_parameters22.rds")) %>% mutate(year=2022) #with SE
# 
# #### use the propagate() function to MonteCarlo estimate the std dev of the home range kernel radius
# params21_list <- split(paramsFULL21, seq(nrow(paramsFULL21)))
# 
# for(i in 1:length(params21_list)) {
#   
#   data <- params21_list[i]
#   
#   #pull the a parameter and its SE, b param and SE
#   a <- c(data[[1]][,2], data[[1]][,4])
#   b <- c(data[[1]][,3], data[[1]][,5])
#   
#   #data for propagate() needs to be one variable per column with rows of estimate, SE, DF (optional)
#   params <- cbind(a,b)
#   
#   #tryCatch() deals with errors in a loop 
#     # --> will continue on if propagate() throws an error instead of killing the loop
#   tryCatch(
#     
#     {out <- propagate::propagate(expression( (log((1/0.01)-1) + a) / (-b) ), params)
#     
#     ## this has estimated the RADIUS of 99% probability of capture -- with std dev
#     
#     #extract standard deviation from MonteCarlo simulation
#     rad_sd <- out$sim[2]*10 #multiply by 10 since radius was in "trap units" and traps are all 10m apart
#     # #extract mean estimate from MonteCarlo simulation
#     # rad_mean <- out$sim[1] #as.numeric() for just the value
#     
#     params21_list[[i]][7] <- rad_sd }, 
#     error=function(e) {params21_list[[i]][7] <- NA} )  #if propagate() throws an error, there will be no 7th column
#   
# }
# 
# #rbind the list back into a df (rows without sd get NA)
# params21_sd <- do.call(bind_rows, params21_list) %>% rename(rad_sd=V7) %>%
#   mutate(rad_sd = ifelse(rad_sd>10, NA, rad_sd)) %>% #get rid of the crazy large estimates for now
#   mutate(rad_sd_sq = rad_sd^2) %>% #std dev of area should be units squared
#   #make columns here look like those in "spaceusedata" dataframe
#   separate_wider_delim(stsb, delim="_", names=c("season", "food", "helm", "sex", "breeder")) %>% 
#   unite("trt", food, helm, sep="_") %>% unite("fxnl", sex, breeder, sep="_") %>%
#   select(year, season, trt, fxnl, rad_sd_sq) #just keep the squared std dev
# 
# 
# #---------repeat for 2022-----------
# 
# 
# params22_list <- split(paramsFULL22, seq(nrow(paramsFULL22)))
# 
# for(i in 1:length(params22_list)) {
#   
#   data <- params22_list[i]
#   
#   a <- c(data[[1]][,2], data[[1]][,4])
#   b <- c(data[[1]][,3], data[[1]][,5])
#   
#   params <- cbind(a,b)
#   
#   tryCatch(
#     
#     {out <- propagate::propagate(expression( (log((1/0.01)-1) + a) / (-b) ), params)
#     
#     ## this has estimated the RADIUS of 90% probability of capture -- with std dev
#     
#     #extract standard deviation from MonteCarlo simulation
#     rad_sd <- out$sim[2]*10 #as.numeric() for just the value
#     # #extract mean estimate from MonteCarlo simulation
#     # rad_mean <- as.numeric(out$sim[1])
#     
#     params22_list[[i]][7] <- rad_sd }, 
#     error=function(e) {params22_list[[i]][7] <- NA} )
#   
# }
# 
# #rbind back into a df (rows without sd get NA)
# params22_sd <- do.call(bind_rows, params22_list) %>% rename(rad_sd=V7) %>%
#   mutate(rad_sd = ifelse(rad_sd>10, NA, rad_sd)) %>% #get rid of the crazy large estimates for now
#   mutate(rad_sd_sq = rad_sd^2) %>% #std dev of area should be units squared
#   #rearchitect columns to match "spaceusedata" df
#   separate_wider_delim(stsb, delim="_", names=c("season", "food", "helm", "sex", "breeder")) %>%
#   unite("trt", food, helm, sep="_") %>% unite("fxnl", sex, breeder, sep="_") %>%
#   select(year, season, trt, fxnl, rad_sd_sq) #just keep sd column
# 
# #-------------------------------
# 
# ## join params_sd data to spaceusedata for plotting
# 
# params2122_sd <- rbind(params21_sd, params22_sd)
# spaceusedata <- spaceusedata %>% left_join(params2122_sd, by=c("year", "season", "trt", "fxnl"))
#
#------- END FOR REALSIES ---------
