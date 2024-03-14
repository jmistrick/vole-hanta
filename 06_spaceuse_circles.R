####### ------------ CONSTRUCT SPACE USE CIRCLES -----------------------

# load packages
library(here)
library(tidyverse)
library(igraph)
library(lubridate)
library(janitor)
library(ggforce) #for geom_circle in ggplot
library(cowplot)

#clear environment
rm(list = ls())

######## HOME RANGE DISTRIBUTION ###########
#negative sigmoidal curve calculated following Wanelik & Farine 2022: DOI 10.1007/s00265-022-03222-5
#where the declining probability (P) of an individual being detected at a distance (d)
#from the centroid of its home range is given by:

P(d) = 1 / (1 + e^(-a-bd))

#where a describes the overall size of the home range
#b describes the steepness of the edge of the home range
#and d is the logarithmic distance from the centroid

# p <- 0.01 #eg probability of detection 1%
# a <- params21[14,2]
# b <- params21[14,3]
#
# (log((1/p)-1) + a) / (-b)

######### ^^ this equation is what I'm using to calculate e.g. the distance from the centroid with 1% probability of finding the animal,
  ##just to plot some figure showing approx HR size to visualize the overlaps



#https://stackoverflow.com/questions/65089775/plotting-custom-functions-in-ggplot-with-variables-from-dataframe



####----------- LOAD DATA -----------------

#params
params21 <- readRDS(here("spaceuse_parameters21.rds"))
params22 <- readRDS(here("spaceuse_parameters22.rds"))







# ### ALL OF THIS goes basically into the figure below
# ### just a lil thing 6/16 ###
# 
# #mean radius in 2021 by season, sex, breed
# y21 <- params21 %>%
#   mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b)) %>%
#   mutate(area = 2*pi*(rad_0.01^2)) %>%
#   separate_wider_delim(stsb, delim="_", names=c("season", "foodtrt", "helmtrt", "sex", "breeder")) %>%
#   group_by(season, sex, breeder) %>%
#   summarize(mean = mean(area),
#             sd = sd(area))
# 
# #visualize M/F breeder/non size per trt from summer to fall
# params21 %>%
#   mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b)) %>%
#   separate_wider_delim(stsb, delim="_", names=c("season", "foodtrt", "helmtrt", "sex", "breeder")) %>%
#   unite(trt, foodtrt, helmtrt) %>%
#   ggplot(aes(x=season, y=rad_0.01, color=sex, shape=breeder)) +
#            geom_point() +
#   facet_wrap(~trt)
# 
# 
# #mean radius in 2021 by season, sex, breed
# y22 <- params22 %>%
#   mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b)) %>%
#   mutate(area = 2*pi*(rad_0.01^2)) %>%
#   separate_wider_delim(stsb, delim="_", names=c("season", "foodtrt", "helmtrt", "sex", "breeder")) %>%
#   group_by(season, sex, breeder) %>%
#   summarize(mean = mean(area))
# 
# #visualize M/F breeder/non size per trt from summer to fall
# params22 %>%
#   mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b)) %>%
#   separate_wider_delim(stsb, delim="_", names=c("season", "foodtrt", "helmtrt", "sex", "breeder")) %>%
#   unite(trt, foodtrt, helmtrt) %>%
#   ggplot(aes(x=season, y=rad_0.01, color=sex, shape=breeder)) +
#   geom_point() +
#   facet_wrap(~trt)
# 
# 
# #mean across both years, all trts
# table <- rbind(params21, params22) %>%
#   mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b)) %>%
#   mutate(area = 2*pi*(rad_0.01^2)) %>%
#   mutate(area = area*10) %>%
#   separate_wider_delim(stsb, delim="_", names=c("season", "foodtrt", "helmtrt", "sex", "breeder")) %>%
#   unite(trt, foodtrt, helmtrt) %>%
#   group_by(season, breeder, sex) %>%
#   summarize(mean = mean(area),
#             sd = sd(area))


### PLOT mean HR SIZE by functional group, trt, year ####
####### NEW FIGURE NOV 2023 jesus creepers and rice, why am I still working on this####
#mean area each year by trt, sex, breed
params21 <- params21 %>% mutate(year=2021)
params22 <- params22 %>% mutate(year=2022)

data <- rbind(params21, params22) %>%
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b)) %>%
  mutate(area = 2*pi*(rad_0.01^2)) %>%
  mutate(area = area*10) %>%
  separate_wider_delim(stsb, delim="_", names=c("season", "foodtrt", "helmtrt", "sex", "breeder")) %>%
  unite(trt, foodtrt, helmtrt) %>%
  unite(fxnl, sex, breeder) %>%
  group_by(year, season, trt, fxnl) %>%
  summarize(mean = mean(area),
            sd = sd(area)) %>%
  mutate(season = factor(season, levels=c("summer", "fall"))) %>%
  mutate(trt = factor(trt, levels=c("unfed_control", "unfed_deworm",
                                    "fed_control", "fed_deworm")))

trt.labs <- as_labeller(c("unfed_control" = "Unfed-Control",
                          "unfed_deworm" = "Unfed-Deworm",
                          "fed_control" = "Fed-Control",
                          "fed_deworm" = "Fed-Deworm"))
year.labs <- as_labeller(c("2021" = "2021",
                          "2022" = "2022"))

png(filename = here("spaceuse21-22_volehantaFig1.png"), height=6, width = 11, units = "in", res=600)
ggplot(aes(x=season, y=mean, color=trt, shape=fxnl), data=data) +
  geom_jitter(size=4, width=0.15) +
  scale_x_discrete(labels=c("summer"="Summer", "fall"="Autumn")) +
  scale_shape_manual(values=c(19, 1, 17, 2),
                     labels=c("Female Breeder", "Female Nonbreeder", "Male Breeder", "Male Nonbreeder")) +
  scale_color_manual(values = c("#B2DF8A", "#33A02C", "#CAB2D6", "#6A3D9A"),
                     name = "Treatment",
                     labels = c("Unfed-Control", "Unfed-Deworm", "Fed-Control",  "Fed-Deworm")) +
  facet_grid(year ~ trt, labeller = labeller(trt=trt.labs, year=year.labs)) +
  labs(y = paste("Mean Space Use", "(m\u00B2)"), x=NULL, shape="Functional Group") +
  guides(color="none") +
  theme_bw() +
  theme(strip.text = element_text(size=16),
        axis.title.y = element_text(size=16, margin=margin(0,20,0,0)),
        axis.text = element_text(size=14, color="#808080"),
        legend.text = element_text(size=13),
        legend.title = element_text(size=14),
        legend.position = "bottom")
dev.off()


#############################################










#MONTHLY centroids
centroids21 <- readRDS(here("monthly_centroids21.rds")) %>% rename(tag = Tag_ID)
centroids22 <- readRDS(here("monthly_centroids22.rds")) %>% rename(tag = Tag_ID)

#load fulltrap data - pull tag, month, site
ft21 <- readRDS(here("fulltrap21_03.04.24.rds"))
trapdat21 <- ft21 %>% select(c(year, season, trt, site, month, tag, sex, season_breeder)) %>%
  group_by(tag, month) %>% slice(1)

ft22 <- readRDS(here("fulltrap22_03.04.24.rds"))
trapdat22 <- ft22 %>% select(c(year, season, trt, site, month, tag, sex, season_breeder)) %>%
  group_by(tag, month) %>% slice(1)



####----------- CREATE 'CIRCLES' df for PLOTTING -----------------

circles21 <- left_join(centroids21, trapdat21, by=c("tag", "month", "site")) %>%
  unite(stsb, season, trt, sex, season_breeder) %>% left_join(params21, by="stsb") %>%
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #remove may from factor
  separate_wider_delim(stsb, delim="_", names=c("season", "food_trt", "helm_trt", "sex", "season_breeder")) %>%
  unite(trt, food_trt, helm_trt) %>%
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b))

circles22 <- left_join(centroids22, trapdat22, by=c("tag", "month", "site")) %>%
  unite(stsb, season, trt, sex, season_breeder) %>% left_join(params22, by="stsb") %>%
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #remove may from factor
  separate_wider_delim(stsb, delim="_", names=c("season", "food_trt", "helm_trt", "sex", "season_breeder")) %>%
  unite(trt, food_trt, helm_trt) %>%
  unite(fxnl_grp, sex, season_breeder, remove=FALSE) %>%
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b))

####### I SHOULD DEFINTELY CONFIRM THIS CALCULATION WITH SOMEONE BECAUSE I'M RULL DUMB AT MATH ################
#matt M-S said it was fine :)



####--------------real quick plots just to see----------------------

# library(ggforce) #for geom_circle in ggplot
#https://ggforce.data-imaginist.com/reference/geom_circle.html

# circles21 %>% filter(site=="vaarinkorpi" & month=="sept") %>%
#   ggplot() +
#   geom_point(aes(x=x, y=y, color=sex)) +
#   geom_circle( aes(x0=x, y0=y, r=rad_0.01, color=sex), alpha=0.5) +
#   geom_rect(aes(xmin = 0, xmax = 11, ymin = 0, ymax = 11),
#             fill=NA, alpha = 0.4, color = "black", linetype=2) +
#   coord_fixed()


####-------------------------------------------





######## PLOT MULTIPLE SITES ACROSS MONTHS #############

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
#2e level is all the months per site (5) excluding May


####----------------------------END----------------------------------


####------------ PLOT [JUST CIRCLES] IN A LOOP (per YEAR)---------------------------

library(gridExtra)

for(i in 1:length(circles21_list)) {

  png(filename = paste("~ERRBODYcircles_monthlycentroids_", "rad0.01_", names(circles21_list)[[i]], "_2021", ".png", sep = ""),
      width=18 , height=5, units="in", res=600)

  p <- list()

  for(j in 1:length(circles21_list[[i]])){

    data <- circles21_list[[i]][[j]]

    #plot
    p[[j]] <- data %>%
      ggplot() +
      geom_point(aes(x=x, y=y, color=sex), show.legend=FALSE) +
      xlim(-1.5,13) + ylim(-1.5,13) +
      geom_circle( aes(x0=x, y0=y, r=rad_0.01, fill=sex, linetype=season_breeder), alpha=0.5) +
      scale_fill_manual(values=c("#f282a780", "#00d0ff80")) +
      geom_rect(aes(xmin = 0, xmax = 11, ymin = 0, ymax = 11),
                fill=NA, alpha = 0.4, color = "#444444", linetype=2) +
      theme_void() +
      theme(legend.position = "bottom",
            axis.ticks = element_blank(),
            axis.text  = element_blank(),
            axis.title = element_blank()) +
      labs(title=paste(names(circles21_list[[i]])[j]),
           fill="Sex") +
      coord_fixed()

  }

  do.call(grid.arrange, c(p, ncol=5, top=paste( c(names(circles21_list)[[i]]), "2021") )) #plot all the plots in list p

  dev.off()

}

####------------------repeat for 2022---------------------------

#define colors for each fxnl group
#because number of fxnl groups per plot varies, want to be sure that each group is always same color
#https://stackoverflow.com/questions/17180115/manually-setting-group-colors-for-ggplot2
fxnl.colors <- c(F_breeder="#c9184a", F_nonbreeder="#ffa9b9", M_breeder="#023e8a", M_nonbreeder="#a3d5ff")

for(i in 1:length(circles22_list)) {

  png(filename = paste("ERRBODYcircles_", "rad0.01_fxnl_", names(circles22_list)[[i]], "_2022", ".png", sep = ""),
      width=18 , height=5, units="in", res=600)

  p <- list()

  for(j in 1:length(circles22_list[[i]])){

    data <- circles22_list[[i]][[j]]

    #plot
    p[[j]] <- data %>%
      ggplot() +
      geom_point(aes(x=x, y=y, color=fxnl_grp)) +
      geom_circle( aes(x0=x, y0=y, r=rad_0.01, fill=fxnl_grp), alpha=0.5) +
      scale_fill_manual(values=fxnl.colors) +
      scale_color_manual(values=fxnl.colors) +
      geom_rect(aes(xmin = 0, xmax = 11, ymin = 0, ymax = 11),
                fill=NA, alpha = 0.4, color = "black", linetype=2) +
      theme(legend.position = "bottom") +
      labs(title=paste(names(circles22_list[[i]])[j])) +
      coord_fixed()

  }

  do.call(grid.arrange, c(p, ncol=5, top=paste( c(names(circles22_list)[[i]]), "2022") )) #plot all the plots in list p

  dev.off()

}

#####----------------------------------------------------------------------------


