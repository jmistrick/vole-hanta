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
library(tidyverse) #v 2.0.0
library(janitor) #v 2.2.0
library(lubridate) #v 1.9.3
library(igraph) #v 1.6.0
library(gridExtra) #v 2.3
library(ggforce) #for geom_circle in ggplot #v 0.4.1
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
params21 <- readRDS(here("spaceuse_parameters21.rds")) %>% mutate(year=2021)
params22 <- readRDS(here("spaceuse_parameters22.rds")) %>% mutate(year=2022)

#MONTHLY centroid locations for each vole
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


### SUMMARIZE mean space use kernel size (1% probability of detection) by functional group, trt, year ####
### using the a and b params for each fxnl group and plotting space use size in area (meters^2) ####

spaceusedata <- rbind(params21, params22) %>%
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b)) %>% #calculate the radius where probability of detection is 1%
  mutate(area = 2*pi*(rad_0.01^2)) %>%
  mutate(area = area*10) %>% #multiply by 10 since a,b parameters are calculated in 'trap units' and traps are 10m apart
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

#plot mean space use area by fxnl group, treatment, and year and save figure as png
png(filename = here("Figure_1_spaceuse_by_fxnl.png"), height=6, width = 12, units = "in", res=600)
ggplot(aes(x=season, y=mean, color=trt, shape=fxnl), data=spaceusedata) +
  geom_jitter(size=5, width=0.15) +
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

#for 2021 data
circles21 <- left_join(centroids21, trapdata21, by=c("tag", "month", "site")) %>%
  unite(stsb, season, trt, sex, season_breeder) %>% left_join(params21, by="stsb") %>%
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #make sure month is a factor
  separate_wider_delim(stsb, delim="_", names=c("season", "food_trt", "helm_trt", "sex", "season_breeder")) %>%
  unite(trt, food_trt, helm_trt) %>%
  unite(fxnl_grp, sex, season_breeder, remove=FALSE) %>%
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b))

#for 2022 data
circles22 <- left_join(centroids22, trapdata22, by=c("tag", "month", "site")) %>%
  unite(stsb, season, trt, sex, season_breeder) %>% left_join(params22, by="stsb") %>%
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #make sure month is a factor
  separate_wider_delim(stsb, delim="_", names=c("season", "food_trt", "helm_trt", "sex", "season_breeder")) %>%
  unite(trt, food_trt, helm_trt) %>%
  unite(fxnl_grp, sex, season_breeder, remove=FALSE) %>%
  mutate(rad_0.01 = (log((1/0.01)-1) + a) / (-b))


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

  png(filename = paste("spaceuse_kernel_", "rad0.01_byfxnl_", names(circles21_list)[[i]], "_2021", ".png", sep = ""),
      width=18 , height=5, units="in", res=600)

  p <- list()

  for(j in 1:length(circles21_list[[i]])){

    data <- circles21_list[[i]][[j]]

    #plot
    p[[j]] <- data %>%
      ggplot() +
      geom_point(aes(x=x, y=y, color=fxnl_grp), show.legend=FALSE) +
      xlim(-1.5,13) + ylim(-1.5,13) +
      geom_circle( aes(x0=x, y0=y, r=rad_0.01, fill=fxnl_grp), alpha=0.5) +
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

  png(filename = paste("spaceuse_kernel_", "rad0.01_byfxnl_", names(circles22_list)[[i]], "_2022", ".png", sep = ""),
      width=18 , height=5, units="in", res=600)

  p <- list()

  for(j in 1:length(circles22_list[[i]])){

    data <- circles22_list[[i]][[j]]

    #plot
    p[[j]] <- data %>%
      ggplot() +
      geom_point(aes(x=x, y=y, color=fxnl_grp), show.legend=FALSE) +
      xlim(-1.5,13) + ylim(-1.5,13) +
      geom_circle( aes(x0=x, y0=y, r=rad_0.01, fill=fxnl_grp), alpha=0.5) +
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

