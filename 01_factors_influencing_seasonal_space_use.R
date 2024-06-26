### 02 - Factors Influencing Seasonal Space Use
### AUTHOR
### 23 March 2024
### this code accompanies the manuscript: "Ecological factors alter how spatial overlap predicts viral 
  # infection dynamics in wild rodent populations"
### Run using R version 4.3.2 (2023-10-31) -- "Eye Holes"

### PURPOSE: 
### THIS CODE tests the how the relationship between the probability of capturing a vole in a trap and the
  # distance from the vole's seasonal centroid (ie space use) differs by sex, reproductive status, and treatment
  # SEPARATELY in summer and fall seasons 2021; summer and fall seasons 2022

#RESULT: Yes, differences in space use exist by sex, repro status, treatment in both seasons
  #thus, it makes sense to calculate space use parameters separately per season/sex/repro/trt

###------------------------------------------------------------------------------

# load packages
library(here) #version 1.0.1
library(tidyverse) #version 2.0.0
library(visreg) #v 2.7.0 #https://pbreheny.github.io/visreg/articles/web/overlay.html
library(cowplot) #v 1.1.2
library(ggtext) #v 0.1.2
library(gtsummary) #v 1.7.2 #https://www.danieldsjoberg.com/gtsummary/articles/tbl_regression.html


#clear environment
rm(list = ls())

##---------------- LOAD THE CAPTURE DATA ----------------------

##OPTION: load either 2021 or 2022 data and run the remainder of the script, 
  #then clear environment, load the other other year and rerun

#load the 2021 fulltrap dataset as an R data file
fulltrap <- readRDS(file = "fulltrap21_03.04.24.rds")

# #load the 2022 fulltrap dataset as an R data file
# fulltrap <- readRDS(file = "fulltrap22_03.04.24.rds")


#pull all the traps and their x, y coordinates, save as df with (trapID, x, y)
traps <- fulltrap %>% group_by(trap) %>% slice(1) %>%
  select(trap, x, y) %>%
  arrange(trap)


##------------------------------------

# Reading in required functions in wanelik_farine_2022_functions.R
source(here("01-1_wanelik_farine_2022_functions.R"))

## The following code was developed by Wanelik & Farine 2022
  # Wanelik, K.M., Farine, D.R. A new method for characterising shared space use networks using animal trapping data.
  # Behav Ecol Sociobiol 76, 127 (2022). https://doi.org/10.1007/s00265-022-03222-5

#with minor additions by __AUTHOR__


##------------ CALCULATE DISTANCES FOR ALL VOLES, RUN GLM BY SEX, BREEDER, TRT ------------------------
##------------------------------------- SUMMER ----------------------------------------------

    data <- fulltrap %>%
      filter(season=="summer") %>%
      select(tag, food_trt, helm_trt, sex, season_breeder, trap, x, y) %>%
      mutate(season_breeder = factor(season_breeder, levels=c("nonbreeder", "breeder"))) %>%
      group_by(tag, trap) %>%
      mutate(Det.count = length(tag)) %>% slice(1) %>% #Det.count is number of times animal was in that trap
      ungroup()

    tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded
    tagdata <- data %>% group_by(tag) %>% slice(1) %>%
      select(tag, food_trt, helm_trt, sex, season_breeder) #all tags with their data

    traps_rep <- do.call("rbind", replicate(length(tags$tag), traps, simplify = FALSE)) #all traps replicated for each tag
    tags_rep <- do.call("rbind", replicate(length(traps$trap), tags, simplify = FALSE)) %>% arrange(tag) #all tags replicated for each trap

    #VERY LARGE df of all tags for all traps
    #new column of Det.obs - whether that animal was observed in that trap
    #new column of Det.count - number of times animal was captured in that trap (0-4 in a trapping occasion)
    fulldata <- cbind(tags_rep, traps_rep) %>%
      left_join(data, by=c("tag", "trap", "x", "y")) %>%
      group_by(tag) %>%
      # fill(c(sex, season_breeder, trt), .direction="downup") %>%
      fill(c(sex, season_breeder, food_trt, helm_trt), .direction="downup") %>%
      ungroup() %>%
      mutate(Det.count = replace_na(Det.count, 0),
             Det.obs = ifelse(Det.count==0, 0, 1)) %>%
      arrange(tag, trap) %>%
      select(tag, sex, season_breeder, food_trt, helm_trt, x, y, Det.obs, Det.count) %>%
      rename(x.trap = x,
             y.trap = y,
             Tag_ID = tag)

    #----------------------- 3. Generating overlap network  ------------------------------#

    # Recalculating (weighted) centroids
    # Weighted because a trap in which an individual was caught multiple times will have a greater influence on its centroid than a trap in which it was caught only once
    matrix_dists_real2 <- fulldata #only one entry per tag/trap - Det.obs is (1-observed,0-not) Det.count is # of times in that trap
    matrix_dists_real2 <- matrix_dists_real2 %>% mutate(Det.count2 = Det.count)
    matrix_dists_real2$Det.count2[matrix_dists_real2$Det.count==0]<-1
    matrix_dists_real2 <- as.data.frame(lapply(matrix_dists_real2, rep, matrix_dists_real2$Det.count2)) #has as many entries (rows) per tag/trap as the animal was captured there
    matrix_dists_real2<-matrix_dists_real2[,-9] #drops the Det.count2 column
    x <- sapply(by(matrix_dists_real2$x.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
    y <- sapply(by(matrix_dists_real2$y.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
    # obs_centroids <- data.frame(Tag_ID=tagdata$tag, sex=tagdata$sex, season_breeder=tagdata$season_breeder, trt=tagdata$trt)
    obs_centroids <- data.frame(Tag_ID=tagdata$tag, sex=tagdata$sex, season_breeder=tagdata$season_breeder, food_trt=tagdata$food_trt, helm_trt=tagdata$helm_trt)
    obs_centroids$x <- x[match(obs_centroids$Tag_ID,names(x))]
    obs_centroids$y <- y[match(obs_centroids$Tag_ID,names(y))]

    # Creating a matrix of distances between each individual's observed centroid and each trap
    matrix_dists_obs <- fulldata #only one entry per individual
    matrix_dists_obs$x.ind <- obs_centroids$x[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)] #x coord of centroid of individual
    matrix_dists_obs$y.ind <- obs_centroids$y[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)] #y coord of centroid of individual
    matrix_dists_obs$Dist <- sqrt((matrix_dists_obs$x.ind-matrix_dists_obs$x.trap)^2 + (matrix_dists_obs$y.ind - matrix_dists_obs$y.trap)^2)
    matrix_dists_obs$Dist.log<-log(matrix_dists_obs$Dist + 1)
    # # Removing unobserved individuals
    # matrix_dists_obs <- matrix_dists_obs[which(!is.na(matrix_dists_obs$Dist)),] #not needed, matrix is only of voles capped that month

    # Estimating home range profiles (negative sigmoidal curves) using this data
    fit.summer <- glm(Det.obs ~ Dist.log + sex + season_breeder +
                        food_trt + helm_trt + food_trt*helm_trt +
                        Dist.log*sex + Dist.log*season_breeder +
                        Dist.log*food_trt + Dist.log*helm_trt +
                        Dist.log*food_trt*helm_trt,
                      data=matrix_dists_obs, family=binomial, control = list(maxit = 50))

    summary(fit.summer)

  ################# SHOULD THERE BE AN INTERACTION of SEX/TRT/BREED and DIST?

## Figures for Sex, Repro, and Food Trt are generated by different code for 2021 vs 2022 - 2022 code commented out
    
    #by sex - 2021 season (pvalue is hard-coded because its so small)
    sex.summer.21 <- visreg(fit.summer, "Dist.log", by="sex", scale='response', rug=FALSE,
                            gg=TRUE, overlay=TRUE) +
      scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
      scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
      scale_color_manual(values=c("#f282a7", "#00d0ff")) +
      scale_fill_manual(values=c("#f282a750", "#00d0ff50")) +
      theme_half_open() +
      theme(legend.position = "bottom",
            axis.title = element_text(size=18),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text = element_text(size=16),
            plot.margin = margin(t = 10, r = 20, b = 20, l = 10, unit = "pt")) +
      labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
      annotate(geom = "text", x=2.4, y=.85, size = 6,
               label = paste("p <0.001")) +
      annotate(geom = "text", x=2.4, y=.9, size = 6,
               label = paste("OR =",
                             round( exp(coef(summary(fit.summer))[8,1]), digits=3) )) +
      #https://github.com/pbreheny/visreg/issues/56 #color points by group separately
      geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=sex),
                 size=3, alpha=0.1, shape=16)
    
  # #by sex - 2022 season (pvalue is hard-coded because it's very small)
  # sex.summer.22 <- visreg(fit.summer, "Dist.log", by="sex", scale='response', rug=FALSE,
  #                 gg=TRUE, overlay=TRUE) +
  #     scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
  #     scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
  #     scale_color_manual(values=c("#f282a7", "#00d0ff")) +
  #     scale_fill_manual(values=c("#f282a750", "#00d0ff50")) +
  #     theme_half_open() +
  #     theme(legend.position = "bottom",
  #           axis.title = element_text(size=18),
  #           axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
  #           axis.text = element_text(size=16),
  #           plot.margin = margin(t = 10, r = 20, b = 20, l = 10, unit = "pt")) +
  #     labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
  #     annotate(geom = "text", x=2.4, y=.85, size = 6,
  #              label = paste("p < 0.001")) +
  #     annotate(geom = "text", x=2.4, y=.9, size = 6,
  #              label = paste("OR =",
  #                            round( exp(coef(summary(fit.summer))[8,1]), digits=3) )) +
  #   #https://github.com/pbreheny/visreg/issues/56 #color points by group separately
  #     geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=sex),
  #                size=3, alpha=0.1, shape=16)


    #by repro - 2021 season (pvalue is hard-coded because its so small)
    repro.summer.21 <- visreg(fit.summer, "Dist.log", by="season_breeder", scale="response", rug=FALSE,
                              gg=TRUE, overlay=TRUE) +
      scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
      scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
      scale_color_manual(values=c("#8b64b9", "#e8ac65")) +
      scale_fill_manual(values=c("#8b64b950", "#e8ac6550")) +
      theme_half_open() +
      theme(legend.position = "bottom",
            axis.title = element_text(size=18),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text = element_text(size=16),
            plot.margin = margin(t = 10, r = 10, b = 20, l = 20, unit = "pt")) +
      labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
      annotate(geom = "text", x=2.4, y=.85, size = 6,
               label = paste("p <0.001")) +
      annotate(geom = "text", x=2.4, y=.9, size = 6,
               label = paste("OR =",
                             round( exp(coef(summary(fit.summer))[9,1]), digits=3) )) +
      geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=season_breeder),
                 size=3, alpha=0.15, shape=16)
    
    # #by repro - 2022 season (pvalue is hard-coded because it's very small)
    # repro.summer.22 <- visreg(fit.summer, "Dist.log", by="season_breeder", scale="response", rug=FALSE,
    #                 gg=TRUE, overlay=TRUE) +
    #   scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
    #   scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
    #   scale_color_manual(values=c("#8b64b9", "#e8ac65")) +
    #   scale_fill_manual(values=c("#8b64b950", "#e8ac6550")) +
    #   theme_half_open() +
    #   theme(legend.position = "bottom",
    #         axis.title = element_text(size=18),
    #         axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    #         axis.text = element_text(size=16),
    #         plot.margin = margin(t = 10, r = 10, b = 20, l = 20, unit = "pt")) +
    #   labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
    #   annotate(geom = "text", x=2.4, y=.85, size = 6,
    #            label = paste("p < 0.001")) +
    #   annotate(geom = "text", x=2.4, y=.9, size = 6,
    #            label = paste("OR =",
    #                          round( exp(coef(summary(fit.summer))[9,1]), digits=3) )) +
    #   geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=season_breeder),
    #              size=3, alpha=0.15, shape=16)

    #by food - 2021 season (pvalue is hard-coded because its so small)
    food.summer.21 <- visreg(fit.summer, "Dist.log", by="food_trt", scale='response', rug=FALSE,
                          gg=TRUE, overlay=TRUE) +
      scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
      scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
      scale_color_manual(values=c("#794624", "#68b63e")) +
      scale_fill_manual(values=c("#79462450", "#68b63e50")) +
      theme_half_open() +
      theme(legend.position = "bottom",
            axis.title = element_text(size=18),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text = element_text(size=16),
            plot.margin = margin(t = 20, r = 20, b = 10, l = 10, unit = "pt")) +
      labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
      annotate(geom = "text", x=2.4, y=.85, size = 6,
               label = paste("p <0.001")) +
      annotate(geom = "text", x=2.4, y=.9, size = 6,
               label = paste("OR =",
                             round( exp(coef(summary(fit.summer))[10,1]), digits=3) )) +
      geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=food_trt),
                 size=3, alpha=0.15, shape=16)

    # #by food - 2022 season
    # food.summer <- visreg(fit.summer, "Dist.log", by="food_trt", scale='response', rug=FALSE,
    #                gg=TRUE, overlay=TRUE) +
    #   scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
    #   scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
    #   scale_color_manual(values=c("#794624", "#68b63e")) +
    #   scale_fill_manual(values=c("#79462450", "#68b63e50")) +
    #   theme_half_open() +
    #   theme(legend.position = "bottom",
    #         axis.title = element_text(size=18),
    #         axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    #         axis.text = element_text(size=16),
    #         plot.margin = margin(t = 20, r = 20, b = 10, l = 10, unit = "pt")) +
    #   labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
    #   annotate(geom = "text", x=2.4, y=.85, size = 6,
    #            label = paste("p =",
    #                          round( coef(summary(fit.summer))[10,4], digits=3) )) +
    #   annotate(geom = "text", x=2.4, y=.9, size = 6,
    #            label = paste("OR =",
    #                          round( exp(coef(summary(fit.summer))[10,1]), digits=3) )) +
    #   geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=food_trt),
    #              size=3, alpha=0.15, shape=16)


    #by worms
    worms.summer <- visreg(fit.summer, "Dist.log", by="helm_trt", scale='response', rug=FALSE,
                    gg=TRUE, overlay=TRUE) +
      scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
      scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
      scale_color_manual(values=c("#808080", "#ffe048")) +
      scale_fill_manual(values=c("#80808070", "#ffe04850")) +
      theme_half_open() +
      theme(legend.position = "bottom",
            axis.title = element_text(size=18),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            axis.text = element_text(size=16),
            plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt")) +
      labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
      annotate(geom = "text", x=2.4, y=.85, size = 6,
               label = paste("p =",
                             round( coef(summary(fit.summer))[11,4], digits=3) )) +
      annotate(geom = "text", x=2.4, y=.9, size = 6,
               label = paste("OR =",
                             round( exp(coef(summary(fit.summer))[11,1]), digits=3) )) +
      geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=helm_trt),
                 size=3, alpha=0.1, shape=16)

    #generate 2021 figure
    png(filename="space_useGLM_summer2021.png", height=12, width=16, units="in", res=600)
    plot_grid(sex.summer.21, repro.summer.21, food.summer.21, worms.summer,
              labels="AUTO", nrow=2)
    dev.off()
    
    # #generate 2022 figure
    # png(filename="space_useGLM_summer2022.png", height=12, width=16, units="in", res=600)
    # plot_grid(sex.summer.22, repro.summer.22, food.summer, worms.summer,
    #           labels="AUTO", nrow=2)
    # dev.off()


##-----------------------------------------------------------------------------------------------------

    ##------------ CALCULATE DISTANCES FOR ALL VOLES, RUN GLM BY SEX,BREEDER,TRT ------------------------
    ##------------------------------------- FALL ----------------------------------------------


    data <- fulltrap %>%
      filter(season=="fall") %>%
      select(tag, food_trt, helm_trt, sex, season_breeder, trap, x, y) %>%
      mutate(season_breeder = factor(season_breeder, levels=c("nonbreeder", "breeder"))) %>%
      group_by(tag, trap) %>%
      mutate(Det.count = length(tag)) %>% slice(1) %>% #Det.count is number of times animal was in that trap
      ungroup()

    tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded
    tagdata <- data %>% group_by(tag) %>% slice(1) %>%
      select(tag, food_trt, helm_trt, sex, season_breeder) #all tags with their data

    traps_rep <- do.call("rbind", replicate(length(tags$tag), traps, simplify = FALSE)) #all traps replicated for each tag
    tags_rep <- do.call("rbind", replicate(length(traps$trap), tags, simplify = FALSE)) %>% arrange(tag) #all tags replicated for each trap

    #VERY LARGE df of all tags for all traps
    #new column of Det.obs - whether that animal was observed in that trap
    #new column of Det.count - number of times animal was captured in that trap (0-4 in a trapping occasion)
    fulldata <- cbind(tags_rep, traps_rep) %>%
      left_join(data, by=c("tag", "trap", "x", "y")) %>%
      group_by(tag) %>%
      fill(c(sex, season_breeder, food_trt, helm_trt), .direction="downup") %>%
      ungroup() %>%
      mutate(Det.count = replace_na(Det.count, 0),
             Det.obs = ifelse(Det.count==0, 0, 1)) %>%
      arrange(tag, trap) %>%
      select(tag, sex, season_breeder, food_trt, helm_trt, x, y, Det.obs, Det.count) %>%
      rename(x.trap = x,
             y.trap = y,
             Tag_ID = tag)

    #----------------------- 3. Generating overlap network  ------------------------------#

    # Recalculating (weighted) centroids
    # Weighted because a trap in which an individual was caught multiple times will have a greater influence on its centroid than a trap in which it was caught only once
    matrix_dists_real2 <- fulldata #only one entry per tag/trap - Det.obs is (1-observed,0-not) Det.count is # of times in that trap
    matrix_dists_real2 <- matrix_dists_real2 %>% mutate(Det.count2 = Det.count)
    matrix_dists_real2$Det.count2[matrix_dists_real2$Det.count==0]<-1
    matrix_dists_real2 <- as.data.frame(lapply(matrix_dists_real2, rep, matrix_dists_real2$Det.count2)) #has as many entries (rows) per tag/trap as the animal was captured there
    matrix_dists_real2<-matrix_dists_real2[,-9] #drops the Det.count2 column
    x <- sapply(by(matrix_dists_real2$x.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
    y <- sapply(by(matrix_dists_real2$y.trap[matrix_dists_real2$Det.obs==1], matrix_dists_real2$Tag_ID[matrix_dists_real2$Det.obs==1],mean),identity)
    # obs_centroids <- data.frame(Tag_ID=tagdata$tag, sex=tagdata$sex, season_breeder=tagdata$season_breeder, trt=tagdata$trt)
    obs_centroids <- data.frame(Tag_ID=tagdata$tag, sex=tagdata$sex, season_breeder=tagdata$season_breeder, food_trt=tagdata$food_trt, helm_trt=tagdata$helm_trt)
    obs_centroids$x <- x[match(obs_centroids$Tag_ID,names(x))]
    obs_centroids$y <- y[match(obs_centroids$Tag_ID,names(y))]

    # Creating a matrix of distances between each individual's observed centroid and each trap
    matrix_dists_obs <- fulldata #only one entry per individual
    matrix_dists_obs$x.ind <- obs_centroids$x[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)] #x coord of centroid of individual
    matrix_dists_obs$y.ind <- obs_centroids$y[match(matrix_dists_obs$Tag_ID,obs_centroids$Tag_ID)] #y coord of centroid of individual
    matrix_dists_obs$Dist <- sqrt((matrix_dists_obs$x.ind-matrix_dists_obs$x.trap)^2 + (matrix_dists_obs$y.ind - matrix_dists_obs$y.trap)^2)
    matrix_dists_obs$Dist.log<-log(matrix_dists_obs$Dist + 1)
    # # Removing unobserved individuals
    # matrix_dists_obs <- matrix_dists_obs[which(!is.na(matrix_dists_obs$Dist)),] #not needed, matrix is only of voles capped that month

    # Estimating home range profiles (negative sigmoidal curves) using this data
    fit.fall <- glm(Det.obs ~ Dist.log + sex + season_breeder +
                      food_trt + helm_trt + food_trt*helm_trt +
                      Dist.log*sex + Dist.log*season_breeder +
                      Dist.log*food_trt + Dist.log*helm_trt +
                      Dist.log*food_trt*helm_trt,
                    data=matrix_dists_obs, family=binomial, control = list(maxit = 50))
    summary(fit.fall)

    ################# SHOULD THERE BE AN INTERACTION of SEX/TRT/BREED and DIST?

### all figures are generated by the same code in 2021 and 2022

#by sex
sex.fall <- visreg(fit.fall, "Dist.log", by="sex", scale='response', rug=FALSE,
       gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
  scale_color_manual(values=c("#f282a7", "#00d0ff")) +
  scale_fill_manual(values=c("#f282a750", "#00d0ff50")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text = element_text(size=16),
        plot.margin = margin(t = 10, r = 20, b = 20, l = 10, unit = "pt")) +
  labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
  annotate(geom = "text", x=2.4, y=.85, size = 6,
           label = paste("p =",
                         round( coef(summary(fit.fall))[8,4], digits=3) )) +
  annotate(geom = "text", x=2.4, y=.9, size = 6,
           label = paste("OR =",
                         round( exp(coef(summary(fit.fall))[8,1]), digits=3) )) +
  geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=sex),
             size=3, alpha=0.15, shape=16)


#by repro
repro.fall <- visreg(fit.fall, "Dist.log", by="season_breeder", scale="response", rug=FALSE,
       gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
  scale_color_manual(values=c("#8b64b9", "#e8ac65")) +
  scale_fill_manual(values=c("#8b64b950", "#e8ac6550")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text = element_text(size=16),
        plot.margin = margin(t = 10, r = 10, b = 20, l = 20, unit = "pt")) +
  labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
  annotate(geom = "text", x=2.4, y=.85, size = 6,
           label = paste("p =",
                         round( coef(summary(fit.fall))[9,4], digits=3) )) +
  annotate(geom = "text", x=2.4, y=.9, size = 6,
           label = paste("OR =",
                         round( exp(coef(summary(fit.fall))[9,1]), digits=3) )) +
  geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=season_breeder),
             size=3, alpha=0.15, shape=16)


#by food
food.fall <- visreg(fit.fall, "Dist.log", by="food_trt", scale='response', rug=FALSE,
       gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
  scale_color_manual(values=c("#794624", "#68b63e")) +
  scale_fill_manual(values=c("#79462450", "#68b63e50")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text = element_text(size=16),
        plot.margin = margin(t = 20, r = 20, b = 10, l = 10, unit = "pt")) +
  labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
  annotate(geom = "text", x=2.4, y=.85, size = 6,
           label = paste("p =",
                         round( coef(summary(fit.fall))[10,4], digits=3) )) +
  annotate(geom = "text", x=2.4, y=.9, size = 6,
           label = paste("OR =",
                         round( exp(coef(summary(fit.fall))[10,1]), digits=3) )) +
  geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=food_trt),
             size=3, alpha=0.15, shape=16)


#by worms
worms.fall <- visreg(fit.fall, "Dist.log", by="helm_trt", scale='response', rug=FALSE,
       gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.01,0.01))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(expand = expansion(mult=c(0.01,0.01))) +
  scale_color_manual(values=c("#808080", "#ffe048")) +
  scale_fill_manual(values=c("#80808070", "#ffe04850")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text = element_text(size=16),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt")) +
  labs(x="Log Distance from Seasonal Centroid", y="Probability of Capture") +
  annotate(geom = "text", x=2.4, y=.85, size = 6,
           label = paste("p =",
                         round( coef(summary(fit.fall))[11,4], digits=3) )) +
  annotate(geom = "text", x=2.4, y=.9, size = 6,
           label = paste("OR =",
                         round( exp(coef(summary(fit.fall))[11,1]), digits=3) )) +
  geom_point(data=matrix_dists_obs, aes(x=Dist.log, y=Det.obs, color=helm_trt),
             size=3, alpha=0.1, shape=16)


#generate 2021 figure
png(filename="space_useGLM_fall2021.png", height=12, width=16, units="in", res=600)
plot_grid(sex.fall, repro.fall, food.fall, worms.fall, labels="AUTO", nrow=2)
dev.off()

# #generate 2022 figure
# png(filename="space_useGLM_fall2022.png", height=12, width=16, units="in", res=600)
# plot_grid(sex.fall, repro.fall, food.fall, worms.fall, labels="AUTO", nrow=2)
# dev.off()





######## pretty OUTPUT MODEL SUMMARY #########

#summary table 2021 summer
fit.summer %>% tbl_regression(exponentiate = TRUE,
                              pvalue_fun = ~ style_pvalue(.x, digits = 2),) %>%
  bold_p(t = 0.10) %>%
  bold_labels() %>%
  italicize_levels() %>%
  gtsummary::as_tibble() %>%
  write.csv(here("space_useGLM_TABLE_summer2021.csv"))
#summary table 2021 fall
fit.fall %>% tbl_regression(exponentiate = TRUE,
                            pvalue_fun = ~ style_pvalue(.x, digits = 2),) %>%
  bold_p(t = 0.10) %>%
  bold_labels() %>%
  italicize_levels() %>% #stop here to get HTML output in RStudio
  gtsummary::as_tibble() %>%
  write.csv(here("space_useGLM_TABLE_fall2021.csv"))

# #summary table 2022 summer
# fit.summer %>% tbl_regression(exponentiate = TRUE,
#                                 pvalue_fun = ~ style_pvalue(.x, digits = 2),) %>%
#       bold_p(t = 0.10) %>%
#       bold_labels() %>%
#       italicize_levels() %>%
#   gtsummary::as_tibble() %>%
#   write.csv(here("space_useGLM_TABLE_summer2022.csv"))
# #summary table 2022 fall
# fit.fall %>% tbl_regression(exponentiate = TRUE,
#                               pvalue_fun = ~ style_pvalue(.x, digits = 2),) %>%
#       bold_p(t = 0.10) %>%
#       bold_labels() %>%
#       italicize_levels() %>% #stop here to get HTML output in RStudio
#   gtsummary::as_tibble() %>%
#   write.csv(here("space_useGLM_TABLE_fall2022.csv"))


##### GLM model diagnostics ######

# #https://rpubs.com/benhorvath/glm_diagnostics
# plot(density(resid(fit.summer, type='pearson')))
# plot(density(resid(fit.fall, type='pearson')))
# 
# plot(density(rstandard(fit.summer, type='pearson')))
# plot(density(rstandard(fit.fall, type='pearson')))
# 
# plot(density(resid(fit.summer, type='deviance')))
# plot(density(resid(fit.fall, type='deviance')))
# 
# plot(density(rstandard(fit.summer, type='deviance')))
# plot(density(rstandard(fit.fall, type='deviance')))
# 
# library(statmod) #v 1.5.0
# plot(density(qresid(fit.summer)))
# plot(density(qresid(fit.fall)))
# 
# par(mfrow=c(1,2))
# plot(density(matrix_dists_obs$Det.obs), main='M_0 y_hat')
# lines(density(predict(fit.fall, type='response')), col='red')
# 
# qqnorm(statmod::qresid(fit.summer)); qqline(statmod::qresid(fit.summer))
# qqnorm(statmod::qresid(fit.fall)); qqline(statmod::qresid(fit.fall))
