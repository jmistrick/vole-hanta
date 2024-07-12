# Run using R version 4.3.2 "Eye Holes"

# ## - Efficacy of deworming treatment
# Author: _______________________
# Code developed with input from __________________________________
# Associated Publication:
# "Ecological factors alter how spatial overlap predicts viral 
    # infection dynamics in wild rodent populations"


#Effect of deworm treatment on helminth presence and helminth infection intensity

#Testing on FULL 2021 & 2022 FEC dataset ( ___ entries)
#This code runs with Fecal Egg Count (FEC) data from __________________
  #These data are still preliminary and are unpublished. Accessed June 5 2024.

#--------------------------------------------------

#load packages
library(here) #v1.0.1
library(tidyverse) #v2.0.0
library(lme4) #v1.1-35.1
library(lmerTest) #v3.1-3 #like lme4 but with pvalues 
library(janitor) #v2.2.0

# #clear environment
rm(list = ls())


############################## LOAD DATA #######################################

#metadata for FEC data
  #this metadata file is unpublished
#load FULL 2021 and 2022 data - INCLUDING May captures, INCLUDING sex/repro=NA
fulltrap21_ALL <- readRDS(file="fulltrap21_ALL_03.04.24.rds")
fulltrap22_ALL <- readRDS(file="fulltrap22_ALL_03.04.24.rds")
fulltrap21.22_ALL <- rbind(fulltrap21_ALL, fulltrap22_ALL)
#one entry per sample ID
metadata <- fulltrap21.22_ALL %>% 
  group_by(samp_id) %>% slice(1) %>% #remove duplicate entries per occasion (ie week recaps) %>%
  dplyr::select(c("year", "site", "trt", "food_trt", "helm_trt",
                  "season", "month", "occasion",
                  "tag", "samp_id", "sex", "season_breeder"))

#preliminary FEC data from _____ v.June 5 2024
  #this data file is unpublished
FECdata <- read_csv(here("FEC data for Janine 6-5-24.csv"), col_names=TRUE, na = c("", "NA", "#N/A")) %>% 
  clean_names()


# #check for voles processed 2x in the same occasion
# processed_twice <- fulltrap21.22_ALL %>% group_by(tag, year, occasion) %>%
#   mutate(dup = n_distinct(samp_id)) %>% filter(dup > 1) %>%
#   select(year, month, occasion, session, tag, samp_id, firstcap, trap) %>%
#   arrange(tag) %>% ungroup()
# #short version of FEC data
# FECshort <- FECdata %>% select(year, sample_id, nematode_epg, nematode_y_n)
# #join FEC data to animals processed twice to decide which sample ID's data to keep
# processed_twiceFEC <- left_join(processed_twice, FECshort, by=c("year", "samp_id"="sample_id"))
# write.csv(processed_twiceFEC, file="processed_twiceFEC.csv")


#pull 2021/22 FEC data, combine with metadata, keep only complete entries
#helm_trt, pre_post, and nematode.y.n are all factors
all21.22FEC <- FECdata %>%
  filter(is.na(exclude_from_fec_analysis)) %>% #remove samples with feces weight less than 0.01, FEC count may be inaccurate
  rename("samp_id" = "sample_id",
         "nematode.epg" = "nematode_epg",
         "nematode.y.n" = "nematode_y_n") %>%
  select(c("samp_id", "nematode.epg", "nematode.y.n")) %>%
  inner_join(metadata, by="samp_id") %>%
  drop_na(nematode.epg) %>% drop_na(nematode.y.n) %>% #remove entries missing FEC data
  mutate(nematode.epg = as.numeric(nematode.epg),
       nematode.y.n = as.factor(nematode.y.n)) %>%
  #remove three duplicate entries (details to follow)
  filter(samp_id != 274) %>% #vole 226170 sampled twice in July 21, remove second sample ID (eggs detected in both samples)
  filter(samp_id != 340) %>% #vole 219980 sampled twice in Aug 21, remove first sample ID (no eggs in first sample, eggs in second)
  filter(samp_id != 712) %>% #vole 219809 was sampled twice in Sept 21, remove second sample ID (eggs detected in both samples)
  filter(samp_id != 1433) %>% #vole 219322 sampled twice in July 22, remove second sample ID (eggs detected in first sample, no eggs in second)
  filter(samp_id != 1999) %>% #vole 219371 sampled twice in Sept 22, remove first sample (count=NA first sample, no eggs in second)
  filter(samp_id != 1612) %>% #vole 219422 sampled twice in Aug 22, remove first sample (no eggs detected in first sample, eggs in second)
  filter(samp_id != 1668) %>% #vole 219525 sampled twice in Aug 22, remove second sample (no eggs in either sample)
  filter(samp_id != 2133) %>% #vole 226624 sampled twice in Oct 22, remove second sample (eggs detected in first sample, no eggs in second)
  filter(samp_id != 996) %>% #vole 226770 sampled twice in Oct 22, remove second sample (no eggs detected in first sample, count=NA in second)
  filter(samp_id != 1562) %>% #vole 227032 sampled twice in Aug 22, remove second sample (no eggs detected in either sample)
  filter(samp_id != 2186) %>% #vole 227049/219527 sampled twice in Oct 22, remove second sample (no eggs detected in either sample)
  filter(samp_id != 1565) %>% #vole 227166 sampled twice in Aug 22, remove second sample (no eggs detected in either sample)
  filter(samp_id != 1745) %>% #vole 227207/219468 sampled twice in Sept 22, remove second sample (eggs in first sample, no eggs in second sample)
  group_by(year, tag) %>% arrange(occasion, .by_group = TRUE) %>%
  mutate(capt_nbr = row_number()) %>% #add column for the capture number (ie the first, second, third capture of given vole - IN THAT YEAR)
  mutate(pre_post = as.factor(case_when(capt_nbr == 1 ~ "pre",
                                        capt_nbr > 1 ~ "post"))) %>%
  # drop_na(sex) %>% drop_na(season_breeder) %>% #model has 17ish voles with sex or repro=NA
  ungroup()

#separate 2021 data
all2021FEC <- all21.22FEC %>% filter(year=="2021")
#separate 2022 data
all2022FEC <- all21.22FEC %>% filter(year=="2022")

################################################################################################
############################ RUN PREVALENCE & INTENSITY MODELS #################################
################################################################################################
############## here we're running two models, prevalence and intensity #########################
########### models have both pre- and post- measurements, interaction term #####################
################################################################################################

#this code runs on the 'all2021FEC' / 'all2022FEC' dataset (models combine pre- and post- data with an 
  #interaction by helm_trt), using 'helm_trt' and 'pre_post' (both FACTORS) as predictors

########### 2021 MODELS ###########

#p = prevalence model
#if needed: adjust the optimizer: https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
p.prepost <- lme4::glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
                   data = all2021FEC, family = binomial) #n=1010
summary(p.prepost) #trt*prepost b=-0.787 (log odds) p=0.008**
exp(coef(summary(p.prepost))[,1]) #trt*prepost b=0.455 (odds ratio)

# ###### GLMM model diagnostics ######
# #https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
# library(DHARMa) #v0.4.6
# #calculate residuals (then run diagnostics on these)
# simulationOutput <- simulateResiduals(fittedModel = p.prepost)
# #qq plot and residual vs fitted
# plot(simulationOutput) #qq look good, scaled residuals look good
# #formal test for overdispersion
# testDispersion(simulationOutput) #looks good, pvalue is large
# #formal test for zero inflation (common type of overdispersion)
# testZeroInflation(simulationOutput) #looks good, pvalue is large
# plotResiduals(simulationOutput, all2021FEC$helm_trt) #not working?
# plotResiduals(simulationOutput, all2021FEC$pre_post) #looks good
# #### OVERALL: DHARMa looks good #####

#------------------------------------------------

#i = intensity (only samples with nematode.epg >0)
i.prepost <- lmerTest::lmer(log(nematode.epg) ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
                  data = subset(all2021FEC, nematode.epg > 0)) #n=434
summary(i.prepost) #trt*prepost b=-0.74 p=0.0039**

# #model diagnostics
# plot(i.prepost) #some patterning
# shapiro.test(resid(i.prepost)) #p=0.004676 some issues of non-normality
# 
# ######## intensity LMM model diagnostics ############
# #LMM model diagnostics > https://goodekat.github.io/redres/
# # devtools::install_github("goodekat/redres")
# library(redres) #v0.0.0.9
# # creates a plot of the conditional studentized residuals versus the fitted values
# plot_redres(i.prepost, type="std_cond")
# plot_redres(i.prepost, type="pearson_cond") #Pearson conditional residuals
# # creates a residual quantile plot for the error term
# plot_resqq(i.prepost) #looks okay
# # creates normal quantile plots for each random effect
# plot_ranef(i.prepost) #looks okay
# ###### OVERALL: model diagnostics look good

#------------------------------------------------

# #summary table of mixed model output
# #https://cran.r-project.org/web/packages/sjPlot/vignettes/tab_mixed.html
# library(sjPlot)
# #save it https://stackoverflow.com/questions/67280933/how-to-save-the-output-of-tab-model
# tab_model(p.prepost, file="SUPP_Table_S3.doc")
# tab_model(i.prepost, file="SUPP_Table S4.doc")



########### 2022 MODELS ###########

#p = prevalence model
#if needed: adjust the optimizer: https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
p.prepost <- lme4::glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
                         data = all2022FEC, family = binomial) #n=1007
summary(p.prepost) #trt*prepost b=-0.703 (log odds) p=0.0468*
exp(coef(summary(p.prepost))[,1]) #trt*prepost b=0.495 (odds ratio)

# ###### GLMM model diagnostics ######
# #https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
# library(DHARMa)
# #calculate residuals (then run diagnostics on these)
# simulationOutput <- simulateResiduals(fittedModel = p.prepost)
# #qq plot and residual vs fitted
# plot(simulationOutput) #qq look good, scaled residuals look good
# #formal test for overdispersion
# testDispersion(simulationOutput) #looks good, pvalue is large
# #formal test for zero inflation (common type of overdispersion)
# testZeroInflation(simulationOutput) #looks good, pvalue is large
# plotResiduals(simulationOutput, all2022FEC$helm_trt) #doesn't work... :/
# plotResiduals(simulationOutput, all2022FEC$pre_post) #looks good
# #### OVERALL: DHARMa looks good #####

#------------------------------------------------

#i = intensity (only samples with nematode.epg >0)
i.prepost <- lmerTest::lmer(log(nematode.epg) ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
                            data = subset(all2022FEC, nematode.epg > 0)) #n=262
summary(i.prepost) #trt*prepost b=-0.036 p=0.915

# #model diagnostics
# plot(i.prepost) #some patterning
# shapiro.test(resid(i.prepost)) #p=7.75e-5 very small, issues of non-normality
# 
# ######## intensity LMM model diagnostics ############
# #LMM model diagnostics > https://goodekat.github.io/redres/
# # devtools::install_github("goodekat/redres")
# library(redres)
# # creates a plot of the conditional studentized residuals versus the fitted values
# plot_redres(i.prepost, type="std_cond")
# plot_redres(i.prepost, type="pearson_cond") #Pearson conditional residuals
# # creates a residual quantile plot for the error term
# plot_resqq(i.prepost) #looks okay, deviates some on the right
# # creates normal quantile plots for each random effect
# plot_ranef(i.prepost) #similar pattern as above
# ###### OVERALL: model diagnostics look okay, not great but not overly alarming

#------------------------------------------------



#################################### VISUALIZE MODEL OUTPUT ############################################

## visualize the interaction effects
## these figures do not appear in the manuscript or the supplement

library(visreg) #v2.7.0 #https://pbreheny.github.io/visreg/articles/web/overlay.html
library(cowplot) #v1.1.2
library(ggtext) #0.1.2


###------- 2021 DATA ------------

### PREVALENCE MODEL ###

#change pre_post to numeric to get a regression line in visreg output
all2021FEC.num <- all2021FEC %>%
  mutate(pre_post = as.numeric(pre_post)) %>%
  mutate(nematode.y.n_num = case_when(nematode.y.n == "0" ~ 0,
                                      nematode.y.n == "1" ~ 1)) #nematode.y.n from factor to numeric to plot points on continous y axis
p.prepost.num <- lme4::glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
                             data = all2021FEC.num, family = binomial)

#visualize infection probability, pretty
visreg(p.prepost.num, "pre_post", by="helm_trt", scale='response', rug=FALSE,
       gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.03,0.05))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(breaks=c(1,2),
                     expand = expansion(mult=c(0.02,0.05))) +
  scale_color_manual(values=c("#808080", "#ffe048")) +
  scale_fill_manual(values=c("#80808070", "#ffe04850")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size=11),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt")) +
  labs(x="Treatment Stage", y="Likelihood of Helminth Infection") +
  annotate(geom = "text", x=1.35, y=1.05, size = 4,
           label = paste("Trt Grp*Trt Stage")) +
  annotate(geom = "text", x=1.35, y=0.97, size = 4,
           label = paste("OR =",
                         round( exp(coef(summary(p.prepost.num))[4,1]), digits=3) )) +
  annotate(geom = "text", x=1.35, y=.90, size = 4,
           label = paste("p =",
                         round( coef(summary(p.prepost.num))[4,4], digits=3) )) +
  geom_jitter(aes(x=pre_post, y=nematode.y.n_num, color=helm_trt),
              data=all2021FEC.num,
              width=0.1, height=0.05,
              size=3, alpha=0.2, shape=16)
  
### INTENSITY MODEL ###
  
#pre_post as numeric to look like a regression line
#visreg doesn't like it when the model has a subset command in it #https://github.com/pbreheny/visreg/issues/99
witheggs2021FEC.num <- subset(all2021FEC.num, nematode.epg > 0)
i.prepost.num <- lmerTest::lmer(log(nematode.epg) ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
                                data = witheggs2021FEC.num)

#visualize infection intensity, pretty
visreg(i.prepost.num, "pre_post", by="helm_trt", scale='response', rug=FALSE,
       gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.03,0.05))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(breaks=c(1,2),
                     expand = expansion(mult=c(0.02,0.05))) +
  scale_color_manual(values=c("#808080", "#ffe048")) +
  scale_fill_manual(values=c("#80808070", "#ffe04850")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt")) +
  labs(x="Treatment Stage", y="ln(eggs per gram of feces)") +
  annotate(geom = "text", x=1.35, y=8.6, size = 4,
           label = paste("Trt Grp*Trt Stage")) +
  annotate(geom = "text", x=1.35, y=8.0, size = 4,
           label = paste("\u03B2 =",
                         round( coef(summary(i.prepost.num))[4,1], digits=3) )) +
  annotate(geom = "text", x=1.35, y=7.5, size = 4,
           label = paste("p =",
                         round( coef(summary(i.prepost.num))[4,5], digits=3) )) +
  geom_jitter(data=witheggs2021FEC.num, aes(x=pre_post, y=log(nematode.epg), color=helm_trt),
              width=0.05, size=3, alpha=0.2, shape=16)




###-------- 2022 DATA ------------

### PREVALENCE MODEL ###

#change pre_post to numeric to get a regression line in visreg output
all2022FEC.num <- all2022FEC %>%
  mutate(pre_post = as.numeric(pre_post)) %>%
  mutate(nematode.y.n_num = case_when(nematode.y.n == "0" ~ 0,
                                      nematode.y.n == "1" ~ 1)) #nematode.y.n from factor to numeric to plot points on continous y axis
p.prepost.num <- lme4::glmer(nematode.y.n ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
                           data = all2022FEC.num, family = binomial)

#visualize infection probability, pretty
visreg(p.prepost.num, "pre_post", by="helm_trt", scale='response', rug=FALSE,
                gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.03,0.05))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(breaks=c(1,2),
                     expand = expansion(mult=c(0.02,0.05))) +
  scale_color_manual(values=c("#808080", "#ffe048")) +
  scale_fill_manual(values=c("#80808070", "#ffe04850")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size=11),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt")) +
  labs(x="Treatment Stage", y="Likelihood of Helminth Infection") +
  annotate(geom = "text", x=1.35, y=1.05, size = 4,
           label = paste("Trt Grp*Trt Stage")) +
  annotate(geom = "text", x=1.35, y=0.97, size = 4,
           label = paste("OR =",
                         round( exp(coef(summary(p.prepost.num))[4,1]), digits=3) )) +
  annotate(geom = "text", x=1.35, y=.90, size = 4,
           label = paste("p =",
                         round( coef(summary(p.prepost.num))[4,4], digits=3) )) +
  geom_jitter(aes(x=pre_post, y=nematode.y.n_num, color=helm_trt),
              data=all2022FEC.num,
              width=0.1, height=0.05,
              size=3, alpha=0.2, shape=16)


### INTENSITY MODEL ###

#pre_post as numeric to look like a regression line
#visreg doesn't like it when the model has a subset command in it #https://github.com/pbreheny/visreg/issues/99
witheggs2022FEC.num <- subset(all2022FEC.num, nematode.epg > 0)
i.prepost.num <- lmerTest::lmer(log(nematode.epg) ~ helm_trt + pre_post + helm_trt:pre_post + (1|tag),
                  data = witheggs2022FEC.num)

#visualize infection intensity, pretty
visreg(i.prepost.num, "pre_post", by="helm_trt", scale='response', rug=FALSE,
       gg=TRUE, overlay=TRUE) +
  scale_y_continuous(expand = expansion(mult=c(0.03,0.05))) + #controls extra white space on axes (cowplot vignette)
  scale_x_continuous(breaks=c(1,2),
                     expand = expansion(mult=c(0.02,0.05))) +
  scale_color_manual(values=c("#808080", "#ffe048")) +
  scale_fill_manual(values=c("#80808070", "#ffe04850")) +
  theme_half_open() +
  theme(legend.position = "bottom",
        axis.title = element_text(size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size=12),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 20, r = 10, b = 10, l = 20, unit = "pt")) +
  labs(x="Treatment Stage", y="ln(eggs per gram of feces)") +
  annotate(geom = "text", x=1.35, y=8.6, size = 4,
           label = paste("Trt Grp*Trt Stage")) +
  annotate(geom = "text", x=1.35, y=8.0, size = 4,
           label = paste("\u03B2 =",
                         round( coef(summary(i.prepost.num))[4,1], digits=3) )) +
  annotate(geom = "text", x=1.35, y=7.5, size = 4,
           label = paste("p =",
                         round( coef(summary(i.prepost.num))[4,5], digits=3) )) +
  geom_jitter(data=witheggs2022FEC.num, aes(x=pre_post, y=log(nematode.epg), color=helm_trt),
              width=0.05, size=3, alpha=0.2, shape=16)
