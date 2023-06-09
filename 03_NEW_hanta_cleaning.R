### AS OF JUNE 8 2023 - this code is running using the "..._STSB" versions of netmets - has breeders/nonbreeders ###

# code run using R version 4.1.2 (2021-11-01) -- "Bird Hippie"

#load libraries
library(here)
library(tidyverse)
library(janitor)
library(lubridate)
library(lme4)
# library(coin) #permutation wilcox test

#clear environment
rm(list = ls())


#### all the steps to generate 'netmets_puuv' have been commented out as of 06.09.23 (post sex.deg debug) - readRDS to load file

# ##########---------- LOAD NETWORK METRICS and VOLE CAPTURE METADATA -----------###########
# 
# # ### Fulltrap and netmets dfs created separately, load 2021 and 2022 data here, clean as needed
# # 
# # ### INDIVIDUAL VOLE METADATA by year
# #     # >> FROM OTHER R PROJECT! << (hence why I'm not using 'here()')
# # fulltrap21 <- readRDS(file="../vole-spatial-v2/fulltrap21_05.10.23.rds") #go up a level from current wd, then down to file
# # fulltrap22 <- readRDS(file="../vole-spatial-v2/fulltrap22_05.10.23.rds")
# # 
# # ###### THESE VERSIONS OF FULLTRAP have all the animals (breeders and nonbreeders)
# # 
# # ## >>NOTE<< overwintered animals may incorrectly have firstcap==1 in 2022
# #     # PIT tag: 21895 Occ2Sess1 at Vaarinkorpi
# #     # PIT tag: 226280 Occ2Sess2 at Kuoppa
# # # FULLtrap <- rbind(fulltrap21, fulltrap22)
# # # OW <- FULLtrap %>%
# # #   group_by(tag) %>% arrange(year, occ.sess, .by_group = TRUE) %>%
# # #   filter(n_distinct(year) > 1)
# # # write.csv(OW, here("overwinter21-22.csv"))
# # 
# # #join 21 and 22 fulltrap data; correct firstcap for 2022
# # fulltrap21.22 <- rbind(fulltrap21, fulltrap22) %>%
# #   group_by(tag) %>% arrange(tag, year, occ.sess, .by_group = TRUE) %>%
# #   mutate(firstcap = ifelse(date_time == min(date_time), 1, 0)) %>%
# #   mutate(firstcap = factor(firstcap)) %>%
# # #remove may; keep one entry per tag,year,month; pull only relevant columns
# #   filter(month!="may") %>%
# #   ungroup() %>% group_by(tag, year, month) %>%
# #   slice(1) %>%
# #   select(year, site, trt, month, tag, samp_id, sex, season_breeder, traps_per_life, caps_per_life)
# # 
# # ### NETWORK METRICS - STSB version - breeders and nonbreeders!
# # netmets21 <- readRDS(file="../vole-spatial-v2/netmets21_STSB.rds") %>%
# #   mutate(year=as.numeric(2021))
# # netmets22 <- readRDS(file="../vole-spatial-v2/netmets22_STSB.rds") %>%
# #   mutate(year=as.numeric(2022))
# # 
# # netmets21.22 <- rbind(netmets21, netmets22)
# # 
# # ### NETWORK METRICS + METADATA
# # 
# # netmets_full <- left_join(netmets21.22, fulltrap21.22, by=c("year", "site", "month", "tag")) %>%
# #   mutate(site=as.factor(site),
# #          year=as.factor(year),
# #          month = as.factor(month),
# #          month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #remove may from factor levels
# #   select(!c(focal_id)) %>% #remove duplicate column
# #   relocate(c(year, trt, site, month, n.node, tag, samp_id, sex), .before = wt.deg)
# # 
# # saveRDS(netmets_full, here("netmets_full_06.09.23.rds"))
# 
# #load the most recent version of netmets_full
# netmets_full <- readRDS(here("netmets_full_06.09.23.rds"))
# 
# 
# 
# 
# # #visualize degree by month for trt and year
# # netmets_full %>%
# #   ggplot(aes(x=month, y=strength, fill=trt)) +
# #   geom_violin() +
# #   facet_wrap(~year, nrow=2)
# 
# 
# #########################################   LOAD & CLEAN PUUV IFA DATA   ########################################
# 
# # #load, clean, format PUUV IFA data
# # ### NEW! Updated! PUUV_IFA data with the goofy samples from 2021/2022 that KWearing re-ran in May 2023 (code updated 6.7.23)
# # puuv_data <- read.csv(here("puuv_ifa_06.07.23.csv")) %>%
# #   clean_names %>%
# #   #populate column of FINAL PUUV status
# #   #kind of a pain now, since samples could be run 1-4x but we want the result of the last run as the 'final' status
# #   #columns are named as 'puuv_run1' 'puuv_run2' 'puuv_run3' 'puuv_run4'
# #   mutate(FINAL_puuv = ifelse(!is.na(puuv_run4), as.character(puuv_run4),
# #                              ifelse(!is.na(puuv_run3), as.character(puuv_run3),
# #                                     ifelse(!is.na(puuv_run2), as.character(puuv_run2), as.character(puuv_run1))))) %>%
# #   mutate(samp_id = as.numeric(id),
# #          # date_run1 = as_date(date_run1, format= "%m/%d/%Y"),
# #          # puuv_run1 = as.factor(puuv_run1),
# #          # date_run2 = as_date(date_run2, format= "%m/%d/%Y"),
# #          # puuv_run2 = as.factor(puuv_run2),
# #          # date_run3 = as_date(date_run3, format= "%m/%d/%Y"),
# #          # puuv_run3 = as.factor(puuv_run3),
# #          # date_run4 = as_date(date_run4, format= "%m/%d/%Y"),
# #          # puuv_run4 = as.factor(puuv_run4),
# #          FINAL_puuv = as.factor(FINAL_puuv)) %>%
# #   drop_na(FINAL_puuv) %>%
# #   dplyr::select(FINAL_puuv, samp_id) %>%
# #   rename(puuv_ifa = FINAL_puuv)
# # #output is a df with two columns, sample ID and PUUV status (0,1)
# # 
# # # Save puuv_data to a rdata file
# # saveRDS(puuv_data, file = here("hantadata_06.09.23.rds"))
# 
# # Restore puuv_data from the rdata file
# puuv_data <- readRDS(file = "hantadata_06.09.23.rds")
# 
# 
# 
# 
# ##################################  COMBINE NETMETS data with PUUV data ######################################
# 
# #join puuv_data to netmets_full (BOTH YEARS!)
# netmets_puuv <- left_join(netmets_full, puuv_data, by="samp_id") %>%
#   #left join on netmets_full because I need to have network data
#   relocate(puuv_ifa, .after="samp_id") %>%
#   select(!samp_id) %>%
#   drop_na(puuv_ifa) #drop any animals without puuv data
# 
# ########## TO DECIDE: KEEP THE ANIMALS WITH NO PUUV STATUS? they're only useful to plot, can't use for models
# ## there are 121 animals that we have network data for but DO NOT know their PUUV status
# 
# 
# ############ REMOVE THE VOLES from netmets_puuv THAT SEROCONVERT POS TO NEG ##################
# 
# #voles that seroconvert PUUV+ to PUUV-
# puuv_pos_neg <- netmets_puuv %>% group_by(tag) %>% arrange(year, month, .by_group = TRUE) %>%
#   dplyr::select(year, month, tag, puuv_ifa) %>%
#   summarise(status_time = toString(puuv_ifa)) %>%
#   filter(str_detect(status_time, "1,\\s0")) #filter for animals that go from pos to neg
# #pull the PIT tags (21 individuals, mostly 1-0, a few 1-0-1)
# problemchildren <- puuv_pos_neg$tag
# #filter netmetsPUUV to remove 'problemchildren' - status is inconclusive or maybe we detected MatAb
# netmets_puuv <- netmets_puuv %>%
#   filter(!tag %in% problemchildren)
# 
# ## NOW NETMETS_PUUV has ONLY the animals that:  ##
#     #       1. Have network data
#     #       2. Have PUUV IFA data
#     #       3. DON'T seroconvert from PUUV+ to PUUV-
# ## NETMETS_PUUV has BREEDERS and NONBREEDERS as of June 8 2023
# 
# ######################### end problem children ################################
# 
# 
# 
# 
# 
# ###############################################################################
# ### exploratoriness of breeders ###
# ## based off of VanderWaal ground squirrel ms ##
# #https://modelr.tidyverse.org/reference/add_predictions.html
# 
# #subset data to one entry per vole per year
# onepertag <- netmets_puuv %>%
#   group_by(year, tag) %>% slice(1)
# 
# #fit power model to ln(traps) ~ ln(caps)
# powmod <- lm(log(traps_per_life) ~ log(caps_per_life), data=onepertag)
# # summary(powmod) #check it
# # # ln(y) = -0.01564 + 0.73605ln(x)
# # # y = 0.9844817x^0.73605
# 
# # #yes, powmod is much better than linreg
# # linmod <- lm(traps_per_life ~ caps_per_life, data=onepertag)
# # summary(linmod)
# #
# # AIC(linmod, powmod)
# 
# #estimate y-hat (predicted values)
# predicted.values <- predict(powmod, type="response")
# 
# #visualize it - colors are observed, dotted line is predicted
# ggplot(aes(x=log(caps_per_life), y=log(traps_per_life), color=sex), data=onepertag) +
#   geom_jitter(alpha=0.5, size=1.5) +
#   geom_line(aes(y=predicted.values), linetype=2, color="black") +
#   labs(title="exploratoriness")
# 
# #dataframe of observed values
# d <- data.frame(onepertag$tag, onepertag$caps_per_life, onepertag$traps_per_life) %>%
#   rename(caps_per_life = onepertag.caps_per_life,
#          traps_per_life = onepertag.traps_per_life,
#          tag = onepertag.tag)
# #add predicted values from power model
# d <- d %>% modelr::add_predictions(powmod) %>%
#   mutate(pred_exp = exp(pred))
# #calculate residuals (both log() and not-ln versions)
# onepertag_pred <- onepertag %>% left_join(d, by=c("tag", "caps_per_life", "traps_per_life")) %>%
#   mutate(explore = log(traps_per_life) - pred,
#          explore_exp = traps_per_life - pred_exp) #residual of obs traps_per_life - expected
# 
# ### EXPLORE is in terms of log(traps_per_life)  [ maybe this is scaled -1 to 1 because it's a ln? ]
# ### EXPLORE_EXP is in terms of traps_per_life (so how many more/fewer traps were you in than expected)
# 
# # #visualize
# # onepertag_pred %>%
# #   ggplot(aes(x=sex, y=explore, fill=sex)) +
# #   geom_violin() +
# #   geom_hline(yintercept = 0, color="black")
# 
# 
# ############ WHAT UNITS SHOULD THE RESIDUALS BE IN ?? ########################
# 
# exploratory <- onepertag_pred %>% select(c("year", "tag", "explore"))
# 
# netmets_puuv <- netmets_puuv %>% left_join(exploratory, by=c("year", "tag"))
# 
# ############################################################################
# 
# 
# 
# 
# 
# # add previous (lagged) degree (degree from previous month influences current PUUV status)
# # add 0,1 for serovert - animals that go 0-0 or 0-1
#     # BUT! the previous month has to be in the same year (don't want 2021 fall to influence 2022 spring)
# netmets_puuv <- netmets_puuv %>% group_by(year, tag) %>%
#   arrange(month, .by_group = TRUE) %>%
#   mutate(prev_wt.deg = lag(wt.deg, n=1),
#          prev_F.deg = lag(F.deg, n=1),
#          prev_M.deg = lag(M.deg, n=1),
#          prev_b.deg = lag(b.deg, n=1),
#          prev_nb.deg = lag(nb.deg, n=1),
#          prev_n.node = lag(n.node, n=1),
#          prev_month = lag(month, n=1)) %>%
#   mutate(prev_curr_puuv = paste(lag(puuv_ifa), puuv_ifa, sep="-")) %>%
#   mutate(serovert = as.factor(case_when(prev_curr_puuv == "0-0" ~ 0,
#                                         prev_curr_puuv == "0-1" ~ 1,
#                                         prev_curr_puuv == "1-1" ~ NA)))
#   # %>% select(!prev_curr_puuv)
# 
# 
# #save all the above code to here
# saveRDS(netmets_puuv, here("netmets_puuv_06.09.23.rds"))

#load netmets_puuv
netmets_puuv <- readRDS(here("netmets_puuv_06.09.23.rds"))

##################################################################################################################











# ############### PUUV PREVALENCE PER SITE #########################
# 
# #what is the prevalence of hanta at each site in each month?
# puuv_prev <- netmets_puuv %>% group_by(year, site, month) %>%
#   summarise(n = length(tag),
#             n_pos = length(which(puuv_ifa == "1")),
#             prev = n_pos/n) %>%
#   mutate(site = factor(site, levels=c("asema", "helmipollo", "hevonen",
#                                       "ketunpesa", "kiirastuli", "mustikka",
#                                       "kuoppa", "radio", "vaarinkorpi",
#                                       "janakkala", "luostari", "puro", "talo")))
# 
# # #plot prevalence each month, looking for sites with multiple prev=0 in a row
# #     #--> these should probably be removed for analysis
# # puuv_prev %>%
# #   ggplot() +
# #   geom_point(aes(x=month, y=prev, color=prev==0)) +
# #   scale_color_manual(name="prevalence = 0",
# #                      values=setNames(c("red","black"),c(T,F)))+
# #   facet_grid(vars(year), vars(site)) +
# #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# ### ERROR BARS based on sample size
# ## essentially something to show if it's 100% prev of 2 animals or 45 animals
# #https://www.rdocumentation.org/packages/epiR/versions/0.9-79/topics/epi.conf
# library(epiR)
# #Method prevalence require a two-column matrix; the first column specifies the number of positives, 
#     #the second column specifies the total number tested. 
# epi.conf.data <- puuv_prev %>% ungroup() %>% select(n_pos, n) %>% as.matrix()
# prevCI <- epi.conf(epi.conf.data, ctype="prevalence", method="exact", conf.level = 0.95)
# puuv_prev <- cbind(puuv_prev, prevCI)
# 
# #plot prevalence each month with error bars
# puuv_prev %>%
#   ggplot(aes(x=month, y=prev, color=prev==0)) +
#   geom_point() +
#   geom_errorbar(aes(y=est, ymin=lower, ymax=upper), width=.2, alpha=0.3) +
#   scale_color_manual(name="prevalence = 0",
#                      values=setNames(c("red","black"),c(T,F))) +
#   facet_grid(vars(year), vars(site)) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# 
# ## visualize infected as a percentage of the total population
# puuv_prev_stack <- netmets_puuv %>% 
#   mutate(month.n = case_when(month=="june" ~ 6,
#                            month=="july" ~ 7,
#                            month=="aug" ~ 8,
#                            month=="sept" ~ 9,
#                            month=="oct" ~ 10),
#          month.n = as.numeric(month.n)) %>%
#   mutate(site = factor(site, levels=c("asema", "helmipollo", "hevonen",
#                                       "ketunpesa", "kiirastuli", "mustikka",
#                                       "kuoppa", "radio", "vaarinkorpi",
#                                       "janakkala", "luostari", "puro", "talo"))) %>%
#   group_by(year, site, month.n) %>%
#   summarise(n = length(tag),
#             n_pos = length(which(puuv_ifa == "1"))) %>%
#   pivot_longer(-c("year", "month.n", "site"), names_to = "group", values_to = "count")
# 
# puuv_prev_stack %>%
#   ggplot(aes(x=month.n, y=count, fill=group)) +
#   geom_area() +
#   facet_grid(vars(year), vars(site))
# 
# ############################## end prevalence ####################################





# #visualize M.deg and F.deg for males and females in each trt
# 
# netmets_puuv %>%
#   ggplot(aes(x=trt, y=F.deg, color=sex)) +
#   geom_violin() +
#   facet_grid(year ~ month)
# 
# netmets_puuv %>%
#   ggplot(aes(x=trt, y=M.deg, fill=sex)) +
#   geom_boxplot() +
#   facet_grid(year ~ month)
# 
# netmets_puuv %>%
#   ggplot(aes(x=trt, y=wt.deg, fill=sex)) +
#   geom_boxplot() +
#   facet_grid(year ~ month)
# 
# ## COMPARE female and male degree of a given animal
# netmets_puuv %>%
#   ggplot(aes(x=F.deg, y=M.deg, color=trt, )) +
#   geom_point(aes(shape=sex)) +
#   scale_shape_manual(values=c(3, 16)) +
#   coord_fixed() +
#   xlim(0,4) + ylim(0,4) +
#   facet_grid(year ~ month)






#########################################################################

#How does CURRENT degree affect CURRENT infection status?

dat <- netmets_puuv

#visualize infection status by weighted deg by sex,trt,year
dat %>%
  ggplot(aes(x=wt.deg, y=puuv_ifa, color=sex, group=sex)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(trt~year, ncol=2)

#infected males have higher deg than infected? females similar effect size, not as significant
mod <- glmer(puuv_ifa ~ wt.deg:sex + sex + season_breeder + explore + trt + month + n.node + year + (1|site),
             family=binomial, data=dat)
summary(mod)

#infected animals have higher degree with same sex? (slightly stronger effect for females, but only marginally signif)
mod <- glmer(puuv_ifa ~ F.deg:sex + M.deg:sex + sex + season_breeder + trt + month + n.node + year + (1|site),
             family=binomial, data=dat)
summary(mod)

#WARNING: Failed to converge
    ## BEN BOLKER says it's okay, the bobyqa ("Bound Optimization BY Quadratic Approximation") isn't actually doing anything
    ## FALSE POSITIVE 'failed to converge'

#fixed with : adding "control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))"
#https://stats.stackexchange.com/questions/164457/r-glmer-warnings-model-fails-to-converge-model-is-nearly-unidentifiable
    #also: https://www.learn-mlms.com/07-module-7.html#learning-objectives-5


#RESULTS - Current infection does seem to make animals have higher degree (more apparent for males than females)

#########################################################################




#########################################################################################################
########### HOW does previous month's network position affect (likelihood of?) seroconversion ###########

## last month's degree, this month's PUUV status
## grouped by year and treatment (so replicate sites are combined)
netmets_puuv_serov <- netmets_puuv %>%
  # unite(yr_trt, year, trt, remove=FALSE) %>%
  # ungroup() %>% group_by(yr_trt) %>%
  drop_na(serovert)

# #summarize number of 0-0 vs number of 0-1 by year_trt
# netmets_puuv_serov %>% group_by(yr_trt) %>%
#   summarise(n_0 = sum(serovert=="0"),
#             n_1 = sum(serovert=="1"))

## RESULT: There are A LOT MORE animals that don't seroconvert, than those that do
## can the model / analysis even be fair/reasonable if we're comparing two things that are so off-weight?


# library(GGally)
# 
# data <- netmets_puuv_serov %>% ungroup() %>% select(!c(avg.wt.deg, wt.deg,
#                                                        puuv_ifa, tag))
# ggpairs(data)


# #previous wt.deg of males/females that serovert vs those that don't
# library(gridExtra)
# 
# wt <- netmets_puuv_serov %>%
#   ggplot(aes(x=sex, y=prev_wt.deg, color=serovert)) +
#   geom_violin() +
#   facet_grid(trt ~ year)
# 
# female <- netmets_puuv_serov %>%
#   ggplot(aes(x=sex, y=prev_F.deg, fill=serovert)) +
#   geom_violin() +
#   facet_grid(trt ~ year)
# 
# male <- netmets_puuv_serov %>%
#   ggplot(aes(x=sex, y=prev_M.deg, fill=serovert)) +
#   geom_violin() +
#   facet_grid(trt ~ year)
# 
# grid.arrange(female, male, ncol=2)


library(kableExtra)

#number pos and neg by year=2021, treatment, sex
netmets_puuv %>% drop_na(prev_wt.deg) %>% 
  filter(year=="2021") %>%
  mutate(trt = case_when(trt=="unfed_control" ~ "Unfed-Control",
                         trt=="unfed_deworm" ~ "Unfed-Deworm",
                         trt=="fed_control" ~"Fed-Control",
                         trt=="fed_deworm" ~ "Fed-Deworm"),
         trt = factor(trt, levels=c("Unfed-Control", "Unfed-Deworm",
                                    "Fed-Control", "Fed-Deworm"))) %>%
  group_by(year, trt, sex, puuv_ifa) %>%
  summarise(n = length(tag)) %>%
  pivot_wider(id_cols = c(year, trt, sex), values_from = n, names_from = puuv_ifa) %>%
  kbl(col.names = c("Year", "Treatment", "Sex", "N Negative", "N Positive"), align = "lllcc") %>%
  kable_styling(full_width=FALSE,
                bootstrap_options = "striped") %>%
  row_spec(0, bold=T, color="black", background="#DAD7D7")

#number pos and neg by year=2022, treatment, sex
netmets_puuv %>% drop_na(prev_wt.deg) %>% 
  filter(year=="2022") %>%
  mutate(trt = case_when(trt=="unfed_control" ~ "Unfed-Control",
                         trt=="unfed_deworm" ~ "Unfed-Deworm",
                         trt=="fed_control" ~"Fed-Control",
                         trt=="fed_deworm" ~ "Fed-Deworm"),
         trt = factor(trt, levels=c("Unfed-Control", "Unfed-Deworm",
                                    "Fed-Control", "Fed-Deworm"))) %>%
  group_by(year, trt, sex, puuv_ifa) %>%
  summarise(n = length(tag)) %>%
  pivot_wider(id_cols = c(year, trt, sex), values_from = n, names_from = puuv_ifa) %>%
  kbl(col.names = c("Year", "Treatment", "Sex", "N Negative", "N Positive"), align = "lllcc") %>%
  kable_styling(full_width=FALSE,
                bootstrap_options = "striped") %>%
  row_spec(0, bold=T, color="black", background="#DAD7D7")


data <- netmets_puuv %>% drop_na(prev_wt.deg) %>% 
  mutate(trt = case_when(trt=="unfed_control" ~ "Unfed-Control",
                         trt=="unfed_deworm" ~ "Unfed-Deworm",
                         trt=="fed_control" ~ "Fed-Control",
                         trt=="fed_deworm" ~ "Fed-Deworm"),
         trt = factor(trt, levels=c("Unfed-Control", "Unfed-Deworm",
                                    "Fed-Control", "Fed-Deworm"))) %>%
  group_by(year, trt, sex, puuv_ifa) %>%
  summarise(n = length(tag)) %>% ungroup() 

# png(here("n_puuv_breeders.png"), width=800, height=600)
# data %>%
#   ggplot(aes(x=sex, y=trt, size=n, color = puuv_ifa)) +
#   scale_size(range=c(2, 36), name = "N Voles") +
#   geom_point(alpha=0.5) +
#   scale_y_discrete(limits=rev) +
#   facet_wrap(~year) +
#   labs(x="Sex", color = "PUUV Serostatus", title="Count of Infected/Uninfected Voles") +
#   theme(axis.title.y = element_blank(),
#         axis.title.x =  element_text(size=20),
#         axis.text = element_text(size=16),
#         strip.text = element_text(size=18),
#         legend.text = element_text(size=14),
#         legend.title = element_text(size=15),
#         legend.position = "bottom")
# dev.off()



# data %>%
#   ggplot(aes(x=puuv_ifa, y=trt, size=n, color = puuv_ifa)) +
#   scale_size(range=c(2, 36), name = "N Voles") +
#   geom_point(alpha=0.5) +
#   scale_y_discrete(limits=rev) +
#   facet_grid(sex~year) +
#   labs(x="Sex", color = "PUUV Serostatus", title="Count of Infected/Uninfected Voles") +
#   theme(axis.title.y = element_blank(),
#         axis.title.x =  element_text(size=20),
#         axis.text = element_text(size=16),
#         strip.text = element_text(size=18),
#         legend.text = element_text(size=14),
#         legend.title = element_text(size=15),
#         legend.position = "bottom")

png(here("n_puuv_breeders_bar.png"), width=700, height=800)
p <- data %>%
ggplot(aes(fill=puuv_ifa, y=n, x=sex)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("#F2AD00", "#B40F20")) +
  facet_grid(trt~year) +
  labs(x="Sex", y="N Voles", fill = "PUUV Serostatus") +
  theme(axis.title.y = element_blank(),
        axis.title.x =  element_text(size=20),
        axis.text = element_text(size=12),
        strip.text.x = element_text(size=16),
        strip.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=15),
        legend.position = "bottom")
#to color the facet labels (it's a whole thing...)
    #https://github.com/tidyverse/ggplot2/issues/2096
g <- ggplot_gtable(ggplot_build(p))
strip_r <- which(grepl('strip-r', g$layout$name))
fills <- c("#B2DF8A", "#33A02C", "#CAB2D6", "#6A3D9A")
k <- 1
for (i in strip_r) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
plot(g)
dev.off()







## current infection status based on previous degree
trt_labs <- as_labeller(c("unfed_control" = "Unfed Control",
                         "unfed_deworm" = "Unfed Deworm",
                          "fed_control" = "Fed Control",
                          "fed_deworm" = "Fed Deworm"))

# wt <- netmets_puuv %>% drop_na(prev_wt.deg) %>%
#   ggplot(aes(x=sex, y=prev_wt.deg, fill=puuv_ifa, color=puuv_ifa)) +
#   geom_violin() +
#   facet_grid(trt ~ year)

library(ggbeeswarm)

p <- netmets_puuv %>% drop_na(prev_wt.deg) %>%
  ggplot(aes(x=sex, y=prev_F.deg, fill=puuv_ifa)) +
  geom_boxplot() +
  geom_beeswarm(dodge.width=0.8, alpha=0.5) +
  # this plots the violin-style point distribution. groupOnX is not really necessary here, 
  # but throws a warning if not included (this bit from Matt Michalska-Smith)
  geom_quasirandom(dodge.width=0.8, groupOnX=TRUE, show.legend=FALSE, alpha=0.5) +
  scale_fill_manual(values=c("#F2AD00", "#B40F20")) +
  scale_color_manual(values=c("#F2AD00", "#B40F20")) +
  facet_grid(year ~ trt, 
             labeller = labeller(trt = trt_labs)) +
  labs(title="Prior Month Female Network Degree", x="Sex") +
  theme(plot.title = element_text(size= 24, hjust = 0.5),
        legend.position = "none",
        strip.text = element_text(size=20),
        axis.text.x = element_text(size=18),
        axis.title.x = element_text(size=20),
        axis.title.y = element_blank(),
        plot.margin = unit(c(1,0.5,1,1), "cm"))
gf <- ggplot_gtable(ggplot_build(p))
strip_t <- which(grepl('strip-t', gf$layout$name))
fills <- c("#B2DF8A", "#33A02C", "#CAB2D6", "#6A3D9A")
k <- 1
for (i in strip_t) {
  j <- which(grepl('rect', gf$grobs[[i]]$grobs[[1]]$childrenOrder))
  gf$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
plot(gf)

p <- netmets_puuv %>% drop_na(prev_wt.deg) %>%
  ggplot(aes(x=sex, y=prev_M.deg, fill=puuv_ifa)) +
  geom_boxplot() +
  geom_beeswarm(dodge.width=0.8, alpha=0.4) +
  # this plots the violin-style point distribution. groupOnX is not really necessary here, 
  # but throws a warning if not included (this bit from Matt Michalska-Smith)
  geom_quasirandom(dodge.width=0.8, groupOnX=TRUE, show.legend=FALSE, alpha=0.4) +
  scale_fill_manual(values=c("#F2AD00", "#B40F20")) +
  scale_color_manual(values=c("#F2AD00", "#B40F20")) +
  facet_grid(year ~ trt, 
             labeller = labeller(trt = trt_labs)) +
  labs(title="Prior Month Male Network Degree", x="Sex", fill="PUUV Serostatus") +
  theme(plot.title = element_text(size= 24, hjust = 0.5),
        legend.position = "none",
        strip.text = element_text(size=20),
        axis.text.x = element_text(size=18),
        axis.title.x = element_text(size=20),
        axis.title.y = element_blank(),
        plot.margin = unit(c(1,1,1,0.5), "cm"))
gm <- ggplot_gtable(ggplot_build(p))
strip_t <- which(grepl('strip-t', gm$layout$name))
fills <- c("#B2DF8A", "#33A02C", "#CAB2D6", "#6A3D9A")
k <- 1
for (i in strip_t) {
  j <- which(grepl('rect', gm$grobs[[i]]$grobs[[1]]$childrenOrder))
  gm$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
plot(gm)


library(cowplot)

png(here("malefemale_prevdeg.png"), width=1400, height=600)
plot_grid(gf, gm, labels="AUTO", label_size = 30)
dev.off()


# ### Current network position and infection status
# 
# wt <- netmets_puuv_serov %>%
#   ggplot(aes(x=sex, y=wt.deg, color=puuv_ifa)) +
#   geom_violin() +
#   facet_grid(trt ~ year)
# 
# female <- netmets_puuv_serov %>%
#   ggplot(aes(x=sex, y=F.deg, color=puuv_ifa)) +
#   geom_violin() +
#   facet_grid(trt ~ year)
# 
# male <- netmets_puuv_serov %>%
#   ggplot(aes(x=sex, y=M.deg, color=puuv_ifa)) +
#   geom_violin() +
#   facet_grid(trt ~ year)
# 
# grid.arrange(wt, female, male, ncol=3)





#likelihood of seroconverting based on previous network position
mod <- glm(serovert ~ prev_M.deg:sex + prev_F.deg:sex + sex + season_breeder + explore + 
             trt + prev_month + prev_n.node + year,
             family=binomial, data=netmets_puuv_serov)
#SINGULARITY IS AN ISSUE to do glmer (1|site) - some sites 0 or 1 seroverter
summary(mod)
#previous degree only affects p|serovert for females:
    #more likely to serovert if female with high female degree
    #less likely to serovert if female with high male degree (... um what)


mod <- glm(serovert ~ prev_M.deg:sex:season_breeder + prev_F.deg:sex:season_breeder + sex + season_breeder + explore + 
             trt + prev_month + prev_n.node + year,
           family=binomial, data=netmets_puuv_serov)
summary(mod)

mod <- glm(serovert ~ prev_b.deg:sex + prev_nb.deg:sex + sex + season_breeder + explore + 
             trt + prev_month + prev_n.node + year,
           family=binomial, data=netmets_puuv_serov)
summary(mod)


mod %>% tbl_regression(exponentiate = TRUE,
                       pvalue_fun = ~ style_pvalue(.x, digits = 2),) %>%
  bold_p(t = 0.10) %>%
  bold_labels() %>%
  italicize_levels()

# #check for singularity (in glmer)
# #from here: https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
# tt <- getME(mod,"theta")
# ll <- getME(mod,"lower")
# min(tt[ll==0]) # if this is less than 0.05 then you have problems

#########################################################################################################






##########################################################################################################
#How does previous degree affect current infection status?

dat <- netmets_puuv %>% drop_na(prev_wt.deg) %>%
  rename(Sex = sex,
         Treatment = trt,
         Previous_Month = prev_month,
         Previous_Network_Size = prev_n.node,
         Year = year,
         Previous_M.degree = prev_M.deg,
         Previous_F.degree = prev_F.deg)

# dat %>%
#   mutate(puuv_ifa = case_when(puuv_ifa == "0" ~ 0,
#                               puuv_ifa == "1" ~ 1)) %>%
#   ggplot(aes(x=prev_wt.deg, y=puuv_ifa, color = Sex)) +
#   geom_point() +
#   geom_smooth() +
#   # geom_smooth(method = "glm", method.args = list(family = "binomial")) +
#   facet_wrap(~Sex)
# 
# dat %>%
#   mutate(puuv_ifa = case_when(puuv_ifa == "0" ~ 0,
#                               puuv_ifa == "1" ~ 1)) %>%
#   ggplot(aes(x=Previous_F.degree, y=puuv_ifa, color = Sex)) +
#   geom_point() +
#   geom_smooth() + 
#   # geom_smooth(method = "glm", method.args = list(family = "binomial")) +
#   facet_wrap(~Sex)
# 
# dat %>%
#   mutate(puuv_ifa = case_when(puuv_ifa == "0" ~ 0,
#                               puuv_ifa == "1" ~ 1)) %>%
#   ggplot(aes(x=Previous_M.degree, y=puuv_ifa, color = Sex)) +
#   geom_point() +
#   geom_smooth() + 
#   # geom_smooth(method = "glm", method.args = list(family = "binomial")) +
#   facet_wrap(~Sex)

# ## at EEID, Megan Tomamichael suggested running model separately by trt to see if effects are washing out other things...
#   #can see piecemealy how different treatment types affect prev degree affecting current infection status, but overall meh
# dat_trt <- dat %>% filter(Treatment=="unfed_control")
# mod <- glmer(puuv_ifa ~ Previous_M.degree:Sex + Previous_F.degree:Sex + Sex + 
#               Previous_Month + Previous_Network_Size + Year + (1|site),
#              family=binomial, data=dat_trt)


############ IS THERE A WAY to have degree be SEX AND BREEDER specific? ie so overlapping with a malebreeder != malenonbreeder
##### I mean, it's easy enough to make a sex-breeder column with 4 categories -- but does that make things more complicated in a model?

mod <- glmer(puuv_ifa ~ prev_wt.deg:Sex:season_breeder + 
                    Sex + season_breeder + explore +
                    Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
                  # control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                  family=binomial, data=dat)
#WARNING: Failed to converge
summary(mod)


mod_sex <- glmer(puuv_ifa ~ Previous_M.degree:Sex:season_breeder + Previous_F.degree:Sex:season_breeder + 
               Sex + season_breeder + explore +
               Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
               # control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=5e5)),
             family=binomial, data=dat)
#WARNING: Failed to converge
summary(mod_sex)

mod_breed <- glmer(puuv_ifa ~ prev_b.deg:Sex:season_breeder + prev_nb.deg:Sex:season_breeder + 
               Sex + season_breeder + explore +
               Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
             family=binomial, data=dat)
summary(mod_breed)


AIC(mod, mod_sex, mod_breed) ## breeder model has a lower AIC than sex model


# #from here: https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
#     #related: https://stats.stackexchange.com/questions/164457/r-glmer-warnings-model-fails-to-converge-model-is-nearly-unidentifiable
#     #also: https://www.learn-mlms.com/07-module-7.html#learning-objectives-5
# #1. check singularity
# tt <- getME(mod,"theta")
# ll <- getME(mod,"lower")
# min(tt[ll==0]) # if this is big, you Gucci (close to 0 is bad)
# #2. Check the gradient calculations (I don't know what this means or what its doing)
# #3. Restart from the previous fit, increase max number of iterations from 10k to 20k
# ss <- getME(mod,c("theta","fixef"))
# m2 <- update(mod, start=ss, control=glmerControl(optCtrl=list(maxfun=2e4)))
# ## YES - this converges
# summary(m2)
# #3. Try a different optimizer - e.g. bobyqa for first AND second phase (default is 1st phase, Nelder-Mead 2nd phase)
# #fixed with : adding "control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))"

mod1 <- glmer(puuv_ifa ~ MODEL_PARAMS_GO_HERE,
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
             family=binomial, data=dat)

## HOWEVER - the model output between mod0 and mod1 is basically the same, so warning was a false positive
    #Ben Bolker says it's okay to use the OG model: 
    #https://stackoverflow.com/questions/33670628/solution-to-the-warning-message-using-glmer

# mod2 <- glmer(puuv_ifa ~ prev_M.deg + prev_F.deg:sex + sex + 
#                trt + prev_month + prev_n.node + year + (1|site),
#              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
#              family=binomial, data=dat)
# 
# mod3 <- glmer(puuv_ifa ~ prev_M.deg + prev_F.deg + sex + 
#                trt + prev_month + prev_n.node + year + (1|site),
#              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
#              family=binomial, data=dat)
# 
# AIC(mod1, mod2, mod3)

summary(mod1)

#GLMM model diagnostics > https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
library(DHARMa)
#calculate residuals (then run diagnostics on these)
simulationOutput <- simulateResiduals(fittedModel = mod_sex)
plot(simulationOutput) #qq plot and residual vs fitted
testDispersion(simulationOutput) #formal test for overdispersion
# testZeroInflation(simulationOutput) #formal test for zero inflation (common type of overdispersion)

plotResiduals(simulationOutput, dat$Sex)
plotResiduals(simulationOutput, dat$season_breeder)
plotResiduals(simulationOutput, dat$explore)
plotResiduals(simulationOutput, dat$Treatment)
plotResiduals(simulationOutput, dat$Month)
plotResiduals(simulationOutput, dat$Previous_Network_Size)
plotResiduals(simulationOutput, dat$Year)
plotResiduals(simulationOutput, dat$Previous_F.degree)
plotResiduals(simulationOutput, dat$Previous_M.degree)

#### OVERALL: model things look pretty good - 
    ## the residuals are fine, it's not overdispersed, singularity is okay etc


######## OUTPUT MODEL SUMMARY #########

library(gtsummary) #https://www.danieldsjoberg.com/gtsummary/articles/tbl_regression.html

gtsumm_mod <- mod %>% tbl_regression(exponentiate = TRUE,
  pvalue_fun = ~ style_pvalue(.x, digits = 2),) %>%
  bold_p(t = 0.10) %>%
  bold_labels() %>%
  italicize_levels()

library(gt)
# Use function from gt package to save table, after converting to 
# gt object using as_gt()
gt::gtsave(as_gt(gtsumm_mod), expand=30, here("regression_model.png"))

####### end output model summary #########


###############################################################################################

