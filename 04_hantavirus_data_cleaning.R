## NEW, tidy version (was previously called "03_NEW_hanta_cleaning.R") created March 8 2024

# code run using R version _______________

#load libraries
library(here)
library(tidyverse)
library(janitor)
library(lubridate)
library(lme4)
library(gtsummary) #pretty regression summary tables

#clear environment
rm(list = ls())


#### all the steps to generate 'netmets_puuv' have been commented out as of 8 March 2024 - readRDS to load file

# ##########---------- LOAD NETWORK METRICS and VOLE CAPTURE METADATA -----------###########
# 
# # ### Fulltrap and netmets dfs created separately, load 2021 and 2022 data here, clean as needed
# # 
# # ### LOAD INDIVIDUAL VOLE METADATA by year (now generated in this R project)
# # ### these versions have all the animals (breeders and nonbreeders) but month=may, sex=NA, and repro=NA removed
# # fulltrap21 <- readRDS(file="fulltrap21_03.04.24.rds") 
# # fulltrap22 <- readRDS(file="fulltrap22_03.04.24.rds")
# # 
# # ## >>NOTE<< overwintered animals may incorrectly have firstcap==1 in 2022
# #     # PIT tag: 219895 Occ2Sess1 at Vaarinkorpi
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
# #   mutate(firstcap = factor(firstcap)) %>% ungroup() %>% 
# # #keep one entry per tag,year,month; pull only relevant columns
# #   group_by(tag, year, month) %>% slice(1) %>%
# #   select(year, site, trt, month, tag, samp_id, sex, season_breeder, traps_per_life, caps_per_life) %>%
# #   ungroup()
# # 
# # ### LOAD NETWORK METRICS data
# # netmets21 <- readRDS(file="network_metrics21.rds") %>% mutate(year=as.numeric(2021))
# # netmets22 <- readRDS(file="network_metrics22.rds") %>% mutate(year=as.numeric(2022))
# # 
# # netmets21.22 <- rbind(netmets21, netmets22)
# # 
# # ### JOIN NETWORK METRICS + METADATA
# # netmets_full <- left_join(netmets21.22, fulltrap21.22, by=c("year", "site", "month", "tag")) %>%
# #   mutate(site=as.factor(site),
# #          year=as.factor(year),
# #          month = as.factor(month),
# #          month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #remove may from factor levels
# #   select(!c(focal_id)) %>% #remove duplicate column for PIT tag number
# #   relocate(c(year, trt, site, month, n.node, tag, samp_id, sex, season_breeder), .before = wt.deg)
# # 
# # saveRDS(netmets_full, here("netmets_full_03.08.24.rds"))
# 
# #load the most recent version of netmets_full
# netmets_full <- readRDS(here("netmets_full_03.08.24.rds"))
# 
# 
# #########################################   LOAD & CLEAN PUUV IFA DATA   ########################################
# 
# # #load, clean, format PUUV IFA data
# # ### NEW! Updated! PUUV_IFA data with the goofy samples from 2021/2022 that KWearing re-ran in May 2023 (puuv_ifa_6.7.23.csv)
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
# ## there are 120 entries that we have network data for a vole but DO NOT know their PUUV status
# 
# 
# ############ REMOVE THE VOLES from netmets_puuv THAT SEROCONVERT POS TO NEG ##################
# 
# #voles that seroconvert PUUV+ to PUUV-
# puuv_pos_neg <- netmets_puuv %>% group_by(tag) %>% arrange(year, month, .by_group = TRUE) %>%
#   dplyr::select(year, month, tag, puuv_ifa) %>%
#   summarise(status_time = toString(puuv_ifa)) %>%
#   filter(str_detect(status_time, "1,\\s0")) #filter for animals that go from pos to neg
# #pull the PIT tags (21 individuals (50 entries), mostly 1-0, a few 1-0-1)
# problemchildren <- puuv_pos_neg$tag
# #filter netmetsPUUV to remove 'problemchildren' - status is inconclusive or maybe we detected MatAb
# netmets_puuv <- netmets_puuv %>%
#   filter(!tag %in% problemchildren)
# 
# ## NOW NETMETS_PUUV has ONLY the animals that:  ##
#     #       1. Were caught in June-October
#     #       2. Have sex and repro data
#     #       3. Have network data
#     #       4. Have PUUV IFA data
#     #       5. DON'T seroconvert from PUUV+ to PUUV-
# 
# ######################### end problem children ################################
# 
# 
# 
# ## Sept 27, 2023
# ## just a note - there are 170 entries from netmets21 and 22 that DO NOT make it in to netmets_PUUV
# ## these are either animals without PUUV data (120 entries) or animals with inconsistent IFA results (50 entries)
# ## HOWEVER - these animals were in the networks when they were built
# ##      so even once they're gone, other voles will have 'degree' measures to them
# ## March 8, 2024
# ## I think this is okay, we know those 170 entries are legit, voles were observed, so they should contribute to the networks
# ## none of the models run with infection status of neighbors
# 
# 
# 
# ###############################################################################
# ### exploratoriness of voles (breeders/non-breeders) ###
# ## based off of VanderWaal ground squirrel ms ##
# #https://modelr.tidyverse.org/reference/add_predictions.html
# 
# #subset data to one entry per vole per year
# onepertag <- netmets_puuv %>%
#   group_by(year, tag) %>% slice(1)
# 
# #fit power model to ln(traps) ~ ln(caps)
# powmod <- lm(log(traps_per_life) ~ log(caps_per_life), data=onepertag)
# summary(powmod) #check it
# # # ln(y) = -0.01609 + 0.73633ln(x)
# # # y = 0.9840388x^0.73633
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
# ######### lagged degree and seroconvert status ###############
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
#          prev_mb.deg = lag(mb.deg, n=1),
#          prev_mnb.deg = lag(mnb.deg, n=1),
#          prev_fb.deg = lag(fb.deg, n=1),
#          prev_fnb.deg = lag(fnb.deg, n=1),
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
# saveRDS(netmets_puuv, here("netmets_puuv_03.08.24.rds"))

#load netmets_puuv
netmets_puuv <- readRDS(here("netmets_puuv_03.08.24.rds"))


#### summarize some quick counts

# #entries per year (in netmets_puuv)
# y1 <- netmets_puuv %>% filter(year=="2021") #1029 in 2021
# n_distinct(y1$tag) #683 unique voles
# y2 <- netmets_puuv %>% filter(year=="2022") #1061 in 2022
# n_distinct(y2$tag) #694 unique voles
# 
# #entries per year (in netmets --> ie those used in the networks that don't have PUUV data)
# nrow(netmets21) #1129 entries in 2021
# n_distinct(netmets21$tag) #742 voles in 2021
# nrow(netmets22) #1131 entries in 2022
# n_distinct(netmets22$tag) #744 voles in 2022


# ## a random thing: Sept 20, 2023 ##
# ## Kris wanted me to report in the VoleHanta ms the number of animals
#     ##captured, recapped etc before detailing the subset of the data used for the study
# 
# ## so that would be the fulltrap_ALL_03.04.24 dataset
# #INCLUDING MAY, INCLUDING sex/repro=NA
# fulltrap21_ALL <- readRDS(file="fulltrap21_ALL_03.04.24.rds") 
# fulltrap22_ALL <- readRDS(file="fulltrap22_ALL_03.04.24.rds")
# fulltrap21.22_ALL <- rbind(fulltrap21_ALL, fulltrap22_ALL) 
# #4286 capture events (including WR) across MAY-Oct 2021 and 2022
# samples <- fulltrap21.22_ALL %>% drop_na(samp_id) %>% group_by(samp_id) %>% slice(1)
# #2309 captures with samples collected
# n_distinct(samples$tag)
# #1487 unique voles were sampled
# samp_recaps <- samples %>% filter(caps_per_life>1) %>% group_by(tag) %>% slice(1)
# #924 unique voles were sampled and recapped



##################################################################################################################




########### a new thing wee hours of 7/7/23 trying to finish this damn ms ###############

# #how many times is the "previous network position" NOT from the immediately prior month?
# 
# df <- netmets_puuv %>% group_by(tag) %>% mutate(month.n = case_when(month=="june" ~ 1,
#                                                               month=="july" ~ 2,
#                                                               month=="aug" ~ 3,
#                                                               month=="sept" ~ 4,
#                                                               month=="oct" ~ 5)) %>%
#   mutate(dif = month.n - lag(month.n)) %>%
#   drop_na(prev_wt.deg)
# #713 data entries with a current PUUV informed by previous network position
# 
# dfsmol <- df %>% filter(dif > 1) #25 times we had to pull from an earlier month that wasn't the previous one
# 
# ### 25/713 = 3.5%

#########################################################################################







############################ DEGREE DISTRIBUTION BY TRT / YEAR ###########################################

#no figures used in manuscript or supplement, summary values reported in results

netmets_puuv %>% group_by(year, trt) %>%
  summarise(mean=mean(wt.deg),
            sd=sd(wt.deg))

# netmets_puuv %>% group_by(year, month, trt) %>%
#   summarise(mean=mean(wt.deg),
#             sd=sd(wt.deg))

netmets_puuv %>% group_by(year) %>%
  summarise(mean=mean(wt.deg),
            sd=sd(wt.deg))

trt_labs <- as_labeller(c(unfed_control="Unfed Control", 
                          unfed_deworm="Unfed Deworm", 
                          fed_control="Fed Control", 
                          fed_deworm="Fed Deworm"))

#increase axis ticks: https://stackoverflow.com/questions/11335836/increase-number-of-axis-ticks
#add mean for each facet: https://stackoverflow.com/questions/44196384/how-to-produce-different-geom-vline-in-different-facets-in-r

meandata <- netmets_puuv %>% group_by(year, trt) %>%
  summarize(mean_x = mean(wt.deg))

#weighted degree by trt, year bars for mean
netmets_puuv %>%
  ggplot(aes(x=wt.deg)) +
  geom_histogram(stat="bin") +
  geom_vline(data=meandata, aes(xintercept=mean_x, color=trt), linewidth=1, 
             show.legend = FALSE) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  facet_grid(trt~year, labeller=labeller(trt=trt_labs)) +
  xlab("Weighted Degree") + ylab("Count")

# #weighted degree by trt, month
# #2021 data
# netmets_puuv %>% filter(year=="2021") %>%
#   ggplot(aes(x=wt.deg)) +
#   geom_histogram(stat="bin") +
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
#   facet_grid(trt~month) +
#   xlab("Weighted Degree") + ylab("Count")
# #2022 data
# netmets_puuv %>% filter(year=="2022") %>%
#   ggplot(aes(x=wt.deg)) +
#   geom_histogram(stat="bin") +
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
#   facet_grid(trt~month) +
#   xlab("Weighted Degree") + ylab("Count")
############################################################################




############### PUUV PREVALENCE PER SITE #########################

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
# #vignette: http://127.0.0.1:18092/library/epiR/doc/epiR_descriptive.html
# library(epiR)
# #Method prevalence require a two-column matrix; the first column specifies the number of positives,
#     #the second column specifies the total number tested.
# epi.conf.data <- puuv_prev %>% ungroup() %>% select(n_pos, n) %>% as.matrix()
# prevCI <- epi.conf(dat=epi.conf.data, ctype="prevalence", method="exact", conf.level = 0.95) 
# #can add N=1000 (population size), design=1 (design effect) but it doesn't change the output
#   #"The design effect is used to adjust the confidence interval around a prevalence or incidence risk
#   #estimate in the presence of clustering. Adjustment for the effect of
#   #clustering can only be made on those prevalence and incidence risk methods that return a standard
#   #error (i.e., method = "wilson" or method = "fleiss")"
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

############################## end prevalence ####################################










######################## How does previous degree affect current infection status? ############################

netmets_puuv_prevdeg <- netmets_puuv %>% drop_na(prev_wt.deg) %>%
  rename(Sex = sex,
         Treatment = trt,
         Previous_Month = prev_month,
         Previous_Network_Size = prev_n.node,
         Year = year,
         Previous_M.degree = prev_M.deg,
         Previous_F.degree = prev_F.deg)

# # number of captures, individual voles (captures is not double # voles because indivs capped 2x only get only entry)
# y1 <- netmets_puuv_prevdeg %>% filter(Year=="2021")
# n_distinct(y1$tag)
# y2 <- netmets_puuv_prevdeg %>% filter(Year=="2022")
# n_distinct(y2$tag)


# ## at EEID, Megan Tomamichael suggested running model separately by trt to see if effects are washing out other things...
#   #can see piecemealy how different treatment types affect prev degree affecting current infection status, but overall meh
# dat_trt <- dat %>% filter(Treatment=="unfed_control")
# mod <- glmer(puuv_ifa ~ Previous_M.degree:Sex + Previous_F.degree:Sex + Sex +
#               Previous_Month + Previous_Network_Size + Year + (1|site),
#              family=binomial, data=dat_trt)


mod_wdegs <- glmer(puuv_ifa ~ prev_wt.deg:Sex + 
                    Sex + season_breeder + explore +
                    Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                  family=binomial, data=netmets_puuv_prevdeg)

mod_wdegb <- glmer(puuv_ifa ~ prev_wt.deg:season_breeder + 
                    Sex + season_breeder + explore +
                    Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                  family=binomial, data=netmets_puuv_prevdeg)

mod_wdegsb <- glmer(puuv_ifa ~ prev_wt.deg:Sex:season_breeder + 
                    Sex + season_breeder + explore +
                    Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                  family=binomial, data=netmets_puuv_prevdeg)
#WARNING: Failed to converge
# summary(mod_wdeg)


mod_sexs <- glmer(puuv_ifa ~ Previous_M.degree:Sex + Previous_F.degree:Sex + 
                   Sex + season_breeder + explore +
                   Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                 family=binomial, data=netmets_puuv_prevdeg)

mod_sexb <- glmer(puuv_ifa ~ Previous_M.degree:season_breeder + Previous_F.degree:season_breeder + 
                   Sex + season_breeder + explore +
                   Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                 family=binomial, data=netmets_puuv_prevdeg)

mod_sexsb <- glmer(puuv_ifa ~ Previous_M.degree:Sex:season_breeder + Previous_F.degree:Sex:season_breeder + 
               Sex + season_breeder + explore +
               Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
               control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
             family=binomial, data=netmets_puuv_prevdeg)
#WARNING: Failed to converge
# summary(mod_sex)

mod_breeds <- glmer(puuv_ifa ~ prev_b.deg:Sex + prev_nb.deg:Sex + 
                     Sex + season_breeder + explore +
                     Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                   family=binomial, data=netmets_puuv_prevdeg)

mod_breedb <- glmer(puuv_ifa ~ prev_b.deg:season_breeder + prev_nb.deg:season_breeder + 
                     Sex + season_breeder + explore +
                     Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                   family=binomial, data=netmets_puuv_prevdeg)

mod_breedsb <- glmer(puuv_ifa ~ prev_b.deg:Sex:season_breeder + prev_nb.deg:Sex:season_breeder + 
               Sex + season_breeder + explore +
               Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
               control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
             family=binomial, data=netmets_puuv_prevdeg)
 #WARNING: Failed to converge
summary(mod_breedsb)

mod_sbs <- glmer(puuv_ifa ~ prev_mb.deg:Sex + prev_mnb.deg:Sex + 
                  prev_fb.deg:Sex + prev_fnb.deg:Sex + 
                  Sex + season_breeder + explore +
                  Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                family=binomial, data=netmets_puuv_prevdeg)

mod_sbb <- glmer(puuv_ifa ~ prev_mb.deg:season_breeder + prev_mnb.deg:season_breeder + 
                  prev_fb.deg:season_breeder + prev_fnb.deg:season_breeder + 
                  Sex + season_breeder + explore +
                  Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                family=binomial, data=netmets_puuv_prevdeg)

mod_sbsb <- glmer(puuv_ifa ~ prev_mb.deg:Sex:season_breeder + prev_mnb.deg:Sex:season_breeder + 
                           prev_fb.deg:Sex:season_breeder + prev_fnb.deg:Sex:season_breeder + 
                     Sex + season_breeder + explore +
                     Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                   family=binomial, data=netmets_puuv_prevdeg)
#WARNING: Failed to converge
# summary(mod_sb)

AIC(mod_wdegs, mod_wdegb, mod_wdegsb, 
    mod_sexs, mod_sexsb, mod_sexb,
    mod_breeds, mod_breedb, mod_breedsb,
    mod_sbs, mod_sbb, mod_sbsb) ## mod_breedsb has lowest AIC score

#generate pretty table of regression coefs
# library(gtsummary)
mod_sbsb %>% tbl_regression(exponentiate = TRUE,
                           pvalue_fun = ~ style_pvalue(.x, digits = 2),) %>%
  bold_p(t = 0.10) %>%
  bold_labels() %>%
  italicize_levels()

# ############## if model fails to converge: ##################
# #from here: https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
#     #related: https://stats.stackexchange.com/questions/164457/r-glmer-warnings-model-fails-to-converge-model-is-nearly-unidentifiable
#     #also: https://www.learn-mlms.com/07-module-7.html#learning-objectives-5
# #1. check singularity (works for glmer, NOT glm)
# #from here: https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
# tt <- getME(mod_breed,"theta")
# ll <- getME(mod_breed,"lower")
# min(tt[ll==0]) # if this is big, you Gucci (close to 0 is bad)
# #2. Check the gradient calculations (I don't know what this means or what its doing)
# #3. Restart from the previous fit, increase max number of iterations from 10k to 20k
# ss <- getME(mod_sex,c("theta","fixef"))
# m2 <- update(mod_sex, start=ss, control=glmerControl(optCtrl=list(maxfun=2e4)))
# ## YES - this converges - but the logLik and param ests are basically the same
# summary(m2)
# #3. Try a different optimizer - e.g. bobyqa for first AND second phase (default is 1st phase, Nelder-Mead 2nd phase)
# #fixed with : adding "control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))"
# 
# mod_try <- glmer(puuv_ifa ~ MODEL_PARAMS_GO_HERE,
#              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
#              family=binomial, data=dat)
# 
# ## HOWEVER - the model output between mod0 and mod1 is basically the same, so warning was a false positive
#     #Ben Bolker says it's okay to use the OG model: 
#     #https://stackoverflow.com/questions/33670628/solution-to-the-warning-message-using-glmer
# ##################################################


####### model diagnostics #########
#AIC best models are mod_breedsb and mod_sbsb

dat <- netmets_puuv_prevdeg

#GLMM model diagnostics > https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
library(DHARMa)
#calculate residuals (then run diagnostics on these)
simulationOutput <- simulateResiduals(fittedModel = mod_sbsb)
plot(simulationOutput) #qq plot and residual vs fitted
testDispersion(simulationOutput) #formal test for overdispersion
# testZeroInflation(simulationOutput) #formal test for zero inflation (common type of overdispersion)

plotResiduals(simulationOutput, dat$Sex)
plotResiduals(simulationOutput, dat$season_breeder)
plotResiduals(simulationOutput, dat$explore)
plotResiduals(simulationOutput, dat$Treatment)
plotResiduals(simulationOutput, dat$Previous_Month)
plotResiduals(simulationOutput, dat$Previous_Network_Size)
plotResiduals(simulationOutput, dat$Year)
plotResiduals(simulationOutput, dat$Previous_F.degree) ## this one is a little off, but honestly not that much so
plotResiduals(simulationOutput, dat$Previous_M.degree)
plotResiduals(simulationOutput, dat$prev_b.deg)
plotResiduals(simulationOutput, dat$prev_nb.deg)

#### OVERALL: model things look pretty good - 
    ## the residuals are fine, it's not overdispersed, singularity is okay etc



######## pretty OUTPUT MODEL SUMMARY #########

library(gtsummary) #https://www.danieldsjoberg.com/gtsummary/articles/tbl_regression.html

gtsumm_mod <- mod_sex %>% tbl_regression(exponentiate = TRUE,
  pvalue_fun = ~ style_pvalue(.x, digits = 2),) %>%
  bold_p(t = 0.10) %>%
  bold_labels() %>%
  italicize_levels()

library(gt)
# Use function from gt package to save table, after converting to 
# gt object using as_gt()
gt::gtsave(as_gt(gtsumm_mod), expand=30, here("regression_model.png"))

####### end output model summary #########




############ visualize pairwise IFA by degree with smooth ##################
# dat %>%
#   mutate(puuv_ifa = case_when(puuv_ifa == "0" ~ 0,
#                               puuv_ifa == "1" ~ 1)) %>%
#   ggplot(aes(x=prev_wt.deg, y=puuv_ifa, color = Sex)) +
#   geom_point() +
#   geom_smooth() +
#   # geom_smooth(method = "glm", method.args = list(family = "binomial")) +
#   facet_wrap(~Sex)
#   # facet_grid(Sex ~ Year)
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
# 
# dat %>%
#   mutate(puuv_ifa = case_when(puuv_ifa == "0" ~ 0,
#                               puuv_ifa == "1" ~ 1)) %>%
#   ggplot(aes(x=prev_b.deg, y=puuv_ifa, color = Sex)) +
#   geom_point() +
#   geom_smooth() +
#   # geom_smooth(method = "glm", method.args = list(family = "binomial")) +
#   facet_wrap(~Sex)
# 
# dat %>%
#   mutate(puuv_ifa = case_when(puuv_ifa == "0" ~ 0,
#                               puuv_ifa == "1" ~ 1)) %>%
#   ggplot(aes(x=prev_nb.deg, y=puuv_ifa, color = Sex)) +
#   geom_point() +
#   geom_smooth() +
#   # geom_smooth(method = "glm", method.args = list(family = "binomial")) +
#   facet_wrap(~Sex)



###############################################################################################









######################### some random visualizations ##################################

# ## COMPARE female and male degree of a given animal, color by trt
# ## not sure what this is useful for, but interesting
# netmets_puuv %>%
#   ggplot(aes(x=F.deg, y=M.deg, color=trt, )) +
#   geom_point(aes(shape=sex)) +
#   scale_shape_manual(values=c(3, 16)) +
#   coord_fixed() +
#   xlim(0,4) + ylim(0,4) +
#   facet_grid(year ~ month)

######################################################################################



