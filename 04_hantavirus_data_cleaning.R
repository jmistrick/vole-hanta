##### NOVEMBER 2024 - blowing it all up to work with new pct overlap data

### 6 - Hantavirus Data Cleaning
### AUTHOR
### 23 March 2024
### this code accompanies the manuscript: "Ecological factors alter how spatial overlap predicts viral 
  # infection dynamics in wild rodent populations"
### Run using R version 4.3.2 (2023-10-31) -- "Eye Holes"

### PURPOSE: 
# THIS CODE loads and cleans the hantavirus infection data (Immunofluorescence assay serology) and links
# the infection data to the spatial data per vole generated in 07_percent_overlap.R

###------------------------------------------------------------------------------


#load libraries
library(here) #v 1.0.1
library(tidyverse) #v 2.0.0
library(janitor) #v 2.2.0
library(lubridate) #v 1.9.3
library(lme4) #v 1.1-35.1
library(gtsummary) #pretty regression summary tables #v 1.7.2

#clear environment
rm(list = ls())


##########---------- LOAD NETWORK METRICS and VOLE CAPTURE METADATA -----------###########

### Fulltrap and netmets dfs created separately, load 2021 and 2022 data here, clean as needed

### LOAD INDIVIDUAL VOLE METADATA by year 
fulltrap21 <- readRDS(file="fulltrap21_11.11.24.rds")
fulltrap22 <- readRDS(file="fulltrap22_11.11.24.rds")
fulltrap23 <- readRDS(file="fulltrap23_11.11.24.rds")

## >>NOTE<< overwintered animals may incorrectly have firstcap==1 in 2022 or 2023
    # PIT tag: 219895 Occ2Sess1 at Vaarinkorpi 2022
    # PIT tag: 226280 Occ2Sess2 at Kuoppa 2022

#join 21, 22, 23 fulltrap data; correct firstcap for 2022 and 2023
fulltrapALL <- rbind(fulltrap21, fulltrap22, fulltrap23) %>%
  group_by(tag) %>% arrange(tag, year, occ.sess, .by_group = TRUE) %>%
  mutate(firstcap = ifelse(date_time == min(date_time), 1, 0)) %>%
  mutate(firstcap = factor(firstcap)) %>% ungroup() %>%
#keep one entry per tag,year,month; pull only relevant columns
  group_by(tag, year, month) %>% slice(1) %>%
  select(year, site, trt, month, tag, samp_id, sex, season_breeder, traps_per_year, caps_per_year) %>%
  ungroup()

### LOAD NETWORK METRICS data
netmets21 <- readRDS(file="pctover_netmets21.rds") 
netmets22 <- readRDS(file="pctover_netmets22.rds") 
netmets23 <- readRDS(file="pctover_netmets23.rds") 
#join 2021 and 2022
netmetsALL <- rbind(netmets21, netmets22, netmets23)

### JOIN NETWORK METRICS + METADATA
netmets_full <- left_join(netmetsALL, fulltrapALL, by=c("year", "site", "month", "tag")) %>%
  mutate(site=as.factor(site),
         year=as.factor(year),
         month = as.factor(month),
         month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #remove may from factor levels
  select(!c(focal_id)) %>% #remove duplicate column for PIT tag number
  relocate(c(year, trt, site, month, n.node, tag, samp_id, sex, season_breeder), .before = wt.deg.in)

saveRDS(netmets_full, here("netmets_full_11.11.24.rds"))

# #load the most recent version of netmets_full
# netmets_full <- readRDS(here("netmets_full_11.11.24.rds"))


#########################################   LOAD & CLEAN PUUV IFA DATA   ########################################

# #load PUUV IFA data
# puuv_data <- read.csv(here("puuv_ifa_11.11.24.csv")) %>%
#   clean_names %>%
#   #populate column of FINAL PUUV status
#   #kind of a pain, since samples could be run 1-4x but we want the result of the last run as the 'final' status
#   #columns are named as 'puuv_run1' 'puuv_run2' 'puuv_run3' 'puuv_run4'
#   mutate(FINAL_puuv = ifelse(!is.na(puuv_run4), as.character(puuv_run4),
#                              ifelse(!is.na(puuv_run3), as.character(puuv_run3),
#                                     ifelse(!is.na(puuv_run2), as.character(puuv_run2), as.character(puuv_run1))))) %>%
#   mutate(samp_id = as.numeric(id),
#          # date_run1 = as_date(date_run1, format= "%m/%d/%Y"),
#          # puuv_run1 = as.factor(puuv_run1),
#          # date_run2 = as_date(date_run2, format= "%m/%d/%Y"),
#          # puuv_run2 = as.factor(puuv_run2),
#          # date_run3 = as_date(date_run3, format= "%m/%d/%Y"),
#          # puuv_run3 = as.factor(puuv_run3),
#          # date_run4 = as_date(date_run4, format= "%m/%d/%Y"),
#          # puuv_run4 = as.factor(puuv_run4),
#          FINAL_puuv = as.factor(FINAL_puuv)) %>%
#   drop_na(FINAL_puuv) %>%
#   dplyr::select(samp_id, FINAL_puuv) %>%
#   rename(puuv_ifa = FINAL_puuv)
# #output is a df with two columns, sample ID and PUUV status (0,1)
# 
# # Save puuv_data to a rdata file
# saveRDS(puuv_data, file = here("hantadata_11.11.24.rds"))

# Restore puuv_data from the rdata file
puuv_data <- readRDS(file = "hantadata_11.11.24.rds")


##################################  COMBINE NETMETS data with PUUV data ######################################

#join puuv_data to netmets_full (BOTH YEARS!)
netmets_puuv <- left_join(netmets_full, puuv_data, by="samp_id") %>%
  #left join on netmets_full because I need to have network data
  relocate(puuv_ifa, .after="samp_id") %>%
  select(!samp_id) %>%
  drop_na(puuv_ifa) #drop any animals without puuv data



### NEW!!! ### NOV 2024 UPDATE - there are 138 that drop out, so 18 new ones with 2023 data

########## NOTE: the step above DROPS voles with no PUUV infection status data
## there are 120 entries that we have network data for a vole but DO NOT know their PUUV status
## these voles ARE in the networks since we know they were there and where they were, 
## but we can't use for models of how spatial overlap impacts infection status



############ REMOVE THE VOLES from netmets_puuv THAT SEROCONVERT POS TO NEG ##################

## PUUV infection is life-long and chronic with no recovery - so voles should not go from 1-0 or 0-1-0 etc.
## this could be an incorrect result from IFA OR a young vole that had maternal antibodies at first capture
## either way, we can't trust these results and will remove ALL entries for voles with inconsistent IFA results
## result: we removed 50 entries with inconsistent hantavirus infection status

#voles that seroconvert PUUV+ to PUUV-
puuv_pos_neg <- netmets_puuv %>% group_by(tag) %>% arrange(year, month, .by_group = TRUE) %>%
  dplyr::select(year, month, tag, puuv_ifa) %>%
  summarise(status_time = toString(puuv_ifa)) %>%
  filter(str_detect(status_time, "1,\\s0")) #filter for animals that go from pos to neg
#pull the PIT tags (29 individuals (72 entries), mostly 1-0, a few 1-0-1)
problemchildren <- puuv_pos_neg$tag
#filter netmetsPUUV to remove 'problemchildren' - status is inconclusive or maybe we detected MatAb
netmets_puuv <- netmets_puuv %>%
  filter(!tag %in% problemchildren)

###################### FOLLOWING HAS NOT BEEN UPDATED WITH 2023 DATA ###########################################
####################### this probably doesn't need to be published #############################################
###### BRIEFLY, who are the voles that convert PUUV + to - ? #######
# #pull the puuv data on "problem children"
# problemfull <- netmets_puuv %>% filter(tag %in% problemchildren) %>%
#   select(year, month, tag, puuv_ifa)
# #grab body mass from original fulltrap file (NOT PUBLISHED)
# mass21 <- readRDS(file="fulltrap21_ALL_03.04.24.rds")
# mass22 <- readRDS(file="fulltrap22_ALL_03.04.24.rds")
# #combine 2021 and 2022 data, grab only the columns we need, filter for only problemchildren voles
# mass21.22 <- rbind(mass21, mass22) %>%
#   group_by(tag, year, month) %>% slice(1) %>%
#   select(tag, year, month, mass) %>%
#   filter(tag %in% problemchildren) %>%
#   mutate(year = as.factor(year))
# #dataframe with the year, month, tag, puuv status, and body mass for all problemchildren
# problemmass <- left_join(problemfull, mass21.22, by=c("tag", "month", "year"))
# problemmass %>% group_by(tag) %>% arrange(year, month) %>% slice(1) %>% ungroup() %>% #one entry per problemchild
#   filter(mass<17) %>% #remove vole that was 22g at first cap and PUUV+ and vole that went 0-1-0 (firstcap mass=17g) - leaving 19 voles
#   summarise(mean = mean(mass),
#             sd = sd(mass),
#             min = min(mass),
#             max = max(mass))
# #write to csv for quick visualizing
# write.csv(problemmass, file="problemchildren_mass.csv")
###############################################################################################################

## NOW NETMETS_PUUV has ONLY the animals that:  ##
    #       1. Were caught in June-October
    #       2. Have sex and repro data
    #       3. Have network data
    #       4. Have PUUV IFA data
    #       5. DON'T seroconvert from PUUV+ to PUUV-

######################### end problem children ################################

############ HAS NOT BEEN UPDATED WITH 2023 DATA #################
## BEFORE GOING ON - please note:
## there are 170 entries from netmets21 and netmets22 that DO NOT make it in to netmets_PUUV
## these are either animals without PUUV IFA data (120 entries) or animals with inconsistent IFA results (50 entries, 21 voles)
## (animals without PUUV IFA data could have been <10g and no blood was taken or...
## ...something else was wrong during processing and for whatever reason, we didn't get a blood sample or didn't have enough serum)
## HOWEVER - these animals were in the networks when they were constructed
##      so even once they're gone, other voles will have 'degree' measures to them
## I have decided this is okay, we know those 170 captures happened, voles were observed, so they should contribute to the networks
## none of the hantavirus infection probability models run with known infection status of neighbors
## and the models are already running on a subset of the full netmets21 and netmets22 datasets since we can
## only model current PUUV infection probability on animals captured at least twice
## THUS, I have convinced myself (and you hopefully) that my decision is fine and we can proceed




############ HAS NOT BEEN UPDATED WITH 2023 DATA #################

# ####################################################################
# ####################################################################
# ##### WHAT WAS THE PREVALENCE OF PUUV in the sampled animals? ######
# 
# ## NOTE! This is only animals with sex and repro data and NETMETS data (JUNE-Oct)
# puuv_prev <- netmets_puuv %>% 
#   group_by(year, tag) %>% #keep animals that were in capped both years in both years
#   arrange(month) %>%
#   slice(n()) %>% #keep the last entry for each tag
#   select(year, month, tag, puuv_ifa) %>%
#   mutate(puuv_num = case_when(puuv_ifa == 1 ~ 1,
#                               puuv_ifa == 0 ~ 0)) %>%
#   mutate(puuv_num = as.numeric(puuv_num)) %>% ungroup()
# 
# n_distinct(puuv_prev$tag)
# #1367 UNIQUE sampled animals 
# 
# #overwintered animals (10 animals were sampled in both 2021 and 2022)
# ow <- puuv_prev %>% group_by(tag) %>% mutate(n=length(tag)) %>% filter(n>1) %>%
#   arrange(tag)
# 
# puuv_prev %>% 
#   group_by(year) %>%
#   summarise(pos = sum(puuv_num),
#             n = length(puuv_num),
#             prev = pos/n)
# 
# #prevalence was 30.3% in 2021, 16.6% in 2022
# #################################################################
# #################################################################



############ HAS NOT BEEN UPDATED WITH 2023 DATA #################

# ###############################################################################
# ### exploratory behavior of voles (breeders/non-breeders) ###
# ## based off of VanderWaal et al 2013 ground squirrel ms: DOI:10.1007/s00265-013-1602-x ##
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
# # linmod <- lm(traps_per_life ~ caps_per_life, data=onepertag)
# # summary(linmod)
# #
# # AIC(linmod, powmod)
# # #yes, powmod is much better than linreg
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
# #pull measures of exploratory behavior
# exploratory <- onepertag_pred %>% select(c("year", "tag", "explore"))
# #join to netmets df
# netmets_puuv <- netmets_puuv %>% left_join(exploratory, by=c("year", "tag"))
# 
# ############################################################################



######### lagged degree and seroconvert status ###############
# add previous (lagged) degree (degree from previous month influences current PUUV status)
# add 0,1 for serovert - animals that go 0-0 or 0-1
    # BUT! the previous month has to be in the same year (eg don't want 2021 fall to influence 2022 spring)
netmets_puuv <- netmets_puuv %>% group_by(year, tag) %>%
  arrange(month, .by_group = TRUE) %>%
  mutate(prev_wt.deg.in = lag(wt.deg.in, n=1),
         prev_bin.in.deg = lag(bin.in.deg, n=1),
         prev_F.deg = lag(F.deg, n=1),
         prev_M.deg = lag(M.deg, n=1),
         prev_b.deg = lag(b.deg, n=1),
         prev_nb.deg = lag(nb.deg, n=1),
         prev_mb.deg = lag(mb.deg, n=1),
         prev_mnb.deg = lag(mnb.deg, n=1),
         prev_fb.deg = lag(fb.deg, n=1),
         prev_fnb.deg = lag(fnb.deg, n=1),
         prev_n.node = lag(n.node, n=1),
         prev_month = lag(month, n=1)) %>%
  mutate(prev_curr_puuv = paste(lag(puuv_ifa), puuv_ifa, sep="-")) %>%
  mutate(serovert = as.factor(case_when(prev_curr_puuv == "0-0" ~ 0,
                                        prev_curr_puuv == "0-1" ~ 1,
                                        prev_curr_puuv == "1-1" ~ NA)))
  # %>% select(!prev_curr_puuv)


#save all the above code to here
saveRDS(netmets_puuv, here("netmets_puuv_11.11.24.rds"))

# #load netmets_puuv
# netmets_puuv <- readRDS(here("netmets_puuv_11.11.24.rds"))



############ HAS NOT BEEN UPDATED WITH 2023 DATA #################

# ##################################################################################################################
# 
# ### Additional Summarized Data to Report in Manuscript ###
# 
# ##TOTAL number of animals captured, recapped etc in 2021 and 2022 (MAY-OCTOBER)
#   # before detailing the subset of the data used for this specific study
# 
# ## so that would be the fulltrap_ALL_03.04.24 dataset (see 01_...data_cleaning.R)
# #INCLUDING MAY captures, INCLUDING sex/repro=NA
# fulltrap21_ALL <- readRDS(file="fulltrap21_ALL_03.04.24.rds")
# fulltrap22_ALL <- readRDS(file="fulltrap22_ALL_03.04.24.rds")
# fulltrap21.22_ALL <- rbind(fulltrap21_ALL, fulltrap22_ALL)
# #4286 capture events (including WR) across MAY-Oct 2021 and 2022
# n_distinct(fulltrap21.22_ALL$tag) #1518 unique voles (including MAY)
# recapped <- fulltrap21.22_ALL %>% filter(caps_per_life>1)
# n_distinct(recapped$tag) #929 voles were recapped at least once
# ## number captured in both years (how many overwintered 2021->2022)
# ow <- fulltrap21.22_ALL %>% group_by(year, tag) %>% slice(1) %>% ungroup() %>%
#   group_by(tag) %>% mutate(n=length(tag)) %>% filter(n>1) %>% arrange(tag)
# n_distinct(ow$tag) #12 voles in both 2021 and 2022
# 
# ##-----------------
# 
# ##DATA USED FOR VARIOUS STEPS OF THE ANALYSIS IN THE MANUSCRIPT:
# 
# #entries per year (in 'netmets' files --> i.e., those used to construct the spatial overlap networks) 
#   # these entries may or may not have PUUV data
# nrow(netmets21) #1129 entries in 2021
# n_distinct(netmets21$tag) #742 voles in 2021
# nrow(netmets22) #1131 entries in 2022
# n_distinct(netmets22$tag) #744 voles in 2022
# 
# n_distinct(netmets21.22$tag) #1476 unique voles in the networks, 10 are in both 2021 and 2022
# 
# ## "Spatial overlap networks were constructed with all voles captured in June-October with recorded sex 
# ## and reproductive status data: 1129 captures (742 unique voles) in 2021 and 1131 captures (744 unique voles) in 2022."
# 
# ## THIS IS NO LONGER INCLUDED IN THE MANUSCRIPT:
# # #entries per year (in netmets_puuv --> i.e., those that do have PUUV data)
# #   # i.e., the above counts minus the 170 entries removed for lacking or inconsistent PUUV data
# # y1 <- netmets_puuv %>% filter(year=="2021") #1029 in 2021
# # n_distinct(y1$tag) #683 unique voles
# # y2 <- netmets_puuv %>% filter(year=="2022") #1061 in 2022
# # n_distinct(y2$tag) #694 unique voles
# # 
# # ## "Of these, we obtained consistent hantavirus serology results for 1,029 captures (683 unique voles) in June-October 
# # ## 2021 and 1,061 captures (694 unique voles) in June-October 2022."
# 
# ##---------------
# 
# #one last thing...
# #how many times is the "previous network position" (used as the explanatory variable in the hanta probability modesl)
#   #NOT from the immediately prior month?
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
# dfsmol <- df %>% filter(dif > 1) 
# #25 times we had to pull from an earlier month that wasn't the previous one (but was still in that calendar year)
# #25/713 = 3.5%
# 
# ############################################ The End #############################################
