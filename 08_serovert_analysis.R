# NEW! analysis of voles that seroconvert, using the NEW! percent overlap measures for network edges



#load libraries
library(here)
library(tidyverse)
library(ggbeeswarm)

#clear environment
rm(list = ls())


#load data

#load full netmets_puuv
netmets_puuv <- readRDS(here("netmets_puuv_11.07.24.rds"))

# #for wanelik/farine overlaps
# netmets_puuv <- readRDS(here("netmets_puuv_03.08.24.rds"))

##-------------------------------------------

# #subset for voles that don't serovert (but do have previous capture)
# novert <- netmets_puuv %>% filter(serovert==0) %>%
#   select(!c(wt.deg.in, bin.in.deg, F.deg, M.deg, b.deg, nb.deg, mnb.deg, fnb.deg, mb.deg, fb.deg))
# 
# 
# 
# #subset to only voles that seroconvert
# serovert <- netmets_puuv %>% filter(serovert==1) %>%
#   select(!c(wt.deg.in, bin.in.deg, F.deg, M.deg, b.deg, nb.deg, mnb.deg, fnb.deg, mb.deg, fb.deg))
# 
# #plot previous degree for seroverts
# serovert %>% ggplot(aes(x=prev_wt.deg.in, fill=trt)) + 
#   geom_histogram() +
#   facet_wrap(~trt)
# 
# #lol there are so many voles with 0 overlaps in the previous month that seroconvert *dies*


dat <- netmets_puuv %>% drop_na(serovert)

#previous weighted in degree
dat %>% 
  filter(prev_wt.deg.in>0) %>% #remove voles with 0 overlaps in previous month
  ggplot(aes(x=trt, y=prev_wt.deg.in, color=serovert)) +
  scale_color_manual(values=c("black", "red")) +
  geom_beeswarm() +
  stat_summary(aes(group = serovert, color=serovert), fun = mean, fun.min = mean, fun.max = mean,
    geom = "crossbar", width = 0.5, lwd = 0.4)

# #for wanelik/farine overlaps
# dat %>% 
#   filter(prev_wt.deg>0) %>% #remove voles with 0 overlaps in previous month
#   ggplot(aes(x=trt, y=prev_wt.deg, color=serovert)) +
#   scale_color_manual(values=c("black", "red")) +
#   geom_beeswarm() +
#   stat_summary(aes(group = serovert, color=serovert), fun = mean, fun.min = mean, fun.max = mean,
#                geom = "crossbar", width = 0.5, lwd = 0.4)

#previous binary degree (number of overlaps)
dat %>% 
  filter(prev_bin.in.deg>0) %>% #remove voles with 0 overlaps in previous month
  ggplot(aes(x=trt, y=prev_bin.in.deg, color=serovert)) +
  scale_color_manual(values=c("black", "red")) +
  geom_beeswarm() +
  stat_summary(aes(group = serovert, color=serovert), fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", width = 0.5, lwd = 0.4)
