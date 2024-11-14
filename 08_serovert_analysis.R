# NEW! analysis of voles that seroconvert, using the NEW! percent overlap measures for network edges



#load libraries
library(here)
library(tidyverse)
library(ggbeeswarm)

#clear environment
rm(list = ls())


#load data

#load full netmets_puuv
netmets_puuv <- readRDS(here("netmets_puuv_11.11.24.rds"))

# #for wanelik/farine overlaps
# netmets_puuv <- readRDS(here("netmets_puuv_03.08.24.rds"))

##-------------------------------------------

#subset for voles that don't serovert (but do have previous capture)
nonvert <- netmets_puuv %>% filter(serovert==0) %>%
  select(!c(wt.deg.in, bin.in.deg, F.deg, M.deg, b.deg, nb.deg, mnb.deg, fnb.deg, mb.deg, fb.deg))
#655 observations 
n_distinct(nonvert$tag) #490 individuals

nonvert %>% group_by(caps_per_year) %>%
  summarise(n=length(tag))


#subset to only voles that seroconvert
serovert <- netmets_puuv %>% filter(serovert==1) %>%
  select(!c(wt.deg.in, bin.in.deg, F.deg, M.deg, b.deg, nb.deg, mnb.deg, fnb.deg, mb.deg, fb.deg))
#82 observations (individuals)

serovert %>% group_by(caps_per_year) %>%
  summarise(n = length(tag))


###----------------- Visualize some stuff ------------------

#plot previous degree for seroverts
serovert %>% ggplot(aes(x=prev_wt.deg.in, fill=trt)) +
  geom_histogram() +
  facet_wrap(~trt)

#lol there are so many voles with 0 overlaps in the previous month that seroconvert *dies*


dat <- netmets_puuv %>% drop_na(serovert)

#previous weighted in-degree
dat %>% 
  filter(prev_wt.deg.in>0) %>% #remove voles with 0 overlaps in previous month
  ggplot(aes(x=prev_month, y=prev_wt.deg.in, color=serovert)) +
  scale_color_manual(values=c("black", "red")) +
  geom_beeswarm() +
  stat_summary(aes(group = serovert, color=serovert), fun = mean, fun.min = mean, fun.max = mean,
    geom = "crossbar", width = 0.5, lwd = 0.4) +
  facet_grid(year~trt)

#previous (wt) breeder-degree
dat %>% 
  filter(prev_b.deg>0) %>% #remove voles with 0 overlaps in previous month
  ggplot(aes(x=prev_month, y=prev_b.deg, color=serovert)) +
  scale_color_manual(values=c("black", "red")) +
  geom_beeswarm() +
  stat_summary(aes(group = serovert, color=serovert), fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", width = 0.5, lwd = 0.4) +
  facet_grid(year~trt)

#previous (wt) male-degree
dat %>% 
  filter(prev_M.deg>0) %>% #remove voles with 0 overlaps in previous month
  ggplot(aes(x=prev_month, y=prev_M.deg, color=serovert)) +
  scale_color_manual(values=c("black", "red")) +
  geom_beeswarm() +
  stat_summary(aes(group = serovert, color=serovert), fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", width = 0.5, lwd = 0.4) +
  facet_grid(year~trt)

#previous (wt) male breeder-degree
dat %>% 
  filter(prev_mb.deg>0) %>% #remove voles with 0 overlaps in previous month
  ggplot(aes(x=prev_month, y=prev_mb.deg, color=serovert)) +
  scale_color_manual(values=c("black", "red")) +
  geom_beeswarm() +
  stat_summary(aes(group = serovert, color=serovert), fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", width = 0.5, lwd = 0.4) +
  facet_grid(year~trt)

#previous (wt) female-degree
dat %>% 
  filter(prev_F.deg>0) %>% #remove voles with 0 overlaps in previous month
  ggplot(aes(x=prev_month, y=prev_F.deg, color=serovert)) +
  scale_color_manual(values=c("black", "red")) +
  geom_beeswarm() +
  stat_summary(aes(group = serovert, color=serovert), fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", width = 0.5, lwd = 0.4) +
  facet_grid(year~trt)

#previous (wt) female breeder-degree
dat %>% 
  filter(prev_fb.deg>0) %>% #remove voles with 0 overlaps in previous month
  ggplot(aes(x=prev_month, y=prev_fb.deg, color=serovert)) +
  scale_color_manual(values=c("black", "red")) +
  geom_beeswarm() +
  stat_summary(aes(group = serovert, color=serovert), fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", width = 0.5, lwd = 0.4) +
  facet_grid(year~trt)


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
  ggplot(aes(x=prev_month, y=prev_bin.in.deg, color=serovert)) +
  scale_color_manual(values=c("black", "red")) +
  geom_beeswarm() +
  stat_summary(aes(group = serovert, color=serovert), fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", width = 0.5, lwd = 0.4) +
  facet_grid(year~trt)


######--------------------------------------------------------------

#standard error (std dev of a population)
SE <- function(x) sd(x)/sqrt(length(x))

#average degree measures (+SE) per year/trt/month of nonverts
nonvert_means.month <- nonvert %>% 
  group_by(year, trt, month) %>%
  summarise(prev.in = mean(prev_wt.deg.in), se.prev.in = SE(prev_wt.deg.in),
            prev.bin.in = mean(prev_bin.in.deg), se.prev.bin.in = SE(prev_bin.in.deg),
            prev.M = mean(prev_M.deg), se.prev.M = SE(prev_M.deg),
            prev.b = mean(prev_b.deg), se.prev.b = SE(prev_b.deg),
            prev.mb = mean(prev_mb.deg), se.prev.mb = SE(prev_mb.deg),
            serovert = factor(0, levels=c("0", "1")))

#join nonvert means to serovert (observed) data
#https://stackoverflow.com/questions/24148423/r-ggplot-with-two-series-points-and-errorbars-with-legends
plotdata <- serovert %>% select(c(year, trt, month, site, tag, sex, season_breeder, 
                                  prev_wt.deg.in, prev_bin.in.deg, prev_M.deg, prev_mb.deg, prev_b.deg, serovert)) %>%
  rename(prev.in = prev_wt.deg.in,
         prev.bin.in = prev_bin.in.deg,
         prev.M = prev_M.deg,
         prev.b = prev_b.deg,
         prev.mb = prev_mb.deg) %>% #rename columns for easier plotting
  rbind(nonvert_means.month) 

## visualize ##

#previous weighted in degree
plotdata %>% ggplot(aes(x=month)) +
  geom_point(aes(y=prev.in, color=serovert, shape=serovert),
             position = position_jitterdodge(0.2, dodge.width = .2)) +
  geom_errorbar(aes( ymin = (prev.in-se.prev.in), 
                     ymax = (prev.in + se.prev.in))) +
  scale_shape_manual(name = "Sero Status", values=c(1,4)) +
  scale_colour_manual(name = "Sero Status", values=c("black", "red")) +
  labs(title="Previous Weighted In Degree") +
  facet_grid(year ~ trt)

#previous binary in degree
plotdata %>% ggplot(aes(x=month)) +
  geom_point(aes(y=prev.bin.in, color=serovert, shape=serovert),
             position = position_jitterdodge(0.2, dodge.width = .2)) +
  geom_errorbar(aes( ymin = (prev.bin.in-se.prev.bin.in), 
                     ymax = (prev.bin.in + se.prev.bin.in))) +
  scale_shape_manual(name = "Sero Status", values=c(1,4)) +
  scale_colour_manual(name = "Sero Status", values=c("black", "red")) +
  labs(title="Previous Binary In Degree") +
  facet_grid(year ~ trt)

#previous male degree
plotdata %>% ggplot(aes(x=month)) +
  geom_point(aes(y=prev.M, color=serovert, shape=serovert),
             position = position_jitterdodge(0.2, dodge.width = .2)) +
  geom_errorbar(aes( ymin = (prev.M-se.prev.M), 
                     ymax = (prev.M + se.prev.M))) +
  scale_shape_manual(name = "Sero Status", values=c(1,4)) +
  scale_colour_manual(name = "Sero Status", values=c("black", "red")) +
  labs(title="Previous Male Degree") +
  facet_grid(year ~ trt)

#previous breeder degree
plotdata %>% ggplot(aes(x=month)) +
  geom_point(aes(y=prev.b, color=serovert, shape=serovert),
             position = position_jitterdodge(0.2, dodge.width = .2)) +
  geom_errorbar(aes( ymin = (prev.b-se.prev.b), 
                     ymax = (prev.b + se.prev.b))) +
  scale_shape_manual(name = "Sero Status", values=c(1,4)) +
  scale_colour_manual(name = "Sero Status", values=c("black", "red")) +
  labs(title="Previous Breeder Degree") +
  facet_grid(year ~ trt)

#previous male breeder degree
plotdata %>% ggplot(aes(x=month)) +
  geom_point(aes(y=prev.mb, color=serovert, shape=serovert),
             position = position_jitterdodge(0.2, dodge.width = .2)) +
  geom_errorbar(aes( ymin = (prev.mb-se.prev.mb), 
                     ymax = (prev.mb + se.prev.mb))) +
  scale_shape_manual(name = "Sero Status", values=c(1,4)) +
  scale_colour_manual(name = "Sero Status", values=c("black", "red")) +
  labs(title="Previous Male Breeder Degree") +
  facet_grid(year ~ trt)
