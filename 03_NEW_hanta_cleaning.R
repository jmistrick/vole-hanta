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



##########---------- LOAD NETWORK METRICS and VOLE CAPTURE METADATA -----------###########

### Fulltrap and netmets dfs created separately, load 2021 and 2022 data here, clean as needed

### INDIVIDUAL VOLE METADATA by year
    # >> FROM OTHER R PROJECT! << (hence why I'm not using 'here()')
fulltrap21 <- readRDS(file="../vole-spatial-v2/fulltrap21_05.10.23.rds") #go up a level from current wd, then down to file
fulltrap22 <- readRDS(file="../vole-spatial-v2/fulltrap22_05.10.23.rds") 

###### THESE VERSIONS OF FULLTRAP do have all the animals (breeders and nonbreeders) but as of 5/11 I only constructed networks for breeders
    ## it doesn't really matter since netmets files only have data on the breeders and I left_join the trapping data to netmets

## >>NOTE<< overwintered animals may incorrectly have firstcap==1 in 2022
    # PIT tag: 21895 Occ2Sess1 at Vaarinkorpi
    # PIT tag: 226280 Occ2Sess2 at Kuoppa
# FULLtrap <- rbind(fulltrap21, fulltrap22)
# OW <- FULLtrap %>% 
#   group_by(tag) %>% arrange(year, occ.sess, .by_group = TRUE) %>%
#   filter(n_distinct(year) > 1)
# write.csv(OW, here("overwinter21-22.csv"))

#join 21 and 22 fulltrap data; correct firstcap for 2022
fulltrap21.22 <- rbind(fulltrap21, fulltrap22) %>%
  group_by(tag) %>% arrange(tag, year, occ.sess, .by_group = TRUE) %>%
  mutate(firstcap = ifelse(date_time == min(date_time), 1, 0)) %>%
  mutate(firstcap = factor(firstcap)) %>%
#remove may; keep one entry per tag,year,month; pull only relevant columns
  filter(month!="may") %>% 
  ungroup() %>% group_by(tag, year, month) %>%
  slice(1) %>%
  select(year, site, trt, month, tag, sex, samp_id, traps_per_life, caps_per_life) 

### NETWORK METRICS
netmets21 <- readRDS(file="../vole-spatial-v2/netmets21.rds") %>%
  mutate(year=as.numeric(2021))
netmets22 <- readRDS(file="../vole-spatial-v2/netmets22.rds") %>%
  mutate(year=as.numeric(2022))

netmets21.22 <- rbind(netmets21, netmets22)

### NETWORK METRICS + METADATA

netmets_full <- left_join(netmets21.22, fulltrap21.22, by=c("tag", "site", "year", "month")) %>%
  mutate(site=as.factor(site),
         year=as.factor(year),
         month = as.factor(month),
         month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>%
  select(!c(focal_id)) %>% #remove duplicate column
  relocate(c(year, site, trt, month, n.node, avg.wt.deg, tag, sex, samp_id), .before = wt.deg)

# saveRDS(netmets_full, here("netmets_full_05.10.23.rds"))

# #visualize degree by month for trt and year
# netmets_full %>%
#   ggplot(aes(x=month, y=wt.deg, fill=trt)) +
#   geom_violin() +
#   facet_wrap(~year, nrow=2)



#########################################   LOAD & CLEAN PUUV IFA DATA   ########################################

#load, clean, format PUUV IFA data
puuv_data <- read.csv(here("puuv_ifa.csv")) %>%
  clean_names %>%
  #populate column of FINAL PUUV status (result of second run if two runs were done, else result of first run)
  mutate(FINAL_puuv = ifelse(is.na(puuv_confirm), as.character(puuv_initial), as.character(puuv_confirm))) %>%
  mutate(samp_id = as.numeric(id),
         date_initial = as_date(date_initial, format= "%m/%d/%Y"),
         puuv_initial = as.factor(puuv_initial),
         date_confirm = as_date(date_confirm, format= "%m/%d/%Y"),
         puuv_confirm = as.factor(puuv_confirm),
         FINAL_puuv = as.factor(FINAL_puuv)) %>%
  drop_na(FINAL_puuv) %>%
  dplyr::select(FINAL_puuv, samp_id) %>%
  rename(puuv_ifa = FINAL_puuv)
#output is a df with two columns, sample ID and PUUV status (0,1)

# Save puuv_data to a rdata file
# saveRDS(puuv_data, file = here("hantadata_02.28.23.rds"))

# Restore puuv_data from the rdata file
# puuv_data <- readRDS(file = "hantadata_02.28.23.rds")

################################## LOAD NETMETS data and COMBINE with PUUV data ##########################################

#join puuv_data to netmets_full (BOTH YEARS!)
netmets_puuv <- left_join(netmets_full, puuv_data, by="samp_id") %>% 
  #left join on netmets_full because I need to have network data
  relocate(puuv_ifa, .after="samp_id") %>%
  select(!samp_id) %>%
  drop_na(puuv_ifa) #drop any animals without puuv data

######## THIS SEEMS EXTREME: there is PUUV data on 2259 samples -- I have network data on 2277 animals
######## But then when we combine these, there are only 2149 animals with both PUUV and network data
    ##### which means we're missing network data for 128 animals that we have PUUV data for
    ##### which seems hard to believe when I only cut out may networks and uusi site


############ REMOVE THE VOLES from net_mets_puuv THAT SEROCONVERT POS TO NEG ##################

#voles that seroconvert PUUV+ to PUUV-
puuv_pos_neg <- netmets_puuv %>% group_by(tag) %>% arrange(year, month, .by_group = TRUE) %>%
  dplyr::select(year, month, tag, puuv_ifa) %>%
  summarise(status_time = toString(puuv_ifa)) %>%
  filter(str_detect(status_time, "1,\\s0")) #filter for animals that go from pos to neg
#pull the PIT tags
problemchildren <- puuv_pos_neg$tag
#filter netmetsPUUV to remove 'problemchildren'
netmets_puuv <- netmets_puuv %>%
  filter(!tag %in% problemchildren)

## NOW NETMETS_PUUV has ONLY the [BREEDING] animals that:  ##
    #       1. Have network data
    #       2. Have PUUV IFA data
    #       3. DON'T seroconvert from PUUV+ to PUUV-


###############################################################################
### exploratoriness of breeders ###
## based off of VanderWaal ground squirrel ms ##
#https://modelr.tidyverse.org/reference/add_predictions.html

onepertag <- netmets_puuv %>%
  group_by(year, tag) %>% slice(1)

onepertag %>% ggplot(aes(x=caps_per_life, y=traps_per_life)) +
  geom_point()

powmod <- lm(log(traps_per_life) ~ log(caps_per_life), data=onepertag)
# summary(powmod) #check it
# ln(y) = -0.001924 + 0.742403ln(x)
# y = 0.9980778x^0.742403

predicted.values <- predict(powmod, type="response")

# #visualize it
# ggplot(aes(x=log(caps_per_life), y=log(traps_per_life), color=sex), data=onepertag) +
#   geom_jitter(alpha=0.5, size=1.5) +
#   geom_line(aes(y=predicted.values), linetype=2, color="black") +
#   labs(title="exploratoriness")

d <- data.frame(onepertag$tag, onepertag$caps_per_life, onepertag$traps_per_life) %>% 
  rename(caps_per_life = onepertag.caps_per_life,
         traps_per_life = onepertag.traps_per_life,
         tag = onepertag.tag)
d <- d %>% modelr::add_predictions(powmod) %>%
  mutate(pred_exp = exp(pred))
onepertag_pred <- onepertag %>% left_join(d, by=c("tag", "caps_per_life", "traps_per_life")) %>%
  mutate(explore = log(traps_per_life) - pred,
         explore_exp = traps_per_life - pred_exp) #residual of obs traps_per_life - expected

### EXPLORE is in terms of log(traps_per_life)  [ maybe this is scaled -1 to 1 because it's a ln? ]
### EXPLORE_EXP is in terms of traps_per_life (so how many more/fewer traps were you in than expected)

# #visualize
# onepertag_pred %>%
#   ggplot(aes(x=sex, y=explore, fill=sex)) +
#   geom_violin() +
#   geom_hline(yintercept = 0, color="black")

# #yes, powmod is much better than linreg
# linmod <- lm(traps_per_life ~ caps_per_life, data=onepertag)
# summary(linmod)
# 
# AIC(linmod, powmod)


##### WHAT UNITS SHOULD THE RESIDUALS BE IN?? ########################

exploratory <- onepertag_pred %>% select(c("year", "tag", "explore"))

netmets_puuv <- netmets_puuv %>% left_join(exploratory, by=c("year", "tag"))

############################################################################








######################### end problem children ################################

# #pull tag and sample ID
# tag_id <- netmets_puuv %>% dplyr::select(tag, samp_id)


# add previous (lagged) degree (degree from previous month influences current PUUV status)
# add 0,1 for serovert - animals that go 0-0 or 0-1 
    # BUT! the previous month has to be in the same year (don't want 2021 fall to influence 2022 spring)
netmets_puuv <- netmets_puuv %>% group_by(year, tag) %>%
  arrange(month, .by_group = TRUE) %>%
  mutate(prev_wt.deg = lag(wt.deg, n=1),
         prev_F.deg = lag(F.deg, n=1),
         prev_M.deg = lag(M.deg, n=1),
         prev_n.node = lag(n.node, n=1),
         prev_month = lag(month, n=1)) %>%
  mutate(prev_curr_puuv = paste(lag(puuv_ifa), puuv_ifa, sep="-")) %>%
  mutate(serovert = as.factor(case_when(prev_curr_puuv == "0-0" ~ 0,
                                        prev_curr_puuv == "0-1" ~ 1,
                                        prev_curr_puuv == "1-1" ~ NA))) 
  # %>% select(!prev_curr_puuv)

##################################################################################################################








############### PUUV PREVALENCE PER SITE #########################

#what is the prevalence of hanta at each site in each month?
puuv_prev <- netmets_puuv %>% group_by(year, site, month) %>%
  summarise(n = length(tag),
            n_pos = length(which(puuv_ifa == "1")),
            prev = n_pos/n) %>%
  mutate(site = factor(site, levels=c("asema", "helmipollo", "hevonen",
                                      "ketunpesa", "kiirastuli", "mustikka",
                                      "kuoppa", "radio", "vaarinkorpi",
                                      "janakkala", "luostari", "puro", "talo")))

#plot prevalence each month, looking for sites with multiple prev=0 in a row
    #--> these should probably be removed for analysis
puuv_prev %>%
  ggplot() +
  geom_point(aes(x=month, y=prev, color=prev==0)) +
  scale_color_manual(name="prevalence = 0",
                     values=setNames(c("red","black"),c(T,F)))+
  facet_grid(vars(year), vars(site)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


### ERROR BARS based on sample size?
## essentially something to show if it's 100% prev of 2 animals or 45 animals

#https://www.rdocumentation.org/packages/epiR/versions/0.9-79/topics/epi.conf
library(epiR)
#Method prevalence require a two-column matrix; the first column specifies the number of positives, 
    #the second column specifies the total number tested. 
epi.conf.data <- puuv_prev %>% ungroup() %>% select(n_pos, n) %>% as.matrix()
prevCI <- epi.conf(epi.conf.data, ctype="prevalence", method="exact", conf.level = 0.95)
puuv_prev <- cbind(puuv_prev, prevCI)


puuv_prev %>%
  ggplot(aes(x=month, y=prev, color=prev==0)) +
  geom_point() +
  geom_errorbar(aes(y=est, ymin=lower, ymax=upper), width=.2, alpha=0.3) +
  scale_color_manual(name="prevalence = 0",
                     values=setNames(c("red","black"),c(T,F))) +
  facet_grid(vars(year), vars(site)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## visualize infected as a percentage of the total population

puuv_prev_stack <- netmets_puuv %>% 
  mutate(month.n = case_when(month=="june" ~ 6,
                           month=="july" ~ 7,
                           month=="aug" ~ 8,
                           month=="sept" ~ 9,
                           month=="oct" ~ 10),
         month.n = as.numeric(month.n)) %>%
  mutate(site = factor(site, levels=c("asema", "helmipollo", "hevonen",
                                      "ketunpesa", "kiirastuli", "mustikka",
                                      "kuoppa", "radio", "vaarinkorpi",
                                      "janakkala", "luostari", "puro", "talo"))) %>%
  group_by(year, site, month.n) %>%
  summarise(n = length(tag),
            n_pos = length(which(puuv_ifa == "1"))) %>%
  pivot_longer(-c("year", "month.n", "site"), names_to = "group", values_to = "count")

puuv_prev_stack %>%
  ggplot(aes(x=month.n, y=count, fill=group)) +
  geom_area() +
  facet_grid(vars(year), vars(site))


############################## end ####################################





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
# 
# ######## the right angle shape is curious to me, and the fact that 
#   #males go vertically and females horizontally
#   #I don't know if this is an artifact or something cool
#   #I also can't wrap my head around what it means... more on this?






#########################################################################

# #How does CURRENT degree affect infection status?
# 
# dat <- netmets_puuv
# 
# mod <- glmer(puuv_ifa ~ wt.deg + sex + trt + month + n.node + year + (1|site),
#              family=binomial, data=dat)
# 
# mod <- glmer(puuv_ifa ~ F.deg:sex + M.deg:sex + sex + trt + month + n.node + year + (1|site),
#              family=binomial, data=dat)
# #WARNING: Failed to converge
#     ## BEN BOLKER says it's okay, the bobyqa ("Bound Optimization BY Quadratic Approximation") isn't actually doing anything
#     ## FALSE POSITIVE 'failed to converge'
# 
# #fixed with : adding "control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))"
# #https://stats.stackexchange.com/questions/164457/r-glmer-warnings-model-fails-to-converge-model-is-nearly-unidentifiable
#     #also: https://www.learn-mlms.com/07-module-7.html#learning-objectives-5
# 
# summary(mod)
# 
# #RESULTS - Current network position does not correlate to infection status

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
mod <- glm(serovert ~ prev_M.deg:sex + prev_F.deg:sex + sex + explore + 
             trt + prev_month + prev_n.node + year,
             family=binomial, data=netmets_puuv_serov)
#SINGULARITY IS AN ISSUE to do glmer (1|site) - some sites 0 or 1 seroverter
summary(mod)
#previous degree only affects p|serovert for females:
    #more likely to serovert if female with high female degree
    #less likely to serovert if female with high male degree (... um what)

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


mod <- glmer(puuv_ifa ~ Previous_M.degree:Sex + Previous_F.degree:Sex + Sex +
               Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
             family=binomial, data=dat)
#WARNING: Failed to converge
summary(mod)

mod_exp <- glmer(puuv_ifa ~ Previous_M.degree:Sex + Previous_F.degree:Sex + Sex + explore +
                Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
              family=binomial, data=dat)
summary(mod_exp)

AIC(mod, mod_exp) #eeee the AIC is nearly 2 smaller without exploratoriness (suggests more complex model isn't better)

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

mod1 <- glmer(puuv_ifa ~ Previous_M.degree:Sex + Previous_F.degree:Sex + Sex +
               Treatment + Previous_Month + Previous_Network_Size + Year + (1|site),
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
simulationOutput <- simulateResiduals(fittedModel = mod1)
plot(simulationOutput) #qq plot and residual vs fitted
testDispersion(simulationOutput) #formal test for overdispersion
# testZeroInflation(simulationOutput) #formal test for zero inflation (common type of overdispersion)

plotResiduals(simulationOutput, dat$Sex)
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






netmets_puuv %>% drop_na(prev_wt.deg) %>%
  ggplot(aes(y=prev_M.deg, x=month, fill=puuv_ifa)) +
  geom_violin() +
  facet_grid(trt ~ year)







#########################################################################################################
######### How does CURRENT network position (wt.deg) affect CURRENT infection status? ##########

# ## last month's degree, this month's PUUV status
#     ## grouped by year and treatment (so replicate sites are combined)
# netmets_puuv_wtdeg <- netmets_puuv %>% unite(yr_trt, year, trt) %>%
#   ungroup() %>% group_by(yr_trt)
# 
# 
# #https://stackoverflow.com/questions/46407320/data-frame-to-nested-list
# wtdeg.list <- lapply(split(netmets_puuv_wtdeg, netmets_puuv_wtdeg$yr_trt, drop = TRUE),
#                function(x) split(x, x[["month"]], drop = TRUE))
# 
# 
# ###### NOW I HAVE A nested list of all the sites in both years, separated to 2e order list elements by month
#   ## I need to loop through all of those, calculate a u-statistic and then permutation and get a pvalue
#   ## ideally 1,000 to 10,000 permutations WITHOUT it taking a millions years
#   ## I have permutation z-score network things that might have useful code
# 
# wtdeg_results.list <- list()
# 
# for(i in 1:length(wtdeg.list)){
# 
#   #for each year-trt combo (there are 8)
#   print(names(wtdeg.list[i]))
# 
#   yr_trt <- list()
# 
#   for(j in 1:length(wtdeg.list[[i]])){
# 
#     #for each month (of puuv status) within a year_trt combo (should be 4)
#     print(names(wtdeg.list[[i]][j]))
# 
# 
#     ### Mann-Whitney U test - done manually to calculate little u (Croft et al. 2008 pg110)
# 
#     # #to run it in R, use this (but you only get one of the U values)
#     # wilcox.test(wt.deg ~ puuv_ifa, data=data)
# 
#     data <- wtdeg.list[[i]][[j]]
# 
#     #To calculate the U statistic, start by sorting the data
#     data.sort <- data[order(data$wt.deg),]
# 
#     #Then assign the ranks and split into two groups
#     data.sort$deg.ranked <- rank(data.sort$wt.deg, ties.method = 'average')
#     pos <- filter(data.sort, puuv_ifa == "1")
#     neg <- filter(data.sort, puuv_ifa == "0")
# 
#     #With the ranks calculated, the sums can be found and the values of U1 & U2 computed
#     R1 <- sum(pos$deg.ranked)
#     # R2 <- sum(neg$deg.ranked)
# 
#     n1 <- length(pos$deg.ranked)
#     n2 <- length(neg$deg.ranked)
# 
#     U1 <- n1 * n2 + n1 * (n1 + 1) / 2 - R1
#     # U2 <- n1 * n2 + n1 * (n2 + 1) / 2 - R2
# 
#     u_obs <- U1 / (n1*n2)
# 
# 
# 
#     # montecarlo - https://www.countbayesie.com/blog/2015/3/3/6-amazing-trick-with-monte-carlo-simulations
# 
#     mc.list <- list()
# 
#     runs = 9999 #remember, your obs data counts a permutation (runs+1 = 1000)
# 
#     for(k in 1:runs){
# 
#       shuffled.data <- transform(data, puuv_ifa = sample(puuv_ifa, replace=FALSE) )
# 
#       #### SOURCE: https://rpubs.com/aaronsc32/measure-cabbages-mann-whitney  #######
# 
#       #To calculate the U statistic, start by sorting the data
#       shuff.sort <- shuffled.data[order(shuffled.data$wt.deg),]
# 
#       #Then assign the ranks and split into two groups
#       shuff.sort$deg.ranked <- rank(shuff.sort$wt.deg, ties.method = 'average')
#       pos <- filter(shuff.sort, puuv_ifa == "1")
#       neg <- filter(shuff.sort, puuv_ifa == "0")
# 
#       #With the ranks calculated, the sums can be found and the values of U1 & U2 computed
#       R1 <- sum(pos$deg.ranked)
#       # R2 <- sum(neg$deg.ranked)
# 
#       n1 <- length(pos$deg.ranked)
#       n2 <- length(neg$deg.ranked)
# 
#       U1 <- n1 * n2 + n1 * (n1 + 1) / 2 - R1
#       # U2 <- n1 * n2 + n1 * (n2 + 1) / 2 - R2
# 
#       u_mc <- U1 / (n1*n2)
# 
#       #### end #####
# 
#       mc.list[[k]] <- u_mc
# 
#     }
# 
#     mc.vec <- unlist(mc.list) #take a list to a vector with unlist()
# 
#     mc.p.value <- (sum(mc.vec > u_obs) + 1) / (runs + 1)
#     mc.p.value.less <- (sum(mc.vec < u_obs) + 1) / (runs + 1)
# 
# 
#     #dataframe to hold results per month
#     month_ID <- as.vector(names(wtdeg.list[[i]][j]))
#     yr_trt[[j]] <- data.frame(month_ID)
# 
#     yr_trt[[j]]$u <- u_obs
#     yr_trt[[j]]$p.value <- mc.p.value
#     yr_trt[[j]]$p.value.less <- mc.p.value.less
# 
#   }
# 
#   wtdeg_results.list[[i]] <- yr_trt
# 
# }
# 
# 
# #name the 8 1st order elements of lagdeg_results.list
# names(wtdeg_results.list) <- names(wtdeg.list)
# 
# #rename the sublist items (months) for each site
# #accounting for the fact that most sites have 5 months of data but some have 4
# for(i in 1:length(wtdeg_results.list)){
# names(wtdeg_results.list[[i]]) <- c("july", "aug", "sept", "oct")
# }
# 
# 
# 
# #collate results
# #this collapses the 2nd order elements down to the 1st order element
# 
# #make a list to store things
# wtdeg_results.summ <- list()
# 
# #loop across all sites and collapse the 2e dfs into 1 df per 1e level
# for(i in 1:length(wtdeg_results.list)){
#   #for all 12 sites
#   summary <- do.call("rbind", wtdeg_results.list[[i]])
#   wtdeg_results.summ[[i]] <- summary
# }
# 
# #remove the row names of the dfs (its the month, but we already have that info)
# for(i in 1:length(wtdeg_results.summ)){
#   row.names(wtdeg_results.summ[[i]]) <- NULL
# }
# 
# #name the 1st order elements
# names(wtdeg_results.summ) <- names(wtdeg_results.list)
# 
# ## make a freiggein huge df
# wtdeg_results <- do.call(rbind.data.frame, wtdeg_results.summ)
# 
# saveRDS(wtdeg_results, here("wtdeg_breeder_results_05.12.23.RDS"))
# 
# saveRDS(wtdeg_results.list, here("wtdeg_breeder_results.list_05.12.23.RDS"))
# 
# saveRDS(wtdeg_results.summ, here("wtdeg_breeder_results.summ_05.12.23.RDS"))
# 
# saveRDS(wtdeg.list, here("wtdeg_breeder.list_05.12.23.RDS"))


##########################################  END  ########################################################
#########################################################################################################








#########################################################################################################
######### How does previous month's network position (wt.deg) affect current infection status? ##########

# ## last month's degree, this month's PUUV status
#     ## grouped by year and treatment (so replicate sites are combined)
# netmets_puuv_lag <- netmets_puuv %>% unite(yr_trt, year, trt) %>%
#   ungroup() %>% group_by(yr_trt) %>%
#   drop_na(prev_wt.deg)
# 
# 
# #https://stackoverflow.com/questions/46407320/data-frame-to-nested-list
# lagdeg.list <- lapply(split(netmets_puuv_lag, netmets_puuv_lag$yr_trt, drop = TRUE),
#                function(x) split(x, x[["month"]], drop = TRUE))
# 
# 
# ###### NOW I HAVE A nested list of all the sites in both years, separated to 2e order list elements by month
#   ## I need to loop through all of those, calculate a u-statistic and then permutation and get a pvalue
#   ## ideally 1,000 to 10,000 permutations WITHOUT it taking a millions years
#   ## I have permutation z-score network things that might have useful code
# 
# lagdeg_results.list <- list()
# 
# for(i in 1:length(lagdeg.list)){
# 
#   #for each year-trt combo (there are 8)
#   print(names(lagdeg.list[i]))
# 
#   yr_trt <- list()
# 
#   for(j in 1:length(lagdeg.list[[i]])){
# 
#     #for each month (of puuv status) within a year_trt combo (should be 4)
#     print(names(lagdeg.list[[i]][j]))
# 
# 
#     ### Mann-Whitney U test - done manually to calculate little u (Croft et al. 2008 pg110)
# 
#     # #to run it in R, use this (but you only get one of the U values)
#     # wilcox.test(wt.deg ~ puuv_ifa, data=data)
# 
#     data <- lagdeg.list[[i]][[j]]
# 
#     #To calculate the U statistic, start by sorting the data
#     data.sort <- data[order(data$prev_wt.deg),]
# 
#     #Then assign the ranks and split into two groups
#     data.sort$prev_deg.ranked <- rank(data.sort$prev_wt.deg, ties.method = 'average')
#     pos <- filter(data.sort, puuv_ifa == "1")
#     neg <- filter(data.sort, puuv_ifa == "0")
# 
#     #With the ranks calculated, the sums can be found and the values of U1 & U2 computed
#     R1 <- sum(pos$prev_deg.ranked)
#     # R2 <- sum(neg$deg.ranked)
# 
#     n1 <- length(pos$prev_deg.ranked)
#     n2 <- length(neg$prev_deg.ranked)
# 
#     U1 <- n1 * n2 + n1 * (n1 + 1) / 2 - R1
#     # U2 <- n1 * n2 + n1 * (n2 + 1) / 2 - R2
# 
#     u_obs <- U1 / (n1*n2)
# 
# 
# 
#     # montecarlo - https://www.countbayesie.com/blog/2015/3/3/6-amazing-trick-with-monte-carlo-simulations
# 
#     mc.list <- list()
# 
#     runs = 9999 #remember, your obs data counts a permutation (runs+1 = 1000)
# 
#     for(k in 1:runs){
# 
#       shuffled.data <- transform(data, puuv_ifa = sample(puuv_ifa, replace=FALSE) )
# 
#       #### SOURCE: https://rpubs.com/aaronsc32/measure-cabbages-mann-whitney  #######
# 
#       #To calculate the U statistic, start by sorting the data
#       shuff.sort <- shuffled.data[order(shuffled.data$prev_wt.deg),]
# 
#       #Then assign the ranks and split into two groups
#       shuff.sort$prev_deg.ranked <- rank(shuff.sort$prev_wt.deg, ties.method = 'average')
#       pos <- filter(shuff.sort, puuv_ifa == "1")
#       neg <- filter(shuff.sort, puuv_ifa == "0")
# 
#       #With the ranks calculated, the sums can be found and the values of U1 & U2 computed
#       R1 <- sum(pos$prev_deg.ranked)
#       # R2 <- sum(neg$deg.ranked)
# 
#       n1 <- length(pos$prev_deg.ranked)
#       n2 <- length(neg$prev_deg.ranked)
# 
#       U1 <- n1 * n2 + n1 * (n1 + 1) / 2 - R1
#       # U2 <- n1 * n2 + n1 * (n2 + 1) / 2 - R2
# 
#       u_mc <- U1 / (n1*n2)
# 
#       #### end #####
# 
#       mc.list[[k]] <- u_mc
# 
#     }
# 
#     mc.vec <- unlist(mc.list) #take a list to a vector with unlist()
# 
#     mc.p.value <- (sum(mc.vec > u_obs) + 1) / (runs + 1)
# 
# 
#     #dataframe to hold results per month
#     month_ID <- as.vector(names(lagdeg.list[[i]][j]))
#     yr_trt[[j]] <- data.frame(month_ID)
# 
#     yr_trt[[j]]$u <- u_obs
#     yr_trt[[j]]$p.value <- mc.p.value
# 
#   }
# 
#   lagdeg_results.list[[i]] <- yr_trt
# 
# }
# 
# 
# #name the 8 1st order elements of lagdeg_results.list
# names(lagdeg_results.list) <- names(lagdeg.list)
# 
# #rename the sublist items (months) for each site
# #accounting for the fact that most sites have 5 months of data but some have 4
# for(i in 1:length(lagdeg_results.list)){
# names(lagdeg_results.list[[i]]) <- c("july", "aug", "sept", "oct")
# }
# 
# 
# 
# #collate results
# #this collapses the 2nd order elements down to the 1st order element
# 
# #make a list to store things
# lagdeg_results.summ <- list()
# 
# #loop across all sites and collapse the 2e dfs into 1 df per 1e level
# for(i in 1:length(lagdeg_results.list)){
#   #for all 12 sites
#   summary <- do.call("rbind", lagdeg_results.list[[i]])
#   lagdeg_results.summ[[i]] <- summary
# }
# 
# #remove the row names of the dfs (its the month, but we already have that info)
# for(i in 1:length(lagdeg_results.summ)){
#   row.names(lagdeg_results.summ[[i]]) <- NULL
# }
# 
# #name the 1st order elements
# names(lagdeg_results.summ) <- names(lagdeg_results.list)
# 
# ## make a freiggein huge df
# lagdeg_results <- do.call(rbind.data.frame, lagdeg_results.summ)
# 
# saveRDS(lagdeg_results, here("lagdeg_breeder_results_05.12.23.RDS"))
# 
# saveRDS(lagdeg_results.list, here("lagdeg_breeder_results.list_05.12.23.RDS"))
# 
# saveRDS(lagdeg_results.summ, here("lagdeg_breeder_results.summ_05.12.23.RDS"))
# 
# saveRDS(lagdeg.list, here("lagdeg_breeder.list_05.12.23.RDS"))
# 
# 
# # ##### RESULTS: previous network position does not impact current infection status

##########################################  END  ########################################################
#########################################################################################################





######### 5.08 versions have breeder/nonbreeders
# Mdeg <- readRDS(here("Mdeg_results2_05.08.23.RDS"))
# Fdeg <- readRDS(here("Fdeg_results_05.08.23.RDS"))
# 
# results2 <- readRDS(here("Mdeg_results2_05.08.23.rds"))



#########################################################################################################
######### How does previous month's [MALE DEGREE] affect current infection status? ##########

# ## last month's MALE degree, this month's PUUV status
#     ## grouped by year and treatment (so replicate sites are combined)
# netmets_puuv_male <- netmets_puuv %>% unite(yr_trt, year, trt) %>%
#   ungroup() %>% group_by(yr_trt) %>%
#   drop_na(prev_M.deg)
# 
# 
# #https://stackoverflow.com/questions/46407320/data-frame-to-nested-list
# Mdeg.list <- lapply(split(netmets_puuv_male, netmets_puuv_male$yr_trt, drop = TRUE),
#                function(x) split(x, x[["month"]], drop = TRUE))
# 
# 
# ###### NOW I HAVE A nested list of all the sites in both years, separated to 2e order list elements by month
#   ## I need to loop through all of those, calculate a u-statistic and then permutation and get a pvalue
#   ## ideally 1,000 to 10,000 permutations WITHOUT it taking a millions years
#   ## I have permutation z-score network things that might have useful code
# 
# Mdeg_results.list <- list()
# 
# for(i in 1:length(Mdeg.list)){
# 
#   #for each year-trt combo (there are 8)
#   print(names(Mdeg.list[i]))
# 
#   yr_trt <- list()
# 
#   for(j in 1:length(Mdeg.list[[i]])){
# 
#     #for each month (of puuv status) within a year_trt combo (should be 4)
#     print(names(Mdeg.list[[i]][j]))
# 
# 
#     ### Mann-Whitney U test - done manually to calculate little u (Croft et al. 2008 pg110)
# 
#     # #to run it in R, use this (but you only get one of the U values)
#     # wilcox.test(wt.deg ~ puuv_ifa, data=data)
# 
#     data <- Mdeg.list[[i]][[j]]
# 
#     #To calculate the U statistic, start by sorting the data
#     data.sort <- data[order(data$prev_M.deg),]
# 
#     #Then assign the ranks and split into two groups
#     data.sort$prev_deg.ranked <- rank(data.sort$prev_M.deg, ties.method = 'average')
#     pos <- filter(data.sort, puuv_ifa == "1")
#     neg <- filter(data.sort, puuv_ifa == "0")
# 
#     #With the ranks calculated, the sums can be found and the values of U1 & U2 computed
#     R1 <- sum(pos$prev_deg.ranked)
#     # R2 <- sum(neg$deg.ranked)
# 
#     n1 <- length(pos$prev_deg.ranked)
#     n2 <- length(neg$prev_deg.ranked)
# 
#     U1 <- n1 * n2 + n1 * (n1 + 1) / 2 - R1
#     # U2 <- n1 * n2 + n1 * (n2 + 1) / 2 - R2
# 
#     u_obs <- U1 / (n1*n2)
# 
# 
# 
#     # montecarlo - https://www.countbayesie.com/blog/2015/3/3/6-amazing-trick-with-monte-carlo-simulations
# 
#     mc.list <- list()
# 
#     runs = 9999 #remember, your obs data counts a permutation (runs+1 = 1000)
# 
#     for(k in 1:runs){
# 
#       shuffled.data <- transform(data, puuv_ifa = sample(puuv_ifa, replace=FALSE) )
# 
#       #### SOURCE: https://rpubs.com/aaronsc32/measure-cabbages-mann-whitney  #######
# 
#       #To calculate the U statistic, start by sorting the data
#       shuff.sort <- shuffled.data[order(shuffled.data$prev_M.deg),]
# 
#       #Then assign the ranks and split into two groups
#       shuff.sort$prev_deg.ranked <- rank(shuff.sort$prev_M.deg, ties.method = 'average')
#       pos <- filter(shuff.sort, puuv_ifa == "1")
#       neg <- filter(shuff.sort, puuv_ifa == "0")
# 
#       #With the ranks calculated, the sums can be found and the values of U1 & U2 computed
#       R1 <- sum(pos$prev_deg.ranked)
#       # R2 <- sum(neg$deg.ranked)
# 
#       n1 <- length(pos$prev_deg.ranked)
#       n2 <- length(neg$prev_deg.ranked)
# 
#       U1 <- n1 * n2 + n1 * (n1 + 1) / 2 - R1
#       # U2 <- n1 * n2 + n1 * (n2 + 1) / 2 - R2
# 
#       u_mc <- U1 / (n1*n2)
# 
#       #### end #####
# 
#       mc.list[[k]] <- u_mc
# 
#     }
# 
#     mc.vec <- unlist(mc.list) #take a list to a vector with unlist()
# 
#     mc.p.value <- (sum(mc.vec > u_obs) + 1) / (runs + 1)
#     mc.p.value.less <- (sum(mc.vec < u_obs) + 1) / (runs + 1)
# 
# 
#     #dataframe to hold results per month
#     month_ID <- as.vector(names(Mdeg.list[[i]][j]))
#     yr_trt[[j]] <- data.frame(month_ID)
# 
#     yr_trt[[j]]$u <- u_obs
#     yr_trt[[j]]$p.value <- mc.p.value
#     yr_trt[[j]]$p.val.less <- mc.p.value.less
# 
#   }
# 
#   Mdeg_results.list[[i]] <- yr_trt
# 
# }
# 
# 
# #name the 8 1st order elements of Mdeg_results.list
# names(Mdeg_results.list) <- names(Mdeg.list)
# 
# #rename the sublist items (months) for each site
# #accounting for the fact that most sites have 5 months of data but some have 4
# for(i in 1:length(Mdeg_results.list)){
# names(Mdeg_results.list[[i]]) <- c("july", "aug", "sept", "oct")
# }
# 
# 
# 
# #collate results
# #this collapses the 2nd order elements down to the 1st order element
# 
# #make a list to store things
# Mdeg_results.summ <- list()
# 
# #loop across all sites and collapse the 2e dfs into 1 df per 1e level
# for(i in 1:length(Mdeg_results.list)){
#   #for all 12 sites
#   summary <- do.call("rbind", Mdeg_results.list[[i]])
#   Mdeg_results.summ[[i]] <- summary
# }
# 
# #remove the row names of the dfs (its the month, but we already have that info)
# for(i in 1:length(Mdeg_results.summ)){
#   row.names(Mdeg_results.summ[[i]]) <- NULL
# }
# 
# #name the 1st order elements
# names(Mdeg_results.summ) <- names(Mdeg_results.list)
# 
# ## make a freiggein huge df
# Mdeg_results <- do.call(rbind.data.frame, Mdeg_results.summ)
# 
# saveRDS(Mdeg_results, here("Mdeg_breeder_results_05.12.23.RDS"))
# 
# saveRDS(Mdeg_results.list, here("Mdeg_breeder_results.list_05.12.23.RDS"))
# 
# saveRDS(Mdeg_results.summ, here("Mdeg_breeder_results.summ_05.12.23.RDS"))
# 
# saveRDS(Mdeg.list, here("Mdeg_breeder.list_05.12.23.RDS"))

##########################################  END  ########################################################
#########################################################################################################











#########################################################################################################
######### How does previous month's [FEMALE DEGREE] affect current infection status? ##########

# ## last month's FEMALE degree, this month's PUUV status
# ## grouped by year and treatment (so replicate sites are combined)
# netmets_puuv_female <- netmets_puuv %>% unite(yr_trt, year, trt) %>%
#   ungroup() %>% group_by(yr_trt) %>%
#   drop_na(prev_F.deg)
# 
# 
# #https://stackoverflow.com/questions/46407320/data-frame-to-nested-list
# Fdeg.list <- lapply(split(netmets_puuv_female, netmets_puuv_female$yr_trt, drop = TRUE),
#                     function(x) split(x, x[["month"]], drop = TRUE))
# 
# 
# ###### NOW I HAVE A nested list of all the sites in both years, separated to 2e order list elements by month
# ## I need to loop through all of those, calculate a u-statistic and then permutation and get a pvalue
# ## ideally 1,000 to 10,000 permutations WITHOUT it taking a millions years
# ## I have permutation z-score network things that might have useful code
# 
# Fdeg_results.list <- list()
# 
# for(i in 1:length(Fdeg.list)){
#   
#   #for each year-trt combo (there are 8)
#   print(names(Fdeg.list[i]))
#   
#   yr_trt <- list()
#   
#   for(j in 1:length(Fdeg.list[[i]])){
#     
#     #for each month (of puuv status) within a year_trt combo (should be 4)
#     print(names(Fdeg.list[[i]][j]))
#     
#     
#     ### Mann-Whitney U test - done manually to calculate little u (Croft et al. 2008 pg110)
#     
#     # #to run it in R, use this (but you only get one of the U values)
#     # wilcox.test(wt.deg ~ puuv_ifa, data=data)
#     
#     data <- Fdeg.list[[i]][[j]]
#     
#     #To calculate the U statistic, start by sorting the data
#     data.sort <- data[order(data$prev_F.deg),]
#     
#     #Then assign the ranks and split into two groups
#     data.sort$prev_deg.ranked <- rank(data.sort$prev_F.deg, ties.method = 'average')
#     pos <- filter(data.sort, puuv_ifa == "1")
#     neg <- filter(data.sort, puuv_ifa == "0")
#     
#     #With the ranks calculated, the sums can be found and the values of U1 & U2 computed
#     R1 <- sum(pos$prev_deg.ranked)
#     # R2 <- sum(neg$deg.ranked)
#     
#     n1 <- length(pos$prev_deg.ranked)
#     n2 <- length(neg$prev_deg.ranked)
#     
#     U1 <- n1 * n2 + n1 * (n1 + 1) / 2 - R1
#     # U2 <- n1 * n2 + n1 * (n2 + 1) / 2 - R2
#     
#     u_obs <- U1 / (n1*n2)
#     
#     
#     
#     # montecarlo - https://www.countbayesie.com/blog/2015/3/3/6-amazing-trick-with-monte-carlo-simulations
#     
#     mc.list <- list()
#     
#     runs = 9999 #remember, your obs data counts a permutation (runs+1 = 1000)
#     
#     for(k in 1:runs){
#       
#       shuffled.data <- transform(data, puuv_ifa = sample(puuv_ifa, replace=FALSE) )
#       
#       #### SOURCE: https://rpubs.com/aaronsc32/measure-cabbages-mann-whitney  #######
#       
#       #To calculate the U statistic, start by sorting the data
#       shuff.sort <- shuffled.data[order(shuffled.data$prev_F.deg),]
#       
#       #Then assign the ranks and split into two groups
#       shuff.sort$prev_deg.ranked <- rank(shuff.sort$prev_F.deg, ties.method = 'average')
#       pos <- filter(shuff.sort, puuv_ifa == "1")
#       neg <- filter(shuff.sort, puuv_ifa == "0")
#       
#       #With the ranks calculated, the sums can be found and the values of U1 & U2 computed
#       R1 <- sum(pos$prev_deg.ranked)
#       # R2 <- sum(neg$deg.ranked)
#       
#       n1 <- length(pos$prev_deg.ranked)
#       n2 <- length(neg$prev_deg.ranked)
#       
#       U1 <- n1 * n2 + n1 * (n1 + 1) / 2 - R1
#       # U2 <- n1 * n2 + n1 * (n2 + 1) / 2 - R2
#       
#       u_mc <- U1 / (n1*n2)
#       
#       #### end #####
#       
#       mc.list[[k]] <- u_mc
#       
#     }
#     
#     mc.vec <- unlist(mc.list) #take a list to a vector with unlist()
#     
#     mc.p.value <- (sum(mc.vec > u_obs) + 1) / (runs + 1)
#     mc.p.value.less <- (sum(mc.vec < u_obs) + 1) / (runs + 1)
#     
#     
#     #dataframe to hold results per month
#     month_ID <- as.vector(names(Fdeg.list[[i]][j]))
#     yr_trt[[j]] <- data.frame(month_ID)
#     
#     yr_trt[[j]]$u <- u_obs
#     yr_trt[[j]]$p.value <- mc.p.value
#     yr_trt[[j]]$p.val.less <- mc.p.value.less
#     
#   }
#   
#   Fdeg_results.list[[i]] <- yr_trt
#   
# }
# 
# 
# #name the 8 1st order elements of Fdeg_results.list
# names(Fdeg_results.list) <- names(Fdeg.list)
# 
# #rename the sublist items (months) for each site
# #accounting for the fact that most sites have 5 months of data but some have 4
# for(i in 1:length(Fdeg_results.list)){
#   names(Fdeg_results.list[[i]]) <- c("july", "aug", "sept", "oct")
# }
# 
# 
# 
# #collate results
# #this collapses the 2nd order elements down to the 1st order element
# 
# #make a list to store things
# Fdeg_results.summ <- list()
# 
# #loop across all sites and collapse the 2e dfs into 1 df per 1e level
# for(i in 1:length(Fdeg_results.list)){
#   #for all 12 sites
#   summary <- do.call("rbind", Fdeg_results.list[[i]])
#   Fdeg_results.summ[[i]] <- summary
# }
# 
# #remove the row names of the dfs (its the month, but we already have that info)
# for(i in 1:length(Fdeg_results.summ)){
#   row.names(Fdeg_results.summ[[i]]) <- NULL
# }
# 
# #name the 1st order elements
# names(Fdeg_results.summ) <- names(Fdeg_results.list)
# 
# ## make a freiggein huge df
# Fdeg_results <- do.call(rbind.data.frame, Fdeg_results.summ)
# 
# saveRDS(Fdeg_results, here("Fdeg_breeder_results_05.12.23.RDS"))
# 
# saveRDS(Fdeg_results.list, here("Fdeg_breeder_results.list_05.12.23.RDS"))
# 
# saveRDS(Fdeg_results.summ, here("Fdeg_breeder_results.summ_05.12.23.RDS"))
# 
# saveRDS(Fdeg.list, here("Fdeg_breeder.list_05.12.23.RDS"))

##########################################  END  ########################################################
#########################################################################################################







###################### FROM MSI 5.11.23 ########################

## results from MSI 5.11.23 ##

wtdeg <- readRDS(here("wtdeg_breeder_results_05.12.23.rds")) %>%
  rownames_to_column("yr_trt")
write.csv(wtdeg, here("wtdeg_breeder_results.csv"))

lagdeg <- readRDS(here("lagdeg_breeder_results_05.12.23.rds")) %>%
  rownames_to_column("yr_trt")
write.csv(lagdeg, here("lagdeg_breeder_results.csv"))

Mdeg <- readRDS(here("Mdeg_breeder_results_05.12.23.rds")) %>%
  rownames_to_column("yr_trt")
write.csv(Mdeg, here("Mdeg_breeder_results.csv"))

Fdeg <- readRDS(here("Fdeg_breeder_results_05.12.23.rds")) %>%
  rownames_to_column("yr_trt")
write.csv(Fdeg, here("Fdeg_breeder_results.csv"))

#################################################################
















########################## old things below ################################







##### compare network metric and puuv status

################## THIS CODE HAS ALREADY BEEN ADDED ABOVE ##################################
# #create a df with lagged degree (degree from previous month influences current PUUV status)
# nm_puuv_lag <- net_mets_puuv %>% group_by(tag) %>%
#   arrange(year, occasion, .by_group = TRUE) %>%
#   mutate(n_ifa = length(tag)) %>%
#   filter(n_ifa > 1) %>% #remove animals with only one IFA 
#   select(!c(n_ifa, wt.deg, net.centralization, net.clust, mod)) %>%
#   mutate(lag_deg = lag(deg, n=1)) %>%
#   mutate(prev_curr_puuv = paste( lag(puuv_ifa), puuv_ifa, sep="-")) %>%
#   drop_na(lag_deg) %>% #drop rows without previous network position
#   filter(prev_curr_puuv == "0-0" | prev_curr_puuv == "0-1") %>% #only want neg-neg or neg-pos
#   mutate(serovert = as.factor(case_when(prev_curr_puuv == "0-0" ~ 0,
#                               prev_curr_puuv == "0-1" ~ 1)))
############################################################################################

# #create year-month-site column for stats
# net_mets_puuv <- net_mets_puuv %>% unite(y_m_site, year, month, site, sep="_", remove=FALSE)
# 
# 
# net_mets_puuv_list <- split(nm_puuv_lag, nm_puuv_lag$y_m_site)
# 
# outlist <- list()
# 
# for(i in 1:length(net_mets_puuv_list)){
#   
#   print(i)
#   
#   outlist[[i]] <- data.frame(y_m_site=character(length=1),
#                              pval=numeric(length=1))
#   
#   outlist[[i]]$y_m_site <- paste(names(net_mets_puuv_list)[[i]])
#   
#   #in case all animals were pos or neg for puuv, can't run wilcox_test
#   if (length(unique(net_mets_puuv_list[[i]]$puuv_ifa))==1) {value <- NA}
#   
#   #if there are no edges in the network, can't run wilcox_test
#   if ( sum(net_mets_puuv_list[[i]]$lag_deg)==0 ) {value <- NA}
#   
#   #if there are animals that are pos and neg for puuv AND the network has edges, do the thing
#   if (length(unique(net_mets_puuv_list[[i]]$puuv_ifa))==2 & sum(net_mets_puuv_list[[i]]$lag_deg)!=0 ) {
#     #run a ?permutation? wilcox test 
#     out <- wilcox_test(lag_deg ~ puuv_ifa, data = net_mets_puuv_list[[i]], distribution = "approximate", conf.int = TRUE)
#     #pull the pvalue
#     value <- pvalue(out)[1]
#   }
#   
#   #write the output to the df
#   outlist[[i]]$pval <- value
#   
# }
# 
# 
# 
# adataframe <- do.call(rbind.data.frame, outlist)
# 
# whocares <- adataframe %>% filter(pval<0.05)
# #there does not appear to be a difference in the degree of infected vs uninfected animals (IN A SINGLE MONTH)
























############# feb 17 2023 ################
## let's just deal with one site in one month  and then we'll go from there

# library(coin)
#this is a permutation-based wilcox test (?)
out <- wilcox_test(deg ~ puuv_ifa, data = talo_sept22, 
            distribution = "approximate", conf.int = TRUE)


value <- pvalue(out)[1]


#take net_mets_puuv down to one entry per month
sitelevel_net_mets_puuv <- net_mets_puuv %>% group_by(year, site, month) %>%
  slice(1) %>%
  dplyr::select(!c(tag, samp_id, puuv_ifa, sex, deg, wt.deg))

sitelevel_net_mets_puuv <- sitelevel_net_mets_puuv %>% left_join(puuv_prev, by=c("year", "site", "month")) %>%
  mutate(prev = as.numeric(prev),
         mod = as.numeric(mod),
         net.centralization = as.numeric(net.centralization),
         net.clust = as.numeric(net.clust))

cor.test(sitelevel_net_mets_puuv$mod, sitelevel_net_mets_puuv$prev, method="spearman", exact=FALSE)

cor.test(sitelevel_net_mets_puuv$net.centralization, sitelevel_net_mets_puuv$prev, method="spearman", exact=FALSE)

cor.test(sitelevel_net_mets_puuv$net.clust, sitelevel_net_mets_puuv$prev, method="spearman", exact=FALSE)








################################################################################################################
#########################################   doin stuff    ######################################################
################################################################################################################
##### lil check to see the voles that seroconvert, and then find the ones that convert from 1 to 0... ##########
################################################################################################################

# net_mets_puuv <- net_mets_puuv %>% group_by(tag) %>% arrange(year, month, .by_group = TRUE)
# net_mets_puuv$puuv_ifa <- as.character(net_mets_puuv$puuv_ifa)
# check <- net_mets_puuv %>% group_by(tag) %>%
#   summarise(status = unique(puuv_ifa)) 
# 
# #save the tag IDs of animals that seroconvert
# duplist <- as.vector(check$tag[duplicated(check$tag)])
# length(duplist) #103 animals seroconvert (2021 and 2022)
# 
# check2 <- check %>% filter(tag %in% duplist) #just the animals that are 0 and 1

# check_longitude <- net_mets_puuv %>% group_by(tag) %>% arrange(year, month, .by_group = TRUE) %>%
#   dplyr::select(year, month, tag, puuv_ifa) %>%
#   summarise(status_time = toString(puuv_ifa)) %>%
#   filter(str_detect(status_time, "0,\\s1|1,\\s0")) #filter for all animals (n=103) that seroconvert (regardless of direction)

# #voles that seroconvert PUUV+ to PUUV-
# puuv_pos_neg <- net_mets_puuv %>% group_by(tag) %>% arrange(year, month, .by_group = TRUE) %>%
#   dplyr::select(year, month, tag, puuv_ifa) %>%
#   summarise(status_time = toString(puuv_ifa)) %>%
#   filter(str_detect(status_time, "1,\\s0")) #filter for animals that go from pos to neg
# 
# problemchildren <- puuv_pos_neg$tag
# write.csv(puuv_pos_neg, here("puuv_pos_neg.csv"))
# 
#grab just the sample ID and notes
# hanta_notes <- read.csv(here("puuv_ifa.csv")) %>%
#   select(id, notes) %>% 
#   rename(samp_id = id)

# #voles that seroconvert PUUV- to PUUV+
# puuv_neg_pos <- net_mets_puuv %>% group_by(tag) %>% arrange(year, month, .by_group = TRUE) %>%
#   dplyr::select(year, month, tag, puuv_ifa) %>%
#   summarise(status_time = toString(puuv_ifa)) %>%
#   filter(str_detect(status_time, "0,\\s1"))

# puuvcheckdata <- hantadata %>% left_join(tag_id, by= c("id" = "samp_id")) %>% filter(tag %in% problemchildren) %>%
#   group_by(tag) %>% arrange(id, .by_group = TRUE) %>% relocate(tag, .after=id)
# 
# write.csv(puuvcheckdata, here("IFAdata_posneg.csv"))

# fulltrap21 <- read_rds(here("fulltrap21_02.16.23.rds")) %>%
#   mutate(year = 2021)
# 
# fulltrap22 <- fulltrap22 %>% select(!new) %>%
#   mutate(year = 2022)
# 
# FULLfulltrap <- rbind(fulltrap21, fulltrap22) %>%
#   select(c(year, site, occ.sess, samp_id, tag, sex, firstcap, mass))
#          
# check <- FULLfulltrap %>% group_by(tag) %>% 
#   summarise(n = n_distinct(year)) %>%
#   filter(n > 1) 
# bothyrs <- check$tag
# overwinter <- FULLfulltrap %>% filter(tag %in% bothyrs) %>% arrange(tag, year, occ.sess)
# write.csv(overwinter, here("overwinters.csv"))

# fullproblemchildren <- inner_join(tag_id, hantadata, by= "samp_id") %>% 
#   inner_join(hanta_notes, by="samp_id") %>%
#   left_join(FULLfulltrap, by=c("samp_id", "tag")) %>%
#   filter(tag %in% problemchildren) %>%
#   select(year, site, occ.sess, tag, samp_id, puuv_ifa, mass, firstcap, sex, notes) %>%
#   arrange(tag, occ.sess)
# 
# write.csv(fullproblemchildren, here("confirm_IFA.csv"))


# #overwintered 2021-22 animals
# OW <- net_mets_puuv %>% group_by(tag) %>% arrange(tag, year, month) %>%
#   mutate(year=as.factor(year),
#          ows = n_distinct(year)) %>% 
#   filter(ows == 2)

################################################### END doin' stuff ####################################################



hantadata <- read.csv(here("puuv_ifa.csv")) %>%
  clean_names %>%
  #populate column of FINAL puuv status
  mutate(FINAL_puuv = ifelse(is.na(puuv_confirm), as.character(puuv_initial), as.character(puuv_confirm))) %>%
  mutate(samp_id = as.numeric(id),
         date_initial = as_date(date_initial, format= "%m/%d/%Y"),
         puuv_initial = as.factor(puuv_initial),
         date_confirm = as_date(date_confirm, format= "%m/%d/%Y"),
         puuv_confirm = as.factor(puuv_confirm),
         FINAL_puuv = as.factor(FINAL_puuv)) %>%
  drop_na(FINAL_puuv) %>%
  select(!id) %>%
  left_join(tag_id, by="samp_id") %>%
  drop_na(tag) %>%
  filter(tag %in% problemchildren) %>%
  arrange(tag, samp_id)







#plot for shiggles
net_mets_puuv22 %>% 
  drop_na(puuv_ifa) %>%
  ggplot(aes(x=puuv_ifa, y=wt.deg, fill=sex)) + geom_boxplot() + facet_grid(trt ~ month)



