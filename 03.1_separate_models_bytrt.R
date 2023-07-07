

# code run using R version 4.1.2 (2021-11-01) -- "Bird Hippie"

#load libraries
library(here)
library(tidyverse)
library(janitor)
library(lubridate)
library(lme4)
library(gtsummary) #pretty regression summary tables

#clear environment
rm(list = ls())

######################### FORMAT DATA ##################################

#load netmets_puuv
netmets_puuv <- readRDS(here("netmets_puuv_06.09.23.rds")) %>%
    rename(Sex = sex,
           Treatment = trt,
           Previous_Month = prev_month,
           Previous_Network_Size = prev_n.node,
           Year = year,
           Previous_M.degree = prev_M.deg,
           Previous_F.degree = prev_F.deg)

#pull the numeric (ie not categorical predictors)
scaled_numeric <- netmets_puuv %>% ungroup() %>% select(explore,
                                                        prev_wt.deg,
                                               Previous_F.degree, Previous_M.degree, 
                                               prev_b.deg, prev_nb.deg,
                                               prev_mb.deg, prev_mnb.deg, prev_fb.deg, prev_fnb.deg, 
                                               Previous_Network_Size)
#scale and center
scaled_numeric <- scale(scaled_numeric, center=TRUE, scale=TRUE)
#pull the nessary categorical predictors
nm_puuv_short <- netmets_puuv %>% select(Year, Treatment, site, Previous_Month, tag, puuv_ifa, Sex, season_breeder) %>%
  mutate(prev_month_num = case_when(Previous_Month == "june" ~ 1,
                                    Previous_Month == "july" ~ 2,
                                    Previous_Month == "aug" ~ 3,
                                    Previous_Month == "sept" ~ 4)) %>%
  mutate(prev_month_num = as.numeric(prev_month_num))
#cbind scaled numeric
nm_puuv_scaled <- cbind(nm_puuv_short, scaled_numeric)  



######################## How does previous degree affect current infection status? ############################


####-------------------------------- UNFED CONTROL SITES, BOTH YEARS ------------------------------------------

nm_UC <- nm_puuv_scaled %>% drop_na(prev_wt.deg) %>%
  filter(Treatment=="unfed_control")

# NOTE: THERE ARE NO NON-BREEDERS that are PUUV+ - have to remove 'season_breeder' to make the models happy


mod_wdegs <- glmer(puuv_ifa ~ prev_wt.deg:Sex + 
                     Sex + explore +
                     Previous_Month + Previous_Network_Size + Year + (1|site),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                   family=binomial, data=nm_UC)

#no other models - there are no infected nonbreeders

########################

mod_sexs <- glmer(puuv_ifa ~ Previous_M.degree:Sex + Previous_F.degree:Sex + 
                    Sex + explore +
                    Previous_Month + Previous_Network_Size + Year + (1|site),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                  family=binomial, data=nm_UC)

#no other models - there are no infected nonbreeders

########################

mod_breeds <- glmer(puuv_ifa ~ prev_b.deg:Sex + prev_nb.deg:Sex + 
                      Sex + explore +
                      Previous_Month + Previous_Network_Size + Year + (1|site),
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                    family=binomial, data=nm_UC)

#no other models - there are no infected nonbreeders

########################

mod_sbs <- glmer(puuv_ifa ~ prev_mb.deg:Sex + prev_mnb.deg:Sex + 
                   prev_fb.deg:Sex + prev_fnb.deg:Sex + 
                   Sex + explore +
                   Previous_Month + Previous_Network_Size + Year + (1|site),
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                 family=binomial, data=nm_UC)

#no other models - there are no infected nonbreeders

########################

AIC(mod_wdegs, 
    mod_sexs,
    mod_breeds, 
    mod_sbs)

#BEST MODEL: mod_breeds
summary(mod_breeds)

########################

#GLMM model diagnostics > https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
library(DHARMa)
#calculate residuals (then run diagnostics on these)
simulationOutput <- simulateResiduals(fittedModel = mod_breeds)
plot(simulationOutput) #qq plot and residual vs fitted

#other specialized goodness-of-fit tests
testDispersion(simulationOutput) #formal test for overdispersion
# testZeroInflation(simulationOutput) #formal test for zero inflation (common type of overdispersion)
testUniformity(simulationOutput)
testOutliers(simulationOutput)
testQuantiles(simulationOutput)
testZeroInflation(simulationOutput)

data <- nm_UC

plotResiduals(simulationOutput, data$Sex)
plotResiduals(simulationOutput, data$season_breeder)
plotResiduals(simulationOutput, data$explore)
# plotResiduals(simulationOutput, data$prev_month_num)
plotResiduals(simulationOutput, data$Previous_Month)
plotResiduals(simulationOutput, data$Previous_Network_Size)
plotResiduals(simulationOutput, data$Year)
# plotResiduals(simulationOutput, data$prev_mb.deg)
# plotResiduals(simulationOutput, data$prev_mnb.deg)
# plotResiduals(simulationOutput, data$prev_fb.deg)
# plotResiduals(simulationOutput, data$prev_fnb.deg)
# plotResiduals(simulationOutput, data$Previous_F.degree) 
# plotResiduals(simulationOutput, data$Previous_M.degree)
plotResiduals(simulationOutput, data$prev_b.deg)
plotResiduals(simulationOutput, data$prev_nb.deg)

### OVERALL: model fit is pretty good
  ### residuals by parameter: only nonbreeder degree is whack


####----------------------------------------------- END ----------------------------------------------------




####-------------------------------- UNFED DEWORM SITES, BOTH YEARS ------------------------------------------


nm_UD <- nm_puuv_scaled %>% drop_na(prev_wt.deg) %>%
  filter(Treatment=="unfed_deworm")


mod_wdegs <- glmer(puuv_ifa ~ prev_wt.deg:Sex + 
                     Sex + season_breeder + explore +
                     Previous_Month + Previous_Network_Size + Year + (1|site),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                   family=binomial, data=nm_UD)

mod_wdegb <- glmer(puuv_ifa ~ prev_wt.deg:season_breeder + 
                     Sex + season_breeder + explore +
                     Previous_Month + Previous_Network_Size + Year + (1|site),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                   family=binomial, data=nm_UD)

mod_wdegsb <- glmer(puuv_ifa ~ prev_wt.deg:Sex:season_breeder + 
                      Sex + season_breeder + explore +
                      Previous_Month + Previous_Network_Size + Year + (1|site),
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                    family=binomial, data=nm_UD)

########################

mod_sexs <- glmer(puuv_ifa ~ Previous_M.degree:Sex + Previous_F.degree:Sex + 
                    Sex + season_breeder + explore +
                    Previous_Month + Previous_Network_Size + Year + (1|site),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                  family=binomial, data=nm_UD)

mod_sexb <- glmer(puuv_ifa ~ Previous_M.degree:season_breeder + Previous_F.degree:season_breeder + 
                    Sex + season_breeder + explore +
                    Previous_Month + Previous_Network_Size + Year + (1|site),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                  family=binomial, data=nm_UD)

mod_sexsb <- glmer(puuv_ifa ~ Previous_M.degree:Sex:season_breeder + Previous_F.degree:Sex:season_breeder + 
                     Sex + season_breeder + explore +
                     Previous_Month + Previous_Network_Size + Year + (1|site),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                   family=binomial, data=nm_UD)

########################

mod_breeds <- glmer(puuv_ifa ~ prev_b.deg:Sex + prev_nb.deg:Sex + 
                      Sex + season_breeder + explore +
                      Previous_Month + Previous_Network_Size + Year + (1|site),
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                    family=binomial, data=nm_UD)

mod_breedb <- glmer(puuv_ifa ~ prev_b.deg:season_breeder + prev_nb.deg:season_breeder + 
                      Sex + season_breeder + explore +
                      Previous_Month + Previous_Network_Size + Year + (1|site),
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                    family=binomial, data=nm_UD)

mod_breedsb <- glmer(puuv_ifa ~ prev_b.deg:Sex:season_breeder + prev_nb.deg:Sex:season_breeder + 
                       Sex + season_breeder + explore +
                       Previous_Month + Previous_Network_Size + Year + (1|site),
                     control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                     family=binomial, data=nm_UD)

########################

mod_sbs <- glmer(puuv_ifa ~ prev_mb.deg:Sex + prev_mnb.deg:Sex + 
                   prev_fb.deg:Sex + prev_fnb.deg:Sex + 
                   Sex + season_breeder + explore +
                   Previous_Month + Previous_Network_Size + Year + (1|site),
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                 family=binomial, data=nm_UD)

mod_sbb <- glmer(puuv_ifa ~ prev_mb.deg:season_breeder + prev_mnb.deg:season_breeder + 
                   prev_fb.deg:season_breeder + prev_fnb.deg:season_breeder + 
                   Sex + season_breeder + explore +
                   Previous_Month + Previous_Network_Size + Year + (1|site),
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                 family=binomial, data=nm_UD)

# NOTE : there is only 1 infected nonbreeder female, can't fit model by functionalgroup:Sex:Breeder - HESSIAN IS SINGULAR

########################

AIC(mod_wdegs, mod_wdegb, mod_wdegsb, 
    mod_sexs, mod_sexsb, mod_sexb,
    mod_breeds, mod_breedb, mod_breedsb,
    mod_sbs, mod_sbb)

#BEST MODEL: mod_breeds
summary(mod_breeds)

########################

#GLMM model diagnostics > https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
library(DHARMa)
#calculate residuals (then run diagnostics on these)
simulationOutput <- simulateResiduals(fittedModel = mod_breeds)
plot(simulationOutput) #qq plot and residual vs fitted

#other specialized goodness-of-fit tests
testDispersion(simulationOutput) #formal test for overdispersion
# testZeroInflation(simulationOutput) #formal test for zero inflation (common type of overdispersion)
testUniformity(simulationOutput)
testOutliers(simulationOutput)
testQuantiles(simulationOutput)
testZeroInflation(simulationOutput)

data <- nm_UD

plotResiduals(simulationOutput, data$Sex)
plotResiduals(simulationOutput, data$season_breeder)
plotResiduals(simulationOutput, data$explore)
# plotResiduals(simulationOutput, data$prev_month_num)
plotResiduals(simulationOutput, data$Previous_Month)
plotResiduals(simulationOutput, data$Previous_Network_Size)
plotResiduals(simulationOutput, data$Year)
# plotResiduals(simulationOutput, data$prev_mb.deg)
# plotResiduals(simulationOutput, data$prev_mnb.deg)
# plotResiduals(simulationOutput, data$prev_fb.deg)
# plotResiduals(simulationOutput, data$prev_fnb.deg)
# plotResiduals(simulationOutput, data$Previous_F.degree) 
# plotResiduals(simulationOutput, data$Previous_M.degree)
plotResiduals(simulationOutput, data$prev_b.deg)
plotResiduals(simulationOutput, data$prev_nb.deg)

### OVERALL: model diagnostics look good, 
    ### residuals by parameter are good except nb degree which is v weird


####----------------------------------------------- END ----------------------------------------------------




####-------------------------------- FED CONTROL SITES, BOTH YEARS ------------------------------------------

nm_FC <- nm_puuv_scaled %>% drop_na(prev_wt.deg) %>%
  filter(Treatment=="fed_control")

### NOTE: there must not be enough infected juvs - breederb and sbb models won't run (Singularity)


mod_wdegs <- glmer(puuv_ifa ~ prev_wt.deg:Sex + 
                     Sex + season_breeder + explore +
                     Previous_Month + Previous_Network_Size + Year + (1|site),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                   family=binomial, data=nm_FC)

mod_wdegb <- glmer(puuv_ifa ~ prev_wt.deg:season_breeder + 
                     Sex + season_breeder + explore +
                     Previous_Month + Previous_Network_Size + Year + (1|site),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                   family=binomial, data=nm_FC)

mod_wdegsb <- glmer(puuv_ifa ~ prev_wt.deg:Sex:season_breeder + 
                      Sex + season_breeder + explore +
                      Previous_Month + Previous_Network_Size + Year + (1|site),
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                    family=binomial, data=nm_FC)

########################

mod_sexs <- glmer(puuv_ifa ~ Previous_M.degree:Sex + Previous_F.degree:Sex + 
                    Sex + season_breeder + explore +
                    Previous_Month + Previous_Network_Size + Year + (1|site),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                  family=binomial, data=nm_FC)

mod_sexb <- glmer(puuv_ifa ~ Previous_M.degree:season_breeder + Previous_F.degree:season_breeder + 
                    Sex + season_breeder + explore +
                    Previous_Month + Previous_Network_Size + Year + (1|site),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                  family=binomial, data=nm_FC)

mod_sexsb <- glmer(puuv_ifa ~ Previous_M.degree:Sex:season_breeder + Previous_F.degree:Sex:season_breeder + 
                     Sex + season_breeder + explore +
                     Previous_Month + Previous_Network_Size + Year + (1|site),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                   family=binomial, data=nm_FC)

########################

mod_breeds <- glmer(puuv_ifa ~ prev_b.deg:Sex + prev_nb.deg:Sex + 
                      Sex + season_breeder + explore +
                      Previous_Month + Previous_Network_Size + Year + (1|site),
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                    family=binomial, data=nm_FC)

# #NO GO - MODEL IS SINGULAR
# mod_breedb <- glmer(puuv_ifa ~ prev_b.deg:season_breeder + prev_nb.deg:season_breeder + 
#                       Sex + season_breeder + explore +
#                       Previous_Month + Previous_Network_Size + Year + (1|site),
#                     control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
#                     family=binomial, data=nm_FC)
# 
# #1. check singularity (works for glmer, NOT glm)
# #from here: https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
# tt <- getME(mod_breedb,"theta")
# ll <- getME(mod_breedb,"lower")
# min(tt[ll==0]) # if this is big, you Gucci (close to 0 is bad)

mod_breedsb <- glmer(puuv_ifa ~ prev_b.deg:Sex:season_breeder + prev_nb.deg:Sex:season_breeder + 
                       Sex + season_breeder + explore +
                       Previous_Month + Previous_Network_Size + Year + (1|site),
                     control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                     family=binomial, data=nm_FC)

########################

mod_sbs <- glmer(puuv_ifa ~ prev_mb.deg:Sex + prev_mnb.deg:Sex + 
                   prev_fb.deg:Sex + prev_fnb.deg:Sex + 
                   Sex + season_breeder + explore +
                   Previous_Month + Previous_Network_Size + Year + (1|site),
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                 family=binomial, data=nm_FC)

# ## NO GOOD - MODEL IS SINGULAR
# mod_sbb <- glmer(puuv_ifa ~ prev_mb.deg:season_breeder + prev_mnb.deg:season_breeder + 
#                    prev_fb.deg:season_breeder + prev_fnb.deg:season_breeder + 
#                    Sex + season_breeder + explore +
#                    Previous_Month + Previous_Network_Size + Year + (1|site),
#                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
#                  family=binomial, data=nm_FC)
# 
# #1. check singularity (works for glmer, NOT glm)
# #from here: https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
# tt <- getME(mod_sbb,"theta")
# ll <- getME(mod_sbb,"lower")
# min(tt[ll==0]) # if this is close to 0, that's a bad thing

mod_sbsb <- glmer(puuv_ifa ~ prev_mb.deg:Sex:season_breeder + prev_mnb.deg:Sex:season_breeder + 
                    prev_fb.deg:Sex:season_breeder + prev_fnb.deg:Sex:season_breeder + 
                    Sex + season_breeder + explore +
                    Previous_Month + Previous_Network_Size + Year + (1|site),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                  family=binomial, data=nm_FC)

########################

AIC(mod_wdegs, mod_wdegb, mod_wdegsb, 
    mod_sexs, mod_sexsb, mod_sexb,
    mod_breeds, mod_breedsb,
    mod_sbs, mod_sbsb)

#BEST MODEL: mod_sbsb
summary(mod_sbsb)

#####################################


#GLMM model diagnostics > https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
library(DHARMa)
#calculate residuals (then run diagnostics on these)
simulationOutput <- simulateResiduals(fittedModel = mod_sbsb)
plot(simulationOutput) #qq plot and residual vs fitted

#other specialized goodness-of-fit tests
testDispersion(simulationOutput) #formal test for overdispersion
# testZeroInflation(simulationOutput) #formal test for zero inflation (common type of overdispersion)
testUniformity(simulationOutput)
testOutliers(simulationOutput)
testQuantiles(simulationOutput)
testZeroInflation(simulationOutput)

data <- nm_FC

plotResiduals(simulationOutput, data$Sex)
plotResiduals(simulationOutput, data$season_breeder)
plotResiduals(simulationOutput, data$explore)
# plotResiduals(simulationOutput, data$prev_month_num)
plotResiduals(simulationOutput, data$Previous_Month)
plotResiduals(simulationOutput, data$Previous_Network_Size)
plotResiduals(simulationOutput, data$Year)
plotResiduals(simulationOutput, data$prev_mb.deg)
plotResiduals(simulationOutput, data$prev_mnb.deg)
plotResiduals(simulationOutput, data$prev_fb.deg)
plotResiduals(simulationOutput, data$prev_fnb.deg)
# plotResiduals(simulationOutput, data$Previous_F.degree) 
# plotResiduals(simulationOutput, data$Previous_M.degree)
# plotResiduals(simulationOutput, data$prev_b.deg)
# plotResiduals(simulationOutput, data$prev_nb.deg)

### OVERALL: everything looks good
  ### all the residuals by parameters look good too

####----------------------------------------------- END ----------------------------------------------------




####-------------------------------- FED DEWORM SITES, BOTH YEARS ------------------------------------------

nm_FD <- nm_puuv_scaled %>% drop_na(prev_wt.deg) %>%
  filter(Treatment=="fed_deworm")


mod_wdegs <- glmer(puuv_ifa ~ prev_wt.deg:Sex + 
                     Sex + season_breeder + explore +
                     Previous_Month + Previous_Network_Size + Year + (1|site),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                   family=binomial, data=nm_FD)

mod_wdegb <- glmer(puuv_ifa ~ prev_wt.deg:season_breeder + 
                     Sex + season_breeder + explore +
                     Previous_Month + Previous_Network_Size + Year + (1|site),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                   family=binomial, data=nm_FD)

mod_wdegsb <- glmer(puuv_ifa ~ prev_wt.deg:Sex:season_breeder + 
                      Sex + season_breeder + explore +
                      Previous_Month + Previous_Network_Size + Year + (1|site),
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                    family=binomial, data=nm_FD)

########################

mod_sexs <- glmer(puuv_ifa ~ Previous_M.degree:Sex + Previous_F.degree:Sex + 
                    Sex + season_breeder + explore +
                    Previous_Month + Previous_Network_Size + Year + (1|site),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                  family=binomial, data=nm_FD)

mod_sexb <- glmer(puuv_ifa ~ Previous_M.degree:season_breeder + Previous_F.degree:season_breeder + 
                    Sex + season_breeder + explore +
                    Previous_Month + Previous_Network_Size + Year + (1|site),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                  family=binomial, data=nm_FD)

mod_sexsb <- glmer(puuv_ifa ~ Previous_M.degree:Sex:season_breeder + Previous_F.degree:Sex:season_breeder + 
                     Sex + season_breeder + explore +
                     Previous_Month + Previous_Network_Size + Year + (1|site),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                   family=binomial, data=nm_FD)

########################

mod_breeds <- glmer(puuv_ifa ~ prev_b.deg:Sex + prev_nb.deg:Sex + 
                      Sex + season_breeder + explore +
                      Previous_Month + Previous_Network_Size + Year + (1|site),
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                    family=binomial, data=nm_FD)

mod_breedb <- glmer(puuv_ifa ~ prev_b.deg:season_breeder + prev_nb.deg:season_breeder + 
                      Sex + season_breeder + explore +
                      Previous_Month + Previous_Network_Size + Year + (1|site),
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                    family=binomial, data=nm_FD)

mod_breedsb <- glmer(puuv_ifa ~ prev_b.deg:Sex:season_breeder + prev_nb.deg:Sex:season_breeder + 
                       Sex + season_breeder + explore +
                       Previous_Month + Previous_Network_Size + Year + (1|site),
                     control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                     family=binomial, data=nm_FD)

########################

mod_sbs <- glmer(puuv_ifa ~ prev_mb.deg:Sex + prev_mnb.deg:Sex + 
                   prev_fb.deg:Sex + prev_fnb.deg:Sex + 
                   Sex + season_breeder + explore +
                   Previous_Month + Previous_Network_Size + Year + (1|site),
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                 family=binomial, data=nm_FD)

mod_sbb <- glmer(puuv_ifa ~ prev_mb.deg:season_breeder + prev_mnb.deg:season_breeder + 
                   prev_fb.deg:season_breeder + prev_fnb.deg:season_breeder + 
                   Sex + season_breeder + explore +
                   Previous_Month + Previous_Network_Size + Year + (1|site),
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                 family=binomial, data=nm_FD)

mod_sbsb <- glmer(puuv_ifa ~ prev_mb.deg:Sex:season_breeder + prev_mnb.deg:Sex:season_breeder + 
                    prev_fb.deg:Sex:season_breeder + prev_fnb.deg:Sex:season_breeder + 
                    Sex + season_breeder + explore +
                    prev_month_num + Previous_Network_Size + Year + (1|site),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                  family=binomial, data=nm_FD)

########################

AIC(mod_wdegs, mod_wdegb, mod_wdegsb, 
    mod_sexs, mod_sexsb, mod_sexb,
    mod_breeds, mod_breedb, mod_breedsb,
    mod_sbs, mod_sbb, mod_sbsb)

## BEST MODEL: mod_sbsb
summary(mod_sbsb)

#####################################


#GLMM model diagnostics > https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
library(DHARMa)
#calculate residuals (then run diagnostics on these)
simulationOutput <- simulateResiduals(fittedModel = mod_sbsb)
plot(simulationOutput) #qq plot and residual vs fitted

#other specialized goodness-of-fit tests
testDispersion(simulationOutput) #formal test for overdispersion
# testZeroInflation(simulationOutput) #formal test for zero inflation (common type of overdispersion)
testUniformity(simulationOutput)
testOutliers(simulationOutput)
testQuantiles(simulationOutput)
testZeroInflation(simulationOutput)

data <- nm_FD

plotResiduals(simulationOutput, data$Sex)
plotResiduals(simulationOutput, data$season_breeder)
plotResiduals(simulationOutput, data$explore)
plotResiduals(simulationOutput, data$prev_month_num)
# plotResiduals(simulationOutput, data$Previous_Month)
plotResiduals(simulationOutput, data$Previous_Network_Size)
plotResiduals(simulationOutput, data$Year)
plotResiduals(simulationOutput, data$prev_mb.deg)
plotResiduals(simulationOutput, data$prev_mnb.deg)
plotResiduals(simulationOutput, data$prev_fb.deg)
plotResiduals(simulationOutput, data$prev_fnb.deg)
# plotResiduals(simulationOutput, data$Previous_F.degree) 
# plotResiduals(simulationOutput, data$Previous_M.degree)
# plotResiduals(simulationOutput, data$prev_b.deg)
# plotResiduals(simulationOutput, data$prev_nb.deg)

### OVERALL: QQplot is good, some quantile deviations
  ## also show up in several predictors - some wibbledy-wobblies and jeremy-beremies
#but this guy says it's fine if your QQplot is okay
#https://stats.stackexchange.com/questions/531601/dharma-quantile-deviations-detected


####--------------------------------------- END ---------------------------------------------


################################################################################################
############################### BEST MODELS for each TREATMENT #################################

#unfed control
nm_UC <- nm_puuv_scaled %>% drop_na(prev_wt.deg) %>%
  filter(Treatment=="unfed_control")

mod_UC <- glmer(puuv_ifa ~ prev_b.deg:Sex + prev_nb.deg:Sex + 
                      Sex + explore +
                      Previous_Month + Previous_Network_Size + Year + (1|site),
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                    family=binomial, data=nm_UC)
## UC RESULTS:
  # Males: previous breeder deg 1.64 increase infection (p=0.071)
  # Males: previous nonbreeder deg 0.05 decrease infection (p=0.071)

plot1 <- nm_UC %>% filter(Sex=="M") %>%
  mutate(puuv_ifa = case_when(puuv_ifa=="0" ~ 0,
                              puuv_ifa=="1" ~ 1)) %>%
  ggplot(aes(x=prev_nb.deg, y=puuv_ifa)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  # labs(title="Males: Previous Non-Reproductive Vole Degree") +
  xlab("Previous Non-Reproductive Vole Degree") +
  ylab("PUUV Infection Status (IFA)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10)))

plot2 <- nm_UC %>% filter(Sex=="M") %>%
  mutate(puuv_ifa = case_when(puuv_ifa=="0" ~ 0,
                              puuv_ifa=="1" ~ 1)) %>%
  ggplot(aes(x=prev_b.deg, y=puuv_ifa)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  # labs(title="Males: Previous Reproductive Vole Degree")+
  xlab("Previous Reproductive Vole Degree") +
  ylab("PUUV Infection Status (IFA)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10)))

library(cowplot)
png(filename = "UC_prev.b-nb.deg_forMales.png",
    width=10, height=5, units="in", res=600)
plot_grid(plot1, plot2, labels="AUTO", label_size = 22)
dev.off()

#unfed deworm
nm_UD <- nm_puuv_scaled %>% drop_na(prev_wt.deg) %>%
  filter(Treatment=="unfed_deworm")

mod_UD <- glmer(puuv_ifa ~ prev_b.deg:Sex + prev_nb.deg:Sex + 
                      Sex + season_breeder + explore +
                      Previous_Month + Previous_Network_Size + Year + (1|site),
                    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                    family=binomial, data=nm_UD)
## RESULTS:
  # (Prev month matters, less likely to be pos in Sept, Oct)
  # Females: previous nonbreeder degree 42.6 (yes 42) increase infection (p=0.021)

## but none of the plots (with just degree) are that compelling
nm_UD %>% filter(Sex=="F") %>%
  mutate(puuv_ifa = case_when(puuv_ifa=="0" ~ 0,
                              puuv_ifa=="1" ~ 1)) %>%
  ggplot(aes(x=prev_nb.deg, y=puuv_ifa)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  labs(title="Infection probability of Females by previous nonbreeder degree")


#fed control
mod_FC <- glmer(puuv_ifa ~ prev_mb.deg:Sex:season_breeder + prev_mnb.deg:Sex:season_breeder + 
                    prev_fb.deg:Sex:season_breeder + prev_fnb.deg:Sex:season_breeder + 
                    Sex + season_breeder + explore +
                    Previous_Month + Previous_Network_Size + Year + (1|site),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                  family=binomial, data=nm_FC)
## RESULTS:
  # Explore: more exploratory animals were less likely 0.58 to be infected (p=0.017)
  # (Year: much less likely to be pos in 2022)
  # (nonbreeder: OR of 0 - very unlikely to be infected?)
  # Female breeder: prev male nb degree 2.85 increase infection (p=0.025)
  # Male nonbreeder: prev male nb degree makes more likely to be infected (with a ridiculously large estimate, marginal pval)
  # ?? M and F nonbreeders with fnb degree ?? Odds Ratio = 0 - does that suggest increasing fnb deg makes you (nearly) definitely PUUV-?

# #EXPLORE: not terribly compelling
# nm_FC %>% 
#   mutate(puuv_ifa = case_when(puuv_ifa=="0" ~ 0,
#                               puuv_ifa=="1" ~ 1)) %>%
#   ggplot(aes(x=explore, y=puuv_ifa)) +
#   geom_point() +
#   geom_smooth(method = "glm", method.args = list(family = "binomial"))

# ## but none of the plots (with just degree) are that compelling
# ## wide 95% CI for Female Breeders and prev male nonbreeder degree
# nm_FC %>% filter(Sex=="F" & season_breeder=="breeder") %>%
#   mutate(puuv_ifa = case_when(puuv_ifa=="0" ~ 0,
#                               puuv_ifa=="1" ~ 1)) %>%
#   ggplot(aes(x=prev_mnb.deg, y=puuv_ifa)) +
#   geom_point() +
#   geom_smooth(method = "glm", method.args = list(family = "binomial"))

# ## ONLY TWO infected male nonbreeders
# nm_FC %>% filter(Sex=="M" & season_breeder=="nonbreeder") %>%
#   mutate(puuv_ifa = case_when(puuv_ifa=="0" ~ 0,
#                               puuv_ifa=="1" ~ 1)) %>%
#   ggplot(aes(x=prev_mnb.deg, y=puuv_ifa)) +
#   geom_point() +
#   geom_smooth(method = "glm", method.args = list(family = "binomial"))
# 
# ## ONLY FOUR infected male nonbreeders
# nm_FC %>% filter(Sex=="F" & season_breeder=="nonbreeder") %>%
#   mutate(puuv_ifa = case_when(puuv_ifa=="0" ~ 0,
#                               puuv_ifa=="1" ~ 1)) %>%
#   ggplot(aes(x=prev_fnb.deg, y=puuv_ifa)) +
#   geom_point() +
#   geom_smooth(method = "glm", method.args = list(family = "binomial"))
# 
# ## ONLY TWO infected male nonbreeders
# nm_FC %>% filter(Sex=="M" & season_breeder=="nonbreeder") %>%
#   mutate(puuv_ifa = case_when(puuv_ifa=="0" ~ 0,
#                               puuv_ifa=="1" ~ 1)) %>%
#   ggplot(aes(x=prev_fnb.deg, y=puuv_ifa)) +
#   geom_point() +
#   geom_smooth(method = "glm", method.args = list(family = "binomial"))


#fed deworm
nm_FD <- nm_puuv_scaled %>% drop_na(prev_wt.deg) %>%
  filter(Treatment=="fed_deworm")

mod_FD <- glmer(puuv_ifa ~ prev_mb.deg:Sex:season_breeder + prev_mnb.deg:Sex:season_breeder + 
                    prev_fb.deg:Sex:season_breeder + prev_fnb.deg:Sex:season_breeder + 
                    Sex + season_breeder + explore +
                    prev_month_num + Previous_Network_Size + Year + (1|site),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                  family=binomial, data=nm_FD)
summary(mod_FD)

## RESULTS:
  # (Male increases)
  # (nonbreeder (strongly) decreases-ie nearly no nonbreeders infected?)
  # (Year: more likely to be positive in 2022)
  # Previous network size had a marginal, negative effect 0.43 on infection (p=0.082)
  # Explore: more exploratory animals were less likely  0.63 to be infected (p=0.030)
  # Female breeder: male breeder deg corr with less likely 0.29 to be infected (p=0.015)
  # Female breeder: female breeder deg corr with 4.26 more likely to be infected (p=0.005)

# #EXPLORE: not terribly compelling
# nm_FD %>% 
#   mutate(puuv_ifa = case_when(puuv_ifa=="0" ~ 0,
#                               puuv_ifa=="1" ~ 1)) %>%
#   ggplot(aes(x=explore, y=puuv_ifa)) +
#   geom_point() +
#   geom_smooth(method = "glm", method.args = list(family = "binomial"))

plot1 <- nm_FD %>% filter(Sex=="F" & season_breeder=="breeder") %>%
  mutate(puuv_ifa = case_when(puuv_ifa=="0" ~ 0,
                              puuv_ifa=="1" ~ 1)) %>%
  ggplot(aes(x=prev_fb.deg, y=puuv_ifa)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  # labs(title="Reproductive Females: Previous Repro. Female Degree") +
  xlab("Previous Reproductive Female Degree") +
  ylab("PUUV Infection Status (IFA)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10)))

# ehh wide CI
plot2 <- nm_FD %>% filter(Sex=="F" & season_breeder=="breeder") %>%
  mutate(puuv_ifa = case_when(puuv_ifa=="0" ~ 0,
                              puuv_ifa=="1" ~ 1)) %>%
  ggplot(aes(x=prev_mb.deg, y=puuv_ifa)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  # labs(title="Reproductive Females: Previous Repro. Male Degree") +
  xlab("Previous Reproductive Male Degree") +
  ylab("PUUV Infection Status (IFA)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10)))

library(cowplot)
png(filename = "FD_prev.fb-mb.deg_forFemaleBreeders.png",
    width=10, height=5, units="in", res=600)
plot_grid(plot2, plot1, labels="AUTO", label_size = 22)
dev.off()


####### pretty OUTPUT MODEL SUMMARY #########

library(gtsummary) #https://www.danieldsjoberg.com/gtsummary/articles/tbl_regression.html

mod_FD %>% tbl_regression(exponentiate = TRUE,
                                         pvalue_fun = ~ style_pvalue(.x, digits = 2),) %>%
  bold_p(t = 0.10) %>%
  bold_labels() %>%
  italicize_levels()

# library(gt)
# # Use function from gt package to save table, after converting to 
# # gt object using as_gt()
# gt::gtsave(as_gt(gtsumm_mod), expand=30, here("regression_model.png"))

####### end output model summary #########
