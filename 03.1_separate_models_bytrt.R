

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


#load netmets_puuv
netmets_puuv <- readRDS(here("netmets_puuv_06.09.23.rds"))


######################## How does previous degree affect current infection status? ############################


# nm_UC <- netmets_puuv %>% drop_na(prev_wt.deg) %>%
#   filter(trt=="unfed_control") %>%
#   rename(Sex = sex,
#          Treatment = trt,
#          Previous_Month = prev_month,
#          Previous_Network_Size = prev_n.node,
#          Year = year,
#          Previous_M.degree = prev_M.deg,
#          Previous_F.degree = prev_F.deg)

# nm_UD <- netmets_puuv %>% drop_na(prev_wt.deg) %>%
#   filter(trt=="unfed_deworm") %>%
#   rename(Sex = sex,
#          Treatment = trt,
#          Previous_Month = prev_month,
#          Previous_Network_Size = prev_n.node,
#          Year = year,
#          Previous_M.degree = prev_M.deg,
#          Previous_F.degree = prev_F.deg)

# nm_FC <- netmets_puuv %>% drop_na(prev_wt.deg) %>%
#   filter(trt=="fed_control") %>%
#   rename(Sex = sex,
#          Treatment = trt,
#          Previous_Month = prev_month,
#          Previous_Network_Size = prev_n.node,
#          Year = year,
#          Previous_M.degree = prev_M.deg,
#          Previous_F.degree = prev_F.deg)

nm_FD <- netmets_puuv %>% drop_na(prev_wt.deg) %>%
  filter(trt=="fed_deworm") %>%
  rename(Sex = sex,
         Treatment = trt,
         Previous_Month = prev_month,
         Previous_Network_Size = prev_n.node,
         Year = year,
         Previous_M.degree = prev_M.deg,
         Previous_F.degree = prev_F.deg)

nm_FD_scaled <- nm_FD %>% ungroup() %>% select(explore, prev_mb.deg, prev_mnb.deg, prev_fb.deg, prev_fnb.deg, Previous_Network_Size)
nm_FD_scaled <- scale(nm_FD_scaled, center=TRUE, scale=TRUE)

nm_FD_small <- nm_FD %>% select(Year, Treatment, site, Previous_Month, tag, puuv_ifa, Sex, season_breeder) %>%
  mutate(prev_month_num = case_when(Previous_Month == "june" ~ 1,
                                    Previous_Month == "july" ~ 2,
                                    Previous_Month == "aug" ~ 3,
                                    Previous_Month == "sept" ~ 4)) %>%
  mutate(prev_month_num = as.numeric(prev_month_num))

nm_FD_new <- cbind(nm_FD_scaled, nm_FD_small)  

#define the data
data <- nm_FD_new

mod_sbsb <- glmer(puuv_ifa ~ prev_mb.deg:Sex:season_breeder + prev_mnb.deg:Sex:season_breeder + 
                    prev_fb.deg:Sex:season_breeder + prev_fnb.deg:Sex:season_breeder + 
                    Sex + season_breeder + explore +
                    prev_month_num + Previous_Network_Size + Year + (1|site),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                  family=binomial, data=data)

############################################################################################################################

#define the data
data <- nm_FD

############################################################################################################################


mod_wdegs <- glmer(puuv_ifa ~ prev_wt.deg:Sex + 
                     Sex + season_breeder + explore +
                     Previous_Month + Previous_Network_Size + Year + (1|site),
                   family=binomial, data=data)

mod_wdegb <- glmer(puuv_ifa ~ prev_wt.deg:season_breeder + 
                     Sex + season_breeder + explore +
                     Previous_Month + Previous_Network_Size + Year + (1|site),
                   family=binomial, data=data)

mod_wdegsb <- glmer(puuv_ifa ~ prev_wt.deg:Sex:season_breeder + 
                      Sex + season_breeder + explore +
                      Previous_Month + Previous_Network_Size + Year + (1|site),
                    family=binomial, data=data)

########################

mod_sexs <- glmer(puuv_ifa ~ Previous_M.degree:Sex + Previous_F.degree:Sex + 
                    Sex + season_breeder + explore +
                    Previous_Month + Previous_Network_Size + Year + (1|site),
                  family=binomial, data=data)

mod_sexb <- glmer(puuv_ifa ~ Previous_M.degree:season_breeder + Previous_F.degree:season_breeder + 
                    Sex + season_breeder + explore +
                    Previous_Month + Previous_Network_Size + Year + (1|site),
                  family=binomial, data=data)

mod_sexsb <- glmer(puuv_ifa ~ Previous_M.degree:Sex:season_breeder + Previous_F.degree:Sex:season_breeder + 
                     Sex + season_breeder + explore +
                     Previous_Month + Previous_Network_Size + Year + (1|site),
                   family=binomial, data=data)

########################

mod_breeds <- glmer(puuv_ifa ~ prev_b.deg:Sex + prev_nb.deg:Sex + 
                      Sex + season_breeder + explore +
                      Previous_Month + Previous_Network_Size + Year + (1|site),
                    family=binomial, data=data)

mod_breedb <- glmer(puuv_ifa ~ prev_b.deg:season_breeder + prev_nb.deg:season_breeder + 
                      Sex + season_breeder + explore +
                      Previous_Month + Previous_Network_Size + Year + (1|site),
                    family=binomial, data=data)

mod_breedsb <- glmer(puuv_ifa ~ prev_b.deg:Sex:season_breeder + prev_nb.deg:Sex:season_breeder + 
                       Sex + season_breeder + explore +
                       Previous_Month + Previous_Network_Size + Year + (1|site),
                     family=binomial, data=data)

########################

mod_sbs <- glmer(puuv_ifa ~ prev_mb.deg:Sex + prev_mnb.deg:Sex + 
                   prev_fb.deg:Sex + prev_fnb.deg:Sex + 
                   Sex + season_breeder + explore +
                   Previous_Month + Previous_Network_Size + Year + (1|site),
                 family=binomial, data=data)

mod_sbb <- glmer(puuv_ifa ~ prev_mb.deg:season_breeder + prev_mnb.deg:season_breeder + 
                   prev_fb.deg:season_breeder + prev_fnb.deg:season_breeder + 
                   Sex + season_breeder + explore +
                   Previous_Month + Previous_Network_Size + Year + (1|site),
                 family=binomial, data=data)

mod_sbsb <- glmer(puuv_ifa ~ prev_mb.deg:Sex:season_breeder + prev_mnb.deg:Sex:season_breeder + 
                    prev_fb.deg:Sex:season_breeder + prev_fnb.deg:Sex:season_breeder + 
                    Sex + season_breeder + explore +
                    Previous_Month + Previous_Network_Size + Year + (1|site),
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
                  family=binomial, data=data)

########################

AIC(mod_wdegs, mod_wdegb, mod_wdegsb, 
    mod_sexs, mod_sexsb, mod_sexb,
    mod_breeds, mod_breedb, mod_breedsb,
    mod_sbs, mod_sbb, mod_sbsb)

#nm_UC - all models are bad, mod_breeds is the least bad
summary(mod_breeds) 
#Male overlap with nonbreeders decreases probability of infection

#nm_UF - some models are bad, some at least didn't throw errors, mod_breeds is the least bad
summary(mod_breeds) 
#Female overlap with nonbreeders increases probability of infection

#nm_FC - some of the breed and sb models had errors
summary(mod_breedsb)
#male breeder overlap with breeders increases prob infection
#male nonbreeder overlap with breeders decreases prob infection
#male nonbreeder overlap with nonbreeders increases prob infection

#nm_FD - some failed to converge, a few errors
summary(mod_sbsb)
#explore has a negative effect on infection
#F breeder overlap malebreeder decreases prob infection
#M breeder overlap malebreeder (marginal) increases prob infection
#F breeder overlap femalebreeder INCREASES prob infection (smallest pvalue)


#####################################


#GLMM model diagnostics > https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
library(DHARMa)
#calculate residuals (then run diagnostics on these)
simulationOutput <- simulateResiduals(fittedModel = mod_sbsb)
plot(simulationOutput) #qq plot and residual vs fitted
testDispersion(simulationOutput) #formal test for overdispersion
# testZeroInflation(simulationOutput) #formal test for zero inflation (common type of overdispersion)

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


