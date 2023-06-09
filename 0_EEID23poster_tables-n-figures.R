#### code from '03_NEW_hanta_cleaning' that was used to generate figures and tables for EEID 2023 poster
## probably would still run off the newer versions of netmets_puuv but I don't actually know
# most of these figures just weren't necessary for the final ms / are too focused on prevalence so I can't publish them in my stuff...


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


############ LOAD NETMETS_PUUV FIRST TO RUN ANY OF THE FOLLOWING CODE #########################


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


############################  MODEL GOES HERE ############################

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


########################## MODEL STUFF GOES HERE ################################

###############################################################################################

