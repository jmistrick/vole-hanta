#load libraries
library(here)
library(tidyverse)

#load data
fulltrap21 <- readRDS(here("fulltrap21_02.16.23.rds"))
fulltrap22 <- readRDS(here("fulltrap22_02.28.23.rds"))
hanta <- readRDS(here("hantadata_02.28.23.rds"))

traphanta_21 <- left_join(fulltrap21, hanta, by="samp_id") %>%
  select(month, occ.sess, site, tag, samp_id, trap, x, y, puuv_ifa) %>%
  filter(month != "may") %>%
  mutate(x_jitter = jitter(x),
         y_jitter = jitter(y))

traphanta_22 <- left_join(fulltrap22, hanta, by="samp_id") %>%
  select(month, occ.sess, site, tag, samp_id, trap, x, y, puuv_ifa) %>%
  filter(month != "may") %>%
  mutate(x_jitter = jitter(x),
         y_jitter = jitter(y))

#plot stat_density with points
#https://ggplot2.tidyverse.org/reference/geom_density_2d.html
##didn't use, but relevant perhaps:
  #https://stackoverflow.com/questions/24078774/overlay-two-ggplot2-stat-density2d-plots-with-alpha-channels

p <- traphanta_21 %>%
  filter(site=="asema") %>%
  ggplot(aes(x = x_jitter, y = y_jitter)) +
  geom_density2d_filled(alpha=0.5) +
  facet_wrap(~month, ncol=5) 

p + geom_point(aes(x=x_jitter, y=y_jitter, color=puuv_ifa)) +
  scale_color_manual(values= c("0" = "black",
                               "1" = "red"))

ggsave(final, file=paste0("_2021_capture_puuv", ".png"), width = 16, height = 4, units = "in")



traphanta21_list <- split(traphanta_21, traphanta_21$site)
traphanta22_list <- split(traphanta_22, traphanta_22$site)


for(i in names(traphanta21_list)){
  
  print(i)
  
  data <- traphanta21_list[[i]]
  
  p <- data %>%
      ggplot(aes(x = x_jitter, y = y_jitter)) +
      geom_density2d_filled(alpha=0.5) +
      facet_wrap(~month, ncol=5)
  
  final <- p + geom_point(aes(x=x_jitter, y=y_jitter, color=puuv_ifa)) +
      scale_color_manual(values= c("0" = "black",
                                   "1" = "red")) +
      theme(legend.position = "none")
  
  ggsave(final, file=paste0(i, "_2021_capture_puuv", ".png"), width = 16, height = 4, units = "in")
    
}


for(i in names(traphanta22_list)){
  
  print(i)
  
  data <- traphanta22_list[[i]]
  
  p <- data %>%
    ggplot(aes(x = x_jitter, y = y_jitter)) +
    geom_density2d_filled(alpha=0.5) +
    facet_wrap(~month, ncol=5)
  
  final <- p + geom_point(aes(x=x_jitter, y=y_jitter, color=puuv_ifa)) +
    scale_color_manual(values= c("0" = "black",
                                 "1" = "red")) +
    theme(legend.position = "none")
  
  ggsave(final, file=paste0(i, "_2022_capture_puuv", ".png"), width = 16, height = 4, units = "in")
  
}

