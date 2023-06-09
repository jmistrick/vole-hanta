####### these are the Mann-Whitney u tests doing pairwise comparison on wt.deg, M.deg, F.deg between infected/uninfected animals
### this code used to be in "03_NEW_hanta_cleaning" until 6.9.23 when Janine decided all those pairwise comparisons probs needed to go elsewhere
## this code takes 5ever to run in R - code was run on MSI to get through all the iterations
## code runs off netmets_puuv file created in '03_NEW_hanta_cleaning'


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
