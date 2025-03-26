#29 October 2024
#the goal of this code is to use the centroids and the radii based on 50%, 95% HR extents...
  #to create networks of percent HR overlap between voles - essentially replacing the networks based on...
  #overlap of HR distributions with something more concrete

#14 Nov 2024
#this code now calculates % overlap (directed) between all voles, each month in all three years
  #output: "pct_overlap_list2#" which is a nested list of adj matrices per month/site
  #output: "edges_summary2#" a long-ass df (edge list) of pct overlap between all voles in a year

#26 March 2025
#this code now creates HRs of variable size:
  #if vole was caught at least twice in a season, their HR_rad is the mean distance between trapped locs and centroid
  #if vole was caught 1x or only in 1 trap, their HR_rad is based on estimated kernel, scaled smaller for the more times
    #they were caught in the same trap

####### all this code was run 3.26.25 so pct_overlap_list21, 22, 23 and edges_summary21, 22, 23 are up to date

# load packages
library(here) #v 1.0.0
library(tidyverse) #v 2.0.0
library(sf)

#clear environment
rm(list = ls())


####----------- LOAD DATA ----------------- (copied from 06_spaceuse_circles.R)

#a and b parameters of negative sigmoidal curve to estimate space use
#params files generated in file: "02_construct_spatial_overlap_networks.R"
params21 <- readRDS(here("spaceuse_parameters21.rds")) %>% mutate(year=2021) 
# %>% select(!c(aSE, bSE)) #remove SE if it was calculated in 02-1_functions_construct_overlap_networks
params22 <- readRDS(here("spaceuse_parameters22.rds")) %>% mutate(year=2022) 
# %>% select(!c(aSE, bSE)) #remove SE if it was calculated in 02-1_functions_construct_overlap_networks
params23 <- readRDS(here("spaceuse_parameters23.rds")) %>% mutate(year=2023) 

#MONTHLY centroid locations for each vole
#centroids files generated in file: "02_construct_spatial_overlap_networks.R"
centroids21 <- readRDS(here("monthly_centroids21.rds")) %>% rename(tag = Tag_ID)
centroids22 <- readRDS(here("monthly_centroids22.rds")) %>% rename(tag = Tag_ID)
centroids23 <- readRDS(here("monthly_centroids23.rds")) %>% rename(tag = Tag_ID)

#load fulltrap data (trim to necessary trapping metadata for each capture)
ft21 <- readRDS(here("fulltrap21_11.11.24.rds"))
trapmetadata21 <- ft21 %>% group_by(season, tag) %>%
  mutate(caps_per_season = length(tag), 
         traps_per_season = n_distinct(trap)) %>% ungroup() %>% #calculate caps / traps per season
  dplyr::select(c(year, season, trt, site, month, tag, sex, season_breeder, caps_per_season, traps_per_season)) %>%
  group_by(tag, month) %>% slice(1) #keep one row per vole/month

ft22 <- readRDS(here("fulltrap22_11.11.24.rds"))
trapmetadata22 <- ft22 %>% group_by(season, tag) %>%
  mutate(caps_per_season = length(tag), 
         traps_per_season = n_distinct(trap)) %>% ungroup() %>% 
  dplyr::select(c(year, season, trt, site, month, tag, sex, season_breeder, caps_per_season, traps_per_season)) %>%
  group_by(tag, month) %>% slice(1) #keep one row per vole/month

ft23 <- readRDS(here("fulltrap23_11.11.24.rds"))
trapmetadata23 <- ft23 %>% group_by(season, tag) %>%
  mutate(caps_per_season = length(tag), 
         traps_per_season = n_distinct(trap)) %>% ungroup() %>% 
  dplyr::select(c(year, season, trt, site, month, tag, sex, season_breeder, caps_per_season, traps_per_season)) %>%
  group_by(tag, month) %>% slice(1) #keep one row per vole/month


####--------------------------------------

### SUMMARIZE mean space use kernel size (50% core and 95% peripheral) by functional group, trt, year ####
### using the a and b params for each fxnl group and calculate radius of space use  ####

#mean space use kernel data - included both core 50% and periphery 95% HR
spaceusekernels <- rbind(params21, params22, params23) %>%
  # mutate(rad_50_trap = (log((1/0.5)-1) + a) / (-b)) %>% #calculate the radius ("trap units") where probability of detection is 50%
  # mutate(rad_50_m = rad_50_trap*10) %>% #multiply by 10 since a,b parameters are calculated in 'trap units' and traps are 10m apart
  mutate(rad_95_trap = (log((1/0.05)-1) + a) / (-b)) %>% #calculate the radius ("trap units") where probability of detection is 95%
  mutate(rad_95_m = rad_95_trap*10) %>% #multiply by 10 since a,b parameters are calculated in 'trap units' and traps are 10m apart
  separate_wider_delim(stsb, delim="_", names=c("season", "foodtrt", "helmtrt", "sex", "season_breeder")) %>%
  unite(trt, foodtrt, helmtrt) %>%
  dplyr::select(-c(a,b,rad_95_m)) %>%
  rename(HR_rad = rad_95_trap) %>% #using the trap units version of 95% peripheral HR
  mutate(season = factor(season, levels=c("summer", "fall"))) %>%
  mutate(trt = factor(trt, levels=c("unfed_control", "unfed_deworm",
                                    "fed_control", "fed_deworm")))

##-----------------------------------------------


################################################################################
## new stuff Nov 21, 2024 - getting some "variance" measure per vole / season ##
############ using SEASONAL centroids and fulltrap dataset #####################
################################################################################

## the goal of the following code is to more explicitly estimate space use for voles with enough captures
  # in a season to have some idea of what they do. Of course, then we have the issue that summer and fall
  # are different lengths, BUT we're just going to go with it
  # Dave Kennedy (and I? March 25, 2025) decided that all voles with at least 2 captures/season get a HR 
      # that reflects their observed space use
  # voles with 1 cap/season OR captures in only 1 trap will have a HR that is a scaled version of the kernel estimate
  # so that HR size estimate gets smaller the more often you're caught in 1 spot


###----------2021-----------

#calculate weighted average trapped location per SEASON per vole (all voles)
season.centroids21 <- ft21 %>% group_by(site, season, tag) %>%
  mutate(s.x = mean(x),
         s.y = mean(y)) %>% #weighted average centroid per SEASON
  select(site, season, tag, s.x, s.y) %>% 
  slice(1) %>% ungroup() #keep one entry per vole

#join fulltrap with SEASONAL centroid locations per vole
#calculate mean trapped dist from centroid per season (not quite a variance... more like std dev)
#https://stats.stackexchange.com/questions/13272/2d-analog-of-standard-deviation
#std deviation lives in the same units as your data, variance is units^2
season.rad21 <- ft21 %>% 
  group_by(season, tag) %>% mutate(caps_per_season = length(tag), 
                                   traps_per_season = n_distinct(trap)) %>% #calc n traps / caps per season
  #do this for all voles with at least 2 captures and caught in different traps
  filter(caps_per_season>1 & traps_per_season>1) %>% ungroup() %>%
  select(site, season, month, occ.sess, tag, x, y) %>% #trim down fulltrap to just what we need
  left_join(season.centroids21, by=c("site", "season", "tag")) %>%
  # mutate(x = x*10, y = y*10,
  #        s.x = s.x*10, s.y = s.y*10) %>% #convert trap units to meters
  mutate(dist = sqrt((x-s.x)^2+(y-s.y)^2)) %>% #calculate dist bw each trapped location and centroid
  group_by(site, season, tag) %>%
  summarise(HR_rad = mean(dist)) #units of dist is CURRENTLY TRAP UNITS (*10 to get meters)

### NOTE ### VOLES with only 1 capture/season or caught x times but only in 1 trap ARE NOT IN THIS DATASET


###-----------2022------------

#calculate weighted average trapped location per SEASON per vole (all voles)
season.centroids22 <- ft22 %>% group_by(site, season, tag) %>%
  mutate(s.x = mean(x),
         s.y = mean(y)) %>% #weighted average centroid per SEASON
  select(site, season, tag, s.x, s.y) %>%
  slice(1) %>% ungroup() #keep one entry per vole

#join fulltrap with SEASONAL centroid locations per vole
#calculate mean trapped dist from centroid per season (not quite a variance... more like std dev)
#https://stats.stackexchange.com/questions/13272/2d-analog-of-standard-deviation
#std deviation lives in the same units as your data, variance is units^2
season.rad22 <- ft22 %>% 
  group_by(season, tag) %>% mutate(caps_per_season = length(tag), 
                                   traps_per_season = n_distinct(trap)) %>% #calc n traps / caps per season
  #do this for all voles with at least 2 captures and caught in different traps
  filter(caps_per_season>1 & traps_per_season>1) %>% ungroup() %>%
  select(site, season, month, occ.sess, tag, x, y) %>% #trim down fulltrap to just what we need
  left_join(season.centroids22, by=c("site", "season", "tag")) %>%
  # mutate(x = x*10, y = y*10,
  #        s.x = s.x*10, s.y = s.y*10) %>% #convert trap units to meters
  mutate(dist = sqrt((x-s.x)^2+(y-s.y)^2)) %>% #calculate dist bw each trapped location and centroid
  group_by(site, season, tag) %>%
  summarise(HR_rad = mean(dist)) #units of dist is CURRENTLY TRAP UNITS (*10 to get meters)

### NOTE ### VOLES with only 1 capture/season or caught x times but only in 1 trap ARE NOT IN THIS DATASET


###----------2023-------------

#calculate weighted average trapped location per SEASON per vole (all voles)
season.centroids23 <- ft23 %>% group_by(site, season, tag) %>%
  mutate(s.x = mean(x),
         s.y = mean(y)) %>% #weighted average centroid per SEASON
  select(site, season, tag, s.x, s.y) %>%
  slice(1) %>% ungroup() #keep one entry per vole

#join fulltrap with SEASONAL centroid locations per vole
#calculate mean trapped dist from centroid per season (not quite a variance... more like std dev)
#https://stats.stackexchange.com/questions/13272/2d-analog-of-standard-deviation
  #std deviation lives in the same units as your data, variance is units^2
season.rad23 <- ft23 %>% 
  group_by(season, tag) %>% mutate(caps_per_season = length(tag), 
                                   traps_per_season = n_distinct(trap)) %>% #calc n traps / caps per season
  #do this for all voles with at least 2 captures and caught in different traps
  filter(caps_per_season>1 & traps_per_season>1) %>% ungroup() %>%
  select(site, season, month, occ.sess, tag, x, y) %>% #trim down fulltrap to just what we need
  left_join(season.centroids23, by=c("site", "season", "tag")) %>%
  # mutate(x = x*10, y = y*10,
  #        s.x = s.x*10, s.y = s.y*10) %>% #convert trap units to meters
  mutate(dist = sqrt((x-s.x)^2+(y-s.y)^2)) %>% #calculate dist bw each trapped location and centroid
  group_by(site, season, tag) %>%
  summarise(HR_rad = mean(dist)) #units of dist is CURRENTLY TRAP UNITS (*10 to get meters)

### NOTE ### VOLES with only 1 capture/season or caught x times but only in 1 trap ARE NOT IN THIS DATASET



##-----------------------------------------------

#2021 join all the data needed for HR circles
homerange21 <- trapmetadata21 %>% #has traps/caps per season
  left_join(season.rad21, by=c("season", "site", "tag")) %>% #add HR_rad for known seasonal HRs (caught at least 2x, in different traps)
  #patch in HR_rad from kernel estimates for voles with only 1 cap / caught in only 1 trap
  rows_patch(spaceusekernels, by=c("year", "season", "trt", "sex", "season_breeder"), unmatched = "ignore") %>%
  #for voles with only 1 cap / caught in only 1 trap, adjust HR to be kernel est / sqrt(number of times observed in 1 trap)
  #ie the more we see a vole in only 1 trap, the smaller their HR estimate will be
  mutate(HR_rad = case_when(traps_per_season == 1 ~ HR_rad/sqrt(caps_per_season),
                            traps_per_season != 1 ~ HR_rad)) %>%
  left_join(centroids21, by=c("tag", "month", "site")) %>% #monthly centroid locations (x,y)
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>%
  unite(fxnl_grp, sex, season_breeder, sep="_")

#2022 join all the data needed for HR circles
homerange22 <- trapmetadata22 %>% #has traps/caps per season
  left_join(season.rad22, by=c("season", "site", "tag")) %>% #add HR_rad for known seasonal HRs (caught at least 2x, in different traps)
  #patch in HR_rad from kernel estimates for voles with only 1 cap / caught in only 1 trap
  rows_patch(spaceusekernels, by=c("year", "season", "trt", "sex", "season_breeder"), unmatched = "ignore") %>%
  #for voles with only 1 cap / caught in only 1 trap, adjust HR to be kernel est / sqrt(number of times observed in 1 trap)
  #ie the more we see a vole in only 1 trap, the smaller their HR estimate will be
  mutate(HR_rad = case_when(traps_per_season == 1 ~ HR_rad/sqrt(caps_per_season),
                            traps_per_season != 1 ~ HR_rad)) %>%
  left_join(centroids22, by=c("tag", "month", "site")) %>% #monthly centroid locations (x,y)
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>%
  unite(fxnl_grp, sex, season_breeder, sep="_")

#2023 join all the data needed for HR circles
homerange23 <- trapmetadata23 %>% #has traps/caps per season
  left_join(season.rad23, by=c("season", "site", "tag")) %>% #add HR_rad for known seasonal HRs (caught at least 2x, in different traps)
  #patch in HR_rad from kernel estimates for voles with only 1 cap / caught in only 1 trap
  rows_patch(spaceusekernels, by=c("year", "season", "trt", "sex", "season_breeder"), unmatched = "ignore") %>%
  #for voles with only 1 cap / caught in only 1 trap, adjust HR to be kernel est / sqrt(number of times observed in 1 trap)
  #ie the more we see a vole in only 1 trap, the smaller their HR estimate will be
  mutate(HR_rad = case_when(traps_per_season == 1 ~ HR_rad/sqrt(caps_per_season),
                                   traps_per_season != 1 ~ HR_rad)) %>%
  left_join(centroids23, by=c("tag", "month", "site")) %>% #monthly centroid locations (x,y)
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>%
  unite(fxnl_grp, sex, season_breeder, sep="_")


##### NOTE - HR_rad is in TRAP UNITS ######


###---------------------------------------------


#so what I want to do is separate homerange2X into a nested list by site/month (like I've done before)
#then loop through each site, and loop through each month
#and each month, make it into a sf (maybe the quicker way)
#and then calculate weighted edges and maybe save those as a nested list for easy plotting?


############# 2021 DATA ####################

#a slick little something from stackoverflow to construct a nested list in one go #blessed
site_month_list <- lapply(split(homerange21, homerange21$site, drop = TRUE),
                              function(x) split(x, x[["month"]], drop = TRUE))

#list for results
pct_overlap_list <- list() 

for(i in 1:length(site_month_list)){
  
  #for each site
  print(i)
  
  #create a list to hold results per site
  pct_overlap_list[[i]] <- list()
  
  for(j in 1:length(site_month_list[[i]])){
    
    print(j)
    
    pct_overlap_list[[i]][[j]] <- list()
    
    #pull the df for that site/month
    data <- site_month_list[[i]][[j]]
    
    #make a matrix of just the coordinates
    coords <- matrix(c(data$x,data$y), ncol = 2)
    #convert to sf object
    coords_sf <- st_as_sf(data, coords=c("x","y")) 
    #buffer each point with the corresponding radius
    buff_sf <- st_buffer(coords_sf, dist=data$HR_rad) #coords_sf is built from data, so the order of the points should be identical so each point gets its correct radius
    #add the area
    buff_sf$area <- st_area(buff_sf)
    
    #calculate overlap among all pairs
    weighted_overlap <- st_intersection(buff_sf, buff_sf) %>% 
      dplyr::mutate(area = st_area(.), 
                    pct_overlap = area / area.1 ) %>% # "area" is the area of overlap, "area.1" is the total area of focal vole's HR 
      tibble::as_tibble() %>%
      dplyr::select(focal = tag.1, 
                    neighbor = tag, 
                    pct_overlap, ) #selecting and changing column names at the same time 
    
    ##at this point, tag.1 is FOCAL and tag is neighbor
    ##answering the question: how much of FOCAL'S HR is shared with NEIGHBOR?
      
    #but we've lost voles that have no overlaps, so we have to add those back
    tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded for site/month
    tags_iter <- do.call("rbind", replicate(length(tags$tag), tags, simplify = FALSE)) #all tags replicated for each trap
    tags_rep <- tags_iter %>% arrange(tag) %>% rename(focal = tag)
    tag_by_tag <- cbind(tags_rep, tags_iter) %>% rename(neighbor = tag) #every vole as focal, paired to every vole (incluing self) as neighbor
    
    #join the overlap to this complete list
    weighted_overlap_full <- tag_by_tag %>% left_join(weighted_overlap, by=c("focal", "neighbor")) %>%
      filter(focal != neighbor) %>% #remove self-overlaps
      replace_na(list(pct_overlap=0)) %>% #replace NAs with 0 for overlaps that didn't occur
      pivot_wider(id_cols = neighbor, names_from = focal, values_from = pct_overlap, values_fill=NA) %>% #convert long format edgelist to wide format adj matrix
      column_to_rownames(var="neighbor") %>% #for tibble
      as.matrix()
    
    ##column in FOCAL, row is NEIGHBOR - "How much of FOCAL's HR do they share with NEIGHBOR?"
    
    #https://stackoverflow.com/questions/77115183/reorder-matrix-rows-and-columns-based-on-alphabetical-order-of-colnames-and-rown
    weighted_overlap_full <- weighted_overlap_full[sort(rownames(weighted_overlap_full)), sort(colnames(weighted_overlap_full))]
    
    # weighted_overlap$site <- names(site_month_list[i])
    # weighted_overlap$month <- names(site_month_list[[i]][j])
    
    pct_overlap_list[[i]][[j]] <- weighted_overlap_full
    
  }
  
  
}


#name the 12 1st order elements of overlap_network_list as the sites
names(pct_overlap_list) <- names(site_month_list)

#rename the sublist items (months) for each site
for(i in 1:length(pct_overlap_list)){
  names(pct_overlap_list[[i]]) <- c("june", "july", "aug", "sept", "oct")
}


saveRDS(pct_overlap_list, "pct_overlap_list21.rds")

#load it
pct_overlap_list21 <- readRDS(here("pct_overlap_list21.rds"))



###---------- 2021 edge list of pct overlap -----------------

#turn all the adjacency matrices in pct_overlap_list into edge lists (so we can rbind next)
edgelist_list <- list()

for(i in 1:length(pct_overlap_list21)){
  
  edgelist_list[[i]] <- list()
  
  for(j in 1:length(pct_overlap_list21[[i]])){
    
    adjmat <- as.data.frame(pct_overlap_list21[[i]][[j]])
    adjmat <- rownames_to_column(adjmat, var="neighbor")
    edges <- pivot_longer(adjmat, -neighbor, names_to="focal", values_to="weight")
    edges$month <- names(pct_overlap_list21[[i]][j])
    edges$site <- names(pct_overlap_list21[i])
    
    edgelist_list[[i]][[j]] <- edges
  }
}


#collate MONTHLY edgelist_list results
#make a list to store things
edges_summary <- list()

#loop across all sites and collapse the dfs per month into one df for the site
for(i in 1:length(edgelist_list)){

  #for all 12 sites
  summary <- do.call("rbind", edgelist_list[[i]])
  edges_summary[[i]] <- summary
}

#name the 12 1st order elements as their sites
names(edges_summary) <- names(site_month_list)

## make edges_summary into freiggein huge df
edges_summary21 <- do.call(rbind.data.frame, edges_summary)
edges_summary21$year <- 2021

edges_summary21 <- edges_summary21 %>% drop_na(weight) %>% filter(weight>0) %>%
  select(year, site, month, focal, neighbor, weight)

saveRDS(edges_summary21, "edges_summary21.rds")

####---------- end ---------------------------------------









###################################################################
###################################################################

############ 2022 ##############

#a slick little something from stackoverflow to construct a nested list in one go #blessed
site_month_list <- lapply(split(homerange22, homerange22$site, drop = TRUE),
                          function(x) split(x, x[["month"]], drop = TRUE))

#list for results
pct_overlap_list <- list() 

for(i in 1:length(site_month_list)){
  
  #for each site
  print(i)
  
  #create a list to hold results per site
  pct_overlap_list[[i]] <- list()
  
  for(j in 1:length(site_month_list[[i]])){
    
    print(j)
    
    pct_overlap_list[[i]][[j]] <- list()
    
    #pull the df for that site/month
    data <- site_month_list[[i]][[j]]
    
    #make a matrix of just the coordinates
    coords <- matrix(c(data$x,data$y), ncol = 2)
    #convert to sf object
    coords_sf <- st_as_sf(data, coords=c("x","y")) 
    #buffer each point with the corresponding radius
    buff_sf <- st_buffer(coords_sf, dist=data$HR_rad) #coords_sf is built from data, so the order of the points should be identical so each point gets its correct radius
    #add the area
    buff_sf$area <- st_area(buff_sf)
    
    #calculate overlap among all pairs
    weighted_overlap <- st_intersection(buff_sf, buff_sf) %>% 
      dplyr::mutate(area = st_area(.), 
                    pct_overlap = area / area.1 ) %>% # "area" is the area of overlap, "area.1" is the total area of focal vole's HR 
      tibble::as_tibble() %>%
      dplyr::select(focal = tag.1, 
                    neighbor = tag, 
                    pct_overlap, ) #selecting and changing column names at the same time 
    
    ##at this point, tag.1 is FOCAL and tag is neighbor
    ##answering the question: how much of FOCAL'S HR is shared with NEIGHBOR?
    
    #but we've lost voles that have no overlaps, so we have to add those back
    tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded for site/month
    tags_iter <- do.call("rbind", replicate(length(tags$tag), tags, simplify = FALSE)) #all tags replicated for each trap
    tags_rep <- tags_iter %>% arrange(tag) %>% rename(focal = tag)
    tag_by_tag <- cbind(tags_rep, tags_iter) %>% rename(neighbor = tag) #every vole as focal, paired to every vole (incluing self) as neighbor
    
    #join the overlap to this complete list
    weighted_overlap_full <- tag_by_tag %>% left_join(weighted_overlap, by=c("focal", "neighbor")) %>%
      filter(focal != neighbor) %>% #remove self-overlaps
      replace_na(list(pct_overlap=0)) %>% #replace NAs with 0 for overlaps that didn't occur
      pivot_wider(id_cols = neighbor, names_from = focal, values_from = pct_overlap, values_fill=NA) %>% #convert long format edgelist to wide format adj matrix
      column_to_rownames(var="neighbor") %>% #for tibble
      as.matrix()
    
    ##column in FOCAL, row is NEIGHBOR - "How much of FOCAL's HR do they share with NEIGHBOR?"
    
    #https://stackoverflow.com/questions/77115183/reorder-matrix-rows-and-columns-based-on-alphabetical-order-of-colnames-and-rown
    weighted_overlap_full <- weighted_overlap_full[sort(rownames(weighted_overlap_full)), sort(colnames(weighted_overlap_full))]
    
    # weighted_overlap$site <- names(site_month_list[i])
    # weighted_overlap$month <- names(site_month_list[[i]][j])
    
    pct_overlap_list[[i]][[j]] <- weighted_overlap_full
    
  }
  
  
}


#name the 12 1st order elements of overlap_network_list as the sites
names(pct_overlap_list) <- names(site_month_list)

#rename the sublist items (months) for each site
for(i in 1:length(pct_overlap_list)){
  ifelse( length(pct_overlap_list[[i]]) == 5,
          names(pct_overlap_list[[i]]) <- c("june", "july", "aug", "sept", "oct"),
          names(pct_overlap_list[[i]]) <- c("july", "aug", "sept", "oct") )
}


saveRDS(pct_overlap_list, "pct_overlap_list22.rds")

#load it
pct_overlap_list22 <- readRDS(here("pct_overlap_list22.rds"))



###---------- 2022 edge list of pct overlap -----------------

#turn all the adjacency matrices in pct_overlap_list into edge lists (so we can rbind next)
edgelist_list <- list()

for(i in 1:length(pct_overlap_list22)){
  
  edgelist_list[[i]] <- list()
  
  for(j in 1:length(pct_overlap_list22[[i]])){
    
    adjmat <- as.data.frame(pct_overlap_list22[[i]][[j]])
    adjmat <- rownames_to_column(adjmat, var="neighbor")
    edges <- pivot_longer(adjmat, -neighbor, names_to="focal", values_to="weight")
    edges$month <- names(pct_overlap_list22[[i]][j])
    edges$site <- names(pct_overlap_list22[i])
    
    edgelist_list[[i]][[j]] <- edges
  }
}


#collate MONTHLY edgelist_list results
#make a list to store things
edges_summary <- list()

#loop across all sites and collapse the dfs per month into one df for the site
for(i in 1:length(edgelist_list)){
  
  #for all 12 sites
  summary <- do.call("rbind", edgelist_list[[i]])
  edges_summary[[i]] <- summary
}

#name the 12 1st order elements as their sites
names(edges_summary) <- names(site_month_list)

## make edges_summary into freiggein huge df
edges_summary22 <- do.call(rbind.data.frame, edges_summary)
edges_summary22$year <- 2022

edges_summary22 <- edges_summary22 %>% drop_na(weight) %>% filter(weight>0) %>%
  select(year, site, month, focal, neighbor, weight)

saveRDS(edges_summary22, "edges_summary22.rds")

####---------- end ---------------------------------------





###################################################################
###################################################################

############ 2023 ##############

#a slick little something from stackoverflow to construct a nested list in one go #blessed
site_month_list <- lapply(split(homerange23, homerange23$site, drop = TRUE),
                          function(x) split(x, x[["month"]], drop = TRUE))

#list for results
pct_overlap_list <- list() 

for(i in 1:length(site_month_list)){
  
  #for each site
  print(i)
  
  #create a list to hold results per site
  pct_overlap_list[[i]] <- list()
  
  for(j in 1:length(site_month_list[[i]])){
    
    print(j)
    
    pct_overlap_list[[i]][[j]] <- list()
    
    #pull the df for that site/month
    data <- site_month_list[[i]][[j]]
    
    #make a matrix of just the coordinates
    coords <- matrix(c(data$x,data$y), ncol = 2)
    #convert to sf object
    coords_sf <- st_as_sf(data, coords=c("x","y")) 
    #buffer each point with the corresponding radius
    buff_sf <- st_buffer(coords_sf, dist=data$HR_rad) #coords_sf is built from data, so the order of the points should be identical so each point gets its correct radius
    #add the area
    buff_sf$area <- st_area(buff_sf)
    
    #calculate overlap among all pairs
    weighted_overlap <- st_intersection(buff_sf, buff_sf) %>% 
      dplyr::mutate(area = st_area(.), 
                    pct_overlap = area / area.1 ) %>% # "area" is the area of overlap, "area.1" is the total area of focal vole's HR 
      tibble::as_tibble() %>%
      dplyr::select(focal = tag.1, 
                    neighbor = tag, 
                    pct_overlap, ) #selecting and changing column names at the same time 
    
    ##at this point, tag.1 is FOCAL and tag is neighbor
    ##answering the question: how much of FOCAL'S HR is shared with NEIGHBOR?
    
    #but we've lost voles that have no overlaps, so we have to add those back
    tags <- data %>% group_by(tag) %>% slice(1) %>% select(tag) #df of all tags recorded for site/month
    tags_iter <- do.call("rbind", replicate(length(tags$tag), tags, simplify = FALSE)) #all tags replicated for each trap
    tags_rep <- tags_iter %>% arrange(tag) %>% rename(focal = tag)
    tag_by_tag <- cbind(tags_rep, tags_iter) %>% rename(neighbor = tag) #every vole as focal, paired to every vole (incluing self) as neighbor
    
    #join the overlap to this complete list
    weighted_overlap_full <- tag_by_tag %>% left_join(weighted_overlap, by=c("focal", "neighbor")) %>%
      filter(focal != neighbor) %>% #remove self-overlaps
      replace_na(list(pct_overlap=0)) %>% #replace NAs with 0 for overlaps that didn't occur
      pivot_wider(id_cols = neighbor, names_from = focal, values_from = pct_overlap, values_fill=NA) %>% #convert long format edgelist to wide format adj matrix
      column_to_rownames(var="neighbor") %>% #for tibble
      as.matrix()
    
    ##column in FOCAL, row is NEIGHBOR - "How much of FOCAL's HR do they share with NEIGHBOR?"
    
    #https://stackoverflow.com/questions/77115183/reorder-matrix-rows-and-columns-based-on-alphabetical-order-of-colnames-and-rown
    weighted_overlap_full <- weighted_overlap_full[sort(rownames(weighted_overlap_full)), sort(colnames(weighted_overlap_full))]
    
    # weighted_overlap$site <- names(site_month_list[i])
    # weighted_overlap$month <- names(site_month_list[[i]][j])
    
    pct_overlap_list[[i]][[j]] <- weighted_overlap_full
    
  }
  
  
}


#name the 12 1st order elements of overlap_network_list as the sites
names(pct_overlap_list) <- names(site_month_list)

#rename the sublist items (months) for each site
#accounting for the fact that most sites have 5 months of data but some have 4
for(i in 1:length(pct_overlap_list)){
  ifelse( length(pct_overlap_list[[i]]) == 5,
          names(pct_overlap_list[[i]]) <- c("june", "july", "aug", "sept", "oct"),
          names(pct_overlap_list[[i]]) <- c("july", "aug", "sept", "oct") )
}


saveRDS(pct_overlap_list, "pct_overlap_list23.rds")

#load it
pct_overlap_list23 <- readRDS(here("pct_overlap_list23.rds"))



###---------- 2023 edge list of pct overlap -----------------

#turn all the adjacency matrices in pct_overlap_list into edge lists (so we can rbind next)
edgelist_list <- list()

for(i in 1:length(pct_overlap_list23)){
  
  edgelist_list[[i]] <- list()
  
  for(j in 1:length(pct_overlap_list23[[i]])){
    
    adjmat <- as.data.frame(pct_overlap_list23[[i]][[j]])
    adjmat <- rownames_to_column(adjmat, var="neighbor")
    edges <- pivot_longer(adjmat, -neighbor, names_to="focal", values_to="weight")
    edges$month <- names(pct_overlap_list23[[i]][j])
    edges$site <- names(pct_overlap_list23[i])
    
    edgelist_list[[i]][[j]] <- edges
  }
}


#collate MONTHLY edgelist_list results
#make a list to store things
edges_summary <- list()

#loop across all sites and collapse the dfs per month into one df for the site
for(i in 1:length(edgelist_list)){
  
  #for all 12 sites
  summary <- do.call("rbind", edgelist_list[[i]])
  edges_summary[[i]] <- summary
}

#name the 12 1st order elements as their sites
names(edges_summary) <- names(site_month_list)

## make edges_summary into freiggein huge df
edges_summary23 <- do.call(rbind.data.frame, edges_summary)
edges_summary23$year <- 2023

edges_summary23 <- edges_summary23 %>% drop_na(weight) %>% filter(weight>0) %>%
  select(year, site, month, focal, neighbor, weight)

saveRDS(edges_summary23, "edges_summary23.rds")

####---------- end ---------------------------------------













###-------------- OLD CODE - NOT USING -----------------------------

# #this was useful, but I don't think I want to make each month a multipoint because it "sees" all the voles at a site as a single entity,
#   #not several different polygons (which is what I did in the loop below)
# #buffer multipoint with different distances (https://stackoverflow.com/questions/47057316/st-buffer-multipoint-with-different-distance)
# 
# #### JUST SOME LEARNING:
# #sfg = st_point / st_polygon (just one GEOMETRY - a point, a polygon, a multipolygon)
# #sfc <- is a simple feature list column - it's the COLUMN of multiple sf objects as a list (ie just the locations of the geometries)
#   #CLASS "sfc" "sfc_point"
# #st_sf <- makes the SIMPLE FEATURE object - which has a geometry column AND the attributes of each geometry (it's basically a special df)
#   #CLASS "sf" "data.frame"
# 
# 
# #subset data for just one site/month (to get code running, will eventually need to get this all into ONE GIANT LOOP)
# #just pull the vole ID, radius, and x-y coordinates of monthly centroid for now
# data <- homerange21 %>% filter(month=="july", site=="hevonen") %>% ungroup() %>%
#   select(x, y, rad_95_trap, tag, sex, breeder)
# 
# ### these steps will make an sf object for a site/month and plot it
# 
# #make a matrix of just the coordinates
# coords <- matrix(c(data$x,data$y), ncol = 2)
# #convert to sf object
# coords_sf <- st_as_sf(data, coords=c("x","y"))
#     #st_as_sf() provide the object to be converted to sf (ie the df) - anything extra in the df becomes an attribute
#     #in the case of point data, coords= give the columns containing the coordinates (these points become the "geometry" column in the sf)
#     #if you don't provide a CRS, sf will treat as Euclidean space
# #buffer each point with the corresponding radius
# buff_sf <- st_buffer(coords_sf, dist=data$rad_95_trap) #coords_sf is built from data, so the order of the points should be identical so each point gets its correct radius
#     #I think now the points are essentially polygons
#     #with no CRS, the numeric dist is assumed to have units of the coordinates (which is fine, all "trap units")
# 
# 
# #plot with ggplot (can ploth with base R too - plot() )
# ggplot(buff_sf) + geom_sf(aes(fill=sex, alpha=0.5)) +
#   geom_sf_label(aes(label=tag)) #include label with vole ID
# 
# 
# ####----------------------
# 
# # #the following code just needs the 'data' file which is currently just one site/month
# 
# #basically here I want to make each vole's HR its own polygon, then combine all those polygons together into a single sf
# #NOT using a multipoint or multipolygon because that acts like all the circles are one layer and doesn't "see" overlaps the same way (I think)
# 
# #empty list to hold results
# poly_list <- list()
# 
# #loop over df rows #https://campus.datacamp.com/courses/intermediate-r-for-finance/loops-3?ex=10
# for(i in 1:nrow(data)){
# 
#   #pull each row (vole) as its own df
#   df  <- as.data.frame(data[i,])
#   coords <- matrix(c(df$x,df$y), ncol = 2) #centroid coordinates
#   tt <- st_as_sf(df, coords=c("x","y"))
#   sf <- st_buffer(tt, df$rad_95_trap) #buffer with radius
#   sf$area = st_area(sf) #add area of polygon IN TRAP UNITS
# 
#   poly_list[[i]] <- sf #save result as one item of the list (list will have one item for each vole that site/month)
# 
# }
# 
# #convert list of sf's to one sf
#   #(https://stackoverflow.com/questions/51312935/convert-a-list-of-sf-objects-into-one-sf)
# single_sf <- do.call(rbind, poly_list)
# 
# #and now some slick code someone else wrote to compute the pairwise overlap between all the polygons, and have directional % overlap
# #https://stackoverflow.com/questions/70009412/how-to-compute-all-pairwise-interaction-between-polygons-and-the-the-percentage
# # nOTE 'st_overlaps' will not capture polygons contained within another, for that you want 'st_intersects'
# 
# #regarding the sf::st_intersection bit: "sf::st_intersection() is vectorized. So it will find & return all the intersections of
#     #the first & second argument for you. In this case, the two arguments are the same set of polygons."
# 
# weighted_overlap <- st_intersection(single_sf, single_sf) %>%
#   dplyr::mutate(area = st_area(.),
#                 pct_overlap = area / area.1 ) %>% # "area" is the area of overlap, "area.1" is the total area of focal vole's HR
#   tibble::as_tibble() %>%
#   dplyr::select(neighbor = tag,
#                 focal = tag.1,
#                 pct_overlap, ) %>% #selecting and changing column names at the same time
#   filter(pct_overlap != 0) %>% #remove pairs with no overlap
#   filter(focal != neighbor) %>% #remove self-overlaps
#   select(focal, neighbor, pct_overlap)
# 
# #the output is a tibble with three columns: focal, neighbor, and pct_overlap
# #pct_overlap answers the question: "how much of focal's HR is shared with neighbor?"


##-----------------------------------

  