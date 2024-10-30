#29 October 2024
#the goal of this code is to use the centroids and the radii based on 50%, 95% HR extents...
  #to create networks of percent HR overlap between voles - essentially replacing the networks based on...
  #overlap of HR distributions with something more concrete


# load packages
library(here) #v 1.0.0
library(tidyverse) #v 2.0.0
library(sf)



####----------- LOAD DATA ----------------- (copied from 06_spaceuse_circles.R)

#a and b parameters of negative sigmoidal curve to estimate space use
#params files generated in file: "02_construct_spatial_overlap_networks.R"
params21 <- readRDS(here("spaceuse_parameters21.rds")) %>% mutate(year=2021) 
# %>% select(!c(aSE, bSE)) #remove SE if it was calculated in 02-1_functions_construct_overlap_networks
params22 <- readRDS(here("spaceuse_parameters22.rds")) %>% mutate(year=2022) 
# %>% select(!c(aSE, bSE)) #remove SE if it was calculated in 02-1_functions_construct_overlap_networks

#MONTHLY centroid locations for each vole
#centroids files generated in file: "02_construct_spatial_overlap_networks.R"
centroids21 <- readRDS(here("monthly_centroids21.rds")) %>% rename(tag = Tag_ID)
centroids22 <- readRDS(here("monthly_centroids22.rds")) %>% rename(tag = Tag_ID)

#load fulltrap data (trapping metadata for each capture) - pull PIT tag, month, site
ft21 <- readRDS(here("fulltrap21_03.04.24.rds"))
trapdata21 <- ft21 %>% select(c(year, season, trt, site, month, tag, sex, season_breeder)) %>%
  group_by(tag, month) %>% slice(1)

ft22 <- readRDS(here("fulltrap22_03.04.24.rds"))
trapdata22 <- ft22 %>% select(c(year, season, trt, site, month, tag, sex, season_breeder)) %>%
  group_by(tag, month) %>% slice(1)


####--------------------------------------

### SUMMARIZE mean space use kernel size (50% core and 95% peripheral) by functional group, trt, year ####
### using the a and b params for each fxnl group and plotting space use size in area (meters^2) ####

#space use data - included both core 50% and periphery 95% HR
spaceusedata <- rbind(params21, params22) %>%
  mutate(rad_50_trap = (log((1/0.5)-1) + a) / (-b)) %>% #calculate the radius ("trap units") where probability of detection is 50%
  mutate(rad_50_m = rad_50_trap*10) %>% #multiply by 10 since a,b parameters are calculated in 'trap units' and traps are 10m apart
  mutate(rad_95_trap = (log((1/0.05)-1) + a) / (-b)) %>% #calculate the radius ("trap units") where probability of detection is 95%
  mutate(rad_95_m = rad_95_trap*10) %>% #multiply by 10 since a,b parameters are calculated in 'trap units' and traps are 10m apart
  separate_wider_delim(stsb, delim="_", names=c("season", "foodtrt", "helmtrt", "sex", "breeder")) %>%
  unite(trt, foodtrt, helmtrt) %>%
  select(-c(a,b,)) %>%
  mutate(season = factor(season, levels=c("summer", "fall"))) %>%
  mutate(trt = factor(trt, levels=c("unfed_control", "unfed_deworm",
                                    "fed_control", "fed_deworm")))

##-----------------------------------------------

#join all the data needed for HR circles
homerange21 <- trapdata21 %>% rename(breeder = season_breeder) %>% 
  left_join(spaceusedata, by=c("year", "season", "trt", "sex", "breeder")) %>%
  left_join(centroids21, by=c("tag", "month", "site")) %>%
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct")))
#join all the data needed for HR circles
homerange22 <- trapdata22 %>% rename(breeder = season_breeder) %>% 
  left_join(spaceusedata, by=c("year", "season", "trt", "sex", "breeder")) %>%
  left_join(centroids22, by=c("tag", "month", "site")) %>%
  mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct")))


###---------------------------------------------

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
# data <- homerange21 %>% filter(month=="june", site=="ketunpesa") %>% ungroup() %>% 
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


#so what I want to do is separate homerange21 into a nested list by site/month (like I've done before)
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
    buff_sf <- st_buffer(coords_sf, dist=data$rad_95_trap) #coords_sf is built from data, so the order of the points should be identical so each point gets its correct radius
    #add the area
    buff_sf$area <- st_area(buff_sf)
    
    #calculate overlap among all pairs
    weighted_overlap <- st_intersection(buff_sf, buff_sf) %>% 
      dplyr::mutate(area = st_area(.), 
                    pct_overlap = area / area.1 ) %>% # "area" is the area of overlap, "area.1" is the total area of focal vole's HR 
      tibble::as_tibble() %>%
      dplyr::select(neighbor = tag, 
                    focal = tag.1, 
                    pct_overlap, ) %>% #selecting and changing column names at the same time 
      filter(pct_overlap != 0) %>% #remove pairs with no overlap
      filter(focal != neighbor) %>% #remove self-overlaps
      select(focal, neighbor, pct_overlap)
    
    weighted_overlap$site <- names(site_month_list[i])
    weighted_overlap$month <- names(site_month_list[[i]][j])
    
    pct_overlap_list[[i]][[j]] <- weighted_overlap
    
  }
  
  
}


#name the 12 1st order elements of overlap_network_list as the sites
names(pct_overlap_list) <- names(site_month_list)

#rename the sublist items (months) for each site
for(i in 1:length(pct_overlap_list)){
  names(pct_overlap_list[[i]]) <- c("june", "july", "aug", "sept", "oct")
}



#collate MONTHLY pct_overlap results
#make a list to store things
pct_overlap_summary <- list()

#loop across all sites and collapse the dfs per month into one df for the site
for(i in 1:length(pct_overlap_list)){
  
  #for all 12 sites
  summary <- do.call("rbind", pct_overlap_list[[i]])
  pct_overlap_summary[[i]] <- summary
}

#name the 12 1st order elements as their sites
names(pct_overlap_summary) <- names(site_month_list)

## make pct_overlap_summary into freiggein huge df
pct_overlap_summary21 <- do.call(rbind.data.frame, pct_overlap_summary)
pct_overlap_summary21$year <- 2021
row.names(pct_overlap_summary21) <- NULL

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
    buff_sf <- st_buffer(coords_sf, dist=data$rad_95_trap) #coords_sf is built from data, so the order of the points should be identical so each point gets its correct radius
    #add the area
    buff_sf$area <- st_area(buff_sf)
    
    #calculate overlap among all pairs
    weighted_overlap <- st_intersection(buff_sf, buff_sf) %>% 
      dplyr::mutate(area = st_area(.), 
                    pct_overlap = area / area.1 ) %>% # "area" is the area of overlap, "area.1" is the total area of focal vole's HR 
      tibble::as_tibble() %>%
      dplyr::select(neighbor = tag, 
                    focal = tag.1, 
                    pct_overlap, ) %>% #selecting and changing column names at the same time 
      filter(pct_overlap != 0) %>% #remove pairs with no overlap
      filter(focal != neighbor) %>% #remove self-overlaps
      select(focal, neighbor, pct_overlap)
    
    weighted_overlap$site <- names(site_month_list[i])
    weighted_overlap$month <- names(site_month_list[[i]][j])
    
    pct_overlap_list[[i]][[j]] <- weighted_overlap
    
  }
  
  
}


#name the 12 1st order elements of overlap_network_list as the sites
names(pct_overlap_list) <- names(site_month_list)

#rename the sublist items (months) for each site
#accounting for the fact that most sites have 5 months of data but some have 4 (1 or 2 sites in 2022)
for(i in 1:length(pct_overlap_list)){
  ifelse( length(pct_overlap_list[[i]]) == 5,
          names(pct_overlap_list[[i]]) <- c("June", "July", "August", "September", "October"),
          names(pct_overlap_list[[i]]) <- c("July", "August", "September", "October") )
}




#collate MONTHLY pct_overlap results
#make a list to store things
pct_overlap_summary <- list()

#loop across all sites and collapse the dfs per month into one df for the site
for(i in 1:length(pct_overlap_list)){
  
  #for all 12 sites
  summary <- do.call("rbind", pct_overlap_list[[i]])
  pct_overlap_summary[[i]] <- summary
}

#name the 12 1st order elements as their sites
names(pct_overlap_summary) <- names(site_month_list)

## make pct_overlap_summary into freiggein huge df
pct_overlap_summary22 <- do.call(rbind.data.frame, pct_overlap_summary)
pct_overlap_summary22$year <- 2022
row.names(pct_overlap_summary22) <- NULL
