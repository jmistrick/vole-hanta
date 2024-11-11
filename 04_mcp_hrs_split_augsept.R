# RESURRECTED! for NEW! VoleHanta analysis Nov 2024

#load libraries
library(here)
library(tidyverse)
# library(rgdal)
library(sp)
library(sf)
library(adehabitatHR)
library(scales) #used for pretty plots in the loop - helps make polygons partly transparent using the alpha argument
library(Hmisc) #has capitalize() function


# library(kableExtra)
#these are supposed to allow you to print pngs of Kables but they don't work
# library(magick)
# library(webshot)
# webshot::install_phantomjs()

#clear environment
rm(list = ls())

#load the fulltrap dataset (make sure it's the most recent version)
fulltrap <- readRDS(here("fulltrap21_11.11.24.rds"))

#################################################################################
################################  prep code  ####################################
#################################################################################

#subset fulltrap into May-Aug & Sept-Oct for use in MCP analysis - only animals captured at least 5 times THAT SEASON
summer_mcp_trap <- fulltrap %>%
  filter(season=="summer") %>% #summer only
  ##something is wrong and caps_per_life isn't always correct... FIX THIS!!
  dplyr::select(!c(caps_per_life, traps_per_life)) %>%
  group_by(tag) %>% mutate(relocs = length(tag)) %>% 
  filter(relocs >= 5) %>% #filter for at least 5 captures
  ungroup()

fall_mcp_trap <- fulltrap %>%
  filter(season == "fall") %>% #fall only
  filter(caps_per_year >= 5) %>% #filter for at least 5 captures (these captures EXCLUDING may)
  ungroup()

# #how many residents per season?
# n_distinct(summer_mcp_trap$tag) #45
# n_distinct(fall_mcp_trap$tag) #55
# #how many voles were resident in both seasons?
# length(intersect(summer_mcp_trap$tag, fall_mcp_trap$tag)) #8
# intersect(summer_mcp_trap$tag, fall_mcp_trap$tag) #who were they?
# #check - how many in summer but not in fall
# length(setdiff(summer_mcp_trap$tag, fall_mcp_trap$tag)) #37
# #check - how many in fall but not in summer
# length(setdiff(fall_mcp_trap$tag, summer_mcp_trap$tag)) #47

#change the x,y coordinates of the trap to easting, northing (required for MCP)
#use A1 at Puro because TREEBRIDGE! is the best - easting: 398049	northing: 6791091
summer_mcp_trap <- summer_mcp_trap %>%
  mutate(x.jitter = jitter(x),
         y.jitter = jitter(y)) %>% #jitter points to keep trap-happy voles
  #easting (x coordinate) should be A1.easting - [(trap.letter.as.number-1)*10] 
  mutate(easting = ( 398049 - ((x.jitter)-1)*10) ) %>% #minus because A1 is on the E side of the grid, so coordinates go W
  relocate(easting, .after = x.jitter) %>%
  #northing (y coordinate) 
  #if remainder of x/2 is 0 (if x is even), multiply y by 2 - if else, multiply by 2 then subtract 1
  mutate(northing = ifelse(
    (x %% 2) == 0, ((y.jitter-1)*10) + 6766143, ((y-1)*10) + 6766143)) %>%
  relocate(northing, .after = y.jitter) %>%
  dplyr::select(-c(x, y, x.jitter, y.jitter)) #remove dummy trap coordinates


fall_mcp_trap <- fall_mcp_trap %>%
  mutate(x.jitter = jitter(x),
         y.jitter = jitter(y)) %>% #jitter points to keep trap-happy voles
  #easting (x coordinate) should be A1.easting - [(trap.letter.as.number-1)*10] 
  mutate(easting = ( 398049 - ((x.jitter)-1)*10) ) %>% #minus because A1 is on the E side of the grid, so coordinates go W
  relocate(easting, .after = x.jitter) %>%
  #northing (y coordinate) 
  #if remainder of x/2 is 0 (if x is even), multiply y by 2 - if else, multiply by 2 then subtract 1
  mutate(northing = ifelse(
    (x %% 2) == 0, ((y.jitter-1)*10) + 6766143, ((y-1)*10) + 6766143)) %>%
  relocate(northing, .after = y.jitter) %>%
  dplyr::select(-c(x, y, x.jitter, y.jitter)) #remove dummy trap coordinates



##############################################################################
############### create MCP_traits df for summer/fall #########################
##############################################################################

#trim down MCP_trap dataset to just one entry per animal, summary traits for that season 
summer_mcp_traits <- summer_mcp_trap %>%
  group_by(tag) %>% arrange(occ.sess) %>% slice(1) %>% #keep only first entry per animal
  dplyr::select(year, season, site, trt, tag, sex, season_breeder) %>%
  ungroup()

fall_mcp_traits <- fall_mcp_trap %>%
  group_by(tag) %>% arrange(occ.sess) %>% slice(1) %>% #keep only first entry per animal
  dplyr::select(year, season, site, trt, tag, sex, season_breeder) %>%
  ungroup()


#pull out tag, sex for summer and fall (for focal/neighbor sex ID)
summer_mcp_sex <- summer_mcp_traits %>%
  dplyr::select(c(tag, sex))

fall_mcp_sex <- fall_mcp_traits %>%
  dplyr::select(c(tag, sex))



###############################################  end prep code   ################################################

######### ATTENTION ATTENTION ATTENTION - the following code will create the CAs and calculate all the metrics
      ### unless you have updated the code, you don't need to run it all and can just skip all this
      ### at the end of this, there is code to pull the rds files for summer_ fall_ and summary dfs


################################################################################################################
######### Loop to Create SpatialPointsDF for each site & polygons for each resident vole in 2021 ###############
### results are: one list of 12 spdf objects, one list of 12 dfs with polygon areas, .png file for each site ###
################################################################################################################

#split() mcp_trap into ... A LIST! by site
summer_mcp_list <- split(summer_mcp_trap, f = summer_mcp_trap$site)
fall_mcp_list <- split(fall_mcp_trap, f = fall_mcp_trap$site)

####### SUMMER #########

#initiate empty dfs for results
summer_sitespdf_list <- list()
summer_CA_list <- list()

for(i in 1:length(summer_mcp_list)){
  
  GRID <- summer_mcp_list[[i]] #pull that month's data

  #Create Spatial Points for all capture locations and assign vole IDs to each location
  #https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
  
  #SPDF needs 3 things: coordinates, data, and coordinate reference system
  coords <- GRID[ , c("easting", "northing")]
  data <- GRID[ ,c("tag", "sex", "season_breeder", "site", "season")]
  crs <- CRS("+proj=utm +zone=35 +units=m +no_defs +ellps=GRS80") #might not be necessary since using UTM coordinates
  #combine all three into a spdf
  spdf <- SpatialPointsDataFrame(coords = coords, data = data, proj4string = crs)

  #create MCPs for our new dataset "merge" by individual animal ID
  cp <- mcp(spdf[,1], percent=100) #(95% is the default) #column 1 has the vole IDs
  #write cp to the sitespdf_list
  summer_sitespdf_list[[i]] <- cp

  ## to get area of the bounding polygon for each MCP (ie 'CA' size)
  #area expressed in hectares if map projection was UTM
  mcpdata <- as.data.frame(cp)
  #write the ID/CA area df to the CAarea_list
  summer_CA_list[[i]] <- mcpdata

  # Plot the home ranges and save to a png file in my working directory
  png(paste("SUMMER_caparea_", names(summer_mcp_list)[i], "_2021", ".png", sep = ""))

  plot(cp, col = cp@data$id, pch = 16, main = paste(capitalize(names(summer_mcp_list)[i]), "SUMMER 21 Capture Areas of Resident Voles", sep = " "))
  plot(cp, col = alpha(1:10, 0.5), add = TRUE)
  #add crosshairs for each recapture location
  plot(spdf, add=TRUE)

  dev.off()

}

#name the 12 1st order elements as their sites
names(summer_sitespdf_list) <- names(summer_mcp_list)
names(summer_CA_list) <- names(summer_mcp_list)

## this has created two lists and created png files for each site
  #first list is basically the 'cp' file for each site - it's the full spatialpointsdf of all polygons
  #second list is the cp@data file that gives ID of each animal and the area of their CA

### for later analysis - make CA_list into df
summer_CA_summary <- do.call(rbind.data.frame, summer_CA_list)

#clean up the df
summer_CA_summary <- summer_CA_summary %>%
  rownames_to_column("name") %>% #row names are the sites, make that a column
  separate(name, c("site", NA)) %>% #separate the site part from the index and get rid of the index
  mutate(site = as.factor(site),
         area.m = area*10000) %>% #make site a factor #convert area from hectares to m^2
  rename(tag = id)

#join mcp_traits to CA_summary
summer_CA_summary <- summer_CA_summary %>%
  left_join(y=summer_mcp_traits, by=c("site", "tag")) %>%
  dplyr::select(year, site, season, trt, tag, sex, season_breeder, area.m) #rearrange and drop area (in hectares)

# SAVE IT
summer_CA_summary21 <- summer_CA_summary
saveRDS(summer_CA_summary21, "summer_CA_summary21.rds")

##########################################################################################
##########################################################################################

####### FALL #########

#initiate empty dfs for results
fall_sitespdf_list <- list()
fall_CA_list <- list()

for(i in 1:length(fall_mcp_list)){
  
  GRID <- fall_mcp_list[[i]] #pull that month's data
  
  #Create Spatial Points for all capture locations and assign vole IDs to each location
  #https://stackoverflow.com/questions/32583606/create-spatialpointsdataframe
  
  #SPDF needs 3 things: coordinates, data, and coordinate reference system
  coords <- GRID[ , c("easting", "northing")]
  data <- GRID[ ,c("tag", "sex", "season_breeder", "site", "season")]
  crs <- CRS("+proj=utm +zone=35 +units=m +no_defs +ellps=GRS80") #might not be necessary since using UTM coordinates
  #combine all three into a spdf
  spdf <- SpatialPointsDataFrame(coords = coords, data = data, proj4string = crs)
  
  #create MCPs for our new dataset "merge" by individual animal ID
  cp <- mcp(spdf[,1], percent=100) #(95% is the default) #column 1 has the vole IDs
  #write cp to the sitespdf_list
  fall_sitespdf_list[[i]] <- cp
  
  ## to get area of the bounding polygon for each MCP (ie 'CA' size)
  #area expressed in hectares if map projection was UTM
  mcpdata <- as.data.frame(cp)
  #write the ID/CA area df to the CAarea_list
  fall_CA_list[[i]] <- mcpdata
  
  # Plot the home ranges and save to a png file in my working directory
  png(paste("FALL_caparea_", names(fall_mcp_list)[i], "_2021", ".png", sep = ""))
  
  plot(cp, col = cp@data$id, pch = 16, main = paste(capitalize(names(fall_mcp_list)[i]), "FALL 21 Capture Areas of Resident Voles", sep = " "))
  plot(cp, col = alpha(1:10, 0.5), add = TRUE)
  #add crosshairs for each recapture location
  plot(spdf, add=TRUE)
  
  dev.off()
  
}

#name the 12 1st order elements as their sites
names(fall_sitespdf_list) <- names(fall_mcp_list)
names(fall_CA_list) <- names(fall_mcp_list)

## this has created two lists and created png files for each site
#first list is basically the 'cp' file for each site - it's the full spatialpointsdf of all polygons
#second list is the cp@data file that gives ID of each animal and the area of their CA

### for later analysis - make CA_list into df
fall_CA_summary <- do.call(rbind.data.frame, fall_CA_list)

#clean up the df
fall_CA_summary <- fall_CA_summary %>%
  rownames_to_column("name") %>% #row names are the sites, make that a column
  separate(name, c("site", NA)) %>% #separate the site part from the index and get rid of the index
  mutate(site = as.factor(site),
         area.m = area*10000) %>% #make site a factor #convert area from hectares to m^2
  rename(tag = id)

#join mcp_traits to CA_summary
fall_CA_summary <- fall_CA_summary %>%
  left_join(y=fall_mcp_traits, by=c("site", "tag")) %>%
  dplyr::select(year, site, season, trt, tag, sex, season_breeder, area.m) #rearrange and drop area (in hectares)

# SAVE IT
fall_CA_summary21 <- fall_CA_summary
saveRDS(fall_CA_summary21, "fall_CA_summary21.rds")

###################################### END SPDF & CA polygons ####################################



#combine summer and fall CA data
CA_summary21 <- rbind(summer_CA_summary21, fall_CA_summary21)

#table of mean capture area (m^2) per functional group
CA_summary21 %>% group_by(trt, sex, season_breeder) %>%
  summarise(avg_CA = mean(area.m),
            n = length(tag))

library(ggridges)

#plot histogram of capture area for each fxnl group by trt, season
CA_summary21 %>% 
  unite("fxnl", sex, season_breeder, sep="_", remove=FALSE) %>%
  ggplot(aes(x=area.m, y=fxnl, fill=fxnl)) +
  xlim(-10, 3000) + #removes one male outlier at 6000m^2 for better visualizing
  geom_density_ridges(stat = "binline", bins = 20, scale = 0.95, draw_baseline = FALSE) +
  facet_grid(season~trt)
























# ##############################################################################################################
# ################################# in a loop for all 12 sites in 2021 #########################################
# ####### Calculate area of overlap & % overlap (how much of each vole's CA overlaps with another's) ###########
# ##############################################################################################################
# 
# summer_pct_overlap_list <- list()
# 
# for(i in 1:length(summer_sitespdf_list)){
# 
#   print(i)
# 
#   cp <- summer_sitespdf_list[[i]]
# 
#   newcp <- st_as_sfc(summer_sitespdf_list[[i]])
# 
#   #calculate amount of overlap between polygons (in hectares)
#   l <- lapply( newcp, function(x) {lapply(newcp, function(y) st_intersection(x,y) %>% st_area() ) })
# 
#   mat <- matrix(unlist(l), ncol = length(newcp), byrow = TRUE)
#   diag(mat) <- NA #set the diagonal to NA so I don't get confused later
#   colnames(mat) <- cp@data$id
#   rownames(mat) <- cp@data$id
#   # mat #shows area of overlap in hectares (because UTM)
# 
#   df <- as.data.frame(mat) #convert matrix to df
#   df <- tibble::rownames_to_column(df, "focal") #row names become column
#   df <- df %>% pivot_longer(!focal, names_to="neighbor", values_to="area_overlap")
#   df <- df %>% filter(area_overlap > 0) #remove 0 and NA values
#   # df #is a df version of the matrix, shows area of overlap in hectares
# 
#   #calculate percent overlap (a directed measure) between polygons
#   l2 <- lapply( newcp, function(x) {lapply(newcp, function(y) st_intersection( x, y ) %>% st_area() / st_area(x)  ) })
# 
#   mat2 <- matrix(unlist(l2), ncol = length(newcp), byrow = TRUE)
#   diag(mat2) <- NA #set the diagonal to NA so I don't get confused later
#   rownames(mat2) <- cp@data$id
#   colnames(mat2) <- cp@data$id
#   # mat2 #shows percent overlap: how much of rowID overlaps with columnID
# 
#   df2 <- as.data.frame(mat2) #convert matrix to df
#   df2 <- tibble::rownames_to_column(df2, "focal") #row names become column
#   df2 <- df2 %>% pivot_longer(!focal, names_to="neighbor", values_to="pct_overlap")
#   df2 <- df2 %>% filter(pct_overlap > 0) #remove 0 and NA values
#   # df2 #is a df version of the matrix, shows percent overlap
# 
#   df_full <- left_join(df, df2) #combine area of overlap and percent overlap into one df
#   # df_full
# 
#   summer_pct_overlap_list[[i]] <- df_full #and write that df for a given site as a list item
# 
# }
# 
# #name the 12 1st order elements as their sites
# names(summer_pct_overlap_list) <- names(summer_mcp_list)
# #this is a list of 12 sites, each item has a df matrix of the amt and percentage overlap between residents
# 
# ## make pct_overlap_list into freiggein huge df
# summer_pct_overlap_summary <- do.call(rbind.data.frame, summer_pct_overlap_list)
# 
# #clean up the df
# summer_pct_overlap_summary <- summer_pct_overlap_summary %>%
#   rownames_to_column("name") %>% #row names are the sites, make that a column
#   separate(name, c("site", NA)) %>% #separate the site part from the index and get rid of the index
#   mutate(site = as.factor(site)) %>% #make site a factor
#   mutate(across(c(area_overlap, pct_overlap), ~round(., digits=2)))
# 
# #join the grid_trts to pct_overlap_summary
# summer_pct_overlap_summary <- summer_pct_overlap_summary %>%
#   left_join(y=grid_trts, by = "site")
# 
# #join the sex of focal animal, neighbor animal
# new <- left_join(summer_pct_overlap_summary, summer_mcp_sex, by = c("focal" = "tag")) %>%
#   rename(focal_sex = sex) %>%
#   relocate(focal_sex, .after=focal)
# summer_pct_overlap_summary <- left_join(new, summer_mcp_sex, by=c("neighbor" = "tag")) %>%
#   rename(neighbor_sex = sex) %>%
#   relocate(neighbor_sex, .after=neighbor)
# #make a new column for focal-neighbor sexes
# summer_pct_overlap_summary <- summer_pct_overlap_summary %>%
#   unite(f_n_sex, focal_sex, neighbor_sex, sep = "_", remove = FALSE) %>% #keep the original columns
#   relocate(f_n_sex, .before=focal) # %>% drop_na() #remove sex=NA (there are none)
# #make f_n_sex a factor, combine M_F and F_M
# summer_pct_overlap_summary$f_n_sex <- as.factor(summer_pct_overlap_summary$f_n_sex)
# levels(summer_pct_overlap_summary$f_n_sex) <- c("F_F", "mixed", "mixed", "M_M")
# 
# ######## NICE - this now has the % overlap and area of overlap for all voles that did overlap in their CA
# ### HOWEVER - we've now lost any data on any individuals that did not overlap, which we'll need to do % overlapping
# ## hence this next bit of code:
# 
# #############################################################################################################
# ############# Calculate percent overlapping (how many of the resident voles have any CA overlap) ############
# #############################################################################################################
# 
# ## pulling from code above where I made cp@data into a big df for all sites (ie the CA_summary df)
# #summarise count of 'resident' voles per site
# n.CAs <- summer_CA_summary %>%
#   group_by(site) %>%
#   summarise(unique(tag)) %>%
#   count() %>%
#   ungroup()
# 
# #summarise number of focal voles from pct_overlap_summary (ie number per site with overlaps)
# n.overlaps <- summer_pct_overlap_summary %>% group_by(site) %>%
#   summarise(unique(focal)) %>%
#   count() %>%
#   rename(n.overlapping = n) %>%
#   ungroup()
# 
# summer_pct_overlapping_summary <- left_join(n.CAs, n.overlaps, by="site") #join the two summaries
# 
# #replace NA with 0 for n.overlapping at Kuoppa, Radio
#   #(there were no overlapping voles there so it didn't have a row in n.overlaps so it got NA in 'pct_ovlpping_summary')
# summer_pct_overlapping_summary <- summer_pct_overlapping_summary %>% replace(is.na(.), 0)
# #add a column for percent overlapping (how many voles have a CA that overlaps with another)
# summer_pct_overlapping_summary <- summer_pct_overlapping_summary %>%
#   mutate(pct.overlapping = round(n.overlapping/n, 2))
# 
# #join the grid_trts to pct_overlapping_summary
# summer_pct_overlapping_summary <- summer_pct_overlapping_summary %>%
#   left_join(y=grid_trts, by = "site")
# 
# 
# ###############################################################################################
# ###############################################################################################
# 
# ####################### CALCULATE PERCENT OVERLAP, PERCENT OVERLAPPING #######################
# 
# ## repeat for fall
# 
# fall_pct_overlap_list <- list()
# 
# for(i in 1:length(fall_sitespdf_list)){
# 
#   print(i)
# 
#   cp <- fall_sitespdf_list[[i]]
# 
#   newcp <- st_as_sfc(fall_sitespdf_list[[i]])
# 
#   #calculate amount of overlap between polygons (in hectares)
#   l <- lapply( newcp, function(x) {lapply(newcp, function(y) st_intersection(x,y) %>% st_area() ) })
# 
#   mat <- matrix(unlist(l), ncol = length(newcp), byrow = TRUE)
#   diag(mat) <- NA #set the diagonal to NA so I don't get confused later
#   colnames(mat) <- cp@data$id
#   rownames(mat) <- cp@data$id
#   # mat #shows area of overlap in hectares (because UTM)
# 
#   df <- as.data.frame(mat) #convert matrix to df
#   df <- tibble::rownames_to_column(df, "focal") #row names become column
#   df <- df %>% pivot_longer(!focal, names_to="neighbor", values_to="area_overlap")
#   df <- df %>% filter(area_overlap > 0) #remove 0 and NA values
#   # df #is a df version of the matrix, shows area of overlap in hectares
# 
#   #calculate percent overlap (a directed measure) between polygons
#   l2 <- lapply( newcp, function(x) {lapply(newcp, function(y) st_intersection( x, y ) %>% st_area() / st_area(x)  ) })
# 
#   mat2 <- matrix(unlist(l2), ncol = length(newcp), byrow = TRUE)
#   diag(mat2) <- NA #set the diagonal to NA so I don't get confused later
#   rownames(mat2) <- cp@data$id
#   colnames(mat2) <- cp@data$id
#   # mat2 #shows percent overlap: how much of rowID overlaps with columnID
# 
#   df2 <- as.data.frame(mat2) #convert matrix to df
#   df2 <- tibble::rownames_to_column(df2, "focal") #row names become column
#   df2 <- df2 %>% pivot_longer(!focal, names_to="neighbor", values_to="pct_overlap")
#   df2 <- df2 %>% filter(pct_overlap > 0) #remove 0 and NA values
#   # df2 #is a df version of the matrix, shows percent overlap
# 
#   df_full <- left_join(df, df2) #combine area of overlap and percent overlap into one df
#   # df_full
# 
#   fall_pct_overlap_list[[i]] <- df_full #and write that df for a given site as a list item
# 
# }
# 
# #name the 12 1st order elements as their sites
# names(fall_pct_overlap_list) <- names(fall_mcp_list)
# #this is a list of 12 sites, each item has a df matrix of the amt and percentage overlap between residents
# 
# ## make pct_overlap_list into freiggein huge df
# fall_pct_overlap_summary <- do.call(rbind.data.frame, fall_pct_overlap_list)
# 
# #clean up the df
# fall_pct_overlap_summary <- fall_pct_overlap_summary %>%
#   rownames_to_column("name") %>% #row names are the sites, make that a column
#   separate(name, c("site", NA)) %>% #separate the site part from the index and get rid of the index
#   mutate(site = as.factor(site)) %>% #make site a factor
#   mutate(across(c(area_overlap, pct_overlap), ~round(., digits=2)))
# 
# #join the grid_trts to pct_overlap_summary
# fall_pct_overlap_summary <- fall_pct_overlap_summary %>%
#   left_join(y=grid_trts, by = "site")
# 
# 
# #join the sex of focal animal, neighbor animal
# new <- left_join(fall_pct_overlap_summary, fall_mcp_sex, by = c("focal" = "tag")) %>%
#   rename(focal_sex = sex) %>%
#   relocate(focal_sex, .after=focal)
# fall_pct_overlap_summary <- left_join(new, fall_mcp_sex, by=c("neighbor" = "tag")) %>%
#   rename(neighbor_sex = sex) %>%
#   relocate(neighbor_sex, .after=neighbor)
# #make a new column for focal-neighbor sexes
# fall_pct_overlap_summary <- fall_pct_overlap_summary %>%
#   unite(f_n_sex, focal_sex, neighbor_sex, sep = "_", remove = FALSE) %>% #keep the original columns
#   relocate(f_n_sex, .before=focal) # %>% drop_na() #remove sex=NA (there are none)
# #make f_n_sex a factor, combine M_F and F_M
# fall_pct_overlap_summary$f_n_sex <- as.factor(fall_pct_overlap_summary$f_n_sex)
# levels(fall_pct_overlap_summary$f_n_sex) <- c("F_F", "mixed", "mixed", "M_M")
# 
# ######## NICE - this now has the % overlap and area of overlap for all voles that did overlap in their CA
# ### HOWEVER - we've now lost any data on any individuals that did not overlap which we'll need to do % overlapping
# ## hence this next bit of code:
# 
# ############################################################################################################
# ############ Calculate percent overlapping (how many of the resident voles have any CA overlap) ############
# ############################################################################################################
# 
# ## pulling from code above where I made cp@data into a big df for all sites (ie the CA_summary df)
# #summarise count of 'resident' voles per site
# n.CAs <- fall_CA_summary %>%
#   group_by(site) %>%
#   summarise(unique(tag)) %>%
#   count() %>%
#   ungroup()
# 
# #summarise number of focal voles from pct_overlap_summary (ie number per site with overlaps)
# n.overlaps <- fall_pct_overlap_summary %>% group_by(site) %>%
#   summarise(unique(focal)) %>%
#   count() %>%
#   rename(n.overlapping = n) %>%
#   ungroup()
# 
# fall_pct_overlapping_summary <- left_join(n.CAs, n.overlaps, by="site") #join the two summaries
# 
# #replace NA with 0 for n.overlapping
# #(there were no overlapping voles there so it didn't have a row in n.overlaps so it got NA in 'pct_ovlpping_summary')
# fall_pct_overlapping_summary <- fall_pct_overlapping_summary %>% replace(is.na(.), 0)
# #add a column for percent overlapping (how many voles have a CA that overlaps with another)
# fall_pct_overlapping_summary <- fall_pct_overlapping_summary %>% mutate(pct.overlapping = round(n.overlapping/n, 2))
# 
# #join the grid_trts to pct_overlapping_summary
# fall_pct_overlapping_summary <- fall_pct_overlapping_summary %>%
#   left_join(y=grid_trts, by = "site")
# 
# 
# #######################################################################################
# ##########################   end percent overlap/lapping  #############################
# #######################################################################################
# 
# 
# 
# #combine summer and fall for a single measure
# 
# #add season identifier to fall_ or summer_ df
# fall_CA_summary <- fall_CA_summary %>%
#   mutate(season = "fall") %>% relocate(season, .after='site')
# fall_pct_overlap_summary <- fall_pct_overlap_summary %>%
#   mutate(season = "fall") %>% relocate(season, .after='site')
# fall_pct_overlapping_summary <- fall_pct_overlapping_summary %>%
#   mutate(season = "fall") %>% relocate(season, .after='site')
# summer_CA_summary <- summer_CA_summary %>%
#   mutate(season = "summer") %>% relocate(season, .after='site')
# summer_pct_overlap_summary <- summer_pct_overlap_summary %>%
#   mutate(season = "summer") %>% relocate(season, .after='site')
# summer_pct_overlapping_summary <- summer_pct_overlapping_summary %>%
#   mutate(season = "summer") %>% relocate(season, .after='site')
# 
# #combine summer and fall dfs, make season a factor with summer first
# CA_summary <- rbind(summer_CA_summary, fall_CA_summary) %>%
#   mutate(season = factor(season, levels=c("summer", "fall")))
# pct_overlap_summary <- rbind(summer_pct_overlap_summary, fall_pct_overlap_summary) %>%
#   mutate(season = factor(season, levels=c("summer", "fall")))
# pct_overlapping_summary <- rbind(summer_pct_overlapping_summary, fall_pct_overlapping_summary) %>%
#   mutate(season = factor(season, levels=c("summer", "fall")))
# 
# #create season_trt variable; sex_trt / sex_trt_season; fnsex_trt / fnsex_trt_season
# CA_summary <- CA_summary %>%
#   unite(season_trt, season, trt, sep="_", remove=FALSE) %>%
#   unite(sex_trt, sex, trt, sep="_", remove=FALSE) %>%
#   unite(season_sex_trt, season, sex_trt, sep="_", remove=FALSE) %>%
#   mutate(season_trt = factor(season_trt),
#          sex_trt = factor(sex_trt),
#          season_sex_trt = factor(season_sex_trt))
# pct_overlap_summary <- pct_overlap_summary %>%
#   unite(season_trt, season, trt, sep="_", remove=FALSE) %>%
#   unite(fnsex_trt, f_n_sex, trt, sep="_", remove=FALSE) %>%
#   unite(season_fnsex_trt, season, fnsex_trt, remove=FALSE) %>%
#   mutate(season_trt = factor(season_trt),
#          fnsex_trt = factor(fnsex_trt),
#          season_fnsex_trt = factor(season_fnsex_trt))
# pct_overlapping_summary <- pct_overlapping_summary %>%
#   unite(season_trt, season, trt, sep="_", remove=FALSE) %>%
#   mutate(season_trt = factor(season_trt))
# 
# 
# #save all these files to RDS so I can pull them for 04_mcp...viz.R file
# saveRDS(summer_CA_summary, file = here("summer_CA_summary_07.01.22.rds"))
# saveRDS(summer_pct_overlap_summary, file = here("summer_pct_overlap_summary_07.01.22.rds"))
# saveRDS(summer_pct_overlapping_summary, file = here("summer_pct_overlapping_summary_07.01.22.rds"))
# 
# saveRDS(fall_CA_summary, file = here("fall_CA_summary_07.01.22.rds"))
# saveRDS(fall_pct_overlap_summary, file = here("fall_pct_overlap_summary_07.01.22.rds"))
# saveRDS(fall_pct_overlapping_summary, file = here("fall_pct_overlapping_summary_07.01.22.rds"))
# 
# saveRDS(CA_summary, file = here("CA_summary_07.01.22.rds"))
# saveRDS(pct_overlap_summary, file = here("pct_overlap_summary_07.01.22.rds"))
# saveRDS(pct_overlapping_summary, file = here("pct_overlapping_summary_07.01.22.rds"))
# 
# ######################## END DF CREATION #############################



########### TO SKIP ALL THE ABOVE AND JUST PULL THE FILES ############
# plotting code for CAs hasn't been updated since before 1 July 2022 #

#load MCP data files from the rdata files
summer_CA_summary <- readRDS(file = here("summer_CA_summary.rds"))
summer_pct_overlap_summary <- readRDS(file = here("summer_pct_overlap_summary.rds"))
summer_pct_overlapping_summary <- readRDS(file = here("summer_pct_overlapping_summary.rds"))

fall_CA_summary <- readRDS(file = here("fall_CA_summary.rds"))
fall_pct_overlap_summary <- readRDS(file = here("fall_pct_overlap_summary.rds"))
fall_pct_overlapping_summary <- readRDS(file = here("fall_pct_overlapping_summary.rds"))

CA_summary <- readRDS(file = here("CA_summary.rds"))
pct_overlap_summary <- readRDS(file = here("pct_overlap_summary.rds"))
pct_overlapping_summary <- readRDS(file = here("pct_overlapping_summary.rds"))



########################## BEGIN ANALYSIS ############################






######################################################################################



plot <- ageclass_season_summary %>%
  drop_na(age_class) %>%
  ggplot(aes(x=age_class, y=mean_n, group=trt, color=trt)) +
  geom_point(size=3, position=position_dodge(.2)) +
  geom_errorbar(aes(ymin=mean_n-sd_n, ymax=mean_n+sd_n), width=.2, size=1,
                position=position_dodge(.2)) +
  labs(y = "Mean Number of Resident Voles per Site",
       x = "Age Class",
       color = "Treatment:") +
  scale_x_discrete(labels=c("sub-adult" = "Sub-adult (\u2264 17g)",
                            "adult" = "Adult (> 17g)")) +
  scale_color_manual(values = c("#B2DF8A", "#33A02C", "#CAB2D6", "#6A3D9A"), 
                     name = "Treatment",
                     labels = c("Unfed-Control", "Unfed-Deworm", "Fed-Control", "Fed-Deworm")) +
  theme(legend.position="bottom",
        axis.title = element_text(size=13),
        axis.text = element_text(size=11),
        legend.text = element_text(size=11),
        strip.text = element_text(size=11)) +
  facet_wrap(~season, labeller = season.labs)

ggsave("ageclasscountseason.png", plot=plot, height=5, width=10, units=c("in"), dpi=600)



#######################################end 12.1.22###############################################


















# ######### FROM JANINE OF OCTOBER 3 2022: PROBABALY STILL NEED TO DEAL WITH THIS(below) #############
# 
# 
# 
# 
# 
# ####### RESIDENTS: HOW MANY overlaps and HOW MUCH do they overlap? ############
# 
# #create a df with every resident and their # of overlaps and avg % overlap each season
# 
# #trim CA to only the necessary columns 
# CA_trim <- CA_summary %>% dplyr::select(site, trt, season, tag, sex, area)
# #join CA and pct_overlap - but left_join pct_overlap to CA so we keep all the animals that had CAs but didn't have overlap
# #(if we left_join CA to pct_overlap we lose entries for animals that had no overlaps)
# #then: fill area_overlap and pct_overlap = 0 for animals will no CA overlap in a given season
# #calculate n_overlaps per focal per season, and avg pct overlap per focal per season
# #keep only one summary entry per focal vole per season (103 total)
# CA_and_overlap_summary <- left_join(CA_trim, pct_overlap_summary, by=c("season"="season", "site"="site", "trt"="trt", 
#                                                                        "tag"="focal", "sex"="focal_sex")) %>%
#   relocate(area, .before="area_overlap") %>%
#   relocate(tag, .after=f_n_sex) %>%
#   relocate(sex, .after=tag) %>%
#   dplyr::select(!c(season_trt, season_fnsex_trt, fnsex_trt, food_trt, helm_trt)) %>%
#   replace_na(list(area_overlap=0, pct_overlap=0)) %>% #if animal had no overlaps, change NA to 0
#   group_by(tag, season) %>%
#   mutate(n_overlaps = n_distinct(neighbor, na.rm=TRUE),
#          avg_pct_overlap = round(mean(pct_overlap),2)) %>% #calculate number of unique neighbors, not counting NAs
#   relocate(c(n_overlaps, avg_pct_overlap), .after=pct_overlap) %>%
#   slice(1) %>% #keep only one entry per focal, season
#   ungroup() %>%
#   group_by(site, season) %>%
#   mutate(n_residents = n_distinct(tag)) %>% #calculate the number of residents per site
#   ungroup() %>%
#   dplyr::select(-c(f_n_sex, neighbor, neighbor_sex, area_overlap, pct_overlap)) %>% #remove overlap-specific columns
#   arrange(season, site, tag)
# 
# #write CA_and_overlap_summary to rds to add to superoverlappers -> 03_degree_distribution_analysis
# saveRDS(CA_and_overlap_summary, file = here("CA_and_overlap_summary.rds"))
# 
# 
# #plot number of overlaps by trt, season
# CA_and_overlap_summary %>% ggplot(aes(x=season, y=n_overlaps, color=trt)) +
#   geom_jitter(height=0.1, width=0.35, size=2, alpha=0.6) +
#   labs(title="Number of Overlaps per Resident Vole")
# CA_and_overlap_summary %>% ggplot(aes(x=season, y=n_overlaps, fill=trt)) +
#   geom_boxplot() +
#   labs(title="Number of Overlaps per Resident Vole")
# #plot number of overlaps by sex, trt, season
# CA_and_overlap_summary %>% 
#   drop_na(sex) %>%
#   ggplot(aes(x=season, y=n_overlaps, color=sex)) +
#   geom_jitter(height=0.1, width=0.35, size=2) +
#   # geom_boxplot() +
#   facet_wrap(~trt) +
#   labs(title="Number of Overlaps per Resident Vole")
# # #plot number of overlaps (remove all the 0s)
# # CA_and_overlap_summary %>% 
# #   filter(n_overlaps != 0) %>%
# #   ggplot(aes(x=season, y=n_overlaps, color=trt)) +
# #   geom_jitter(height=0.1, width=0.35, size=2, alpha=0.6) +
# #   # geom_boxplot() +
# #   # geom_violin() +
# #   labs(title="Number of Overlaps per Resident Vole")
# 
# #summary tables for number of overlaps
# sumstats_noverlaps <- CA_and_overlap_summary %>%
#   group_by(season, trt, sex) %>%
#   summarise(Mean = round(mean(n_overlaps),2),
#             SD = round(sd(n_overlaps),2),
#             Min = min(n_overlaps),
#             Max = max(n_overlaps)) %>%
#   mutate(season=case_when(season == "summer" ~ "Summer",
#                           season == "fall" ~ "Autumn")) %>%
#   mutate(trt=case_when(trt == "unfed_control" ~ "Unfed-Control",
#                        trt == "unfed_deworm" ~ "Unfed-Deworm",
#                        trt == "fed_control" ~ "Fed-Control",
#                        trt == "fed_deworm" ~ "Fed-Deworm")) %>%
#   rename(c(Season=season, Treatment=trt, Sex=sex))
# kbl(sumstats_noverlaps) %>%
#   kable_styling(bootstrap_options = c("striped", "condensed"), full_width = F)
# 
# #number overlaps by season, sex
# CA_and_overlap_summary %>%
#   group_by(season, sex) %>%
#   summarise(mean = round(mean(n_overlaps),2),
#             sd = round(sd(n_overlaps),2),
#             min = min(n_overlaps),
#             max = max(n_overlaps))
# 
# 
# #plot mean pct overlap by trt, season
# CA_and_overlap_summary %>% ggplot(aes(x=season, y=avg_pct_overlap, color=trt)) +
#   geom_jitter(height=0, width=0.35, size=2, alpha=0.75) +
#   labs(title="Average Percent Overlap per Vole")
# CA_and_overlap_summary %>% ggplot(aes(x=season, y=avg_pct_overlap, fill=trt)) +
#   geom_boxplot() +
#   labs(title="Average Percent Overlap per Vole")
# #plot mean pct overlap by sex, trt, season
# CA_and_overlap_summary %>% 
#   drop_na(sex) %>%
#   ggplot(aes(x=season, y=avg_pct_overlap, color=sex)) +
#   geom_jitter(height=0, width=0.35, size=2) +
#   # geom_boxplot() +
#   facet_wrap(~trt) +
#   labs(title="Average Percent Overlap per Vole")
# # #plot mean pct overlap (remove all the 0s)
# # CA_and_overlap_summary %>% 
# #   filter(n_overlaps != 0) %>%
# #   ggplot(aes(x=season, y=avg_pct_overlap, fill=trt)) +
# #   geom_boxplot() +
# #   labs(title="Average Percent Overlap per Vole")
# 
# #summary tables for pct overlap
# sumstats_pctoverlap <- pct_overlap_summary %>%
#   group_by(season, trt) %>%
#   summarise(Mean = round(mean(pct_overlap),2),
#             SD = round(sd(pct_overlap),2),
#             Min = min(pct_overlap),
#             Max = max(pct_overlap)) %>%
#   mutate(season=case_when(season == "summer" ~ "Summer",
#                           season == "fall" ~ "Autumn")) %>%
#   mutate(trt=case_when(trt == "unfed_control" ~ "Unfed-Control",
#                        trt == "unfed_deworm" ~ "Unfed-Deworm",
#                        trt == "fed_control" ~ "Fed-Control",
#                        trt == "fed_deworm" ~ "Fed-Deworm")) %>%
#   rename(c(Season=season, Treatment=trt))
# kbl(sumstats_pctoverlap) %>%
#   kable_styling(bootstrap_options = c("striped", "condensed"), full_width = F)
# 
# #summary tables for pct overlap with fnsex
# sumstats_pctoverlap <- pct_overlap_summary %>%
#   group_by(season, trt, f_n_sex) %>%
#   summarise(Mean = round(mean(pct_overlap),2),
#             SD = round(sd(pct_overlap),2),
#             Min = min(pct_overlap),
#             Max = max(pct_overlap)) %>%
#   mutate(season=case_when(season == "summer" ~ "Summer",
#                           season == "fall" ~ "Autumn")) %>%
#   mutate(trt=case_when(trt == "unfed_control" ~ "Unfed-Control",
#                        trt == "unfed_deworm" ~ "Unfed-Deworm",
#                        trt == "fed_control" ~ "Fed-Control",
#                        trt == "fed_deworm" ~ "Fed-Deworm")) %>%
#   mutate(f_n_sex=case_when(f_n_sex == "M_M" ~ "Male/Male",
#                            f_n_sex == "F_F" ~ "Female/Female",
#                            f_n_sex == "mixed" ~ "Male/Female")) %>%
#   rename(c(Season=season, Treatment=trt,  Sex_Pairing = f_n_sex))
# kbl(sumstats_pctoverlap) %>%
#   kable_styling(bootstrap_options = c("striped", "condensed"), full_width = F) 
# 
# ###### this doesn't have any of the 0 overlaps in it ####
# #plot pct overlap by sex, trt, season
# pct_overlap_summary %>% 
#   drop_na(focal_sex) %>%
#   ggplot(aes(x=season, y=pct_overlap, color=focal_sex)) +
#   geom_jitter(height=0, width=0.35, size=2) +
#   # geom_boxplot() +
#   facet_wrap(~trt) +
#   labs(title="Percent Overlap per Focal Vole")






