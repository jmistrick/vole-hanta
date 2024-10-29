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


#join all the data needed for HR circles
homerange21 <- trapdata21 %>% rename(breeder = season_breeder) %>% 
  left_join(spaceusedata, by=c("year", "season", "trt", "sex", "breeder")) %>%
  left_join(centroids21, by=c("tag", "month", "site"))
#join all the data needed for HR circles
homerange22 <- trapdata22 %>% rename(breeder = season_breeder) %>% 
  left_join(spaceusedata, by=c("year", "season", "trt", "sex", "breeder")) %>%
  left_join(centroids22, by=c("tag", "month", "site"))


###---------------------------------------------


#buffer multipoint with different distances (https://stackoverflow.com/questions/47057316/st-buffer-multipoint-with-different-distance)

#subset data for just one site/month (to get code running, will eventually need to get this all into ONE GIANT LOOP)
#just pull the vole ID, radius, and x-y coordinates of monthly centroid for now
data <- homerange21 %>% filter(month=="june", site=="ketunpesa") %>% ungroup() %>% 
  mutate(id = tag) %>%
  select(x, y, rad_95_trap, id)
#make a matrix of just the coordinates
coords <- matrix(c(data$x,data$y), ncol = 2)
#convert to sf object
tt <- st_as_sf(data, coords=c("x","y"))
#buffer each point with the corresponding radius
tt_buff <- sf::st_buffer(tt, data$rad_95_trap)

#plot with ggplot (can ploth with base R too - plot() )
ggplot(tt_buff) + geom_sf(aes(fill=rad_95_trap, alpha=0.5)) +
  geom_sf_label(aes(label=id)) #include label with vole ID


####----------------------

# #the following code just needs the 'data' file which is currently just one site/month
# data <- homerange21 %>% filter(month=="june", site=="ketunpesa") %>% ungroup() %>% 
#   mutate(id = tag) %>%
#   select(x, y, rad_95_trap, id)

#basically here I want to make each vole's HR its own polygon, then combine all those polygons together into a single sf
#NOT using a multipoint or multipolygon because that acts like all the circles are one layer and doesn't "see" overlaps the same way (I think)

#empty list to hold results
poly_list <- list()

#loop over df rows #https://campus.datacamp.com/courses/intermediate-r-for-finance/loops-3?ex=10
for(i in 1:nrow(data)){
  
  #pull each row (vole) as its own df
  df  <- as.data.frame(data[i,])
  coords <- matrix(c(df$x,df$y), ncol = 2) #centroid coordinates
  tt <- st_as_sf(df, coords=c("x","y"))
  sf <- st_buffer(tt, df$rad_95_trap) #buffer with radius
  sf$area = st_area(sf) #add area of polygon IN TRAP UNITS
  
  poly_list[[i]] <- sf #save result as one item of the list (list will have one item for each vole that site/month)
  
}

#convert list of sf's to one sf 
  #(https://stackoverflow.com/questions/51312935/convert-a-list-of-sf-objects-into-one-sf)
single_sf <- do.call(rbind, poly_list)

#and now some slick code someone else wrote to compute the pairwise overlap between all the polygons, and have directional % overlap
#https://stackoverflow.com/questions/70009412/how-to-compute-all-pairwise-interaction-between-polygons-and-the-the-percentage
# nOTE 'st_overlaps' will not capture polygons contained within another, for that you want 'st_intersects'
#output is a matrix which is NOT SYMMETRICAL
  #basically it's calling the sf twice and pairing all the circles against e/o or something
pctoverlap_mat <- sf::st_intersection(single_sf, single_sf) %>% 
  dplyr::mutate(
    area = sf::st_area(.),
    proportion = area / area.1 ) %>% #the area.1 is because both things called single_sf have an 'area' column ;)
  tibble::as_tibble() %>%
  dplyr::select(
    id_1 = id,
    id_2 = id.1,
    proportion, ) %>% 
  tidyr::pivot_wider(
    names_from = id_1,
    values_from = proportion,
    values_fill = 0 )

#pivot the matrix wider to be a long format table with focal vole, neighbor vole and...
  #the amount of the focal's HR that is overlapped by the neighbor
## I THINK I did things correctly so the right voles are focal vs neighbor, might be good to triple check though
#this code also drops voles that do not overlap and removes self overlaps
    # depending on how I make the edge list to make this into a network, might want the 0's back but definitely not the self-loops
pctoverlap_mat %>% rename(focal = id_2) %>%
  pivot_longer(!focal, names_to = "neighbor", values_to="pct_overlap") %>%
  filter(pct_overlap>0) %>%
  mutate(neighbor=as.numeric(neighbor)) %>% filter(neighbor != focal) %>%
  select(focal, neighbor, pct_overlap)



