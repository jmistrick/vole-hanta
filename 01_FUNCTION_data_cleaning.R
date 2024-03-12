#code to clean data for "Vole Hanta" manuscript
#run using R version 4.3.2 "Eye Holes"

#now running on 03.04.24 versions of vole_capture_data and week_recap_data
  #these versions are the most up-to-date with Katy Wearing's corrections to all 3 years of data (2021-2023)
  #all updates were discussed and confirmed by Katy W, Jasmine Veitch, Janine M, and Dyess Harp in Jan/Feb 2024


# load packages
library(here)
library(tidyverse)
library(igraph)
library(lubridate)
library(janitor) 

#clear environment
rm(list = ls())


########### TO RUN THIS CODE IN THIS FILE ################
# #clean 2021 data
# processingdata = "vole_capture_data_03.04.24.csv"
# WRdata= "week_recap_data_03.04.24.csv"
# yr=2021
# fulltrap_output = "fulltrap21_03.04.24.rds"

#clean 2022 data
processingdata = "vole_capture_data_03.04.24.csv"
WRdata= "week_recap_data_03.04.24.csv"
yr=2022
# fulltrap_output = "fulltrap22_03.04.24.rds"
#########################################################


## This function cleans and combines the vole processing and the week recap data (input as .csv files)
## Inputs: processing data = "filename.csv" - (fed into here()) THIS SHOULD BE THE FULL 2021, 2022 CLEANED DATA
##          WRdata = "filename.csv" - same as above
##          year = number (2021 or 2022) to get the year you want
##          fulltrap_output = "filename.rds" - where to save the final 'fulltrap' file
## Output: fulltrap file with all the capture/WR data combined, plus new stuff calculated for desired year

# clean_data <- function(processingdata, WRdata, yr, fulltrap_output) {

  ###############################   ENTRY & CLEANING VOLE CAPTURE DATA   ##################################

  #load data
  voledata <- read.csv(here(processingdata))
  #march 03 2024 version - has cleaned/updated PIT tags, vole sexes (2021 and 2022 season)
    #cleaning was done on full 2021,2022,2023 capture data

  #clean names
  voledata <- voledata %>%
    clean_names %>%
    #rename column
    rename(session = sess, samp_id = id) %>%
    #mutate columns to their proper formats
    mutate(year = as.numeric(year),
           occasion = as.numeric(occasion), #as.factor(occasion),
           site = as.character(site),
           food_trt = as.factor(food_trt),
           helm_trt = as.factor(helm_trt),
           date = as_date(date, format= "%m/%d/%Y"),
           session = as.numeric(session),
           trap = as.character(trap),
           tag = as.character(tag),
           samp_id = as.numeric(samp_id), #necessary to make sure leading 0s don't mess anything up
           sex = as.factor(sex),
           ow = as.factor(ow),
           per = as.factor(per),
           nip = as.factor(nip),
           preg = as.factor(preg),
           test = as.factor(test),
           head = as.numeric(head),
           mass = as.numeric(mass),
           ticks = as.numeric(ticks),
           ear_ed = as.factor(ear_ed),
           saliva_sr = as.factor(saliva_sr),
           smear_bs = as.factor(smear_bs),
           bloodspin_bc = as.factor(bloodspin_bc),
           bloodrna_br = as.factor(bloodrna_br),
           fecalegg_fe = as.factor(fecalegg_fe),
           fecalrna_fr = as.factor(fecalrna_fr),
           deworm = as.factor(deworm),
           new = as.factor(new),
           fate = as.factor(fate),
           handler = as.factor(handler),
           notes = as.character(notes)) %>%
    #filter for just one year
    filter(year==yr) %>%
    #remove uusi 2022 data
    filter(site!="uusi") %>% #only relevant for 2022
    #remove any animals without a trap location
    filter(!is.na(trap)) %>%
    #remove any animals without a tag id
    filter(!is.na(tag)) %>%
    #remove animals found dead when setting/supplementing
    filter(!session == "0") %>%
    filter(!is.na(session))

  ########## THE ABOVE CODE removes all animals that have missing data or were found dead NOT during a trapping occasion ##################
  # voledata %>% filter(is.na(trap))
  # #219825 - euthanized 9.1.21 during field course (caught off the grid)
  # #219917 - euthanized 9.1.21 during field course (caught off the grid)
  # voledata %>% filter(is.na(tag)) #these are fine to remove, nothing important here
  # voledata %>% filter(session=="0")
  # #these are animals found dead when setting - but I don't know when they went into the trap (assumed right after we left)
  # #they aren't technically session 4 captures, though so I don't know what I'd do with them if I left them in
  # voledata %>% filter(is.na(session)) #these animals were found DT when supplementing

  ########## RIGHT NOW, all of the DP, DT, and S animals (except those that were part of ^ above) are still in the dataset ################
  # #remove animals euthanized for terminal sampling
  # filter(!fate == "S")
  # #remove animals DP or DT
  # filter(!fate == "DP") %>%
  # filter(!fate == "DT")

  #check for spelling errors, extra groups, weird data
  # df$column[df$column == "old_value"] <- new_value #find and replace
  # unique(voledata$site)
  # unique(voledata$year)
  # unique(voledata$occasion)
  # unique(voledata$session)
  # unique(voledata$food_trt)
  # unique(voledata$helm_trt)
  # unique(voledata$sex)
  # unique(voledata$fate)
  # unique(voledata$ow)
  # unique(voledata$per)
  # unique(voledata$nip)
  # unique(voledata$preg)
  # unique(voledata$test)
  # range(voledata$head, na.rm = TRUE)
  # range(voledata$mass, na.rm = TRUE)
  # range(voledata$ticks, na.rm = TRUE)


  ################################ create a time column to replace session (for CMRnet) ##################################
  voledata <-
    voledata %>%
    mutate(time = case_when(
      session == "1" ~ "06:00:00",
      session == "2" ~ "18:00:00",
      session == "3" ~ "06:00:00",
      session == "4" ~ "18:00:00",)) %>%
    #so we didn't actually check traps 12hrs apart (more like 6am, 4pm) but this is cleaner for the distance moved analysis
    relocate(time, .after = session)

  #turn the time into a lubridate time
  voledata$time <- hms(voledata$time)
  #combine date and time columns into a lubridate time
  voledata$date_time <- ymd_hms(paste(voledata$date, voledata$time))
  #move date and time around
  voledata <-
    voledata %>%
    relocate(date, .after = year) %>%
    relocate(session, .after = occasion) %>%
    relocate(date_time, .after = date) %>%
    dplyr::select(-time)
  ######################################################### END ############################################################

  ################################ convert trap number to a grid coordinate ####################################
  voledata<-
    voledata %>%
    #separate trap ID (letter/number) into column of the letter (x) and number (y)
    separate(trap, into = c("x", "y"), sep = "(?<=[A-Z])(?=[0-9])", remove=FALSE) %>%
    #recode each letter as a number
    #this makes a new column with the letter turned into its corresponding number
    mutate(new_x = case_when(
      x == "A" ~ 1,
      x == "B" ~ 2,
      x == "C" ~ 3,
      x == "D" ~ 4,
      x == "E" ~ 5,
      x == "F" ~ 6,
      x == "G" ~ 7,
      x == "H" ~ 8,
      x == "I" ~ 9,
      x == "J" ~ 10,
      x == "K" ~ 11)) %>%
    #move that new column right after the original letter column
    relocate(new_x, .after = x) %>%
    #remove that original letter column
    dplyr::select(-x) %>%
    #rename the new number x column
    rename(x = new_x) %>%
    #make sure y and x are both numeric, not characters
    mutate(x = as.numeric(x), y = as.numeric(y))


  voledata<-
    voledata %>%
    #adjust the trap number (y) so it corresponds to a grid number
    #if remainder of x/2 is 0 (if x is even), multiply y by 2 - if else multiply by 2 then subtract 1
    mutate(y_new = ifelse((x %% 2) == 0, ((voledata$y)*2), (((voledata$y)*2)-1))) %>%
    relocate(y_new, .after = y) %>%
    #remove the old y column and rename the new one
    dplyr::select(-y) %>%
    rename(y = y_new)
  ######################################################### END ##################################################################

  #remove un-needed columns from datatable - make a smaller version for easier analysis
  trap <- voledata %>%
    dplyr::select(-id_as_number, -date, -ow, -ticks, -ear_ed, -saliva_sr, -smear_bs,
                  -bloodspin_bc, -bloodrna_br, -fecalegg_fe, -fecalrna_fr, -deworm, -new, -handler, -notes)

  ##### SAMPLE_ID and YEAR are back in the data (May 2023)


  #############################  LOAD AND CLEAN WEEK RECAP DATA   ###############################################

  #load the data
  wr_data <- read.csv(here(WRdata))
  #march 03 2024 version - has cleaned/updated PIT tags, vole sexes (2021 and 2022 season)
  #cleaning was done on full 2021,2022,2023 capture data

  #several columns (sex, per, nip, preg, test, fate, handler) have "" or "not noted" --> change these to NA
  wr_data$sex[wr_data$sex == "not noted"] <- NA
  wr_data$sex[wr_data$sex == ""] <- NA
  wr_data$per[wr_data$per == "not noted"] <- NA
  wr_data$per[wr_data$per == ""] <- NA
  wr_data$nip[wr_data$nip == "not noted"] <- NA
  wr_data$nip[wr_data$nip == ""] <- NA
  wr_data$preg[wr_data$preg == "not noted"] <- NA
  wr_data$preg[wr_data$preg == ""] <- NA
  wr_data$test[wr_data$test == "not noted"] <- NA
  wr_data$test[wr_data$test == ""] <- NA
  wr_data$fate[wr_data$fate == ""] <- NA
  wr_data$handler[wr_data$handler == ""] <- NA

  #clean names
  wr_data <- wr_data %>%
    remove_empty(which="rows") %>%
    clean_names %>%
    #rename column
    rename(session = sess) %>%
    #mutate columns to their proper formats
    mutate(year = as.numeric(year),
           occasion = as.numeric(occasion), #as.factor(occasion),
           site = as.character(site),
           food_trt = as.factor(food_trt),
           helm_trt = as.factor(helm_trt),
           date = as_date(date, format= "%m/%d/%Y"),
           session = as.numeric(session),
           trap = as.character(trap),
           tag = as.character(tag),
           sex = as.factor(sex),
           per = as.factor(per),
           nip = as.factor(nip),
           preg = as.factor(preg),
           test = as.factor(test),
           fate = as.factor(fate),
           handler = as.factor(handler),
           notes = as.character(notes)) %>%
    #filter for desired year
    filter(year==yr) %>%
    #remove uusi data
    filter(site!="uusi") %>% #only relevant for 2022
    filter(!is.na(trap)) %>% #remove NA trap
    filter(!is.na(tag)) #remove NA tag

  # #check for missing data
  # wr_data %>% filter(is.na(session)) #no missing sessions
  # wr_data %>% filter(is.na(date)) #no missing dates

  #check for spelling errors, extra groups, weird data
  # unique(wr_data$year)
  # unique(wr_data$occasion)
  # unique(wr_data$site)
  # unique(wr_data$food_trt)
  # unique(wr_data$helm_trt)
  # unique(wr_data$session)
  # unique(wr_data$sex)
  # unique(wr_data$per)
  # unique(wr_data$nip)
  # unique(wr_data$preg)
  # unique(wr_data$test)
  # unique(wr_data$fate)

  ############ DEAD ANIMALS ARE STILL IN THE DATASET ###################
  # #remove animals euthanized for terminal sampling
  # filter(!fate == "S") %>%
  # #remove animals DP or DT
  # filter(!fate == "DP") %>%
  # filter(!fate == "DT")

  ################################ create a time column to replace session (for CMRnet) ##################################
  wr_data <-
    wr_data %>%
    mutate(time = case_when(
      session == "1" ~ "06:00:00",
      session == "2" ~ "18:00:00",
      session == "3" ~ "06:00:00",
      session == "4" ~ "18:00:00",)) %>%
    relocate(time, .after = session)

  #turn the time into a lubridate time
  wr_data$time <- hms(wr_data$time)
  #combine date and time columns into a lubridate time
  wr_data$date_time <- ymd_hms(paste(wr_data$date, wr_data$time))
  #move date and time around
  wr_data <-
    wr_data %>%
    relocate(date, .after = year) %>%
    relocate(session, .after = occasion) %>%
    relocate(date_time, .after = date) %>%
    dplyr::select(-time)
  ######################################################### END ############################################################

  ################################ convert trap number to a grid coordinate ####################################
  wr_data<-
    wr_data %>%
    separate(trap, into = c("x", "y"), sep = "(?<=[A-Z])(?=[0-9])", remove=FALSE) %>%
    #recode each letter as a number
    #this makes a new column with the letter turned into its corresponding number
    mutate(new_x = case_when(
      x == "A" ~ 1,
      x == "B" ~ 2,
      x == "C" ~ 3,
      x == "D" ~ 4,
      x == "E" ~ 5,
      x == "F" ~ 6,
      x == "G" ~ 7,
      x == "H" ~ 8,
      x == "I" ~ 9,
      x == "J" ~ 10,
      x == "K" ~ 11)) %>%
    #move that new column right after the original letter column
    relocate(new_x, .after = x) %>%
    #remove that original letter column
    dplyr::select(-x) %>%
    #rename the new number x column
    rename(x = new_x) %>%
    #make sure y and x are both numeric, not characters
    mutate(x = as.numeric(x), y = as.numeric(y))


  wr_data<-
    wr_data %>%
    #adjust the trap number (y) so it corresponds to a grid number
    #if remainder of x/2 is 0 (if x is even), multiply y by 2 - if else multiply by 2 then subtract 1
    #not sure why, but this code doesn't work if you string it together with the above code... the values for y are all wrong
    mutate(y_new = ifelse((x %% 2) == 0, ((wr_data$y)*2), (((wr_data$y)*2)-1))) %>%
    relocate(y_new, .after = y) %>%
    #remove the old y column and rename the new one
    dplyr::select(-y) %>%
    rename(y = y_new)
  ######################################################### END ############################################################

  #remove un-needed columns from data table - make a smaller version for easier analysis
  recap <- wr_data %>%
    dplyr::select(-date, -handler, -notes)


  #############################  COMBINE AND CLEAN FULLTRAP DATA   ###############################################

  #compare columns and their classes - confirm everything matches
  # compare_df_cols(trap, recap)

  #combining trap and recap dataframes
  fulltrap <- bind_rows(trap, recap)

  #change food_trt to "fed" "unfed"
  fulltrap$food_trt <- fct_recode(fulltrap$food_trt, "fed"="supplement", "unfed"="control")
  fulltrap$food_trt <-fct_relevel(fulltrap$food_trt, c("unfed", "fed"))
  #create trt column
  fulltrap <- fulltrap %>%
    unite(trt, food_trt, helm_trt, sep = "_", remove = FALSE)
  #and then relevel it
  fulltrap$trt <-fct_relevel(fulltrap$trt, c("unfed_control", "unfed_deworm", "fed_control", "fed_deworm"))

  #add month column to replace occasion - keep occasion, session since later code might need numbers
  fulltrap <- fulltrap %>% mutate(month = case_when(occasion == "1" ~ "may",
                                                    occasion == "2" ~ "june",
                                                    occasion == "3" ~ "july",
                                                    occasion == "4" ~ "aug",
                                                    occasion == "5" ~ "sept",
                                                    occasion == "6" ~ "oct")) %>%
    mutate(month = factor(month, levels=c("may", "june", "july", "aug", "sept", "oct"))) %>%
    relocate(month, .after=date_time)

  #FIRSTCAP column: give a 1 if the capture is the first occurrence of that tag, else 0
  fulltrap <- fulltrap %>%
    unite(occ.sess, occasion, session, sep = ".", remove = FALSE) %>% #make a new occ.sess column so I can do things in order
    group_by(tag) %>%
    mutate(firstcap = ifelse(occ.sess == min(occ.sess), 1, 0)) %>%
    mutate(firstcap = factor(firstcap)) %>%
    ungroup()
  #this will record all the WRs as 0 as well
  ################ END FIRSTCAP/NEW ###########################

  # #in case igraph is already running, it masks "%--%" and the lubridate code won't run
  # # library(needs)
  # # prioritize(lubridate)
  # #create a measure of TIME KNOWN ALIVE (also date of first cap, last capture)
  # fulltrap <- fulltrap %>%
  #   group_by(tag) %>%
  #   arrange(date_time) %>%
  #   mutate(first_seen = min(date_time),
  #          last_seen = max(date_time)) %>%
  #   mutate(days_known_alive = round( as.duration(first_seen %--% last_seen) / ddays(1) , digits=2) ) %>%
  #   ungroup()

  ## KEEPING traps_per_life and caps_per_life for exploratoriness measurement later
  
  #create caps_per_life column - number of captures of that individual
  fulltrap <- fulltrap %>%
    group_by(tag) %>%
    mutate(caps_per_life = length(tag)) %>%
    ungroup()
  #count of number of unique traps per animal in lifetime
  fulltrap <- fulltrap %>%
    group_by(tag) %>%
    mutate(traps_per_life = length(unique(trap))) %>%
    relocate(traps_per_life, .after=caps_per_life) %>%
    ungroup()
  # #create a RECAPPED column (binary - were you recapped or not)
  # fulltrap <- fulltrap %>%
  #   group_by(tag) %>%
  #   mutate(recapped = case_when(
  #     caps_per_life == "1" ~ 0,
  #     caps_per_life > "1" ~ 1)) %>%
  #   mutate(recapped = factor(recapped)) %>%
  #   relocate(recapped, .before=caps_per_life) %>%
  #   ungroup()

  #season column
  fulltrap <- fulltrap %>%
    mutate(season = ifelse(month=="sept" | month=="oct", "fall", "summer")) %>%
    mutate(season = factor(season, levels=c("summer", "fall"))) %>%
    relocate(season, .after=month)
  
  # #caps_per_season
  # fulltrap <- fulltrap %>%
  #   group_by(tag, season) %>%
  #   mutate(caps_per_season = length(tag)) %>%
  #   mutate(res = ifelse(caps_per_season >= 5, "resident", "nonresident")) %>% #add resident status
  #   mutate(res = factor(res)) %>%
  #   relocate(res, .after = traps_per_life) %>%
  #   ungroup()
  # ################ NOTE: RESIDENT = 5 caps per SEASON -- NOT 5 caps per LIFE #############
  # #traps_per_season
  # fulltrap <- fulltrap %>%
  #   group_by(tag, season) %>%
  #   mutate(traps_per_season = length(unique(trap))) %>%
  #   ungroup()
  # 
  # #caps_per_occ column - number of captures per occasion
  # fulltrap <- fulltrap %>%
  #   group_by(occasion, tag) %>%
  #   mutate(caps_per_occ = length(tag)) %>%
  #   ungroup()
  # #count of number of unique traps per animal per occasion
  # ##n_distinct() is a dplyr wrapper for length(unique())
  # fulltrap <- fulltrap %>%
  #   group_by(tag, occasion) %>%
  #   mutate(traps_per_occ  = length(unique(trap))) %>%
  #   ungroup()

  #make a column (current_breeder) of "breeder"/"nonbreeder" based on 1 for per/nip/preg or test IN THAT CAPTURE
  #NAs persist
  #then make a new column (ever_breeder) based on current_breeder column
  #where: if ever the animal was breeder during any capture, it gets "breeder"
  #use NA.rm in this to make sure NAs are ignored and get replaced with the right word
  #HOWEVER - if animal NEVER had breeding condition recorded in any capture, ever_breeder will be NA
  fulltrap <- fulltrap %>%
    #current_breeder in this occasion
    mutate(current_breeder = case_when(per=="1" | nip=="1" | preg == "1" ~ "breeder",
                                       test == "1" ~ "breeder",
                                       per=="0" & nip=="0" & preg == "0" ~ "nonbreeder",
                                       test == "0" ~ "nonbreeder")) %>%
    relocate(current_breeder, .after=sex) %>%
    #season_breeder if ever in the season (June-Aug, Sept-Oct)
    group_by(season, tag) %>%
    mutate(season_breeder_init = ifelse(all(is.na(current_breeder)), NA,
                                 ifelse(any(current_breeder == "breeder", na.rm = TRUE), "breeder", "nonbreeder"))) %>%
    relocate(season_breeder_init, .after = current_breeder) %>%
    #ever_breeder if ever in June-Oct
    group_by(tag) %>%
    mutate(ever_breeder = ifelse(all(is.na(current_breeder)), NA,
                                 ifelse(any(current_breeder == "breeder", na.rm = TRUE), "breeder", "nonbreeder"))) %>%
    relocate(ever_breeder, .after = season_breeder_init) %>%
    #adjust fall non-breeders to make sure they're Fall breeders if they were in the Summer
    mutate(season_breeder = ifelse(season=="fall" & season_breeder_init=="nonbreeder" & ever_breeder=="breeder",
                                    "breeder", season_breeder_init)) %>%
    mutate(season_breeder = factor(season_breeder, levels=c("breeder", "nonbreeder"))) %>%
    relocate(season_breeder, .after=season_breeder_init) %>%
    select(!c(current_breeder, ever_breeder, season_breeder_init))

  #reorganize this unholy mess
  fulltrap <- fulltrap %>%
    relocate(c(tag,samp_id,firstcap,trap,x,y), .after=helm_trt)

  #aight, I only want NAs if we DO NOT have that information at all in any way, so...
  #fill within an occasion for data not collected on WRs but collected during processing
  #ie: sex, breeding condition (per,nip,preg,test), current breeder, head, mass
  fulltrap <- fulltrap %>%
    group_by(tag, occasion) %>%
    arrange(session, .by_group = TRUE) %>%
    fill(c(samp_id, sex, per, nip, preg, test,
           season_breeder, mass, head), .direction="downup") %>% #fill missing data within an occasion
    ungroup()

  #### 2021 VOLES WITH SEX=NA (DUE TO FUCKUPS and other things and it's okay, I'm not mad about it)
  # 226128 # 226211 # 226769 # 219682
  ## some 2022 voles (at least 2) have sex=NA because we couldn't determine the correct sex


  ####################################################################################################

  # ## FIRST! for data purposes, save fulltrap_ALL_03.04.24.rds to get counts of all captured voles
  #   ## including May captures, captures with sex=NA or repro=NA
  # 
  # saveRDS(fulltrap, file=here("fulltrap22_ALL_03.04.24.rds"))
  
  
  ##### FOR PUBLICATION: skim down the data file to as little data as necessary for analyses

  #remove data that is not used for vole-hanta analysis
  fulltrap <- fulltrap %>%
    select(!c(per, nip, preg, test, head, mass, fate)) %>% #remove data columns not needed for analysis
    filter(month != "may") %>% #drop may data since not all sites had captures (may data not used in analysis)
    mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #adjust levels
    drop_na(sex) %>% #remove animals with sex=NA (since we can't assign them a HR)
    drop_na(season_breeder) #remove animals without season_breeder data (since we can't assign them a HR)

  
  #save fulltrap to an rdata file so I can pull it for other scripts
  saveRDS(fulltrap, file = here(fulltrap_output))

  # Restore fulltrap from the rdata file
  # fulltrap <- readRDS(file = OUTPUT FILE NAME GOES HERE)

  #################### END ########################


# } #not currently running as a function
