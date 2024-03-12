# ARCHIVED on 3.12.24 - no longer needed in analysis

### code to check things in '03_NEW_hanta_cleaning'
#looking at the voles that seroconvert 1-0 and overwintered animals
#JM removed 6.9.23


################################################################################################################
#########################################   doin stuff    ######################################################
################################################################################################################
##### lil check to see the voles that seroconvert, and then find the ones that convert from 1 to 0... ##########
################################################################################################################

# net_mets_puuv <- net_mets_puuv %>% group_by(tag) %>% arrange(year, month, .by_group = TRUE)
# net_mets_puuv$puuv_ifa <- as.character(net_mets_puuv$puuv_ifa)
# check <- net_mets_puuv %>% group_by(tag) %>%
#   summarise(status = unique(puuv_ifa)) 
# 
# #save the tag IDs of animals that seroconvert
# duplist <- as.vector(check$tag[duplicated(check$tag)])
# length(duplist) #103 animals seroconvert (2021 and 2022)
# 
# check2 <- check %>% filter(tag %in% duplist) #just the animals that are 0 and 1

# check_longitude <- net_mets_puuv %>% group_by(tag) %>% arrange(year, month, .by_group = TRUE) %>%
#   dplyr::select(year, month, tag, puuv_ifa) %>%
#   summarise(status_time = toString(puuv_ifa)) %>%
#   filter(str_detect(status_time, "0,\\s1|1,\\s0")) #filter for all animals (n=103) that seroconvert (regardless of direction)

# #voles that seroconvert PUUV+ to PUUV-
# puuv_pos_neg <- net_mets_puuv %>% group_by(tag) %>% arrange(year, month, .by_group = TRUE) %>%
#   dplyr::select(year, month, tag, puuv_ifa) %>%
#   summarise(status_time = toString(puuv_ifa)) %>%
#   filter(str_detect(status_time, "1,\\s0")) #filter for animals that go from pos to neg
# 
# problemchildren <- puuv_pos_neg$tag
# write.csv(puuv_pos_neg, here("puuv_pos_neg.csv"))
# 
#grab just the sample ID and notes
# hanta_notes <- read.csv(here("puuv_ifa.csv")) %>%
#   select(id, notes) %>% 
#   rename(samp_id = id)

# #voles that seroconvert PUUV- to PUUV+
# puuv_neg_pos <- net_mets_puuv %>% group_by(tag) %>% arrange(year, month, .by_group = TRUE) %>%
#   dplyr::select(year, month, tag, puuv_ifa) %>%
#   summarise(status_time = toString(puuv_ifa)) %>%
#   filter(str_detect(status_time, "0,\\s1"))

# puuvcheckdata <- hantadata %>% left_join(tag_id, by= c("id" = "samp_id")) %>% filter(tag %in% problemchildren) %>%
#   group_by(tag) %>% arrange(id, .by_group = TRUE) %>% relocate(tag, .after=id)
# 
# write.csv(puuvcheckdata, here("IFAdata_posneg.csv"))

# fulltrap21 <- read_rds(here("fulltrap21_02.16.23.rds")) %>%
#   mutate(year = 2021)
# 
# fulltrap22 <- fulltrap22 %>% select(!new) %>%
#   mutate(year = 2022)
# 
# FULLfulltrap <- rbind(fulltrap21, fulltrap22) %>%
#   select(c(year, site, occ.sess, samp_id, tag, sex, firstcap, mass))
#          
# check <- FULLfulltrap %>% group_by(tag) %>% 
#   summarise(n = n_distinct(year)) %>%
#   filter(n > 1) 
# bothyrs <- check$tag
# overwinter <- FULLfulltrap %>% filter(tag %in% bothyrs) %>% arrange(tag, year, occ.sess)
# write.csv(overwinter, here("overwinters.csv"))

# fullproblemchildren <- inner_join(tag_id, hantadata, by= "samp_id") %>% 
#   inner_join(hanta_notes, by="samp_id") %>%
#   left_join(FULLfulltrap, by=c("samp_id", "tag")) %>%
#   filter(tag %in% problemchildren) %>%
#   select(year, site, occ.sess, tag, samp_id, puuv_ifa, mass, firstcap, sex, notes) %>%
#   arrange(tag, occ.sess)
# 
# write.csv(fullproblemchildren, here("confirm_IFA.csv"))


# #overwintered 2021-22 animals
# OW <- net_mets_puuv %>% group_by(tag) %>% arrange(tag, year, month) %>%
#   mutate(year=as.factor(year),
#          ows = n_distinct(year)) %>% 
#   filter(ows == 2)

################################################### END doin' stuff ####################################################