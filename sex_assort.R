## FILE ARCHIVED - NO LONGER USED IN ANALYSIS (but still cool to see)

# load packages
library(here)
library(tidyverse)
library(igraph)
library(lubridate)
library(janitor)

#clear environment
rm(list = ls())

ft21 <- readRDS(here("fulltrap21_03.04.24.rds"))
data = ft21
networks_file = "overlap_networks21.rds"
sex_assort_file = "sex_assort21.rds"

# ft22 <- readRDS(here("fulltrap22_03.04.24.rds"))
# data = ft22
# networks_file = "overlap_networks22.rds"
# sex_assort_file = "sex_assort22.rds"


sex_assort <- function(data, networks_file, sex_assort_file){

  ##---------------- LOAD THE DATA ----------------------

  #running using assortnet (D.Farine) because it can work with weighted networks
  library(assortnet)

  #load, clean the fulltrap dataset (make sure it's the most recent version)
  fulltrap <- data %>%
    filter(month != "may") %>% #drop may data since not all sites had captures
    mutate(month = factor(month, levels=c("june", "july", "aug", "sept", "oct"))) %>% #adjust levels, remove may
    drop_na(sex) %>% #remove animals with sex=NA (since we can't assign then a HR)
    unite("sts", season, trt, sex, remove=FALSE) #add sts column to match params_summary

  # load the network data
  overlap_network_list <- readRDS(here(networks_file))

  #create tag_sex df for assortativity by sex
  tag_sex <- fulltrap %>% group_by(tag) %>% slice(1) %>%
    select(tag, sex) %>%
    arrange(tag)

  #create tag_sb df for assortativity by sex-breeder
  tag_sb <- fulltrap %>% unite(sb, sex, season_breeder) %>%
    group_by(tag) %>% slice(1) %>%
    select(tag, sb) %>%
    arrange(tag)


  ##--------------------------------------------

  # Create list to store results
  sex_assort_list <- list()

  # Calculate network metrics
  for(i in 1:length(overlap_network_list)){

    #length(overlap_network_list) == for each of the 12 sites

    #for each site
    print(names(overlap_network_list[i]))
    site.id <- names(overlap_network_list[i])

    site <- list()

    for(j in 1:length(overlap_network_list[[i]])){

      #(j in 1:length(overlap_network_list[[i]])) == for the 4-5 trapping occasions
      #would be (j in 1:5) except that Kiirastuli 2022 only has 4 :P *thhhbbbbt*

      #for each month
      print(names(overlap_network_list[[i]][j]))
      month <- names(overlap_network_list[[i]][j])

      #dataframe to hold results per month
      site[[j]] <- data.frame(month)

      adjmat <- overlap_network_list[[i]][[j]]
      diag(adjmat) <- 0 #matrix diagonal is NA - assortnet needs it to be 0

      #create WEIGHTED NETWORK from adjacency matrix
      inet <- graph_from_adjacency_matrix(adjmat, mode="undirected", weighted = TRUE, diag = FALSE)

      ids <- get.vertex.attribute(inet, "name") #tag ids for all the animals on the grid

      ### FOR ASSORTATIVITY BY SEX
      #filter tag_sex for only ids caught this site/occ
      ids_sex <- tag_sex %>% filter(tag %in% ids)
      #calculate assortativity
      out <- assortment.discrete(adjmat, as.vector(ids_sex$sex), weighted=TRUE, SE=FALSE, M=1, na.rm=FALSE)
        #out is a list, $r has the assort coef across all individuals, $mixing_matrix has ppn of edges by sex

      check <- n_distinct(ids_sex$sex) #how many sexes are represented? (to catch site/month when only M or F are present)

      mat <- out$mixing_matrix #save just the mixing matrix
      #pull the percent of M/F overlap, F/F, and M/M
        #check is to make sure NA is input if there was only one sex of breeders in that month
      site[[j]]$fm <- ifelse(check==1, NA, mat["M","F"]*2) #double the fm overlaps since network is undirected
      site[[j]]$ff <- ifelse(check==1, NA, mat["F","F"])
      site[[j]]$mm <- ifelse(check==1, NA, mat["M","M"])

      ### FOR ASSORTATIVITY BY SEX-BREEDER
      #filter tag_sex for only ids caught this site/occ
      ids_sb <- tag_sb %>% filter(tag %in% ids)
      #calculate assortativity
      out_sb <- assortment.discrete(adjmat, as.vector(ids_sb$sb), weighted=TRUE, SE=FALSE, M=1, na.rm=FALSE)
      #out is a list, $r has the assort coef across all individuals, $mixing_matrix has ppn of edges by sex-breeder

      #if sex-breeder combo isn't present, value saved==0, else it's a number
      mb <- nrow(ids_sb %>% filter(sb=="M_breeder"))
      mnb <- nrow(ids_sb %>% filter(sb=="M_nonbreeder"))
      fb <- nrow(ids_sb %>% filter(sb=="F_breeder"))
      fnb <- nrow(ids_sb %>% filter(sb=="F_nonbreeder"))

      mat_sb <- out_sb$mixing_matrix #save just the mixing matrix
      #pull the percent of all the overlaps (hold onto your pants, this is going to get cra')
      site[[j]]$fbmb <- ifelse(fb=="0" | mb=="0", NA, mat_sb["M_breeder","F_breeder"]*2) #double since network is undirected
      site[[j]]$fbfb <- ifelse(fb=="0", NA, mat_sb["F_breeder","F_breeder"])
      site[[j]]$mbmb <- ifelse(mb=="0", NA, mat_sb["M_breeder","M_breeder"])
      site[[j]]$fbmnb <- ifelse(fb=="0" | mnb=="0", NA, mat_sb["M_nonbreeder","F_breeder"]*2)
      site[[j]]$fbfnb <- ifelse(fb=="0" | fnb=="0", NA, mat_sb["F_nonbreeder","F_breeder"]*2)
      site[[j]]$mbmnb <- ifelse(mnb=="0" | mb=="0", NA, mat_sb["M_nonbreeder","M_breeder"]*2)
      site[[j]]$mbfnb <- ifelse(fnb=="0" | mb=="0", NA, mat_sb["M_breeder","F_nonbreeder"]*2)
      site[[j]]$fnbfnb <- ifelse(fnb=="0", NA, mat_sb["F_nonbreeder","F_nonbreeder"])
      site[[j]]$mnbmnb <- ifelse(mnb=="0", NA, mat_sb["M_nonbreeder","M_nonbreeder"])
      site[[j]]$fnbmnb <- ifelse(mnb=="0" | fnb=="0", NA, mat_sb["M_nonbreeder","F_nonbreeder"]*2)

    }

    #write the list 'site' as a 1st-order item in wt_net_mets_list
    sex_assort_list[[i]] <- site

  }

  #name the 12 1st order elements of wt_nets_list as the sites
  names(sex_assort_list) <- names(overlap_network_list)

  #rename the sublist items (months) for each site
  #accounting for the fact that most sites have 5 months of data but some have 4
  for(i in 1:length(sex_assort_list)){
    ifelse( length(sex_assort_list[[i]]) == 5,
            names(sex_assort_list[[i]]) <- c("june", "july", "aug", "sept", "oct"),
            names(sex_assort_list[[i]]) <- c("july", "aug", "sept", "oct") )
  }


  ############## condense sex_assort_list to make it easier to use for analysis #######################

  #collate results
  #this collapses the 2nd order elements (network metrics for a single month) down to the 1st order element (site)
  #so now wt_net_mets_list_summary is a list of 12 dfs, each is all the network metrics for a site across all the months

  #make a list to store things
  sex_assort_summary <- list()

  #loop across all sites and collapse the dfs per occasion into one df for the site
  for(i in 1:length(sex_assort_list)){

    #for all 12 sites
    summary <- do.call("rbind", sex_assort_list[[i]])
    sex_assort_summary[[i]] <- summary
  }

  #remove the row names of the dfs (its the month, but we already have that info)
  for(i in 1:length(sex_assort_summary)){
    row.names(sex_assort_summary[[i]]) <- NULL
  }

  #name the 12 1st order elements as their sites
  names(sex_assort_summary) <- names(overlap_network_list)

  ## make sex_assort_summary into freiggein huge df
  sex_assort_summary <- do.call(rbind.data.frame, sex_assort_summary)

  #clean up the df
  sex_assort_summary <- sex_assort_summary %>%
    rownames_to_column("name") %>% #row names are the sites, make that a column
    separate(name, c("site", NA)) %>% #separate the site part from the index and get rid of the index
    mutate(site = as.factor(site)) #make site a factor

  #save it
  saveRDS(sex_assort_summary, here(sex_assort_file))


}



# #visualize one year
# sex_assort_summary %>%
#   left_join(read.csv(here("grid_trts.csv")), by="site") %>% unite(trt, food_trt, helm_trt) %>%
#   select(!c(ff, fm, mm)) %>%
#   pivot_longer(-c(site, trt, month), names_to = "group", values_to = "pct") %>%
#   mutate(month = factor(month, levels=c('june', "july", "aug", "sept", "oct"))) %>%
#   mutate(group = factor(group, levels=c("fbfb", "fbfnb", "fnbfnb", "fbmnb", "fbmb",
#                                         "fnbmnb", "mbfnb", "mbmb", "mbmnb", "mnbmnb"))) %>%
#   ggplot(aes(fill=group, y=pct, x=month)) +
#   geom_bar(position="fill", stat="identity") +
#   scale_fill_manual(values=c("#a00000", "#ff4b4b", "#ff7c7c",
#                              "#ffbf00",
#                              "#800080", "#efbbff",
#                              "#2a940a",
#                              "#0200b9", "#007dff", "#91e7ff")) +
#   facet_wrap(~trt, nrow=1)



###### both years ########

#combo 2021 and 2022
sex_assort21 <- readRDS(here("sex_assort21.rds")) %>%
  mutate(year="2021")
sex_assort22 <- readRDS(here("sex_assort22.rds")) %>%
  mutate(year="2022")

sex_assort21.22 <- rbind(sex_assort21, sex_assort22)



#visualize
png(filename = "sex-breed_assort.png", width=12 , height=6, units="in", res=600)

sex_assort21.22 %>%
  left_join(read.csv(here("grid_trts.csv")), by="site") %>% unite(trt, food_trt, helm_trt) %>%
  select(!c(ff, fm, mm)) %>%
  pivot_longer(-c(year, site, trt, month), names_to = "group", values_to = "pct") %>%
  mutate(trt = factor(trt, levels=c("unfed_control", "unfed_deworm", "fed_control", "fed_deworm"))) %>%
  mutate(month = factor(month, levels=c('june', "july", "aug", "sept", "oct"))) %>%
  mutate(group = factor(group, levels=c("fbfb", "fbfnb", "fnbfnb", "fbmnb", "fbmb",
                                        "fnbmnb", "mbfnb", "mbmb", "mbmnb", "mnbmnb"))) %>%
  ggplot(aes(fill=group, y=pct, x=month)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=c("#a00000", "#ff4b4b", "#ff7c7c",
                             "#ffbf00",
                             "#800080", "#efbbff",
                             "#2a940a",
                             "#0200b9", "#007dff", "#91e7ff")) +
  facet_grid(year~trt) +
  labs(main="Assortativity by Sex, Breeding Status", y="Percent of Overlaps", x="Month")

dev.off()


#2021 - just sex assort (ignoring breeder)
sex_assort21 %>%
  left_join(read.csv(here("grid_trts.csv")), by="site") %>% unite(trt, food_trt, helm_trt) %>%
  select(c(year, site, trt, month, ff, fm, mm)) %>%
  pivot_longer(-c(year, site, trt, month), names_to = "group", values_to = "pct") %>%
  mutate(trt = factor(trt, levels=c("unfed_control", "unfed_deworm", "fed_control", "fed_deworm"))) %>%
  mutate(month = factor(month, levels=c('june', "july", "aug", "sept", "oct"))) %>%
  mutate(group = factor(group, levels=c("ff", "fm", "mm"))) %>%
  ggplot(aes(fill=group, y=pct, x=month)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values=c("#a00000",
                             "#efbbff",
                             "#007dff")) +
  facet_wrap(~trt, ncol=4) +
  labs(main="Assortativity by Sex, Breeding Status", y="Percent of Overlaps", x="Month")

