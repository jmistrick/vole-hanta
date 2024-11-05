#Construct Spatial Overlap Networks
#Calculate Network Metrics


###-------------------------------------------------

#load packages
library(here)
library(tidyverse)
library(tidygraph)


#clear environment
rm(list = ls())


######------------------------------------------------------------------------


## Function for creating spatial overlap networks for individual voles for all the sites, months in a given year

#previously here: function: "create_overlap_networks"
#the output of this function was the file "overlap_networks"

# overlapnets21 <- readRDS(here("overlap_networks21.rds"))

#overlapnets is a nested list of 12 1e items - 1 per site
  #and 5 2e items per site - one per month
    #under each 2e item is a square matrix, of all the tags (numeric L-R, T-B - same order) and the weighted overlap measure between
    #the two voles. NAs are populated on the diagonal (and presumably could be for any overlaps that don't occur)
#this file (ie adjacency matrix per site/month) is what's fed into "calculate_network_metrics" to build the networks


##this file is essentially the file created in "07_percent_overlap.R"
  #creating the same weighted (NOW DIRECTED) adjacency matrices, but using circular HRs at the 95% capture probability from the
  #Wanelik/Farine space use kernels

#load the 2021 overlap list
pctoverlap_list21 <- readRDS(here("pct_overlap_list21.rds"))

# #load the 2022 overlap list
pctoverlap_list22 <- readRDS(here("pct_overlap_list22.rds"))


#load fulltrap data (trapping metadata for each capture) 
ft21 <- readRDS(here("fulltrap21_03.04.24.rds"))

ft22 <- readRDS(here("fulltrap22_03.04.24.rds"))



#####-------------------------------------------------------------------------------------

#### THIS CODE calculates the network metrics (mostly just weighted degree and its relatives)
# for all the voles at all the times
# output is "wt_net_mets_summary" file which can be used for downstream analysis

## Input: data = the FULL fulltrap df for that year
##        networks_file = networks file generated in create_overlap_networks
##        netmets_file = file to be generated, df of network metrics
## Output: netmets_file = dataframe of network metrics for every vole in every occasion it was captured


library(igraph)

calculate_network_metrics <- function(data, networks_file, netmets_file){

  ##---------------- LOAD THE DATA ----------------------
  
  #load, clean the fulltrap dataset - add columns for sts and sb
  fulltrap <- ft21 %>%
    unite("sts", season, trt, sex, remove=FALSE) %>% #add sts (season,treatment,sex) column to match params_summary
    unite("sb", sex, season_breeder, remove=FALSE) #add sb (sex, breeding status) column for degree by functional group

  # load the network data
  overlap_network_list <- pctoverlap_list22


  ##--------------------------------------------

  # Create list to store results
  wt_net_mets_list <- list()

  # Calculate network metrics
  for(i in 1:length(overlap_network_list)){
    #using 'length(overlap_network_list)' to loop through each of the 12 sites

    #for each site
    print(names(overlap_network_list[i]))
    site.id <- names(overlap_network_list[i])

    site <- list()

    for(j in 1:length(overlap_network_list[[i]])){
      #using 'length(overlap_network_list[[i]]))' to loop through all trapping occasions (ie months) with capture data per site

      #for each month
      print(names(overlap_network_list[[i]][j]))
      month.id <- names(overlap_network_list[[i]][j])

      adjmat <- overlap_network_list[[4]][[3]]

      #create WEIGHTED NETWORK from adjacency matrix
      inet <- graph_from_adjacency_matrix(adjmat, weighted = TRUE, mode="directed", diag = FALSE)

      #create UNWEIGHTED NETWORK (just for binary degree)
      adjmat_bin <-ifelse(adjmat>0,1,0)
      inet_bin <- graph_from_adjacency_matrix(adjmat_bin, weighted=NULL, mode="directed", diag=FALSE)

      ids <- get.vertex.attribute(inet, "name") #tag ids for all the animals on the grid
      month <- rep(names(overlap_network_list[[i]])[j],length(ids)) #capture month

      #dataframe to hold results per month
      site[[j]] <- data.frame(ids, month)

      #network metrics to calculate
      site[[j]]$wt.deg.in <- igraph::strength(inet, mode="in") #this is the sum of all IN degree weights for a node
      ## as I have it coded, IN DEGREE is the % overlap a focal vole has with each neighbor
        ## ie what percent of focal's HR is shared with neighbor x

      # site[[j]]$norm.wt.deg <- (igraph::strength(inet))/((igraph::gorder(inet))-1) #your strength/(total nodes-you)

      #binary degree (number of overlaps)
      site[[j]]$bin.in.deg <- igraph::degree(inet_bin, mode="in")

      #number of nodes (voles) in the network
      site[[j]]$n.node <- rep(igraph::gorder(inet), length(ids))
      
      
      #####-------------- calculating sex- and repro-specific degree measures -----------------------
      
      # ###### Calculate Male-degree, Female-degree #########
      # 
      # ###### Calculate Breeder(Reproductive)-degree, Nonbreeder(Nonreproductive)-degree #########
      # # yes, I know I use different language in the manuscript vs here in the code, I'm indecisive and spent
      # # too much time going back and forth. Just accept that anywhere I say "Breeder" here I mean "Reproductive"
      # # and likewise with "nonbreeder" and "nonreproductive"  ¯\_(ツ)_/¯ 
      # 
      # #make a (new) network from the adj matrix
      # g <- graph_from_adjacency_matrix(adjmat, mode="directed", weighted=TRUE, diag = FALSE)
      # #technically, our adj matrix is symmetrical so the network is undirected,
      # #BUT because I want count the number of ties going in/out from a given node,
      # #So I need to have all the 'out' ties listed in one column of the edgelist (so effectively a directed edgelist)
      # #So for these purposes, I built a directed network (but edge weights are symmetrical)
      # 
      # #subset metadata
      # metadata <- fulltrap %>% #starts with fulltrap, need to get down to a 'traits' version
      #   filter(site==site.id & month==month.id) %>%
      #   group_by(tag) %>% slice(1)  #one entry per vole per month
      # 
      # #set "sex" as a vertex attribute
      # g <- set.vertex.attribute(graph=g, name="sex", index=V(g), value=metadata$sex)
      # #set "breeder" as a vertex attribute
      # g <- set.vertex.attribute(graph=g, name="breeder", index=V(g), value=metadata$season_breeder)
      # #set "sb" as vertex attribute (sb=season, breeder)
      # g <- set.vertex.attribute(graph=g, name="sb", index=V(g), value=metadata$sb)
      # #save vectors
      # # focal_sex <- get.vertex.attribute(g, "sex") #female is 1, male is 2
      # # focal_breed <- get.vertex.attribute(g, "breeder") #breeder is 1, nonbreeder is 2
      # focal_id <- get.vertex.attribute(g, "name") #this isn't strictly necessary, just to make sure the focal vole is correct
      # 
      # ## code from Matt Michalska-Smith for male strength/female strength
      # sex_to <- get.vertex.attribute(g, "sex")[get.edgelist(g, names=FALSE)[,2]]
      # #since vertices already have meaningful names, call names=FALSE to return vertex indices instead
      # degree_to_M <- strength(g, mode="out", weights=((sex_to == "M")*get.edge.attribute(g, "weight"))) #this need to be mode="OUT" (the "to" individual)
      # degree_to_F <- strength(g, mode="out", weights=((sex_to == "F")*get.edge.attribute(g, "weight"))) #this need to be mode="OUT"
      # 
      # ## code from Matt M-S for breeder/nonbreeder strength
      # breed_to <- get.vertex.attribute(g, "breeder")[get.edgelist(g, names=FALSE)[,2]]
      # #since vertices already have meaningful names, call names=FALSE to return vertex indices instead
      # degree_to_b <- strength(g, mode="out", weights=((breed_to == "breeder")*get.edge.attribute(g, "weight"))) #this need to be mode="OUT" (the "to" individual)
      # degree_to_nb <- strength(g, mode="out", weights=((breed_to == "nonbreeder")*get.edge.attribute(g, "weight"))) #this need to be mode="OUT"
      # 
      # ## code from Matt M-S for sb strength across 4 pairings
      # sb_to <- get.vertex.attribute(g, "sb")[get.edgelist(g, names=FALSE)[,2]]
      # #since vertices already have meaningful names, call names=FALSE to return vertex indices instead
      # degree_to_mb <- strength(g, mode="out", weights=((sb_to == "M_breeder")*get.edge.attribute(g, "weight"))) #this need to be mode="OUT" (the "to" individual)
      # degree_to_mnb <- strength(g, mode="out", weights=((sb_to == "M_nonbreeder")*get.edge.attribute(g, "weight"))) #this need to be mode="OUT"
      # degree_to_fb <- strength(g, mode="out", weights=((sb_to == "F_breeder")*get.edge.attribute(g, "weight"))) #this need to be mode="OUT"
      # degree_to_fnb <- strength(g, mode="out", weights=((sb_to == "F_nonbreeder")*get.edge.attribute(g, "weight"))) #this need to be mode="OUT"
      # 
      # #network metrics to calculate
      # site[[j]]$focal_id <- focal_id #not necessary, but a good check to make sure I'm pulling info for the right vole
      # # site[[j]]$focal_sex <- focal_sex #not necessary, but a good check to make sure I'm pulling info for the right vole
      # #degree by sex
      # site[[j]]$F.deg <- degree_to_F
      # site[[j]]$M.deg <- degree_to_M
      # #degree by breeding status
      # site[[j]]$b.deg <- degree_to_b
      # site[[j]]$nb.deg <- degree_to_nb
      # #degree by sex-breedingstatus
      # site[[j]]$mb.deg <- degree_to_mb
      # site[[j]]$mnb.deg <- degree_to_mnb
      # site[[j]]$fb.deg <- degree_to_fb
      # site[[j]]$fnb.deg <- degree_to_fnb
      
      #####------------------ end sex- and repro- degree measures -------------------
      

    }

    #write the list 'site' as a 1st-order item in wt_net_mets_list
    wt_net_mets_list[[i]] <- site

  }

  #name the 12 1st order elements of wt_nets_list as the sites
  names(wt_net_mets_list) <- names(overlap_network_list)

  #rename the sublist items (months) for each site
  #accounting for the fact that most sites have 5 months of data but some have 4 (1 or 2 sites in 2022)
  for(i in 1:length(wt_net_mets_list)){
    ifelse( length(wt_net_mets_list[[i]]) == 5,
            names(wt_net_mets_list[[i]]) <- c("June", "July", "August", "September", "October"),
            names(wt_net_mets_list[[i]]) <- c("July", "August", "September", "October") )
  }


  ################################### ABOUT wt_net_mets_list #############################################
  #the output of wt_net_mets_list is a list of 12 1st-order items, 1 per site
  #under each site, there are 5 second-order items, 1 per trapping occasion (named by month)
  #each of those 2nd-order items (months) is a df with tag ID, month, and all the network metrics...
  #for ONLY THE ANIMALS captured on that grid, during that occasion
  ########################################################################################################


  ############## condense wt_net_mets_list to make it easier to use for analysis #######################

  #collate results
  #this collapses the 2nd order elements (network metrics for a single month) down to the 1st order element (site)
  #so now wt_net_mets_list_summary is a list of 12 dfs, each is all the network metrics for a site across all the months

  #make a list to store things
  wt_net_mets_list_summary <- list()

  #loop across all sites and collapse the dfs per occasion into one df for the site
  for(i in 1:length(wt_net_mets_list)){

    #for all 12 sites
    summary <- do.call("rbind", wt_net_mets_list[[i]])
    wt_net_mets_list_summary[[i]] <- summary
  }

  #remove the row names of the dfs (its the month, but we already have that info)
  for(i in 1:length(wt_net_mets_list_summary)){
    row.names(wt_net_mets_list_summary[[i]]) <- NULL
  }

  #name the 12 1st order elements as their sites
  names(wt_net_mets_list_summary) <- names(overlap_network_list)

  ## make net_mets_list_summary into freiggein huge df
  wt_net_mets_summary <- do.call(rbind.data.frame, wt_net_mets_list_summary)

  #clean up the df
  wt_net_mets_summary <- wt_net_mets_summary %>%
    rownames_to_column("name") %>% #row names are the sites, make that a column
    separate(name, c("site", NA)) %>% #separate the site part from the index and get rid of the index
    mutate(site = as.factor(site)) %>% #make site a factor
    rename(tag=ids)

  #save it
  saveRDS(wt_net_mets_summary, here(netmets_file))


}


