#ARCHIVED THIS FILE on 3.12.24 - no longer used in analysis

###### MAY 5 2023 - THIS CODE HAS BEEN ADDED TO "02_FUNCTIONS_construct_overlap_networks"  #############


## vertex attribute fuckery
#running off of 'nets_list'

library(here)
library(tidyverse)
library(igraph)

########-------------------------------- LOAD THE DATA ------------------------------------#########

# load the network data
overlap_network_list <- readRDS(file="../vole-spatial-v2/overlapnets22.rds")

#load fulltrap_traits (metadata)
fulltrap22_traits <- readRDS(file = "../vole-spatial-v2/fulltrap22_05.03.23.rds") %>%
  group_by(tag, month) %>% arrange(tag, month, occ.sess, .by_group = TRUE) %>% slice(1) %>% #one entry per vole per month
  select(year, month, season, occ.sess, site, trt, tag, samp_id, sex) %>%
  drop_na(sex) #didn't build networks with animals with no sex
  

########------------------------- LOOP THROUGH NETWORK LIST, CALC M/F DEGREE --------------------------#########


sex_degree.list <- list()

for(i in 1:length(overlap_network_list)){
  
  #for each site
  print(names(overlap_network_list[i])) 
  site.id <- names(overlap_network_list[i])
  
  site <- list()
  
  for(j in 1:length(overlap_network_list[[i]])){
    
    #for each month
    print(names(overlap_network_list[[i]][j])) 
    month.id <- names(overlap_network_list[[i]][j])
    
    ###### Calculate Male-degree, Female-degree #########
    
    adjmat <- overlap_network_list[[i]][[j]] #pull the adjmatrix for the given occasion
    
    #make a network from the adj matrix
    g <- graph.adjacency(adjmat, mode="directed", weighted=TRUE, diag = FALSE) 
        #technically, our adj matrix is symmetrical so the network is undirected, 
        #BUT because I want count the number of ties going in/out from a given node, 
        #So I need to have all the 'out' ties listed in one column of the edgelist (so effectively a directed edgelist)
        #So for these purposes, I built a directed network (but edge weights are symmetrical)
    
    #remove edges with weight less than threshold
    g2 <- delete.edges(g, which(E(g)$weight<=0.05))
    
    #subset metadata
    metadata <- fulltrap22_traits %>% filter(month==month.id & site==site.id)
    
    #set "sex" as a vertex attribute
    g2 <- set.vertex.attribute(graph=g2, name="sex", index=V(g2), value=metadata$sex)
    #save vectors
    focal_sex <- get.vertex.attribute(g2, "sex") #female is 1, male is 2
    focal_id <- get.vertex.attribute(g2, "name")
    # #visualize it
    # plot(g2, vertex.color=V(g2)$sex,
    #      edge.arrow.size = 0.8)
    
    ## code from Matt M-S for male strength/female strength
    sex_to <- get.vertex.attribute(g2, "sex")[get.edgelist(g2, names=FALSE)[,2]] 
        #since vertices already have meaningful names, call names=FALSE to return vertex indices instead
    weight_to <- get.edge.attribute(g2, "weight")[get.edgelist(g2, names=FALSE)[,2]]
    degree_to_M <- strength(g2, mode="out", weights=((sex_to == "M")*weight_to)) #this need to be mode="OUT" (the "to" individual)
    degree_to_F <- strength(g2, mode="out", weights=((sex_to == "F")*weight_to)) #this need to be mode="OUT"
    
    output <- as.data.frame(cbind(focal_id, focal_sex, degree_to_F, degree_to_M))
    
    site[[j]] <- output
    
  }
  
  sex_degree.list[[i]] <- site
  
}


#name the 1st order elements
names(sex_degree.list) <- names(overlap_network_list)

#rename the sublist items (months) for each site
#accounting for the fact that most sites have 5 months of data but some have 4
for(i in 1:length(sex_degree.list)){
  ifelse( length(sex_degree.list[[i]]) == 5,
          names(sex_degree.list[[i]]) <- c("june", "july", "aug", "sept", "oct"),
          names(sex_degree.list[[i]]) <- c("july", "aug", "sept", "oct") )
}

#collate results
#this collapses the 2nd order elements down to the 1st order element

#make a list to store things
sex_degree.summ <- list()

#loop across all sites and collapse the 2e dfs into 1 df per 1e level
for(i in 1:length(sex_degree.list)){
  #for all 12 sites
  summary <- do.call("rbind", sex_degree.list[[i]])
  #row names to column (month)
  summary <- summary %>% rownames_to_column(var="rowname") %>% separate_wider_delim(rowname, delim=".", names=c("month", NA))
  sex_degree.summ[[i]] <- summary
}

#name the 1st order elements
names(sex_degree.summ) <- names(sex_degree.list)

## make a freiggein huge df
sex_degree.results <- do.call(rbind.data.frame, sex_degree.summ) %>%
  rownames_to_column(var="rowname") %>% separate_wider_delim(rowname, delim=".", names=c("site", NA)) %>%
  mutate(focal_sex = case_when(focal_sex == "1" ~ "F",
                               focal_sex == "2" ~ "M"))


























######## Code from Matt M-S for male strength/female strength in toy network ############
# library(igraph)
# library(tidyverse)
# set.seed(0)
# g <- tidygraph::play_erdos_renyi(100, 0.3, directed=FALSE) %>%
#   mutate(type = rbernoulli(n(), 0.4) %>% as.integer() %>% factor(labels=c("M", "F"))) %>%
#   tidygraph::as.igraph()
# sex_to <- get.vertex.attribute(g, "type")[get.edgelist(g)[,2]]
# degree_to_M <- strength(g, mode="all", weights=sex_to == "M")
# degree_to_F <- strength(g, mode="all", weights=sex_to == "F")
# plot(g)
#########################################################################################










# ### NOW LOOP IT
# #could be occasion j+1 since May is coded as 1
# #and I could probably pull the site name pretty easily
# 
# 
# # Calculate network metrics - using the nets_list list (from dynamicnetcreate() function)
# for(i in 1:length(nets_list)){
#   
#   #length(nets_list) == for each of the 12 sites
#   
#   print(i)
#   site <- list()
#   
#   for(j in 1:5){
#     
#     #(j in 1:5) is for the 5 trapping occasions (NO MAY NETWORKS NOW), this is an easy enough thing to write in...
#     #and there is no easy way to pull the number of occasions (5) from nets_list
#     #it used to be: length(igraph_list[[1]][[1]]) == for each of the 5 trapping occasions - but that's really complicated
#     
#     ongrid <- nets_list[[i]][[3]][,c(1,(j+1))] #this pulls tag ids & the correct column of the netwindows df
#     ongrid <- ongrid %>% filter(.[,2] == 1) #filter for only animals on the grid
#     ids <- as.vector(ongrid$ids)
#     
#     adjmat <- nets_list[[i]][[2]][,,j] #pull the adjmatrix for the given occasion (includes all animals)
#     
#     adjmat <- adjmat[ids,ids] #subset the adjmatrix for only the animals on the grid that month
#     
#     #make a network from the subset adj matrix for a given occasion (j)
#     inet <- graph.adjacency(adjmat >= 1, mode="undirected", weighted=NULL) 
#     #this code will make a binary network where there is only one edge between 2 voles, even if they overlapped multiple times in 48hr
#     
#     tag <- ids #pull the tag numbers for all the animals on that grid
#     month <- rep(j,length(ids)) #this puts j==occasion # in a column for all animals
#     
#   }
# }
