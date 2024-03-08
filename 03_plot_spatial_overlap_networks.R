# load packages
library(here)
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(colorBlindness)
library(patchwork) #effectively cowplot for combining ggraph plots


#clear environment
rm(list = ls())


#######----------------- DATA CLEANING & NETWORK CONSTRUCTION -----------------------###############

#the following code will run off the overlap_nets files generated in "02_construct_spatial_overlap_networks.R"


#######----------------- NETWORKS PLOTTING -----------------------###############

## LOAD 2021 DATA ##

#network data
overlap_network_list <- readRDS(here("overlap_networks21.rds"))

#fulltrap files, for metadata
metadata <- readRDS(here("fulltrap21_03.04.24.rds")) %>% ##breeders and nonbreeders are here
  unite(sb, sex, season_breeder, remove = FALSE)

#MONTHLY centroids = - to align the points in space to monthly centroid location
#jitter to points to avoid overlapping nodes in network if two voles have centroids at the same location
centroids <- readRDS(here("monthly_centroids21.rds")) %>%
  rename(tag = Tag_ID) %>%
  mutate(jitter_x = jitter(x, 20),
         jitter_y = jitter(y, 20))

############################################



## LOAD 2022 DATA ##

# #network data
# overlap_network_list <- readRDS(here("overlap_networks22.rds"))
# 
# #fulltrap files, for metadata
# metadata <- readRDS(here("fulltrap22_03.04.24.rds")) %>% ##breeders and nonbreeders are here
#   unite(sb, sex, season_breeder, remove = FALSE)
# 
# #MONTHLY centroids = - to align the points in space to monthly centroid location
# #jitter to points to avoid overlapping nodes in network if two voles have centroids at the same location
# centroids <- readRDS(here("monthly_centroids22.rds")) %>%
#   rename(tag = Tag_ID) %>%
#   mutate(jitter_x = jitter(x, 20),
#          jitter_y = jitter(y, 20))

##############################################


######## PLOT MULTIPLE SITES ACROSS MONTHS #############

#create custom color scale
mycolors <- c("F_breeder"="#c9184a", "F_nonbreeder"="#ff758f", "M_breeder"="#023e8a", "M_nonbreeder"="#a3d5ff")


#create list to store generated graphs (1 per month per site)
graph_list <- list()

for(i in 1:length(overlap_network_list)) {
  
  #1st level is one list per site
  graph_list[[i]] <- list()
  
  for(j in 1:length(overlap_network_list[[i]])){
    
    #pull the correct site, month of network data
    data <- overlap_network_list[[i]][[j]]
    
    site.id <- names(overlap_network_list[i])
    month.id <- names(overlap_network_list[[i]][j])
    
    ## PLOTTING FROM TIDYGRAPH OBJECT ##
    
    #create metadata to get functional group for node color
    netmeta <- metadata %>% #starts with fulltrap, need to get down to a 'traits' version
      filter(site==site.id & month==month.id) %>% #subset metadata for current site, month
      group_by(tag) %>% slice(1) %>% #one entry per vole per month
      select(tag, sex, season_breeder, sb) %>% #tag should be first column to match to adj mat
      arrange(tag) #make sure they're in numeric order to match adj mat
    
    #matrix of centroid data - for node x,y location in graph
    centroidmat <- centroids %>%
      filter(site==site.id & month==month.id) %>%
      arrange(tag) %>%
      select(jitter_x,jitter_y) %>%
      as.matrix()
    #assign the row names as voleIDS just so no one gets confused
    rownames(centroidmat) <- netmeta$tag
    
    #create tidygraph object from adj matrix
    tidyg <- as_tbl_graph(data, directed=FALSE) %>% left_join(netmeta, by=c('name'='tag'))
    tg <- tidyg %>% activate(edges) %>% arrange(weight) %>% filter(weight>0)
    
    #NOTE: when plotting with ggraph:
    #ggraph(layout=...) a matrix with the node position - make sure it's in the same order as the nodes
    #format as numeric matrix two column x, y for all the nodes in the graph
    
    monthgraph <- ggraph(tg, layout=centroidmat) + #set 'layout' to matrix of monthly centroids
      geom_edge_link(aes(colour=weight, width=weight, group=I(1))) + # add edges to the plot (colored by weight)
      ## even with sorting the edges by weight, they won't necessarily plot with thickest on top
      ## unless you either 1) use 'geom_edge_link0() or 2) add 'group=I(1)' to the aes call of geom_edge_link()
      ## Matt Michalska-Smith is a lifesaver <3
      # geom_node_label(aes(label=name)) + # optional, add node labels to the plot
      geom_node_point(aes(fill=sb), color="black", pch=21, size=6) +
      scale_fill_manual(values=mycolors, 
                        labels=c("Female Breeder", "Female Nonbreeder", "Male Breeder", "Male Nonbreeder")) +
      # scale_fill_manual(values=c("#c9184a", "#ff758f", "#023e8a", "#a3d5ff"),
      #                   labels=c("Female Breeder", "Female Nonbreeder", "Male Breeder", "Male Nonbreeder")) +
      scale_edge_width(range=c(0,3), guide="none") + #scale edge width by weight
      scale_edge_colour_gradient(low="#F0F0F0", high="#000000", guide="none") + # set the (gray)scale, remove legend
      theme_void() +
      labs(title=paste(names(overlap_network_list[[i]])[j]),
           fill="Functional Group:",
           edge_colour="Spatial Overlap") +
      theme(plot.title = element_text(size= 16, hjust = 0.5),
            legend.title = element_text(size=12),
            legend.text = element_text(size=11))
    
    #the output of ggraph() is a nested list - (class ggraph, gg, ggplot)
    #needs to be added to graph_list using [[i]][[j]] <- can store a list as item [[j]] and NOT [[i]][j] <- only stores the first item of (list) g2 as [j]
    graph_list[[i]][[j]] <- monthgraph
    
    #check colors for colorblind friendly
    # cvdPlot(plot)
    
  }
  
  
}



## combine networks per site, print as .png files

## PRINT 2021 GRAPHS ##

for(i in 1:length(graph_list)){
  
  #site level (1st level of nested list graph_list)
  site <- graph_list[[i]]
  
  #use patchwork package to concatenate ggraph plots
  patchwork::wrap_plots(site, nrow=1) + plot_annotation(title=paste(names(overlap_network_list)[[i]], "2021")) +
    plot_layout(guides="collect") & theme(legend.position = "bottom",
                                          plot.margin=margin(c(0,30,0,30)),
                                          legend.box.margin=margin(20,0,0,0))
  
  #save the composite figure as .png
  ggsave(filename = paste("spatial_overlap_", "A_fxnlgrp_", names(overlap_network_list)[[i]], "_2021", ".png", sep = ""),
         plot=last_plot(),
         width=20, height=4.5, units="in",
         dpi=300)
  
}

### NOTE! The labels for the colors in the legend will be wrong if all four fxnl groups are not represented in the graph

######################################

# ## PRINT 2022 GRAPHS ##
# 
# for(i in 1:length(graph_list)){
#   
#   #site level (1st level of nested list graph_list)
#   site <- graph_list[[i]]
#   
#   #use patchwork package to concatenate ggraph plots
#   patchwork::wrap_plots(site, nrow=1) + plot_annotation(title=paste(names(overlap_network_list)[[i]], "2022")) +
#     plot_layout(guides="collect") & theme(legend.position = "bottom",
#                                           plot.margin=margin(c(0,30,0,30)),
#                                           legend.box.margin=margin(20,0,0,0))
#   #save the composite figure as .png
#   ggsave(filename = paste("spatial_overlap_", "fxnlgrp_", names(overlap_network_list)[[i]], "_2022", ".png", sep = ""),
#          plot=last_plot(),
#          width=20, height=4.5, units="in",
#          dpi=300)
#   
# }

### NOTE! The labels for the colors in the legend will be wrong if all four fxnl groups are not represented in the graph


######################################