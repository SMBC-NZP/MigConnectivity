/dontrun{
library(raster)
library(sf)
library(ks)
library(spdep)
library(mapview)

# This script simulates the type of surface observed with genoscape data
# First, we need a blank canvas to start from with real spatial information
# we'll start with a map of the United States. The simulated genoscape is 
# created using the following steps - 
# 1) Aggregate states into N 'populations'
# 2) Generate random points within each 'population' 
# 3) Generate a kernel density to simulate the probability within 
#    each N population (note - an empty raster needs to be made that
#                       spans all N populations)
# 4) Stack the populations together into a single surface 

#1) Empty template then aggregate  
States <- raster::getData("GADM", level = 1, country = "United States")

# convert to sf
States_sf <- st_as_sf(States)

# Make clusters based on geography 
States_sf$Groups <- 0

# South east group
Group1 <- c("Virginia","North Carolina","South Carolina")

# south 
Group2 <- c("Georgia","Florida","Alabama","Mississippi","Louisiana")

# South west
Group3 <- c("Texas","Oklahoma","New Mexico","Arizona")

# northern 
Group4 <- c("Kansas","Arkansas","Missouri","Tennessee","Kentucky","West Virginia","Maryland")

# Set groups # 
States_sf$Groups[States_sf$NAME_1 %in% Group1] <- 1
States_sf$Groups[States_sf$NAME_1 %in% Group2] <- 2
States_sf$Groups[States_sf$NAME_1 %in% Group3] <- 3
States_sf$Groups[States_sf$NAME_1 %in% Group4] <- 4

# merge those groups 
polyG1 <- st_union(States_sf[States_sf$Groups==1,])
polyG2 <- st_union(States_sf[States_sf$Groups==2,])
polyG3 <- st_union(States_sf[States_sf$Groups==3,])
polyG4 <- st_union(States_sf[States_sf$Groups==4,])

# create the genetic surfaces # 
genPops <- simGeneticPops(popBoundaries = list(polyG1,polyG2,polyG3,polyG4),
                          popNames = c("South_east","South","South_west","North"))
						  
# plot the output #
plot(genPops)
}