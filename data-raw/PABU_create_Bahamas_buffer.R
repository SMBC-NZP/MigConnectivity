# This script was written to add a buffer around the Bahamas #

library(raster)
library(sf)
library(MigConnectivity)


origTargetSites <- sf::st_read("data-raw/PABU_nonbreeding_regions.shp")

# there are 7 geometries in the file - need to determine which one
# the bahamas belongs to

# cuba and bahamas are in the second geometry
# plot(origTargetSites[2,]$geometry)

# this splits EVERY polygon into it's own row #
exploded_polys <- st_cast(origTargetSites[2,],"POLYGON")

# we need a polygon to intersect just the bahamas so we can
# id those

bahama_square <- 'POLYGON((-79.651 27.228,-74.972 27.290,-74.596 23.371,-79.207 23.559,-79.651 27.228))'

bahama_sq <- sf::st_as_sfc(bahama_square, crs = 4326)

mat_of_in <- sf::st_within(exploded_polys, bahama_sq, sparse = FALSE)

exploded_polys$merge <- apply(mat_of_in,1,sum)

bahamas <- exploded_polys[exploded_polys$merge==1,]

bahamas <- st_union(bahamas)

bahama_buff <- bahamas %>%
               st_transform('ESRI:102010') %>%
               st_buffer(dist=25000) %>%
               st_transform(4326)

# quick plot to make sure they don't overlap with Cuba
# plot(exploded_polys$geometry)
# plot(bahama_buff, add = TRUE)
# plot(bahamas,add = TRUE)

# now add bahamas buffer back in #

zone2 <- st_union(exploded_polys[exploded_polys$merge==0,], bahama_buff) %>%
         st_combine()

zone2df <- st_set_geometry(origTargetSites[1,],NULL)

newZone <- st_as_sf(zone2df,geometry=zone2)

# now add back into larger

newTargetRegions <- rbind(origTargetSites[1,],
                          newZone,
                          origTargetSites[3:7,])

#sf::st_write(newTargetRegions,
#             "data-raw/PABU_nonbreeding_regions_buffer.shp")
