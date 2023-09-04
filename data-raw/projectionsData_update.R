
projections <- vector('list',11)


names(projections) <- c('EquidistConic','WGS84','Lambert','NorthAmerica','SouthAmerica',
                      'Europe','Africa','AsiaNorth','AsiaSouth','Arctic','Antarctic')

# World Equidistant Conic #
# the short code doesn't work so need to provide full WKT #
projections[[1]] <- 'PROJCS["World_Equidistant_Conic",GEOGCS["GCS_WGS_1984",DATUM["WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Equidistant_Conic"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["Longitude_Of_Center",0],PARAMETER["Standard_Parallel_1",60],PARAMETER["Standard_Parallel_2",60],PARAMETER["Latitude_Of_Center",0],UNIT["Meter",1],AUTHORITY["EPSG","54027"]]'

# WGS84
projections[[2]] <- '+init=epsg:4326'

# Lambert
projections[[3]] <- '+init=epsg:9820'

# North America Equidistant
projections[[4]] <- 'PROJCS["North_America_Equidistant_Conic",GEOGCS["GCS_North_American_1983",DATUM["North_American_Datum_1983",SPHEROID["GRS_1980",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Equidistant_Conic"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["Longitude_Of_Center",-96],PARAMETER["Standard_Parallel_1",20],PARAMETER["Standard_Parallel_2",60],PARAMETER["Latitude_Of_Center",40],UNIT["Meter",1],AUTHORITY["EPSG","102010"]]'

# South America Equidistant
projections[[5]] <- 'PROJCS["South_America_Equidistant_Conic",GEOGCS["GCS_South_American_1969",DATUM["South_American_Datum_1969",SPHEROID["GRS_1967_Truncated",6378160,298.25]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Equidistant_Conic"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["Longitude_Of_Center",-60],PARAMETER["Standard_Parallel_1",-5],PARAMETER["Standard_Parallel_2",-42],PARAMETER["Latitude_Of_Center",-32],UNIT["Meter",1],AUTHORITY["EPSG","102032"]]'

# Europe
projections[[6]] <- 'PROJCS["Europe_Equidistant_Conic",GEOGCS["GCS_European_1950",DATUM["European_Datum_1950",SPHEROID["International_1924",6378388,297]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Equidistant_Conic"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["Longitude_Of_Center",10],PARAMETER["Standard_Parallel_1",43],PARAMETER["Standard_Parallel_2",62],PARAMETER["Latitude_Of_Center",30],UNIT["Meter",1],AUTHORITY["EPSG","102031"]]'

# Africa
projections[[7]] <- 'PROJCS["Africa_Equidistant_Conic",GEOGCS["GCS_WGS_1984",DATUM["WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Equidistant_Conic"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["Longitude_Of_Center",25],PARAMETER["Standard_Parallel_1",20],PARAMETER["Standard_Parallel_2",-23],PARAMETER["Latitude_Of_Center",0],UNIT["Meter",1],AUTHORITY["EPSG","102023"]]'

# AsiaNorth
projections[[8]] <- 'PROJCS["Asia_North_Equidistant_Conic",GEOGCS["GCS_WGS_1984",DATUM["WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Equidistant_Conic"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["Longitude_Of_Center",95],PARAMETER["Standard_Parallel_1",15],PARAMETER["Standard_Parallel_2",65],PARAMETER["Latitude_Of_Center",30],UNIT["Meter",1],AUTHORITY["EPSG","102026"]]'

# AsiaSouth
projections[[9]] <- 'PROJCS["Asia_South_Equidistant_Conic",GEOGCS["GCS_WGS_1984",DATUM["WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Equidistant_Conic"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["Longitude_Of_Center",125],PARAMETER["Standard_Parallel_1",7],PARAMETER["Standard_Parallel_2",-32],PARAMETER["Latitude_Of_Center",-15],UNIT["Meter",1],AUTHORITY["EPSG","102029"]]'

# Arctic
projections[[10]] <- 'PROJCS["North_Pole_Azimuthal_Equidistant",GEOGCS["GCS_WGS_1984",DATUM["WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Azimuthal_Equidistant"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["Central_Meridian",0],PARAMETER["Latitude_Of_Origin",90],UNIT["Meter",1],AUTHORITY["EPSG","102016"]]'

# Antarctic
projections[[11]] <- 'PROJCS["South_Pole_Azimuthal_Equidistant",GEOGCS["GCS_WGS_1984",DATUM["WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Azimuthal_Equidistant"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["Central_Meridian",0],PARAMETER["Latitude_Of_Origin",-90],UNIT["Meter",1],AUTHORITY["EPSG","102019"]]'

usethis::use_data(projections, internal = FALSE, overwrite = TRUE)

