# Generate projection strings list

#########################
# Eight Psi scenarios
# Breeding as rows, Non-breeding as columns
#########################
WGS84<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

Lambert<-"+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96
  +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

EquidistConic<-"+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007
  +b=6371007 +units=m +no_defs"

NorthAmerica <- "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0
  +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

SouthAmerica <- "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=-5 +lat_2=-42 +x_0=0
  +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

Europe <- "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=43 +lat_2=62 +x_0=0 +y_0=0
  +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

AsiaNorth <- "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=15 +lat_2=65 +x_0=0 +y_0=0
  +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

AsiaSouth <- "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=7 +lat_2=-32 +x_0=0 +y_0=0
  +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

Africa <- "+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=20 +lat_2=-23 +x_0=0 +y_0=0
  +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

Arctic <- "+proj=aeqd +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84
  +datum=WGS84 +units=m +no_defs"

Antarctic <- "+proj=aeqd +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84
  +datum=WGS84 +units=m +no_defs"

projections <- list(EquidistConic = EquidistConic, WGS84 = WGS84,
                    Lambert = Lambert, NorthAmerica = NorthAmerica,
                    SouthAmerica = SouthAmerica, Europe = Europe,
                    Africa = Africa, AsiaNorth = AsiaNorth,
                    AsiaSouth = AsiaSouth, Arctic = Arctic,
                    Antarctic = Antarctic)
devtools::use_data(projections, overwrite = T)
