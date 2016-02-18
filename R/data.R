#' Example transition probabilities (psis) between origin and target sites
#'
#' A dataset containing example psi matrices used in Cohen et al. (in prep.).
#'
#' @format A named list with 8 transition probability matrices in it. The
#' direction is from origin site (rows) to target sites (columns), so each
#' row of each matrix sums to 1.  The psi matrices are:
#' \itemize{
#'   \item Full Mix: Full mixing from all origin sites to all target sites
#'   \item Avoid One Site: All origin sites have the same transition
#'      probabilities, mostly avoiding target site 4
#'   \item Full Connectivity: Each origin site transitions to only one target
#'      site
#'   \item Half Mix: Origin sites A and B mix fully between target sites 1 and
#'      2, but don't move to target sites 3 or 4, while origin sites C and D
#'      mix fully between target sites 3 and 4, but don't move to target sites
#'      1 or 2
#'   \item Low: Simulation scenario labelled "Moderate Connectivity" in Cohen
#'      et al. (2014)
#'   \item Medium: Simulation scenario labelled "Strong Connectivity" in Cohen
#'      et al. (2014)
#'   \item One Site Preference: Three origin sites have full mixing, but origin
#'      site D only goes to target site 4
#'   \item Negative: Artificial transition probability scenario developed to
#'      produce a negative MC value under some circumstances
#' }
"samplePsis"

#' Example origin and target site positions and distances on a 2-D plane
#'
#' \code{sampleOriginPos} is a dataset containing example origin site positions
#' from 12 scenarios used in Cohen et al. (in prep.).  For the same 12
#' scenarios, \code{sampleOriginDist} contains the origin site distances,
#' \code{sampleTargetPos} contains the target site positions, and
#' \code{sampleTargetDist} contains the target site distances.
#'
#' @format Each dataset is a named list with 12 matrices in it, representing 12
#' scenarios. The position matrices each have 2 columns (x and y position) and
#' 4 rows (for each origin or target site).  The distance matrices are
#' symmetrical and 4 x 4. The 12 scenarios are:
#' \itemize{
#'   \item Linear: Both origin and target sites arranged in horizontal linear
#'      fashion, with equal distances between each adjacent site
#'   \item B Dist BC*2: Linear, but the central origin sites are twice as far
#'      from each other as the edge sites are from the adjacent origin sites
#'   \item B Dist BC/2: Linear, but the central origin sites are half as far
#'      from each other as the edge sites are from the adjacent origin sites
#'   \item B Dist CD*2: Linear, but the last two origin sites are twice as far
#'      from each other as the other adjacent origin sites
#'   \item B Dist CD/2: Linear, but the last two origin sites are half as far
#'      from each other as the other adjacent origin sites
#'   \item B Grid: Origin sites arranged on a grid, target sites arranged
#'      linearly, both with all adjacent sites (excluding diagonals)
#'      equidistant
#'   \item NB Dist 23*2: Linear, but the central target sites are twice as far
#'      from each other as the edge sites are from the adjacent target sites
#'   \item NB Dist 23/2: Linear, but the central target sites are half as far
#'      from each other as the edge sites are from the adjacent target sites
#'   \item NB Dist 34*2: Linear, but the last two target sites are twice as far
#'      from each other as the other adjacent target sites
#'   \item NB Dist 34/2: Linear, but the last two target sites are half as far
#'      from each other as the other adjacent target sites
#'   \item NB Grid: Target sites arranged on a grid, origin sites arranged
#'      linearly, both with all adjacent sites (excluding diagonals)
#'      equidistant
#'   \item B/NB Grid: Origin and target sites each arranged on a grid, both
#'      with all adjacent sites (excluding diagonals) equidistant
#' }
"sampleOriginPos"

#' @rdname sampleOriginPos
"sampleOriginDist"

#' @rdname sampleOriginPos
"sampleTargetPos"

#' @rdname sampleOriginPos
"sampleTargetDist"

#' Example origin site abundances and relative abundances
#'
#' \code{sampleOriginN} is a dataset containing example origin site abundances
#' from 5 scenarios used in Cohen et al. (in prep.).  For the same 5
#' scenarios, \code{sampleOriginRelN} contains the relative abundances.
#'
#' @format Each dataset is a named list with 5 vectors in it. Each vector has
#' 4 elements (for the 4 origin sites).  The relative abundance vectors each
#' sum to 1. The 5 scenarios are:
#' \itemize{
#'   \item Base: Equal abundance at each origin site
#'   \item B Doub: The second origin site has twice the abundance of the other
#'      three sites
#'   \item B Half: The second origin site has half the abundance of the other
#'      three sites
#'   \item D Doub: The last origin site has twice the abundance of the other
#'      three sites
#'   \item D Half: The last origin site has half the abundance of the other
#'      three sites
#' }
"sampleOriginN"

#' @rdname sampleOriginN
"sampleOriginRelN"

#' Example RMark transition probability estimates from simulated data
#'
#' A dataset containing RMark transition probability estimates from simulated
#' mark-recapture-recovery data from Cohen et al. (2014).  These all represent
#' the intermediate scenario for all settings (moderate connectivity, low
#' re-encounter, 100,000 banded in each breeding area).  Each estimate can
#' be used in \code{estMC} function to estimate MC with uncertainty.
#'
#' @format A list with 100 RMark estimates in it.
"cmrExamples"

#' Example relative abundance estimates from simulated data
#'
#' A dataset containing mcmc relative abundance estimates from simulated
#' BBS-type data from Cohen et al. (in prep). Each estimate can
#' be used in \code{estMC} function to estimate MC with uncertainty.
#'
#' @format A list with 100 mcmc (coda) estimates in it.
"abundExamples"

#' Relative abundance estimates of Ovenbird within two breeding regions
#'
#' A dataset containing relative abundance estimates from BBS data
#' reported in Cohen et al. (in prep). These estimates can
#' be used in \code{estMC} function as \code{originRelAbund} in conjuction
#' with archival light-level geolocator and GPS locations.
#'
#' @format A vector of length two with relative abundance estimates.
"originRelAbund"

#' Archival light-level geolocator bias estimates
#'
#' Location bias estimates in light-level geolocator estimates calculated
#' using birds captured at known locations in Florida, Jamaica and Puerto Rico.
#' Location bias is reported in meters.
#'
#' @format A vector of length two with bias esimates in geolocator locations.
"geo.bias"

#' Archival light-level geolocator covariance matrix
#'
#' Covariance estimates in light-level geolocator estimates calculated
#' using birds captured at known locations in Florida, Jamaica and Puerto Rico.
#' Covariance is reported in meters.
#'
#' @format A 2x2 matrix of covariance estimates.
"geo.vcov"

#' Ovenbird capture locations
#'
#' Capture locations for 39 Ovenbirds caught during the breeding season who carried
#' either a light-level geolocator or PinPoint-10 GPS tag. Ovenbirds were captured
#' at Hubbard Brook Experimental Forest, NH and Jug Bay Wetland Sanctuary, MD. These
#' data are used as \code{originPoints} in the \code{estMC} function. \code{coords.x1}
#' and \code{coords.x2} represent the longitude and latitude of the capture sites,
#' respectively. The data are projected in Lambert Conformal Conic.
#'
#' "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80
#' +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"
#'
#' @format SpatialPoints in Lambert Conformal Conic
"originPoints"

#' Ovenbird non-breeding locations
#'
#' Non-breeding locations for 39 Ovenbirds caught during the breeding season who carried
#' either a light-level geolocator or PinPoint-10 GPS tag. Ovenbirds were captured
#' at Hubbard Brook Experimental Forest, NH and Jug Bay Wetland Sanctuary, MD. These
#' data are used as \code{originPoints} in the \code{estMC} function. \code{coords.x1}
#' and \code{coords.x2} represent the longitude and latitude of the capture sites,
#' respectively. The data are projected in Lambert Conformal Conic.
#'
#' "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80
#' +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"
#'
#' @format \code{SpatialPoints} in Lambert Conformal Conic
"targetPoints"

#' Archival light-level geolocator or PinPoint-10 GPS tag
#'
#' \code{logical} vector indicating whether location estimates were obtained with a
#' light-level geolocator (\code{TRUE}) or PinPoint-10 GPS tag (\code{FALSE})
#'
#' @format \code{logical} of length 39
"isGL"

#' Ovenbird non-breeding distribution target sites
#'
#'Non-breeding distribution target sites used in Cohen et al. (in prep) to estimate
#'MC of Ovenbirds tracked with light-level geolocators and PinPoint-10 GPS tags. There
#'are three non-breeding target sites 1) Florida, United States, 2) Cuba, and 3) Hispaniola
#'(Dominican Republic and Haiti).
#'
#'
#' @format \code{SpatialPolygons}
"targetSites"

#' Distance between Ovenbird non-breeding distribution target sites
#'
#'The pairwise Great Circle Distance between the center of the polygons contained within
#'\code{targetSites}. See "Ovenbird non-breeding distribution target sites" or \code{targetSites}
#'
#'
#' @format square distance matrix
"targetDist"

#' Ovenbird breeding distribution origin sites
#'
#'Breeding distribution prigin sites used in Cohen et al. (in prep) to estimate
#'MC of Ovenbirds tracked with light-level geolocators and PinPoint-10 GPS tags. There
#'are two breeding origin sites, one that encompasses NH and another that encompasses MD capture
#'deployment locations.
#'
#'
#' @format \code{SpatialPolygons}
"originSites"

#' Distance between Ovenbird breeding distribution origin sites
#'
#'The pairwise Great Circle Distance between the center of the polygons contained within
#'\code{originSites}. See "Ovenbird breeding distribution origin sites" or \code{originSites}
#'
#'
#' @format square distance matrix
"originDist"
