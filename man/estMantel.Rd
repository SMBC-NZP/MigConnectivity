% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estConnectivity.R
\name{estMantel}
\alias{estMantel}
\alias{estCorr}
\title{Estimate Mantel correlation (rM) from geolocator, GPS, and/or raster data.}
\usage{
estMantel(
  targetPoints = NULL,
  originPoints = NULL,
  isGL,
  geoBias = NULL,
  geoVCov = NULL,
  targetSites = NULL,
  nBoot = 1000,
  nSim = ifelse(any(isRaster & isGL), 5000, ifelse(any(isGL), 1000, ifelse(any(isRaster),
    10, 1))),
  verbose = 0,
  alpha = 0.05,
  resampleProjection = "ESRI:102010",
  maxTries = 300,
  maintainLegacyOutput = FALSE,
  originSites = NULL,
  isTelemetry = !isGL,
  isRaster = FALSE,
  captured = "origin",
  geoBiasOrigin = geoBias,
  geoVCovOrigin = geoVCov,
  targetRaster = NULL,
  originRaster = NULL,
  dataOverlapSetting = c("dummy", "none", "named"),
  originRelAbund = NULL,
  targetRelAbund = NULL
)

estCorr(
  targetPoints = NULL,
  originPoints = NULL,
  isGL,
  geoBias = NULL,
  geoVCov = NULL,
  targetSites = NULL,
  nBoot = 1000,
  nSim = ifelse(any(isRaster & isGL), 5000, ifelse(any(isGL), 1000, ifelse(any(isRaster),
    10, 1))),
  verbose = 0,
  alpha = 0.05,
  resampleProjection = "ESRI:102010",
  maxTries = 300,
  maintainLegacyOutput = FALSE,
  originSites = NULL,
  isTelemetry = !isGL,
  isRaster = FALSE,
  captured = "origin",
  geoBiasOrigin = geoBias,
  geoVCovOrigin = geoVCov,
  targetRaster = NULL,
  originRaster = NULL,
  dataOverlapSetting = c("dummy", "none", "named"),
  originRelAbund = NULL,
  targetRelAbund = NULL
)
}
\arguments{
\item{targetPoints}{A \code{POINTS} from sf
object, with length number of animals tracked.  Each point indicates the
point estimate location in the non-release season.}

\item{originPoints}{A \code{POINTS} from sf
object, with length number of animals tracked.  Each point indicates the
release location of an animal.}

\item{isGL}{Indicates whether or which animals were tracked with geolocators
Should be either single TRUE or FALSE value, or vector with length of
number of animals tracked, with TRUE for animals in  \code{targetPoints}
with geolocators and FALSE for animals without.}

\item{geoBias}{For GL data, vector of length 2 indicating expected bias
in longitude and latitude of \code{targetPoints}, in
\code{resampleProjection} units (default meters).}

\item{geoVCov}{For GL data, 2x2 matrix with expected variance/covariance
in longitude and latitude of \code{targetPoints}, in
\code{resampleProjection} units (default meters).}

\item{targetSites}{A \code{SpatialPolygons}, \code{SpatialPolygonsDataFrame},
or \code{POLYGONS} sf object indicating valid target location(s). Not
needed unless you want to mask out certain areas (e.g. water) and
\code{captured} is "origin" or you want to use a weighted bootstrap based on
\code{targetRelAbund} for animals captured on the target side.}

\item{nBoot}{Number of bootstrap runs. Animals are sampled with replacement
for each, to estimate sampling uncertainty.}

\item{nSim}{Tuning parameter for GL or raster data. Affects only the speed;
1000 seems to work well with our GL data.  Should be integer > 0.}

\item{verbose}{0 (default) to 3. 0 prints no output during run. 1 prints
a line every 100 bootstraps.  2 prints a line every bootstrap.
3 also prints the number of draws (for tuning nSim only).}

\item{alpha}{Level for confidence/credible intervals provided.}

\item{resampleProjection}{Projection when sampling from geolocator
bias/error. This projection needs units = m. Default is Equidistant
Conic. The default setting preserves distances around latitude = 0 and
longitude = 0. Other projections may work well, depending on the location
of \code{targetPoints}.}

\item{maxTries}{Maximum number of times to run a single GL bootstrap before
exiting with an error.  Default is 300.  Set to NULL to never stop.  This
parameter was added to prevent GL setups where some sample points never
land on target sites from running indefinitely.}

\item{maintainLegacyOutput}{version 0.4.0 of \code{MigConnectivity}
updated the structure of the estimates. If you have legacy code that refers
to elements within a \code{estMigConnectivity} object, you can set this
to TRUE to also keep the old structure. Defaults to FALSE.}

\item{originSites}{A \code{SpatialPolygons}, \code{SpatialPolygonsDataFrame},
or \code{POLYGONS} sf object indicating valid origin location(s). Not
needed unless you want to mask out certain areas (e.g. water) and
\code{captured} is "target" or you want to use a weighted bootstrap based on
\code{originRelAbund} for animals captured on the origin side.}

\item{isTelemetry}{Indicates whether or which animals were tracked with
telemetry/GPS (no location uncertainty on either end).
Should be either single TRUE or FALSE value, or vector with length of
number of animals tracked, with TRUE or FALSE for each animal in data.}

\item{isRaster}{Indicates whether or which animals were tracked with
intrinsic markers (e.g., genetics or isotopes), with location uncertainty
expressed as a raster of probabilities by grid cells, either in
\code{targetRaster} or \code{originRaster}. Should be either single TRUE or
FALSE value, or vector with length of number of animals tracked, with TRUE
or FALSE for each animal in data.}

\item{captured}{Indicates whether or which animals were captured in the
origin sites, the target sites, or neither (another phase of the annual
cycle). Location uncertainty will only be applied where the animal was not
captured. So this doesn't matter for telemetry data. Should be either single
"origin" (default), "target", or "neither" value, or a character vector with
length of number of animals tracked, with "origin", "target", or "neither"
for each animal.}

\item{geoBiasOrigin}{For GL data where \code{captured}!="origin", vector of
length 2 indicating expected bias in longitude and latitude of
\code{originPoints}, in \code{resampleProjection} units (default meters).}

\item{geoVCovOrigin}{For GL data where \code{captured}!="origin", 2x2 matrix
with expected variance/covariance in longitude and latitude of
\code{targetPoints}, in \code{resampleProjection} units (default meters).}

\item{targetRaster}{For intrinsic tracking data, the results of
\code{isoAssign} or a similar function of class \code{intrinsicAssign} or
class \code{RasterBrick}/\code{RasterStack}, for example from the package
\code{assignR}. In any case, it expresses location uncertainty on target
range, through a raster of probabilities by grid cells.}

\item{originRaster}{For intrinsic tracking data, the results of
\code{isoAssign} or a similar function of class \code{intrinsicAssign} or
class \code{RasterBrick}/\code{RasterStack}, for example from the package
\code{assignR}. In any case, it expresses location uncertainty on origin
range, through a raster of probabilities by grid cells.}

\item{dataOverlapSetting}{When there is more than one type of data, this
setting allows the user some flexibility for clarifying which type(s) of
data apply to which animals. Setting "dummy" (the default) indicates that
there are dummy values within each dataset for the animals that isGL,
isTelemetry, etc. don't have that data type (FALSE values). If no animals
have a data type, no dummy values are required. If no animals have more than
one type of data, the user can simplify processing their data by choosing
setting "none" here. In this case, there should be no dummy values, and only
the animals with a type of data should be included in that dataset. The
third setting ("named") is not yet implemented, but will eventually allow
another way to allow animals with more than one type of data with named
animals linking records. When there is only one type of data, it is fastest
to leave this on the default.}

\item{originRelAbund}{the proportion of the total abundance in each of B
\code{originSites}. Used to set up the bootstrap to be weighted by relative
abundance (for animals captured on the origin side). Either a numeric vector
of length B that sums to 1, or an mcmc object (such as is produced by
\code{\link{modelCountDataJAGS}}) or matrix with at least B columns.
If there are more than B columns, the relevant columns should be
labeled "relN[1]" through "relN[B]". Optional, but if you don't set it and
at least some animals are captured on the origin side, there's potential for
rM to be biased (if sampling isn't proportional to abundance).}

\item{targetRelAbund}{the proportion of the total abundance in each of W
\code{targetSites}. Used to set up the bootstrap to be weighted by relative
abundance (for animals captured on the target side). Either a numeric vector
of length W that sums to 1, or an mcmc object (such as is produced by
\code{\link{modelCountDataJAGS}}) or matrix with at least W columns.
If there are more than W columns, the relevant columns should be
labeled "relN[1]" through "relN[W]". Optional, but if you don't set it and
at least some animals are captured on the target side, there's potential for
rM to be biased (if sampling isn't proportional to abundance).}
}
\value{
\code{estMantel} returns a list with elements:
\describe{
  \item{\code{corr}}{List containing estimates of rM:
   \itemize{
    \item{\code{sample}} \code{nBoot} sampled values for Mantel
     correlation. Provided to allow the user to compute own summary
     statistics.
    \item{\code{mean, se, simpleCI, bcCI, median, point}} Summary
     statistics for Mantel correlation bootstraps.
   }
  }
  \item{\code{input}}{List containing the inputs to \code{estMantel}}
}
}
\description{
Resampling of uncertainty for migratory connectivity strength, as quantified
by Mantel correlation (rM), from geolocators, GPS, and/or raster (e.g.,
genoscape or isotope) data.
}
\examples{
\donttest{
data('OVENdata')
rM1 <- estMantel(isGL=OVENdata$isGL,#Logical vector: light-level GL(T)/GPS(F)
                 geoBias = OVENdata$geo.bias, # Geolocator location bias
                 geoVCov = OVENdata$geo.vcov, # Location covariance matrix
                 targetSites = OVENdata$targetSites,#Nonbreeding/target sites
                 originPoints = OVENdata$originPoints, # Capture Locations
                 targetPoints = OVENdata$targetPoints, # Target locations
                 verbose = 1,   # output options
                 nBoot = 10, # This is set low for example
                 resampleProjection = sf::st_crs(OVENdata$targetSites))
rM1
str(rM1, max.level = 2)
}
}
\references{
Cohen, E. B., J. A. Hostetler, M. T. Hallworth, C. S. Rushing, T. S. Sillett,
and P. P. Marra. 2018. Quantifying the strength of migratory connectivity.
Methods in Ecology and Evolution 9: 513 - 524.
\doi{10.1111/2041-210X.12916}
}
\seealso{
\code{\link{estMC}}
}
