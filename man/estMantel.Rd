% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estConnectivity.R
\name{estMantel}
\alias{estMantel}
\title{Estimate Mantel correlation (rM) from geolocator and/or GPS data.}
\usage{
estMantel(
  targetPoints,
  originPoints,
  isGL,
  geoBias = NULL,
  geoVCov = NULL,
  targetSites = NULL,
  nBoot = 1000,
  nSim = 1000,
  verbose = 0,
  alpha = 0.05,
  resampleProjection = "ESRI:102010",
  maxTries = 300,
  maintainLegacyOutput = FALSE
)
}
\arguments{
\item{targetPoints}{A \code{SpatialPoints} from sp or \code{POINTS} from sf object, with length number of
animals tracked.  Each point indicates the point estimate location in
the non-release season}

\item{originPoints}{A \code{SpatialPoints} from sp or \code{POINTS} from sf object, with length number of
animals tracked.  Each point indicates the release location of an animal}

\item{isGL}{Indicates whether or which animals were tracked with geolocators
Should be either single TRUE or FALSE value, or vector with length of
number of animals tracked, with TRUE for animals in
\code{targetPoints} with geolocators and FALSE for animals with GPS.}

\item{geoBias}{For GL data, vector of length 2 indicating expected bias
in longitude and latitude of \code{targetPoints}, in
\code{resampleProjection} units (default meters)}

\item{geoVCov}{For GL data, 2x2 matrix with expected variance/covariance
in longitude and latitude of \code{targetPoints}, in
\code{resampleProjection} units (default meters)}

\item{targetSites}{A \code{SpatialPolygons}, \code{SpatialPolygonsDataFrame},
 or \code{POLYGONS} sf object indicating valid target location(s).  Not needed unless you want
to mask out certain areas (e.g. water)}

\item{nBoot}{Number of bootstrap runs. Animals are sampled with replacement for each,
to estimate sampling uncertainty}

\item{nSim}{Tuning parameter for GL data. Affects only the speed; 1000 seems
to work well with our data.  Should be integer > 0}

\item{verbose}{0 (default) to 3. 0 prints no output during run. 1 prints
a line every 100 bootstraps.  2 prints a line every bootstrap.
3 also prints the number of draws (for tuning nSim for GL data only)}

\item{alpha}{Level for confidence/credible intervals provided.}

\item{resampleProjection}{Projection when sampling from geolocator
bias/error. This projection needs units = m. Default is Equidistant
Conic. The default setting preserves distances around latitude = 0 and
longitude = 0. Other projections may work well, depending on the location
of \code{targetSites}}

\item{maxTries}{Maximum number of times to run a single GL bootstrap before
exiting with an error.  Default is 300.  Set to NULL to never stop.  This
parameter was added to prevent GL setups where some sample points never
land on target sites from running indefinitely}

\item{maintainLegacyOutput}{version 0.4.0 of \code{MigConnectivity}
updated the structure of the estimates. If you have legacy code that refers
to elements within a \code{estMigConnectivity} object, you can set this
to TRUE to also keep the old structure. Defaults to FALSE}
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
Resampling of uncertainty for rM from SpatialPoints geolocators and/or GPS
data.
}
\examples{
data('OVENdata')
rM1 <- estMantel(isGL=OVENdata$isGL,#Logical vector: light-level GL(T)/GPS(F)
                 geoBias = OVENdata$geo.bias, # Geolocator location bias
                 geoVCov = OVENdata$geo.vcov, # Location covariance matrix
                 targetSites = OVENdata$targetSites, # Non-breeding target sites
                 originPoints = OVENdata$originPoints, # Capture Locations
                 targetPoints = OVENdata$targetPoints, # Target locations
                 verbose = 1,   # output options
                 nBoot = 100, # This is set low for example
                 resampleProjection = raster::projection(OVENdata$targetSites))
rM1
str(rM1, max.level = 2)
}
\references{
Cohen, E. B., J. A. Hostetler, M. T. Hallworth, C. S. Rushing, T. S. Sillett,
and P. P. Marra. 2018. Quantifying the strength of migratory connectivity.
Methods in Ecology and Evolution 9: 513 - 524.
\href{http://doi.org/10.1111/2041-210X.12916}{doi:10.1111/2041-210X.12916}
}
\seealso{
\code{\link{estMC}}
}
