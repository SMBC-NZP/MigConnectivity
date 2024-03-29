% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simConnectivity.R
\name{simTelemetryData}
\alias{simTelemetryData}
\title{Simulate telemetry/GPS data}
\usage{
simTelemetryData(
  psi,
  sampleSize,
  originRelAbund = NULL,
  originSites = NULL,
  targetSites = NULL,
  captured = "origin",
  S = 1,
  p = 1,
  requireEveryOrigin = FALSE
)
}
\arguments{
\item{psi}{Transition probabilities between B origin and W target sites.
A matrix with B rows and W columns where rows sum to 1.}

\item{sampleSize}{If captured is "origin", either a vector of
length B with the number of simulated animals to release with geolocators
at each of the B origin sites or a single integer with the total number of
simulated animals to release with GPS at origin sites (in which
case, the origin sites will be sampled according to the relative
abundance). If captured is "target", either a vector of length W with the
number of simulated animals to release with GPS at each of the W target
sites or a single integer with the total number of simulated animals to
release at target sites (in which case, the target sites will be
sampled according to their relative abundance).}

\item{originRelAbund}{Relative abundances at B origin sites. Numeric vector
of length B that sums to 1. Optional unless providing target data and/or
sample size of length 1.}

\item{originSites}{A polygon spatial layer (sf - MULTIPOLYGON) defining the
geographic representation of sites in the origin season. If left NULL, the
simulation won't provide origin points.}

\item{targetSites}{A polygon spatial layer (sf - MULTIPOLYGON) defining the
geographic representation of sites in the target season. If left NULL, the
simulation won't provide target points.}

\item{captured}{Either "origin" (the default) or "target".}

\item{S}{Survival probabilities of released animals. Probably only
relevant for simulating archival tags. Either a matrix with B rows and W
columns (if survival depends on both origin site and target site), a vector
of length W (if survival depends only on target site), or a single number
(if survival is the same for all animals). Default 1 (all tagged animals
survive a year).}

\item{p}{Recapture probabilities of released animals. Only relevant for
simulating archival tags. Either a vector of length B (if captured on origin
and recapture depends on origin site), a vector of length W (if captured on
target and recapture depends on target site), or a single number (if
recapture is the same for all animals). Default 1 (all animals that survive
are recaptured).}

\item{requireEveryOrigin}{If TRUE, the function will throw an error if it
looks like at least one origin site has no animals released in or
migrating to it, or if it can, keep simulating until representation is
met. This helps estTransition not throw an error. Default FALSE.}
}
\value{
\code{simTelemetryData} returns a list with the elements:
\describe{
  \item{\code{originAssignment}}{Vector with true origin site of each animal}
  \item{\code{targetAssignment}}{Vector with true target site of each animal}
  \item{\code{originPointsTrue}}{True origin location of each animal, type sf,
   same projection as originSites}
  \item{\code{targetPointsTrue}}{True target location of each animal, type sf,
   same projection as targetSites}
  \item{\code{originPointsObs}}{Observed origin location of each animal that
   survived and was recaptured, type sf, same projection as originSites. Same
   as originPointsTrue when S and p==1}
  \item{\code{targetPointsObs}}{Observed target location of each animal that
   survived and was recaptured, type sf, same projection as targetSites. Same
   as targetPointsTrue when S and p==1}
  \item{\code{lived}}{0/1 vector for each animal, indicating which survived}
  \item{\code{recaptured}}{0/1 vector for each animal, indicating which were
   recaptured}
  \item{\code{input}}{List containing the inputs to function}
}
}
\description{
Simulate telemetry/GPS data
}
