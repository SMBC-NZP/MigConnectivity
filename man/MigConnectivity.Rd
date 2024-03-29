% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MigConnectivity.R
\docType{package}
\name{MigConnectivity}
\alias{MigConnectivity}
\title{MigConnectivity: A package for quantifying migratory connectivity pattern and
strength for migratory animals}
\description{
The MigConnectivity package allows the user to estimate or calculate
transition probabilities for migratory animals between any two phases of the
annual cycle, using a variety of different data types, with the function
\code{\link{estTransition}}. The user can also estimate or calculate the
strength of migratory connectivity (MC), a standardized metric to quantify
the extent to which populations co-occur between two phases of the annual
cycle. MC is independent of data type and accounts for the relative abundance
of populations distributed across a seasonal range. The package includes
functions to estimate MC (\code{\link{estStrength}}) and the more traditional
metric of migratory connectivity strength (Mantel correlation; rM;
\code{\link{estMantel}}) incorporating uncertainty from multiple
sources of sampling error. Description of the MC metric can be found in Cohen
et al. (2018).
}
\section{Key MigConnectivity Functions}{


\code{\link{estTransition}}: Estimate psi (transition probabilities between
 locations in two phases of the annual cycle)

\code{\link{estStrength}}: Estimate MC, migratory connectivity strength
}

