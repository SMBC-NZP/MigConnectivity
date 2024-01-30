# MigConnectivity 0.4.5.9000

* Fixed fixing zero in estTransition (MCMC mode)

# MigConnectivity 0.4.5

* Responding to CRAN comments
* Switched relative https links to absolute (on web)
* Added references (with DOI links) to DESCRIPTION
* Added/changed return value for plot functions
* Removed messages to console (unless verbose is on) for estTransition(). Switched cat() calls to message() in getIsoMap(). Added verbose argument to weightAssign().
* Made sure par, options, and working directory revert to user settings on function exits.

# MigConnectivity 0.4.4

* Changes for submitting to CRAN
* Several edits to examples
* Function getIsoMap() (and those that call it) no longer automatically save files to the current directory. New argument mapDirectory allows the user to set the directory map files are saved it/read from. The default (NULL) saves to a temporary directory, and we recommend changing this
* Argument surface to function getIsoMap() (and those that call it) has been deprecated and will disappear in a future release

# MigConnectivity 0.4.3

* We edited the package dependencies to remove use of the packages maptools, rgdal, rgeos, and sp, which are retiring
* We sped up the analysis of raster data (usually intrinsic markers such as isotopes) in estTransition (most cases) by converting it to probability table data
* We extended estMantel to allow resampling to be weighted by abundance estimates, which should reduce bias in many cases
* We updated simCountData to reflect current BBS analysis models. modelCountDataJAGS has not been updated yet
* We reduced the length of abundExamples to 10 to reduce package size
* We made a few minor bug fixes

# MigConnectivity 0.4.2

* New function calcTransition: Calculate maximum-likelihood Ψ (transition probabilities between regions in two phases of the annual cycle) without estimating uncertainty or accounting for potential assignment error (e.g., from geolocators). If there are no CMR data, Ψ is proportional to the number of observed animals in each combination of origin and target regions. If there are CMR data, the function accounts for possible differences in reencounter probability
* New function reverseTransition: Reverse Ψ and R (origin relative abundances) estimates to calculate or estimate target site to origin site transition probabilities (Γ), target region relative abundances (W), and origin/target site combination probabilities (Π)
* New function simCMRData: Simulate capture-mark-reencounter (CMR) migratory movement data
* New function simGLData: Simulate geolocator (GL) migratory movement data
* New function simProbData: Simulate Dirichlet-based probability table data
* New function simTelemetryData: Simulate telemetry/GPS migratory movement data
* Significant edits to estTransition, including letting the user choose estimator ("bootstrap" or "MCMC", the former is the default) and allowing for more data types with each. You can now use any combination of CMR, GL, telemetry, probability table (genoscape), and raster (genoscape or isotope) data with the bootstrap estimator and CMR and telemetry data with the MCMC estimator
* New and updated vignettes that demonstrate new functionality
* Several bug fixes

# MigConnectivity 0.4.1

* Fixed bug with estMantel() choosing points off the globe if not given targetSites

# MigConnectivity 0.4.0

* New function estTransition estimates transition probabilities (psi) from banding/ringing data or any combination of geolocator, telemetry, genetics, and isotope data, with location error on either side of transition (or even both)
* New function estStrength estimates MC from previously estimated transition probabilities (from any source) and relative abundances. The existing function estMC still does some of these things, but we will cease active development on it now
* We also added basic plotting of MC and psi estimates and clearer output from printing estimates
* We updated the package to work with sf spatial data as well as sp

# MigConnectivity 0.3.1

* We added the function calcMantel, for completeness. This function calculates Mantel's correlation (rM) from either SpatialPoints objects or (if distances are already calculated) distance matrices
* We also sped up the calculation of distances used in calculating and estimating rM (in calcMantel, estMantel, and optionally in estMC)
* estMC no longer requires geoBias or geoVCov for GPS data

# MigConnectivity 0.3.0

* Added a function that allows users to generate probabilistic assignments from intrinsic markers (only stable-hydrogen isotopes are supported for the time being). The isoAssign function provides users with a variety of assignment options. These include using relative abundance as prior information, and using weighted isotope and relative abundance assignments to minimize assignment area and assignment error. This function returns a number of objects that users may be interested in such as a) individual probabilistic assignments, b) likely-vs-unlikely assignment regions, c) population level assignments and d) assignment locations that can be used to estimate migratory connectivity. The estMC function, when called using intrinsic markers uses the targetSites generated within the isoAssign function by default (can be set by user but not recommended at this time). The targetSites are generated from isotope bands equivalent to the standard deviation used to generate the isotope assignments. Preliminary simulations suggest targetSites created in this manner provide the least biased estimate of MC from intrinsic data
* Added a function to allow users to estimate the relative weighting for isotope and relative abundance that minimizes assignment area and assignment error. The resulting weights can be used in the isoAssign function
* Added support for estimating migratory connectivity strength from isotope data. This update allows users to take the output from isoAssign and estimate MC from isotope assignments while accounting for assignment uncertainty
* Added functions diffMC and diffMantel, for estimating differences (with uncertainty) between two or more INDEPENDENT estimates of MC or rM, respectively. These functions work well for comparing species
* We also added S3 types and summary/print functions for estMC, estMantel, and isoAssign outputs. The print and summary functions will only work with estimates of migratory connectivity strength from v0.3.0 or later
* Basic plot functions were added for visualizing weightAssign results and isoAssign objects
* Finally, we added reference sections to many of the help files and streamlined vignettes to increase build speed

# MigConnectivity 0.2.5

* We sped up calcMC, especially for large numbers of sites
* We made sure at least two target sites are included in each bootstrap sample in estMC, thus ensuring that inference is based on the desired number of bootstrap samples (you can't calculate MC without at least two target sites)
* Added argument maxTries to estMC and estMantel, which provides a limit to how many times the functions will try to resample points within a single bootstrap run, to prevent them from running for days without ever producing results. This issue is usually caused by combinations of targetSites, geoBias, and geoVCov that generate most sample points outside targetSites
* And one minor bug fix: nSim wasn't being passed along in estMC, so setting it away from the default did nothing

# MigConnectivity 0.2.4

* Minor extension: added list "projections" to package data store some possible map projection strings
* Minor bug fix: in the estMC function if targetSites == SpatialPolygonDataFrame - you get an error crs(not identical). However, if the targetSites ==SpatialPolygons - you don't get the error. Status: Fixed
* Minor bug fix: estMC and estMantel summary statistics don't come out if there are any NAs or NaNs in the sample. Status: Fixed
* Other minor changes to documentation and verbose settings

# MigConnectivity 0.2.3

* Added stats & utils packages to import list
* Updated code for installation so that vignettes get built locally
* Editorial changes to vignette so users can see all the code in the vignette without running off the page

# MigConnectivity 0.2.2

* Added p-value option to estMC
* Changed arguments to calcMC and estMC, including changing order and adding sampleSize argument
* Sped up functions a bit
* Expanded and clarified vignette, including new examples
* Added estMantel function for estimation of Mantel's correlation (rM) without estimating MC
* Added resampleProjection argument to estMC, allows users to set projections for bootstrapping location uncertainty

# MigConnectivity 0.1.10

* Incorporated fix for small abundances into main functions
