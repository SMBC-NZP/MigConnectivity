## ----setup, include = FALSE------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo = FALSE, message = FALSE, warning = FALSE------------------------------------------
oo <- getOption("rmarkdown.html_vignette.check_title")
on.exit(options(rmarkdown.html_vignette.check_title = oo))
options(rmarkdown.html_vignette.check_title = FALSE)

library(MigConnectivity)

## ----message=FALSE, warning = FALSE, error=FALSE, echo = FALSE-------------------------------
set.seed(1234)
data(OVENdata)

newDir <- tempdir()
baseURL <- 'https://github.com/SMBC-NZP/MigConnectivity/blob/devpsi2/data-raw/'
file.name <- "OVENpsi.rds"
url1 <- paste0(baseURL, file.name, '?raw=true')
temp <- paste(newDir, file.name, sep = '/')
utils::download.file(url1, temp, mode = 'wb')
OVENpsi <- readRDS(temp)
unlink(temp)

OVEN_mc <- estStrength(originDist = OVENdata$originDist,
                      targetDist = OVENdata$targetDist,
                      originRelAbund = OVENdata$originRelAbund,
                      psi = OVENpsi,
                      nSamples = 5000,
                      verbose = 0)


## ----eval = FALSE----------------------------------------------------------------------------
#  # Read in the estStrength results for OVEN
#  OVEN_mc <- readRDS("OVENmc.rds")
#  
#  # Read in the estStrength results for YEWA
#  YEWA_mc <- readRDS("YEWAmc.rds")

## ----echo = FALSE----------------------------------------------------------------------------
# Read in the estStrength results for YEWA
newDir <- tempdir()
baseURL <- 'https://github.com/SMBC-NZP/MigConnectivity/blob/devpsi2/data-raw/YEWA/'
file.name <- "YEWAmc.rds"
url1 <- paste0(baseURL, file.name, '?raw=true')
temp <- paste(newDir, file.name, sep = '/')
utils::download.file(url1, temp, mode = 'wb')
YEWA_mc <- readRDS(temp)
unlink(temp)

## --------------------------------------------------------------------------------------------
diffStrength(estimates = list(Ovenbird = OVEN_mc,
                              Yellow_Warbler = YEWA_mc))


