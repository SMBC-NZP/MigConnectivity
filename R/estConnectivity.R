###############################################################################
#' Estimate MC, migratory connectivity strength
#'
#' Resampling of uncertainty for MC (migratory connectivity strength)
#' from estimates of psi (transition probabilities) and/or relative abundance.
#' Psi estimates can come from an estMigConnectivity object, an RMark psi
#' matrix, MCMC samples, or other samples expressed in array form.
#' Abundance estimates for each origin site can be
#' either just point estimates (no uncertainty propagated) or MCMC samples.
#' Other inputs include distances between origin sites, distances between target
#' sites, and sample size used to estimate psi.
#'
#'
#' @param originDist Distances between the B origin sites. Symmetric B by B
#'  matrix
#' @param targetDist Distances between the W target sites. Symmetric W by W
#'  matrix
#' @param originRelAbund Relative abundance estimates at B origin sites. Either
#'  a numeric vector of length B that sums to 1, or an mcmc object (such as is
#'  produced by \code{\link{modelCountDataJAGS}}) or matrix with at least
#'  \code{nSamples} rows. If there are more than B columns, the relevant columns
#'  should be labeled "relN[1]" through "relN[B]"
#' @param psi Transition probabilities between B origin and W target sites.
#'  Either a matrix with B rows and W columns where rows sum to 1, an array with
#'  dimensions x, B, and W (with x samples of the transition probability matrix
#'  from another model), an 'estPsi' object (result of calling estTransition),
#'  or a MARK object with estimates of transition probabilities
#' @param sampleSize Total sample size of animals that psi will be estimated
#'  from. Should be the number of animals released in one of the origin sites
#'  and observed in one of the target sites (or vice-versa). Optional, but
#'  recommended, unless psi is an estPsi object (in which case this function can
#'  pull it from there)
#' @param originSites If \code{psi} is a MARK object, this must be a numeric
#'  vector indicating which sites are origin
#' @param targetSites If \code{psi} is a MARK object, this must be a numeric
#'  vector indicating which sites are target
#' @param originNames Optional. Vector of names for the origin sites. Mostly for
#'  internal use
#' @param targetNames Optional. Vector of names for the target sites. Mostly for
#'  internal use
#' @param nSamples Number of times to resample \code{psi} and/or
#'  \code{originRelAbund}. The purpose is to estimate sampling uncertainty;
#'  higher values here will do so with more precision
#' @param row0 If \code{originRelAbund} is an mcmc object or array, this can be
#'  set to 0 (default) or any greater integer to specify where to stop ignoring
#'  samples ("burn-in")
#' @param verbose 0 (default) to 2. 0 prints no output during run. 1 prints
#'  a progress update and summary every 100 samples. 2 prints a
#'  progress update and summary every sample
#' @param alpha Level for confidence/credible intervals provided. Default (0.05)
#'  gives 95 percent CI
#' @param approxSigTest Should function compute approximate one-sided
#'  significance tests (p-values) for MC from the resampling? Default is
#'  FALSE
#' @param sigConst Value to compare MC to in significance test. Default is 0
#' @param maintainLegacyOutput version 0.4.0 of \code{MigConnectivity}
#'  updated the structure of the estimates. If you have legacy code that refers
#'  to elements within an \code{estMigConnectivity} object (results of
#'  \code{estMC}), you can set this to TRUE to also keep the old structure.
#'  Defaults to FALSE
#' @param returnAllInput if TRUE (the default) the output includes all of the
#'  inputs. If FALSE, only the inputs currently used by another MigConnectivity
#'  function are included in the output.
#'
#' @return \code{estStrength} returns a list with the elements:
#' \describe{
#'   \item{\code{MC}}{List containing estimates of migratory connectivity
#'    strength:
#'    \itemize{
#'      \item{\code{sample}} \code{nSamples} sampled values for
#'       MC. Provided to allow the user to compute own summary statistics.
#'      \item{\code{mean}} Mean of \code{MC$sample}. Main estimate of MC,
#'       incorporating parametric uncertainty.
#'      \item{\code{se}} Standard error of MC, estimated from SD of
#'       \code{MC$sample}.
#'      \item{\code{simpleCI}} Default\code{1 - alpha} confidence interval for
#'       MC, estimated as \code{alpha/2} and \code{1 - alpha/2} quantiles of
#'       \code{MC$sample}.
#'      \item{\code{bcCI}} Bias-corrected \code{1 - alpha} confidence interval
#'       for MC. May be preferable to \code{MC$simpleCI} when \code{MC$mean} is
#'       the best estimate of MC. \code{MC$simpleCI} is preferred when
#'       \code{MC$median} is a better estimator. When \code{MC$mean==MC$median},
#'       these should be identical.  Estimated as the
#'       \code{pnorm(2 * z0 + qnorm(alpha / 2))} and
#'       \code{pnorm(2 * z0 + qnorm(1 - alpha / 2))} quantiles of \code{MC$sample},
#'       where z0 is the proportion of \code{MC$sample < MC$mean}.
#'      \item{\code{hpdCI}} \code{1 - alpha} credible interval for MC,
#'       estimated using the highest posterior density (HPD) method.
#'      \item{\code{median}} Median of MC, alternate point estimate also
#'       including parametric uncertainty.
#'      \item{\code{point}} Simple point estimate of MC, using the point
#'      estimates of \code{psi} and \code{originRelAbund} (usually the mean
#'      values), not accounting for sampling error.
#'      \item{\code{simpleP}} Approximate p-value for MC, estimated as the
#'      proportion of bootstrap iterations where MC < \code{sigConst} (or MC >
#'      \code{sigConst} if \code{pointMC < sigConst}).  Note that if the
#'      proportion is 0, a default value of 0.5 / \code{nSamples} is provided,
#'      but this is best interpreted as p < 1 / \code{nSamples}.  NULL when
#'      \code{approxSigTest==FALSE}.
#'      \item{\code{bcP}} Approximate bias-corrected p-value for MC, estimated as
#'      \code{pnorm(qnorm(simpleP) - 2 * z0)}, where z0 is the proportion of
#'      \code{sampleMC < meanMC}.  May be a better approximation of the p-value
#'      than \code{simpleP}, but many of the same limitations apply.  NULL when
#'      \code{approxSigTest==FALSE}.
#'    }
#'   }
#'   \item{\code{input}}{List containing the inputs to \code{estStrength}.}
#' }
#'
#' @export
#'
#' @seealso \code{\link{calcMC}}, \code{\link{estTransition}},
#'   \code{\link{estMC}}, \code{\link{estMantel}},
#'   \code{\link{plot.estMigConnectivity}}
#' @example inst/examples/estStrengthExamples.R
estStrength <- function(originDist, targetDist, originRelAbund, psi,
                        sampleSize = NULL,
                        originSites=NULL, targetSites=NULL,
                        originNames = NULL, targetNames = NULL,
                        nSamples = 1000, row0 = 0, verbose=0,
                        alpha = 0.05, approxSigTest = FALSE, sigConst = 0,
                        maintainLegacyOutput = FALSE,
                        returnAllInput = TRUE) {
  nOriginSites <- nrow(originDist)
  nTargetSites <- nrow(targetDist)
  absAbund <- !is.null(sampleSize)
  if (coda::is.mcmc(originRelAbund) || coda::is.mcmc.list(originRelAbund)) {
    originRelAbund <- as.matrix(originRelAbund)
  }
  if (is.matrix(originRelAbund) && all(dim(originRelAbund)>1)) {
    abundFixed <- FALSE
    if (dim(originRelAbund)[2]>nOriginSites)
      abundParams <- paste('relN[', 1:nOriginSites, ']', sep='')
    else if (dim(originRelAbund)[2]==nOriginSites)
      abundParams <- 1:nOriginSites
    else
      stop('Number of origin sites must be constant between distance matrix and abundance')
    if (dim(originRelAbund)[1] >= nSamples)
      abundRows <- round(seq(from = row0 + 1, to = dim(originRelAbund)[1],
                             length.out = nSamples))
    else
      stop("You need at least nSamples rows to originRelAbund")
    originRelAbund <- as.matrix(originRelAbund[abundRows, abundParams])
    abundBase <- colMeans(originRelAbund)
  }
  else {
    abundFixed <- TRUE
    if (length(originRelAbund)!=nOriginSites)
      stop('Number of origin sites must be constant between distance matrix and abundance')
    abundBase <- originRelAbund
  }
  if (is.matrix(psi)) {
    psiFixed <- TRUE
    psiVCV <- NULL
    if (nrow(psi)!=nOriginSites || ncol(psi)!=nTargetSites)
      stop('Size of psi matrix must be consistant with distance matrices')
    psiBase <- psi
    psiIn <- psi
  }
  else if (inherits(psi, "mark")) {
    psiFixed <- FALSE
    if (!is.numeric(originSites) || !is.numeric(targetSites))
      stop('Must specify which RMark Psi parameters represent origin and target sites')
    psiVCV <- psi$results$beta.vcv
    psiBase <- RMark::TransitionMatrix(RMark::get.real(psi, "Psi",
                                                       se=TRUE))[originSites,
                                                              targetSites]
    if (any(diag(psi$results$beta.vcv) < 0))
      stop("Can't sample model, negative beta variances")
    psiIn <- psi
  }
  else if (is.array(psi)) {
    if (length(dim(psi))!=3)
      stop('Psi should either be 2-(for fixed transition probabilities) or 3-dimensional array')
    psiFixed <- FALSE
    if (dim(psi)[2]!=nOriginSites || dim(psi)[3]!=nTargetSites)
      stop('Size of psi array must be consistent with distance matrices')
    psiBase <- apply(psi, 2:3, mean)
    psiVCV <- NULL
    if (dim(psi)[1]>=nSamples)
      psiSamples <- round(seq(from = 1, to = dim(psi)[1],
                              length.out = nSamples))
    else
      psiSamples <- sample.int(dim(psi)[1], nSamples, replace = TRUE)
    psiIn <- psi
  }
  else if (inherits(psi, "estPsi") || inherits(psi, "estMC")) {
    psiFixed <- FALSE
    psiIn <- psi
    psi <- psi$psi$sample
    if (dim(psi)[2]!=nOriginSites || dim(psi)[3]!=nTargetSites)
      stop('Size of psi array must be consistant with distance matrices')
    psiBase <- apply(psi, 2:3, mean)
    psiVCV <- NULL
    if (dim(psi)[1]>=nSamples)
      psiSamples <- 1:nSamples
    else
      psiSamples <- sample.int(dim(psi)[1], nSamples, replace = TRUE)
    if (is.null(sampleSize))
      sampleSize <- psiIn$input$sampleSize # or wherever we put it
  }
  pointMC <- ifelse(absAbund,
                    calcMC(originDist, targetDist, originRelAbund = abundBase,
                           psi = psiBase, sampleSize = sampleSize),
                    calcMC(originDist, targetDist, originRelAbund = abundBase,
                           psi = psiBase))
  sampleMC <- rep(NA, nSamples)
  psi.array <- array(NA, c(nSamples, nOriginSites, nTargetSites),
                     dimnames = list(NULL, originNames, targetNames))
  for (i in 1:nSamples) {
    if (verbose > 1 || verbose == 1 && i %% 100 == 0)
      cat("\tSample", i, "of", nSamples, "at", date(), "\n")
    # Generate random transition probability matrices
    if (psiFixed)
      psiNew <- psiBase
    else if (is.null(psiVCV))
      psiNew <- psi[psiSamples[i],,]
    else
      psiNew <- makePsiRand(psi, originSites, targetSites)
    psi.array[i, , ] <- psiNew
    # Point estimates of breeding densities
    if (abundFixed)
      abundNew <- abundBase
    else
      abundNew <- originRelAbund[i, abundParams]
    # Calculate MC for new psis and relative breeding densities
    sampleMC[i] <- ifelse(absAbund,
                          calcMC(originDist, targetDist, originRelAbund = abundNew,
                                 psi = psiNew, sampleSize = sampleSize),
                          calcMC(originDist, targetDist, originRelAbund = abundNew,
                                 psi = psiNew))
    if (verbose > 1 || verbose == 1 && i %% 100 == 0)
      cat(" MC mean:", mean(sampleMC, na.rm=TRUE),
          "SD:", sd(sampleMC, na.rm=TRUE),
          "low quantile:", quantile(sampleMC, alpha/2, na.rm=TRUE),
          "high quantile:", quantile(sampleMC, 1-alpha/2, na.rm=TRUE), "\n")
  }
  meanMC <- mean(sampleMC, na.rm=TRUE)
  medianMC <- median(sampleMC, na.rm=TRUE)
  seMC <- sd(sampleMC, na.rm=TRUE)
  # Calculate confidence intervals using quantiles of sampled MC
  simpleCI <- quantile(sampleMC, c(alpha/2, 1-alpha/2), na.rm=TRUE, type = 8,
                       names = FALSE)
  z0 <- qnorm(sum((sampleMC)<meanMC)/nSamples)
  bcCI <- quantile(sampleMC, pnorm(2*z0+qnorm(c(alpha/2, 1-alpha/2))),
                       na.rm=TRUE, names = FALSE)
  MC.mcmc <- coda::as.mcmc(sampleMC) # Ha!
  hpdCI <- as.vector(coda::HPDinterval(MC.mcmc, 1-alpha))
  if (!approxSigTest)
    simpleP <- bcP <- NULL
  else {
    if (pointMC > sigConst)
      simpleP <- sum(sampleMC < sigConst) / nSamples
    else
      simpleP <- sum(sampleMC > sigConst) / nSamples
    if (simpleP == 0)
      simpleP <- 0.5 / nSamples
    bcP <- pnorm(qnorm(simpleP) - 2 * z0)
    if (pointMC < sigConst)
      bcP <- 1 - bcP
  }
  if (returnAllInput) {
    input <- list(originDist = originDist, targetDist = targetDist,
                  originRelAbund = originRelAbund, psi = psiIn,
                  sampleSize = sampleSize, originSites = originSites,
                  targetSites = targetSites,
                  originNames = originNames,
                  targetNames = targetNames,
                  nSamples = nSamples, row0 = row0,
                  verbose = verbose, alpha = alpha,
                  approxSigTest = approxSigTest, sigConst = sigConst,
                  maintainLegacyOutput = maintainLegacyOutput,
                  returnAllInput = TRUE)
  }
  else {
    input <- list(sampleSize = sampleSize,
                  originNames = originNames, targetNames = targetNames,
                  alpha = alpha, returnAllInput = FALSE)
  }
  if (maintainLegacyOutput) {
    mc <- list(sampleMC=sampleMC, samplePsi = psi.array, pointPsi = psiBase,
                pointMC=pointMC, meanMC=meanMC,
                medianMC=medianMC, seMC=seMC, simpleCI=simpleCI,
                bcCI=bcCI, hpdCI=hpdCI, simpleP = simpleP, bcP = bcP,
                sampleCorr = NULL, pointCorr = NULL,
                meanCorr = NULL, medianCorr = NULL, seCorr=NULL,
                simpleCICorr=NULL, bcCICorr=NULL, inputSampleSize = sampleSize,
                alpha = alpha, sigConst = sigConst,
                MC = list(sample = sampleMC, mean = meanMC, se = seMC,
                          simpleCI = simpleCI, bcCI = bcCI, hpdCI = hpdCI,
                          median = medianMC, point = pointMC,
                          simpleP = simpleP, bcP = bcP),
                input = input)
  }
  else {
    mc <- list(MC = list(sample = sampleMC, mean = meanMC, se = seMC,
                          simpleCI = simpleCI, bcCI = bcCI, hpdCI = hpdCI,
                          median = medianMC, point = pointMC,
                          simpleP = simpleP, bcP = bcP),
               input = input)
  }
  class(mc) <- c("estMC", "estMigConnectivity")
  return(mc)
}


estTransitionJAGS <- function (banded, reencountered,
                               originAssignment = NULL, targetAssignment = NULL,
                               alpha = 0.05, nSamples = 1000, verbose=0,
                               originNames = NULL, targetNames = NULL,
                               nThin = 1, nBurnin = 5000, nChains = 3,
                               fixedZero = NULL, psiPrior = NULL,
                               returnAllInput = TRUE) {
  jags.inits <- vector("list", nChains)
  if (is.null(banded)) {
    if (is.null(originAssignment) || is.null(targetAssignment)) {
      stop("If running estTransition without bootstrap, need to provide banding and/or telemetry (through originAssignment and targetAssignment) data")
    }
    nDim <- 0
    nTargetSites <- max(length(unique(targetAssignment)),
                        length(targetNames))
    nOriginSites <- max(length(unique(originAssignment)),
                        length(originNames))
    if (is.null(originNames))
      originNames <- LETTERS[1:nOriginSites]
    if (is.null(targetNames))
      targetNames <- as.character(1:nTargetSites)

    jags.data <- list(npop = nOriginSites, ndest = nTargetSites)
    sampleSize <- 0
    # Parameters to monitor
    params <- c("psi")
    for (i in 1:nChains)
      jags.inits[[i]] <- list(m0 = matrix(runif(nOriginSites * nTargetSites),
                                          nOriginSites, nTargetSites))
  }
  else {
    nDim <- length(dim(reencountered))
    if (is.null(originNames))
      originNames <- dimnames(reencountered)[[1]]
    if (is.null(targetNames))
      targetNames <- dimnames(reencountered)[[nDim]]
    nTargetSites <- dim(reencountered)[nDim]
    nOriginSites <- dim(reencountered)[1]
    if (nDim == 2) {
      nAges <- 1
      if (length(banded)!=nOriginSites)
        stop('Number of origin sites must be constant between reencountered and banded')
      nfound <- apply(reencountered, 1, sum)
      sampleSize <- sum(nfound)
      tmat <- cbind(reencountered, banded - nfound)
      # Data
      jags.data <- list(recmat = tmat, npop = nOriginSites,
                        ndest = nTargetSites, nreleased = banded)
      for (i in 1:nChains)
        jags.inits[[i]] <- list(m0 =  matrix(runif(nOriginSites * nTargetSites),
                                             nOriginSites, nTargetSites),
                                r = rbeta(nTargetSites, 1, 1))
    }
    else {
      nAges <- dim(banded)[2]
      if (dim(banded)[1]!=nOriginSites)
        stop('Number of origin sites must be consistant between reencountered and banded')
      if (dim(reencountered)[2]!=nAges)
        stop('Number of ages must be consistant between banded and reencountered')
      nfound <- apply(reencountered, 1:2, sum)
      sampleSize <- sum(nfound)
      tmat <- array(NA, c(nOriginSites, nAges, nTargetSites + 1))
      tmat[ , , 1:nTargetSites] <- reencountered
      tmat[ , , 1 + nTargetSites] <- banded - nfound
      # Data
      jags.data <- list(recmat = tmat, npop = nOriginSites, nages = nAges,
                        ndest = nTargetSites, nreleased = banded)
      for (i in 1:nChains)
        jags.inits[[i]] <- list(m0 = matrix(runif(nOriginSites * nTargetSites),
                                            nOriginSites, nTargetSites),
                                r = matrix(rbeta(nTargetSites * nAges, 1, 1),
                                           nAges, nTargetSites))
    }
    # Parameters to monitor
    params <- c("psi", "r")
  }
  if (is.null(psiPrior)) {
    psiPrior <- matrix(1, nOriginSites, nTargetSites)
  }
  jags.data$psiPrior <- psiPrior
  if (!is.null(originAssignment)) {
    telmat <- table(factor(originAssignment, levels = 1:nOriginSites),
                    factor(targetAssignment, levels = 1:nTargetSites))
    jags.data$telmat <- telmat
    jags.data$ntel <- as.vector(table(factor(originAssignment,
                                             levels = 1:nOriginSites)))
    sampleSize <- sampleSize + sum(jags.data$ntel)
  }

  if (!is.null(fixedZero)) {
    psiFixed <- matrix(NA, nOriginSites, nTargetSites)
    for (i in 1:nrow(fixedZero)) {
      psiFixed[fixedZero[i, 1], fixedZero[i, 2]] <- 0
    }
    jags.data$m0 <- psiFixed
  }
  file <- tempfile(fileext = ".txt")
  sink(file = file)
  cat("
model{
  # psi prior
  for(i in 1:npop){
    for(k in 1:ndest){
      m0[i, k] ~ dbeta(psiPrior[i, k], 1)
      psi[i, k] <- m0[i, k] / sum(m0[i, 1:ndest])
    } #k
  }#i")
  if (!is.null(banded)){
    if (nAges==1) {
      cat("
  # model for recoveries with known number of ringed
  for (i in 1:npop){
    for(k in 1:ndest){
      p[i, k] <- psi[i, k] * r[k]
    }
    p[i, (ndest+1)] <- 1 - sum(p[i, 1:ndest])
  }

  for(k in 1:ndest) {
    r[k] ~ dunif(0, 1)
  }
  for(i in 1:npop){
    recmat[i, 1:(ndest+1)] ~ dmulti(p[i, 1:(ndest+1)], nreleased[i])
  }")
    }
    else {
      cat("
  for (j in 1:nages) {
        for(k in 1:ndest) {
          r[j, k] ~ dunif(0, 1)
        }
  }
  for (i in 1:npop){
    for (j in 1:nages) {
      for(k in 1:ndest){
        p[i, j, k] <- psi[i, k] * r[j, k]
      }
      p[i, j, (ndest+1)] <- 1 - sum(p[i, j, 1:ndest])
    }
  }
  for(i in 1:npop){
    for (j in 1:nages) {
      recmat[i, j, 1:(ndest+1)] ~ dmulti(p[i, j, 1:(ndest+1)], nreleased[i, j])
    }
  }")
    }
  }
  if (!is.null(originAssignment)) {
    cat("
  for(i in 1:npop){
    telmat[i, 1:ndest] ~ dmulti(psi[i, 1:ndest], ntel[i])
  }
")
  }
  cat("}")
  sink()
  # print(nAges)
  # print(file)
  # print(jags.data)
  # print(rowSums(jags.data$recmat))
  # print(jags.inits)
  out <- R2jags::jags(data = jags.data, inits = jags.inits, params,
                      file,
                      n.chains = nChains, n.thin = nThin,
                      n.iter = nBurnin + ceiling(nSamples * nThin / nChains),
                      n.burnin = nBurnin, DIC = FALSE,
                      progress.bar = ifelse(verbose==0, 'none', 'text'))

  if (nChains > 1) {
    maxRhat <- max(out$BUGSoutput$summary[, "Rhat"], na.rm = TRUE)
    if (maxRhat < 1.1)
      cat("Successful convergence based on Rhat values (all < 1.1).\n")
    else
      cat("**WARNING** Rhat values indicate convergence failure.\n")
  }
  file.remove(file)
  psi <- out$BUGSoutput$sims.list$psi
  dimnames(psi) <- list(NULL, originNames, targetNames)
  meanPsi <- out$BUGSoutput$mean$psi
  simpleCIPsi <- apply(psi, 2:3, quantile, probs = c(alpha/2, 1-alpha/2),
                       na.rm=TRUE, names = FALSE)
  psi.matrix <- array(c(psi), c(dim(psi)[[1]], nOriginSites * nTargetSites),
                      list(NULL, paste(rep(originNames, nTargetSites),
                                       rep(targetNames, each = nOriginSites),
                                       sep = "#")))
  psi.mcmc <- coda::as.mcmc(psi.matrix)
  hpdCI <- coda::HPDinterval(psi.mcmc, 1-alpha)
  hpdCI <- array(hpdCI, c(nOriginSites, nTargetSites, 2),
                 list(originNames, targetNames, c("lower", "upper")))
  hpdCI <- aperm(hpdCI, c(3, 1, 2))
  bcCIPsi <- array(NA, dim = c(2, nOriginSites, nTargetSites),
                   dimnames = list(NULL, originNames, targetNames))
  for (i in 1:nOriginSites) {
    for (j in 1:nTargetSites) {
      psi.z0 <- qnorm(sum(psi[, i, j] < meanPsi[i, j], na.rm = TRUE) /
                        length(which(!is.na(psi[, i, j]))))
      bcCIPsi[ , i, j] <- quantile(psi[, i, j],
                                   pnorm(2 * psi.z0 + qnorm(c(alpha/2, 1-alpha/2))),
                                   na.rm=TRUE, names = FALSE)
    }
  }
  results <- list(psi = list(sample = psi, mean = meanPsi,
                             se = out$BUGSoutput$sd$psi,
                             simpleCI = simpleCIPsi, bcCI = bcCIPsi,
                             hpdCI = hpdCI,
                             median = out$BUGSoutput$median$psi))
  if (!is.null(out$BUGSoutput$sims.list$r)) {
    if (nAges == 1) {
      bcCIr <- array(NA, dim = c(2, nTargetSites),
                     dimnames = list(NULL, targetNames))
      for (j in 1:nTargetSites) {
        psi.z0 <- qnorm(sum(out$BUGSoutput$sims.list$r[ ,j] <
                              out$BUGSoutput$mean$r[j], na.rm = TRUE) /
                          length(which(!is.na(out$BUGSoutput$sims.list$r[ ,j]))))
        bcCIr[ , j] <- quantile(out$BUGSoutput$sims.list$r[ ,j],
                                pnorm(2 * psi.z0 + qnorm(c(alpha/2, 1-alpha/2))),
                                na.rm=TRUE, names = FALSE)
      }
      colnames(out$BUGSoutput$sims.list$r) <- names(out$BUGSoutput$mean$r) <-
        names(out$BUGSoutput$sd$r) <- names(out$BUGSoutput$median$r) <-
        targetNames
      results$r <- list(sample = out$BUGSoutput$sims.list$r,
                       mean = out$BUGSoutput$mean$r,
                       se = out$BUGSoutput$sd$r,
                       simpleCI = apply(out$BUGSoutput$sims.list$r, 2, quantile,
                                        probs = c(alpha/2, 1-alpha/2), type = 8),
                       bcCI = bcCIr, median = out$BUGSoutput$median$r)
    }
    else {
      bcCIr <- array(NA, dim = c(2, nAges, nTargetSites),
                     dimnames = list(NULL, NULL, targetNames))
      for (i in 1:nAges) {
        for (j in 1:nTargetSites) {
          psi.z0 <- qnorm(sum(out$BUGSoutput$sims.list$r[ , i, j] <
                              out$BUGSoutput$mean$r[i, j], na.rm = TRUE) /
                          length(which(!is.na(out$BUGSoutput$sims.list$r[,i,j]))))
          bcCIr[,i,j] <- quantile(out$BUGSoutput$sims.list$r[ , i, j],
                                  pnorm(2 * psi.z0 + qnorm(c(alpha/2,1-alpha/2))),
                                  na.rm=TRUE, names = FALSE)
        }
      }
      dimnames(out$BUGSoutput$sims.list$r)[[3]]<-colnames(out$BUGSoutput$mean$r) <-
        colnames(out$BUGSoutput$sd$r) <- colnames(out$BUGSoutput$median$r) <-
        targetNames
      results$r <- list(sample = out$BUGSoutput$sims.list$r,
                       mean = out$BUGSoutput$mean$r,
                       se = out$BUGSoutput$sd$r,
                       simpleCI = apply(out$BUGSoutput$sims.list$r, 2:3, quantile,
                                        probs = c(alpha/2, 1-alpha/2), type = 8),
                       bcCI = bcCIr, median = out$BUGSoutput$median$r)
    }
  }
  if (returnAllInput) {
    results$input <- list(banded = banded, reencountered = reencountered,
                          originAssignment = originAssignment,
                          targetAssignment = targetAssignment,
                          sampleSize = sampleSize, alpha = alpha,
                          nSamples = nSamples, verbose=verbose,
                          originNames = originNames, targetNames = targetNames,
                          nThin = nThin, nBurnin = nBurnin, nChains = nChains,
                          fixedZero = fixedZero, psiPrior = psiPrior,
                          method = "MCMC", returnAllInput = TRUE)
  }
  else {
    results$input <- list(sampleSize = sampleSize, alpha = alpha,
                          originNames = originNames, targetNames = targetNames,
                          method = "MCMC", returnAllInput = FALSE)
  }
  results$BUGSoutput <- out$BUGSoutput
  return(results)
}

estTransitionBoot <- function(originSites = NULL,
                              targetSites = NULL,
                              originPoints = NULL,
                              targetPoints = NULL,
                              originAssignment = NULL,
                              targetAssignment = NULL,
                              originNames = NULL,
                              targetNames = NULL,
                              nBoot = 1000,
                              isGL = FALSE,
                              isTelemetry = FALSE,
                              isRaster = FALSE,
                              isProb = FALSE,
                              captured = "origin",
                              geoBias = NULL,
                              geoVCov = NULL,
                              geoBiasOrigin = geoBias,
                              geoVCovOrigin = geoVCov,
                              targetRaster = NULL,
                              originRaster = NULL,
                              verbose = 0,
                              alpha = 0.05,
                              resampleProjection = 'ESRI:102010',#MigConnectivity::projections$EquidistConic,
                              nSim = ifelse(any(isRaster), 10, 1000),
                              maxTries = 300,
                              dataOverlapSetting = "dummy",
                              fixedZero = NULL,
                              targetRelAbund = NULL,
                              banded = NULL,
                              reencountered = NULL,
                              method = "bootstrap",
                              m = NULL,
                              returnAllInput = TRUE) {
  # Input checking and assignment
  if (any(captured != "origin" & captured != "target" & captured != "neither")){
    stop("captured should be 'origin', 'target', 'neither', or a vector of those options")}
  if (!(verbose %in% 0:3)){
    stop("verbose should be integer 0-3 for level of output during bootstrap: 0 = none, 1 = every 10, 2 = every run, 3 = number of draws")}
  if (length(geoBias)!=2 && any(isGL & (captured == "origin" | captured == "neither"))){
    stop("geoBias should be vector of length 2 (expected bias in longitude and latitude of targetPoints, in resampleProjection units, default meters)")}
  if (!isTRUE(all.equal(dim(geoVCov), c(2, 2), check.attributes = FALSE)) &&
      any(isGL & (captured == "origin" | captured == "neither"))){
    stop("geoVCov should be 2x2 matrix (expected variance/covariance in longitude and latitude of targetPoints, in resampleProjection units, default meters)")}
  if ((is.null(originPoints) && is.null(originRaster) && is.null(originSites)) &&
      is.null(originAssignment) && is.null(banded)){
    stop("Need to define either originAssignment, originSites, originRaster, originPoints, or banded")}
  if ((is.null(targetPoints) && is.null(targetRaster) &&
       is.null(targetSites)) && is.null(targetAssignment) && is.null(reencountered)){
    stop("Need to define either targetAssignment, targetSites, targetRaster, targetPoints, or reencountered")}
  if ((is.null(banded) && !is.null(reencountered) ||
       !is.null(banded)) && is.null(reencountered)){
    stop("Need to define both banded and reencountered")}
  if(inherits(originSites,"SpatialPolygonsDataFrame")){
    originSites <- sf::st_as_sf(originSites)}
  if(inherits(targetSites,"SpatialPolygonsDataFrame")){
    targetSites <- sf::st_as_sf(targetSites)}

  targetStats <- assignRasterStats(targetRaster)
  targetPointsAssigned <- targetStats$PointsAssigned
  targetSingleCell <- targetStats$SingleCell
  targetRasterXYZ <- targetStats$RasterXYZ
  targetRasterXYZcrs <- targetStats$RasterXYZcrs
  targetRaster <- targetStats$Raster

  originStats <- assignRasterStats(originRaster)
  originPointsAssigned <- originStats$PointsAssigned
  originSingleCell <- originStats$SingleCell
  originRasterXYZ <- originStats$RasterXYZ
  originRasterXYZcrs <- originStats$RasterXYZcrs
  originRaster <- originStats$Raster

  if (dataOverlapSetting != "dummy") {
    if (verbose > 0)
      cat("Configuring data overlap settings\n")
    temp <- reassignInds(dataOverlapSetting = dataOverlapSetting,
                         originPoints = originPoints,
                         targetPoints = targetPoints,
                         originAssignment = originAssignment,
                         targetAssignment = targetAssignment,
                         isGL = isGL, isTelemetry = isTelemetry,
                         isRaster = isRaster, isProb = isProb,
                         captured = captured,
                         originRasterXYZ = originRasterXYZ,
                         originSingleCell = originSingleCell,
                         targetRasterXYZ = targetRasterXYZ,
                         targetSingleCell = targetSingleCell,
                         targetSites = targetSites, originSites = originSites)
    originPoints <- temp$originPoints; targetPoints <- temp$targetPoints
    originAssignment <- temp$originAssignment
    targetAssignment <- temp$targetAssignment
    isGL <- temp$isGL; isTelemetry <- temp$isTelemetry
    isRaster <- temp$isRaster; isProb <- temp$isProb
    originRasterXYZ <- temp$originRasterXYZ
    if (!is.null(originRasterXYZ)){
      colnames(originRasterXYZ) <- c("x", "y", paste0("lyr.", 1:length(isRaster)))
      originRaster <- terra::rast(originRasterXYZ, crs = originRasterXYZcrs,
                                  extent = terra::ext(originRaster), type = "xyz")
    }
    originSingleCell <- temp$originSingleCell
    targetRasterXYZ <- temp$targetRasterXYZ
    if (!is.null(targetRasterXYZ)){
      colnames(targetRasterXYZ) <- c("x", "y", paste0("lyr.", 1:length(isRaster)))
      targetRaster <- terra::rast(targetRasterXYZ, crs = targetRasterXYZcrs,
                                  extent = terra::ext(targetRaster), type = "xyz")
    }
    targetSingleCell <- temp$targetSingleCell
  }
  if (any(isProb & (captured != "target")) && (is.null(targetAssignment) || length(dim(targetAssignment))!=2)){
    stop("With probability assignment (isProb==TRUE) animals captured at origin, targetAssignment must be a [number of animals] by [number of target sites] matrix")}
  if (any(isProb & captured != "origin") && (is.null(originAssignment) || length(dim(originAssignment))!=2)){
    stop("With probability assignment (isProb==TRUE) animals captured at target, originAssignment must be a [number of animals] by [number of origin sites] matrix")}

  if (is.null(targetPoints) && is.null(originPoints) &&
      is.null(targetAssignment) && is.null(originAssignment) &&
      is.null(targetRaster) && is.null(originRaster))
    nAnimals <- 0
  else
    nAnimals <- max(nrow(targetPoints), nrow(originPoints), length(isGL),
                  length(isTelemetry), length(isRaster), length(isProb),
                  min(length(targetAssignment), dim(targetAssignment)[1]),
                  min(length(originAssignment), dim(originAssignment)[1]),
                  ifelse(is.null(targetRaster), 0,
                         ifelse(targetPointsAssigned, dim(targetSingleCell)[3],
                         dim(targetRasterXYZ)[2] - 2)),
                  ifelse(is.null(originRaster), 0,
                         ifelse(originPointsAssigned, dim(originSingleCell)[3],
                                dim(originRasterXYZ)[2] - 2)),
                  length(captured))
  nAnimalsTotal <- nAnimals + sum(banded) #+ sum(reencountered)#
  isCMR <- c(rep(FALSE, nAnimals), rep(TRUE, nAnimalsTotal - nAnimals))
  if (length(isGL)==1){
    isGL <- c(rep(isGL, nAnimals), rep(FALSE, nAnimalsTotal - nAnimals))
  }
  else
    isGL <- c(isGL, rep(FALSE, nAnimalsTotal - nAnimals))
  if (length(isTelemetry)==1){
    isTelemetry <- c(rep(isTelemetry, nAnimals),
                     rep(FALSE, nAnimalsTotal - nAnimals))
  }
  else
    isTelemetry <- c(isTelemetry, rep(FALSE, nAnimalsTotal - nAnimals))
  if (length(isRaster)==1){
    isRaster <- c(rep(isRaster, nAnimals), rep(FALSE, nAnimalsTotal - nAnimals))
  }
  else
    isRaster <- c(isRaster, rep(FALSE, nAnimalsTotal - nAnimals))
  if (length(isProb)==1){
    isProb <- c(rep(isProb, nAnimals), rep(FALSE, nAnimalsTotal - nAnimals))
  }
  else
    isProb <- c(isProb, rep(FALSE, nAnimalsTotal - nAnimals))
  if (length(captured)==1){captured <- rep(captured, nAnimals)}

  isCMR <- c(rep(FALSE, nAnimals), rep(TRUE, nAnimalsTotal - nAnimals))
  if (!is.null(banded)) {
    captured <- c(captured, rep("origin", nAnimalsTotal - nAnimals)) #sum(banded)
  }
  if (nAnimals > 0)
    if (any(!isGL[1:nAnimals] & !isTelemetry[1:nAnimals] & !isRaster[1:nAnimals] &
            !isProb[1:nAnimals]))
      stop("For each individual animal (not in banded) one of the following must be set to TRUE:
           isGL, isTelemetry, isRaster, or isProb")
  if (method=="m-out-of-n-bootstrap" && is.null(m))
    m <- ceiling(nAnimalsTotal / 4) # don't know if this is a good default or not!
  else if (method == "bootstrap")
    m <- nAnimalsTotal

  if (any(isRaster & !isProb)) {
    conversion <- rasterToProb(originSites = originSites,
                               targetSites = targetSites,
                               originAssignment = originAssignment,
                               targetAssignment = targetAssignment,
                               isRaster = isRaster, isProb = isProb,
                               captured = captured, targetRaster = targetRaster,
                               originRaster = originRaster)
    if (!conversion$failTarget)
      targetAssignment <- conversion$targetAssignment
    if (!conversion$failOrigin)
      originAssignment <- conversion$originAssignment
    isRaster <- conversion$isRaster
    isProb <- conversion$isProb
  }

  if (!is.null(originPoints) && !is.null(originSites))
    if(!identical(sf::st_crs(originPoints), sf::st_crs(originSites)))
      # project if needed
      originPoints <- sf::st_transform(originPoints, sf::st_crs(originSites))
  if (!is.null(targetPoints) && !is.null(targetSites))
    if(!identical(sf::st_crs(targetPoints), sf::st_crs(targetSites)))
      # project if needed
      targetPoints <- sf::st_transform(targetPoints, sf::st_crs(targetSites))

  # IF originAssignment is NULL - we need to generate originAssignments from
  # the data provided
  if (is.null(originAssignment)){
    if (verbose > 0)
      cat("Creating originAssignment\n")
    # if geolocator, telemetry and captured in origin then simply get the origin site
    if (all(isGL | isTelemetry | captured != "target") && !is.null(originPoints)){
      originAssignment <- suppressMessages(unclass(sf::st_intersects(x = originPoints,
                                                          y = originSites,
                                                          sparse = TRUE)))
    # if raster and not captured in origin sites then determine the origin site
    }
    else if (all(isRaster & captured != "origin")) {
      # if isRaster == TRUE and captured != origin
      # WEIGHTED XY COORDIANTES FROM THE RASTER
      # get geographically weighted median value
      xyOriginRast <- apply(originRasterXYZ[,3:ncol(originRasterXYZ)],
                            MARGIN = 2,
                            FUN = function(x){
                        # xy <-cbind(weighted.mean(originRasterXYZ[,1], w = x, na.rm = TRUE),
                        #           weighted.mean(originRasterXYZ[,2], w = x, na.rm = TRUE))
                  #select the cell with the highest posterior probability #
                              xy <- cbind(originRasterXYZ[which.max(x)[1],1],
                                          originRasterXYZ[which.max(x)[1],2])
                        return(xy)})
      # returns a point estimate for each bird - turn it into a sf object
      xyOriginRast <- t(xyOriginRast)
      colnames(xyOriginRast) <- c("x","y")
      # right now the assignment CRS is WGS84 - should be the same as the origin raster

      #cat("--originRasterXYZcrs -- \n")
      originAssignRast <- sf::st_as_sf(data.frame(xyOriginRast),
                                       coords = c("x","y"),
                                       crs = originRasterXYZcrs)
                                       # crs = sf::st_crs(originRasterXYZcrs))
                                      # crs = 4326)
      # transform to match originSites
      originAssignRast <- sf::st_transform(originAssignRast, sf::st_crs(originSites))
      originAssignment <- suppressMessages(unclass(sf::st_intersects(x = originAssignRast,
                                            y = originSites,
                                            sparse = TRUE)))
    }   # originAssignment <- what # need point assignment for raster (mean location?)
    else if (!is.null(originPoints))
      # originAssignment <- what # points over where we have them, raster assignment otherwise
      originAssignment <- suppressMessages(unclass(sf::st_intersects(x = originPoints,
                                                                 y = originSites,
                                                                sparse = TRUE)))
    else
      originAssignment <- NULL
    if (!is.null(originAssignment)) {
      originAssignment[lengths(originAssignment)==0] <- NA
      if (any(lengths(originAssignment)>1)){
        warning("Overlapping originSites may cause issues\n")
        originAssignment <- lapply(originAssignment, function (x) x[1])
      }
      originAssignment <- array(unlist(originAssignment))
      duds <- is.na(originAssignment) & captured[1:nAnimals] == "origin"
      if (any(duds)){
        if (verbose > 0)
          cat("Not all origin capture locations are within originSites. Assigning to closest site\n")
        warning("Not all origin capture locations are within originSites. Assigning to closest site.\n",
                "Affects animals: ", paste(which(duds), collapse = ","))
        originAssignment[duds] <-
          sf::st_nearest_feature(x = originPoints[duds,],
                                 y = originSites)

      }
    }
    if (!is.null(reencountered)) {
      originAssignment <- array(c(originAssignment,
                                  rep(1:length(banded), banded)))
    }
  }
  else if (!is.null(reencountered)) {
    if (is.array(originAssignment) && length(dim(originAssignment))>1){
      originAssignment <- rbind(originAssignment,
                                array(0, c(nAnimalsTotal - nAnimals,
                                           dim(originAssignment)[2])))
      place <- nAnimals
      for (i in 1:length(banded)) {
        originAssignment[place + 1:banded[i], i] <- 1
        place <- place + banded[i]
      }
    }
    else {
      nOriginSites <- length(banded)
      originAssignment <- array(c(originAssignment,
                                rep(1:nOriginSites, banded)))
    }
  }

  if (is.null(targetAssignment)){
    if (verbose > 0)
      cat("Creating targetAssignment\n")
    if (all(isGL | isTelemetry | captured != "origin")) {
      targetAssignment <- suppressMessages(unclass(sf::st_intersects(x = targetPoints,
                                                    y = targetSites,
                                                    sparse = TRUE)))
    }
    else if (all(isRaster & captured != "target")){
      #targetAssignment <- what # need point assignment for raster (mean location?)
      xyTargetRast <- apply(targetRasterXYZ[,3:ncol(targetRasterXYZ)],
                          MARGIN = 2,
                          FUN = function(x){
                            #xy <-cbind(Hmisc::wtd.quantile(targetRasterXYZ[,1],probs = 0.5, weight = x, na.rm = TRUE),
                            #           Hmisc::wtd.quantile(targetRasterXYZ[,2],probs = 0.5, weight = x, na.rm = TRUE))
                            # xy <- cbind(weighted.mean(targetRasterXYZ[,1], w = x, na.rm = TRUE),
                            #            weighted.mean(targetRasterXYZ[,2], w = x, na.rm = TRUE))
                  # Select cell with the maximum posterior probability #
                            xy <- cbind(targetRasterXYZ[which.max(x)[1],1],
                                        targetRasterXYZ[which.max(x)[1],2])
                            return(xy)})
      # returns a point estimate for each bird - turn it into a sf object
      xyTargetRast <- t(xyTargetRast)
      colnames(xyTargetRast) <- c("x","y")
      # right now the assignment CRS is WGS84 - should be the same as the origin raster
      targetAssignRast <- sf::st_as_sf(data.frame(xyTargetRast),
                                       coords = c("x","y"),
                                       crs = targetRasterXYZcrs)
                                      # crs = sf::st_crs(targetRasterXYZcrs))
                                       #crs = 4326)
      # transform to match originSites
      targetSites_wgs <- sf::st_transform(targetSites, 4326)
      #targetAssignRast <- sf::st_transform(targetAssignRast, sf::st_crs(targetSites))
      #targetAssignment <- suppressMessages(array(unlist(unclass(sf::st_intersects(x = targetAssignRast,
      #                                                           y = targetSites_wgs,
      #                                                           sparse = TRUE)))))
      # Check which points are in target sites
      targetAssignment <- suppressMessages(unclass(sf::st_intersects(x = targetAssignRast,
                                                      y = targetSites_wgs,
                                                      sparse = TRUE)))

    # NEED TO ADD A CATCH HERE TO ASSIGN THE MAX_prob to CLOSEST TARGET REGION
    # IF INTERSECTS IS (empty)
    }
    else if (!is.null(targetPoints))
   #   targetAssignment <- what # points over where we have them, raster assignment otherwise
      targetAssignment <- suppressMessages(unclass(sf::st_intersects(x = targetPoints,
                                                        y = targetSites,
                                                        sparse = TRUE)))
    else
      targetAssignment <- NULL
   if (!is.null(targetAssignment)) {
     targetAssignment[lengths(targetAssignment)==0] <- NA
     if (any(lengths(targetAssignment)>1)){
       warning("Overlapping targetSites may cause issues\n")
       targetAssignment <- lapply(targetAssignment, function(x) x[1])
     }
     targetAssignment <- array(unlist(targetAssignment))
     duds <- is.na(targetAssignment) & captured[1:nAnimals] == "target"
     if (any(duds)){
       if (verbose > 0)
         cat("Not all target capture locations are within targetSites. Assigning to closest site\n")
       warning("Not all target capture locations are within targetSites. Assigning to closest site.\n",
               "Affects animals: ", paste(which(duds), collapse = ","))
       targetAssignment[duds] <-
         sf::st_nearest_feature(x = targetPoints[duds,],
                                y = targetSites)

     }
   }
   if (!is.null(reencountered)) {
     nTargetSites <- dim(reencountered)[2]
     nOriginSites <- dim(reencountered)[1]
     for (j in 1:nOriginSites)
       targetAssignment <- array(c(targetAssignment,
                                   rep(c(1:nTargetSites, NA),
                                       c(reencountered[j, ],
                                         banded[j] - sum(reencountered[j, ])))))
   }
  }
  else if (!is.null(reencountered)) {
    if (is.array(targetAssignment) && length(dim(targetAssignment))>1){
      nTargetSites <- dim(reencountered)[2]
      nOriginSites <- dim(reencountered)[1]
      targetAssignment <- rbind(targetAssignment,
                                array(0,
                                      c(nAnimalsTotal - nAnimals, nTargetSites)))
      place <- nAnimals
      for (j in 1:nOriginSites) {
        for (i in 1:nTargetSites) {
          targetAssignment[place + 1:reencountered[j, i], i] <- 1
          place <- place + reencountered[j, i]
        }
        targetAssignment[place + 1:(banded[j] - sum(reencountered[j, ])), ] <- NA
        place <- place + banded[j] - sum(reencountered[j, ])
      }
    }
    else {
      nTargetSites <- dim(reencountered)[2]
      nOriginSites <- dim(reencountered)[1]
      for (j in 1:nOriginSites)
        targetAssignment <- array(c(targetAssignment,
                                    rep(c(1:nTargetSites, NA),
                                        c(reencountered[j, ],
                                          banded[j] - sum(reencountered[j, ])))))
    }
  }

  if(is.null(originSites)){
    if (is.array(originAssignment) && length(dim(originAssignment))>1){
      nOriginSites <- ncol(originAssignment)
    }
    else {
      nOriginSites <- length(unique(originAssignment))
    }
  }
  else{
    nOriginSites <- nrow(originSites)
  }

  if(is.null(targetSites)){
    if (is.array(targetAssignment) && length(dim(targetAssignment))>1){
      nTargetSites <- ncol(targetAssignment)
    }
    else {
      nTargetSites <- max(length(unique(targetAssignment[!is.na(targetAssignment)])),
                          length(targetNames),
                          dim(reencountered)[2])
    }
  }
  else {
    nTargetSites <- nrow(targetSites)
  }
  # if (length(targetPoints)!=nAnimals &&
  #     dim(targetAssignment)[length(dim(targetAssignment))]!=nAnimals ||
  #     nrow(originAssignment)!=nAnimals)
  #   stop("isGL should be the same length as originAssignment/originPoints and targetPoints/targetAssignment (number of animals)")
  # if (any(is.na(originAssignment)))
  #   stop("NAs in origin sites (make sure all points fall within polygons)")
  if(!is.null(originPoints) && is.na(sf::st_crs(originPoints)))
    stop('Coordinate system definition needed for originPoints')
  if(!is.null(originSites) && is.na(sf::st_crs(originSites)))
    stop('Coordinate system definition needed for originSites')
  if(!is.null(targetPoints) && is.na(sf::st_crs(targetPoints))){
    stop('Coordinate system definition needed for targetPoints')
  }
  if(!is.null(targetSites) && is.na(sf::st_crs(targetSites))){
    stop('Coordinate system definition needed for targetSites')
  }
  if(!is.null(targetPoints)){
    targetPoints <- sf::st_transform(targetPoints, crs = resampleProjection)}
  if(!is.null(targetSites)){
    targetSites <- sf::st_transform(targetSites, crs = resampleProjection)}
  if(!is.null(originPoints)){
    originPoints <- sf::st_transform(originPoints, crs = resampleProjection)}
  if(!is.null(originSites)){
    originSites <- sf::st_transform(originSites, crs = resampleProjection)}
  # Names - this isn't going to work as currently written for sf objects #
  if (is.null(targetNames)){
    # if (is.null(targetSites[[1]]) || anyDuplicated(targetSites[[1]])){
    #   targetNames <- as.character(1:nTargetSites)}else{
    targetNames <- as.character(1:nTargetSites)
  }

  if (is.null(originNames)){
#    if (is.null(originSites[[1]]) || anyDuplicated(originSites[[1]])){
    if (nOriginSites > 26){
      originNames <- as.character(1:nOriginSites)
    }else{
      originNames <- LETTERS[1:nOriginSites]}
    # }else{
    #   originNames <- 1:nOriginSites
    # }
  }

  targetPointsInSites <- FALSE

  if (targetPointsAssigned && !is.null(targetSites) && any(isRaster)) {
    if (verbose > 0){
      cat('Checking if single cell target points in targetSites, may take a moment\n')}
    targetPointSample2 <- apply(targetSingleCell,
                                FUN = function(x){sf::st_as_sf(data.frame(x),
                                                               coords = c("Longitude", "Latitude"),
                                                               crs = 4326)},
                                MARGIN = 3)
    if(!sf::st_crs(targetSites)==sf::st_crs(targetPointSample2[[1]])){
      targetPointSample2 <- sapply(targetPointSample2, sf::st_transform, crs = resampleProjection)
    }

    targetCon <- sapply(targetPointSample2, FUN = function(z){
      suppressMessages(as.numeric(unclass(sf::st_intersects(x = z, y = targetSites,
                                           sparse = TRUE))))})

    if (!any(is.na(targetCon)))
      targetPointsInSites <- TRUE
    else if (verbose > 0)
      cat('Single cell target points supplied, but some points (proportion',
          format(sum(is.na(targetCon))/length(targetCon), digits = 2), ') not in targetSites\n')
  }else {targetCon <- NULL}

  originPointsInSites <- FALSE
  if (originPointsAssigned && !is.null(originSites) && any(isRaster)) {
    if (verbose > 0)
      cat('Checking if single cell origin points in originSites, may take a moment\n')
    nSamples <- dim(originSingleCell)[1]
    originPointSample2 <- apply(originSingleCell,
                                FUN = function(x){sf::st_as_sf(data.frame(x),
                                                               coords = c("Longitude", "Latitude"),
                                                               crs = 4326)},
                                MARGIN = 3)
    if(!sf::st_crs(originSites)==sf::st_crs(originPointSample2[[1]])){
      originPointSample2 <- sapply(originPointSample2, sf::st_transform,
                                   crs = resampleProjection)
    }

    originCon <- sapply(originPointSample2, FUN = function(z){
      suppressMessages(as.numeric(unclass(sf::st_intersects(x = z, y = originSites,
                                           sparse = TRUE))))})
    if (!any(is.na(originCon)))
      originPointsInSites <- TRUE
    else if (verbose > 0){
      cat('Single cell origin points supplied, but some points (proportion',
          sum(is.na(originCon))/length(originCon), ') not in originSites\n')
    }
  }else{
    originCon <- NULL
  }

  # if (!is.null(targetPoints) && nAnimals < nAnimalsTotal) {
  #   if (verbose > 0)
  #     cat("Filling in dummy target points for CMR data")
  #   dummyPoint <- targetPoints[1, ]
  #   for (i in 1:(nAnimalsTotal - nAnimals)) {
  #     targetPoints <- rbind(targetPoints,
  #                           dummyPoint)
  #   }
  # }
  # if (!is.null(originPoints) && nAnimals < nAnimalsTotal) {
  #   if (verbose > 0)
  #     cat("Filling in dummy origin points for CMR data")
  #   dummyPoint <- originPoints[1, ]
  #   for (i in 1:(nAnimalsTotal - nAnimals)) {
  #     originPoints <- rbind(originPoints,
  #                           dummyPoint)
  #   }
  # }

  if (!is.null(targetRelAbund) && any(captured=="target")) {
    if (coda::is.mcmc(targetRelAbund) || coda::is.mcmc.list(targetRelAbund)) {
      targetRelAbund <- as.matrix(targetRelAbund)
    }
    if (is.matrix(targetRelAbund) && dim(targetRelAbund)[1]>1) {
      abundFixed <- FALSE
      if (dim(targetRelAbund)[2]>nTargetSites)
        abundParams <- paste('relN[', 1:nTargetSites, ']', sep='')
      else if (dim(targetRelAbund)[2]==nTargetSites)
        abundParams <- 1:nTargetSites
      else
        stop('Number of target sites must be constant between sites and abundance')
      if (dim(targetRelAbund)[1] >= nBoot)
        abundRows <- round(seq(from = 1, to = dim(targetRelAbund)[1],
                                      length.out = nBoot))
      else
        stop("You need at least nSamples rows to targetRelAbund")
      targetRelAbund <- as.matrix(targetRelAbund[abundRows, abundParams])
      abundBase <- colMeans(targetRelAbund)
    }
    else {
      abundFixed <- TRUE
      if (length(targetRelAbund)!=nTargetSites)
        stop('Number of target sites must be constant between sites and abundance')
      abundBase <- targetRelAbund
      targetRelAbund <- matrix(targetRelAbund, nBoot, nTargetSites, TRUE)
    }
    weights <- array(0, c(nBoot, nAnimalsTotal))
    if (length(dim(targetAssignment))==2) {
      ta <- apply(targetAssignment, 1, which.max)
      if (is.list(ta)) {
        ta[lengths(ta)==0] <- NA
        ta <- unlist(ta)
      }
    }
    else
      ta <- targetAssignment
    nOriginAnimals <- sum(captured != "target")
    nTargetAnimals <- rep(NA, nTargetSites)
    for (i in 1:nTargetSites) {
      nTargetAnimals[i] <- sum(captured=="target" & ta==i)
      if (nTargetAnimals[i] > 0)
        weights[ , captured=="target" & ta==i] <- targetRelAbund[ , i] /
          nTargetAnimals[i] * (nAnimalsTotal - nOriginAnimals) / nAnimalsTotal
    }
    if (nOriginAnimals>0) {
      if (all(nTargetAnimals>0))
        weights[ , captured!="target"] <- 1/nAnimalsTotal
      else {
        t0 <- which(nTargetAnimals>0)
        weights[ , captured!="target"] <- 1/nAnimalsTotal +
          rowSums(targetRelAbund[ , t0]) *
          sum(nTargetAnimals) / nAnimalsTotal / nOriginAnimals
        # weights[captured!="target" & ta %in% t0] <- (1 - sum(weights)) /
        #   sum(captured!="target" & ta %in% t0)
        if (sum(captured!="target" & ta %in% t0) == 0)
          warning("Not all target sites have likely data. Estimates will be biased.")
      }
    }
    else if (any(nTargetAnimals==0)) {
      warning("Not all target sites have data. Estimates will be biased.")
    }
  }
  else {
    weights <- NULL
    if (any(captured!="origin"))
      warning("Unless target site data were collected in proportion to abundance ",
              "at those sites, estimates will be biased. We recommend including ",
              "an estimate of target site relative abundance with the ",
              "targetRelAbund argument, to allow resampling in proportion to ",
              "abundance.")
  }

  sites.array <- psi.array <- array(0, c(nBoot, nOriginSites, nTargetSites),
                                    dimnames = list(1:nBoot, originNames,
                                                    targetNames))
  if (!is.null(banded))
    r.array <- array(0, c(nBoot, nTargetSites),
                     dimnames = list(1:nBoot, targetNames))
  else
    r.array <- NULL
  if (is.null(dim(originAssignment))){
    originAssignment <- array(originAssignment)}
  if (is.null(dim(targetAssignment))){
    targetAssignment <- array(targetAssignment)}

  if (is.null(fixedZero)) {
    nFixedZero <- 0
  }
  else {
    nFixedZero <- nrow(fixedZero)
  }
  countFailed <- 0
  failed <- FALSE
  if (length(dim(originAssignment))==2){
    pointOriginAssignment <- apply(originAssignment, 1, which.max)
  }
  else{
    pointOriginAssignment <- as.vector(originAssignment)
  }
  if (length(pointOriginAssignment) > nAnimals){
    if (nAnimals==0)
      pointOriginAssignment <- NULL
    else
      pointOriginAssignment <- pointOriginAssignment[1:nAnimals]
  }
  if (length(dim(targetAssignment))==2){
    pointTargetAssignment <- apply(targetAssignment, 1, which.max)
  }
  else{
    pointTargetAssignment <- as.vector(targetAssignment)
  }
  if (length(pointTargetAssignment) > nAnimals) {
    if (nAnimals==0)
      pointTargetAssignment <- NULL
    else
      pointTargetAssignment <- pointTargetAssignment[1:nAnimals]
  }
  if (length(pointTargetAssignment) == length(pointOriginAssignment)) {
   #pointSites <- table(pointOriginAssignment, pointTargetAssignment)
   psi_r <- calcTransition(banded, reencountered,
                           originAssignment = pointOriginAssignment,
                           targetAssignment = pointTargetAssignment,
                           originNames = originNames,
                           targetNames = targetNames,
                           method = "BFGS")
   pointPsi <- psi_r$psi
   point_r <- psi_r$r
  }
  else {
    pointPsi <- NULL
    point_r <- NULL
  }
  boot <- 1
  if (verbose > 0)
    cat("Starting bootstrap\n")
  while (boot <= nBoot) {
    if (verbose > 1 || verbose == 1 && boot %% 100 == 0)
      cat("Bootstrap Run", boot, "of", nBoot, "at", date(), "\n")
    # Make sure have animals from every origin site
    origin.sample <- c() # Start with zero origin sites
    while (length(unique(origin.sample)) < nOriginSites) { #2
      # Sample individual animals with replacement
      animal.sample <- sample.int(m, replace=TRUE, prob = weights[boot,])
      if (any(captured[animal.sample]!='origin')) {
        if (length(dim(originAssignment))==2)
          assignment <- originAssignment[animal.sample, , drop = FALSE]
        else
          assignment <- originAssignment[animal.sample, drop = FALSE]
        oSamp <- locSample(isGL = (isGL[animal.sample] & captured[animal.sample]!='origin'),
                           isRaster = (isRaster[animal.sample] & captured[animal.sample]!='origin'),
                           isProb = (isProb[animal.sample] & captured[animal.sample]!='origin'),
                           isTelemetry = (isTelemetry[animal.sample] |
                                            isCMR[animal.sample] |
                                            captured[animal.sample]=='origin'),
                           geoBias = geoBiasOrigin,
                           geoVCov = geoVCovOrigin,
                           points = originPoints[animal.sample, ],
                           matvals = originRasterXYZ[, c(1:2, animal.sample + 2)],
                           matvals_crs = originRasterXYZcrs,
                           singleCell = originSingleCell[,,animal.sample],
                           overlap1 = originCon[,animal.sample],
                           pointsInSites = originPointsInSites,
                           assignment = assignment,
                           sites = originSites,
                           resampleProjection = resampleProjection,
                           nSim = nSim,
                           maxTries = maxTries)
        if (!is.null(oSamp$notfind)) {
          oSamp$notfind$Animal <- animal.sample[oSamp$notfind$Animal]
          notfind <- unique(oSamp$notfind)
          stop('maxTries (',maxTries,') reached during origin location sampling, exiting. ',
               'Animal(s) where location sampling failed to fall in sites:\n',
               paste(utils::capture.output(print(notfind, row.names = FALSE)), collapse = "\n"),
               '\nExamine originSites',
               ifelse(any(notfind$isGL),
                      ', geoBiasOrigin, geoVcovOrigin, originPoints', ''),
               ifelse(any(notfind$isRaster), ', originRaster', ''),
               ifelse(any(notfind$isTelemetry), ', originPoints/captured', ''),
               ', and resampleProjection to determine why sampled points fell outside sites.')
        }
        origin.sample <- oSamp$site.sample
        if (verbose > 2)
          cat(' ', oSamp$draws, 'origin draw(s) (of length', nSim, 'and of', maxTries, 'possible).\n')
      }
      else {
        # Get origin population for each animal sampled
        if (length(dim(originAssignment))==2){
          origin.sample <- apply(originAssignment[animal.sample, ], 1, which.max)
          if (is.list(origin.sample)) {
            origin.sample[lengths(origin.sample)==0] <- NA
            origin.sample <- unlist(origin.sample)
          }
        }
        else
          origin.sample <- originAssignment[animal.sample]
      }
    }
    if (any(captured[animal.sample]!="target")) {
      if (length(dim(targetAssignment))==2)
        assignment <- targetAssignment[animal.sample, , drop = FALSE]
      else
        assignment <- targetAssignment[animal.sample, drop = FALSE]
      tSamp <- locSample(isGL = (isGL[animal.sample] & captured[animal.sample] != "target"),
                         isRaster = (isRaster[animal.sample] & captured[animal.sample] != "target"),
                         isProb = (isProb[animal.sample] & captured[animal.sample] != "target"),
                         isTelemetry = (isTelemetry[animal.sample] |
                                          isCMR[animal.sample] |
                                          captured[animal.sample] == "target"),
                         geoBias = geoBias, geoVCov = geoVCov,
                         points = targetPoints[animal.sample, ],
                         matvals = targetRasterXYZ[, c(1:2, animal.sample + 2)],
                         matvals_crs = targetRasterXYZcrs,
                         singleCell = targetSingleCell[,,animal.sample],
                         pointsInSites = targetPointsInSites,
                         overlap1 = targetCon[, animal.sample],
                         sites = targetSites,
                         assignment = assignment,
                         resampleProjection = resampleProjection, nSim = nSim,
                         maxTries = maxTries)
      if (!is.null(tSamp$notfind)) {
        tSamp$notfind$Animal <- animal.sample[tSamp$notfind$Animal]
        notfind <- unique(tSamp$notfind)
        stop('maxTries (',maxTries,') reached during target location sampling, exiting. ',
             'Animal(s) where location sampling failed to fall in sites:\n',
             paste(utils::capture.output(print(notfind, row.names = FALSE)), collapse = "\n"),
             '\nExamine targetSites',
             ifelse(any(notfind$isGL),
                    ', geoBiasOrigin, geoVcovOrigin, targetPoints', ''),
             ifelse(any(notfind$isRaster), ', targetRaster', ''),
             ifelse(any(notfind$isTelemetry), ', targetPoints/captured', ''),
             ', and resampleProjection to determine why sampled points fell outside sites.')
      }
      target.sample <- tSamp$site.sample
      target.sample[target.sample==0] <- NA
      #target.point.sample <- tSamp$target.point.sample
      if (verbose > 2)
        cat(' ', tSamp$draws, 'target draw(s) (of length', nSim, 'and of', maxTries, 'possible).\n')
    }
    else {
      # Get target population for each animal sampled
      if (length(dim(targetAssignment))==2){
        target.sample <- apply(targetAssignment[animal.sample], 1, which.max)
        if (is.list(target.sample)) {
          target.sample[lengths(target.sample)==0] <- NA
          target.sample <- unlist(target.sample)
        }
      }
      else
        target.sample <- targetAssignment[animal.sample]
    }
    # Now that we have breeding and non-breeding site for point, add to transition count matrix
    sites <- table(origin.sample,
                   target.sample,
                   useNA = "no")
    sites.array[boot, as.integer(rownames(sites)),
                as.integer(colnames(sites))] <- sites
    if (nFixedZero > 0) {
      for (i in 1:nFixedZero) {
        if (sites.array[boot, fixedZero[i, 1], fixedZero[i, 2]] > 0) {
          failed <- TRUE
          countFailed <- countFailed + 1
          if (countFailed > nBoot * 100)
            stop("estTransition stopped because getting very high number of bootstrap runs:\n",
                 countFailed, " where animals use transitions fixed to zero.\n",
                 "You should examine fixedZero and the data to make sure those ",
                 "transition probabilities are really zero")
          sites.array[boot,,] <- 0
          break
        }
      }
      if (failed) {
        failed <- FALSE
        next
      }
    }
    if (is.null(banded)) {
      banded.sample <- NULL
      reencountered.sample <- NULL
    }
    else {
      reencountered.sample <- table(factor(origin.sample[isCMR[animal.sample]],
                                           levels = 1:nOriginSites),
                                    factor(target.sample[isCMR[animal.sample]],
                                           levels = 1:nTargetSites),
                                    useNA = "no")
      banded.sample <- table(factor(origin.sample[isCMR[animal.sample]],
                                    levels = 1:nOriginSites))
    }
    # Use new function that allows for CMR data
    psi_r <- calcTransition(banded.sample, reencountered.sample,
                            originAssignment = origin.sample[!isCMR[animal.sample]],
                            targetAssignment = target.sample[!isCMR[animal.sample]],
                            originNames = originNames,
                            targetNames = targetNames,
                            method = "BFGS")
    if (any(is.na(psi_r$psi)) || any(psi_r$psi < 0) || any(psi_r$psi > 1)) {
      if (verbose > 2)
        cat(' Bootstrap estimate producing nonsense psi; drawing again\n')
      next
    }
    psi.array[boot, , ] <- psi_r$psi
    if (!is.null(banded))
      r.array[boot, ] <- psi_r$r
    boot <- boot + 1
  }
  if (countFailed > 0)
    warning(countFailed, " bootstrap ",ifelse(countFailed>1, "runs", "run"),
            " failed due to animals using transitions fixed at zero with ",
            nBoot, " successful. If this ratio is high, you should examine ",
            "fixedZero and the data to make sure those transition ",
            "probabilities really are zero\n")
  if (method=="bootstrap") {
    meanPsi <- apply(psi.array, 2:3, mean)
    medianPsi <- apply(psi.array, 2:3, median)
    sePsi <- apply(psi.array, 2:3, sd)
    simpleCIPsi <- apply(psi.array, 2:3, quantile, probs = c(alpha/2, 1-alpha/2),
                         na.rm=TRUE, names = FALSE)
    psi.matrix <- array(c(psi.array), c(nBoot, nOriginSites * nTargetSites),
                        list(NULL, paste(rep(originNames, nTargetSites),
                                         rep(targetNames, each = nOriginSites),
                                         sep = "#")))
    psi.mcmc <- coda::as.mcmc(psi.matrix)
    hpdCI <- coda::HPDinterval(psi.mcmc, 1-alpha)
    hpdCI <- array(hpdCI, c(nOriginSites, nTargetSites, 2),
                   list(originNames, targetNames, c("lower", "upper")))
    hpdCI <- aperm(hpdCI, c(3, 1, 2))
    bcCIPsi <- array(NA, dim = c(2, nOriginSites, nTargetSites),
                     dimnames = list(NULL, originNames, targetNames))
    for (i in 1:nOriginSites) {
      for (j in 1:nTargetSites) {
        psi.z0 <- qnorm(sum(psi.array[, i, j] < meanPsi[i, j], na.rm = TRUE) /
                          length(which(!is.na(psi.array[, i, j]))))
        bcCIPsi[ , i, j] <- quantile(psi.array[, i, j],
                                     pnorm(2 * psi.z0 + qnorm(c(alpha/2, 1-alpha/2))),
                                     na.rm=TRUE, names = FALSE)
      }
    }
    if (!is.null(r.array)){
      mean.r <- apply(r.array, 2, mean)
      median.r <- apply(r.array, 2, median)
      se.r <- apply(r.array, 2, sd)
      simpleCIr <- apply(r.array, 2, quantile, probs = c(alpha/2, 1-alpha/2),
                         na.rm=TRUE, names = FALSE)
      r.mcmc <- coda::as.mcmc(r.array)
      hpdCIr <- coda::HPDinterval(r.mcmc, 1-alpha)
      hpdCIr <- aperm(hpdCIr, c(2, 1))
      bcCIr <- array(NA, dim = c(2, nTargetSites),
                     dimnames = list(c("lower", "upper"), targetNames))
      for (j in 1:nTargetSites) {
        r.z0 <- qnorm(sum(r.array[ ,j] < mean.r[j], na.rm = TRUE) /
                          length(which(!is.na(r.array[ ,j]))))
        bcCIr[ , j] <- quantile(r.array[ ,j],
                                pnorm(2 * r.z0 + qnorm(c(alpha/2, 1-alpha/2))),
                                na.rm=TRUE, names = FALSE)
      }
      r <- list(sample = r.array, mean = mean.r,
                se = se.r, simpleCI = simpleCIr,
                bcCI = bcCIr, hpdCI = hpdCIr, median = median.r)
    }
    else
      r <- NULL
  }
  else {
    meanPsi <- apply(psi.array, 2:3, mean) * m / nAnimals
    medianPsi <- apply(psi.array, 2:3, median)
    sePsi <- apply(psi.array, 2:3, sd) * sqrt(m / nAnimals)
    simpleCIPsi <- apply(psi.array, 2:3, quantile, probs = c(alpha/2, 1-alpha/2),
                         na.rm=TRUE, names = FALSE)
  }
  if (returnAllInput) {
    input <-list(sampleSize = nAnimals, originSites = originSites,
                 targetSites = targetSites,
                 originPoints = originPoints,
                 targetPoints = targetPoints,
                 originAssignment = originAssignment,
                 targetAssignment = targetAssignment,
                 originNames = originNames,
                 targetNames = targetNames,
                 nSamples = nBoot,
                 isGL=isGL, isTelemetry = isTelemetry,
                 isRaster = isRaster, isProb = isProb,
                 captured = captured,
                 geoBias=geoBias, geoVCov=geoVCov,
                 geoBiasOrigin = geoBiasOrigin,
                 geoVCovOrigin = geoVCovOrigin,
                 targetRaster = targetRaster,
                 originRaster = originRaster,
                 verbose = verbose,
                 alpha = alpha,
                 resampleProjection = resampleProjection,
                 nSim = nSim, maxTries = maxTries,
                 dataOverlapSetting = dataOverlapSetting,
                 fixedZero = fixedZero,
                 targetRelAbund = targetRelAbund,
                 banded = banded,
                 reencountered = reencountered,
                 method = method,
                 m = m,
                 returnAllInput = TRUE)
  }
  else {
    input <-list(sampleSize = nAnimals,
                 originNames = originNames,
                 targetNames = targetNames,
                 alpha = alpha,
                 method = method,
                 returnAllInput = FALSE)
  }
  return(list(psi = list(sample = psi.array, mean = meanPsi, se = sePsi,
                         simpleCI = simpleCIPsi, bcCI = bcCIPsi, hpdCI = hpdCI,
                         median = medianPsi, point = pointPsi),
              r = r,
              input = input,
              BUGSoutput = NULL))
}

#' Estimate psi (transition probabilities between locations in two phases of
#' the annual cycle)
#'
#' Estimation and resampling of uncertainty for psi (transition probabilities
#' between origin sites in one phase of the annual cycle and target sites in
#' another for migratory animals). Data can be from any combination of
#' geolocators (GL), telemetry/GPS, intrinsic markers such as isotopes and
#' genetics, and band/ring reencounter data.
#'
#' @param originSites A polygon spatial layer (sf - MULTIPOLYGON)
#'  defining the geographic representation of sites in the
#'  origin season.
#' @param targetSites A polygon spatial layer (sf - MULTIPOLYGON)
#'  defining the geographic representation of sites in the
#'  target season.
#' @param originPoints A \code{sf} or \code{SpatialPoints} object, with number
#'  of rows or length being the number of animals tracked. Each point indicates
#'  the origin location of an animal (or point estimate of same, for GL animals
#'  released on target sites). Note that to simplify input of multiple
#'  data types both between and for the same animal, if origin points are
#'  provided for any animal, they must be provided for all except banding data
#'  (can be dummy values), unless \code{dataOverlapSetting} is set to "none".
#' @param targetPoints For GL or telemetry data, a \code{sf} or
#'  \code{SpatialPoints} object, with length or number of rows number of animals
#'  tracked. Each point indicates the point estimate location of an animal in
#'  the target season. Note that to simplify input of multiple
#'  data types both between and for the same animal, if target points are
#'  provided for any animal, they must be provided for all except banding data
#'  (can be dummy values), unless \code{dataOverlapSetting} is set to "none".
#' @param originAssignment Assignment of animals to origin season sites. Either
#'  an integer vector with length number of animals tracked or a matrix of
#'  probabilities with number of animals tracked rows and number of origin sites
#'  columns (and rows summing to 1). The latter only applies to animals released
#'  in the target sites where there is uncertainty about their origin site, for
#'  example from genetic population estimates from the rubias package.
#'  Optional, but some combination of these inputs should be defined. Note that
#'  if \code{originAssignment} is a probability table, animals with known origin
#'  sites can have 1 in that column and 0s in all others. Also note that if
#'  \code{method} is "MCMC", anything in \code{originAssignment} and
#'  \code{targetAssignment} will be assumed to represent animals tracked via
#'  telemetry, with known origin and target sites.
#' @param targetAssignment Assignment of animals to target season sites. Either
#'  an integer vector with length number of animals tracked or a matrix of
#'  probabilities with number of animals tracked rows and number of target sites
#'  columns (and rows summing to 1). The latter only applies to animals released
#'  in the origin sites where there is uncertainty about their target site, for
#'  example from genetic population estimates from the rubias package.
#'  Optional, but some combination of these inputs needs to be defined. Note
#'  that if \code{targetAssignment} is a probability table, animals with known
#'  target sites can have 1 in that column and 0s in all others.
#' @param originNames Optional, but recommended to keep track. Vector of names
#'  for the origin sites. If not provided, the function will either try to get
#'  these from another input or provide default names (capital letters).
#' @param targetNames Optional, but recommended to keep track. Vector of names
#'  for the target sites. If not provided, the function will either try to get
#'  these from another input or provide default names (numbers).
#' @param nSamples Number of post-burn-in MCMC samples to store (\code{method}
#'  == "MCMC") OR number of bootstrap runs for \code{method}
#'  == "bootstrap". In the latter case, animals are sampled with replacement
#'  for each. For all, the purpose is to estimate sampling uncertainty.
#' @param isGL Indicates whether or which animals were tracked with geolocators.
#'  Should be either single TRUE or FALSE value, or vector with length of
#'  number of animals tracked, with TRUE or FALSE for each animal in data
#'  (except those in \code{banded}, which are handled separately). For
#'  TRUE animals, the model applies \code{geoBias} and \code{geoVCov} to
#'  \code{targetPoints} where \code{captured} == "origin" or "neither" and
#'  \code{geoBiasOrigin} and \code{geoVCovOrigin} to
#'  \code{originPoints} where \code{captured} == "target" or "neither".
#'  Geolocator data should be entered as \code{originPoints} and
#'  \code{targetPoints}.
#' @param isTelemetry Indicates whether or which animals were tracked with
#'  telemetry/GPS (no location uncertainty on either end).
#'  Should be either single TRUE or FALSE value, or vector with length of
#'  number of animals tracked, with TRUE or FALSE for each animal in data
#'  (except those in \code{banded}, which are handled separately).
#'  Telemetry data can be entered as points or using the \code{targetAssignment}
#'  and \code{originAssignment} arguments.
#' @param isRaster Indicates whether or which animals were tracked with
#'  intrinsic markers (e.g., genetics or isotopes), with location uncertainty
#'  expressed as a raster of probabilities by grid cells, either in
#'  \code{targetRaster} or \code{originRaster}. Should be either single TRUE or
#'  FALSE value, or vector with length of number of animals tracked, with TRUE
#'  or FALSE for each animal in data (except those in \code{banded}, which are
#'  handled separately).
#' @param isProb Indicates whether or which animals were tracked with
#'  intrinsic markers (e.g., genetics or isotopes), with location uncertainty
#'  expressed as a probability table, either in \code{targetAssignment} or
#'  \code{originAssignment}. Should be either single TRUE or FALSE value, or
#'  vector with length of number of animals tracked, with TRUE or FALSE for each
#'  animal in data (except those in \code{banded}, which are handled separately).
#' @param captured Indicates whether or which animals were captured in the
#'  origin sites, the target sites, or neither (another phase of the annual
#'  cycle). Location uncertainty will only be applied where the animal was not
#'  captured. So this doesn't matter for telemetry data, and is assumed to be
#'  "origin" for band return data. Should be either single "origin" (default),
#'  "target", or "neither" value, or a character vector with length of number of
#'  animals tracked, with "origin", "target", or "neither" for each animal.
#' @param geoBias For GL data, vector of length 2 indicating expected bias
#'  in longitude and latitude of \code{targetPoints}, in
#'  \code{resampleProjection} units (default meters).
#' @param geoVCov For GL data, 2x2 matrix with expected variance/covariance
#'    in longitude and latitude of \code{targetPoints}, in
#'    \code{resampleProjection} units (default meters).
#' @param geoBiasOrigin For GL data where \code{captured}!="origin", vector of
#'  length 2 indicating expected bias in longitude and latitude of
#'  \code{originPoints}, in \code{resampleProjection} units (default meters).
#' @param geoVCovOrigin For GL data where \code{captured}!="origin", 2x2 matrix
#'  with expected variance/covariance in longitude and latitude of
#'  \code{targetPoints}, in \code{resampleProjection} units (default meters).
#' @param targetRaster For intrinsic tracking data, the results of
#'  \code{isoAssign} or a similar function of class \code{intrinsicAssign} or
#'  class \code{RasterBrick}/\code{RasterStack}, for example from the package
#'  \code{assignR}. In any case, it expresses location uncertainty on target
#'  range, through a raster of probabilities by grid cells.
#' @param originRaster For intrinsic tracking data, the results of
#'  \code{isoAssign} or a similar function of class \code{intrinsicAssign} or
#'  class \code{RasterBrick}/\code{RasterStack}, for example from the package
#'  \code{assignR}. In any case, it expresses location uncertainty on origin
#'  range, through a raster of probabilities by grid cells.
#' @param banded For band return data, a vector or matrix of the number of
#'  released animals from each origin site (including those never reencountered
#'  in a target site). If a matrix, the second dimension is taken as the number
#'  of age classes of released animals; the model estimates reencounter
#'  probability by age class but assumes transition probabilities are the same.
#'  Note that this age model is currently implemented only for \code{method}
#'  set to "MCMC", and only when banding data is analyzed alone (no telemetry
#'  data).
#' @param reencountered For band return data, either a matrix with B rows and W
#'  columns or a B x [number of ages] x W array. Number of animals reencountered
#'  on each target site (by age class banded as) by origin site they came from.
#' @param verbose 0 (default) to 3. 0 prints no output during run (except on
#'  convergence for \code{method} set to "MCMC"). 1 prints an update every 100
#'  samples or bootstraps (or a status bar for "MCMC").  2 prints an update
#'  every sample or bootstrap. 3 also prints the number of draws (for
#'  tuning \code{nSim}).
#' @param alpha Level for confidence/credible intervals provided. Default (0.05)
#'  gives 95 percent CI.
#' @param resampleProjection Projection when sampling from location uncertainty.
#'  Default is Equidistant Conic. The default setting preserves distances
#'  around latitude = 0 and longitude = 0. Other projections may work well,
#'  depending on the location of sites. Ignored unless data are entered using
#'  sites and points and/or rasters.
#' @param nSim Tuning parameter for GL or intrinsic data. Affects only the
#'  speed; 1000 seems to work well with our GL data and 10 for our intrinsic
#'  data, but your results may vary. For data combinations, we put the default
#'  higher (5000) to allow for more data conflicts. Should be integer > 0.
#'  Ignored when \code{method} is "MCMC".
#' @param maxTries Maximum number of times to run a single GL/intrinsic
#'  bootstrap before exiting with an error. Default is 300; you may want to make
#'  a little higher if your \code{nSim} is low and \code{nSamples} is high. Set
#'  to NULL to never exit. This parameter was added to prevent setups where some
#'  sample points never land on target sites from running indefinitely.
#' @param nBurnin For \code{method} set to "MCMC", \code{estTransition} runs a
#'  \code{JAGS} multinomial non-Markovian transitions model, for which it needs
#'  the number of burn-in samples before beginning to store results. Default
#'  5000.
#' @param nChains For \code{method} set to "MCMC", \code{estTransition} runs a
#'  \code{JAGS} multinomial non-Markovian transitions model, for which it needs
#'  the number of MCMC chains (to test for convergence). Default 3.
#' @param nThin For \code{method} set to "MCMC", \code{estTransition} runs a
#'  \code{JAGS} multinomial non-Markovian transitions model, for which it needs
#'  the thinning rate. Default 1.
#' @param dataOverlapSetting When there is more than one type of data, this
#'  setting allows the user some flexibility for clarifying which type(s) of
#'  data apply to which animals. Setting "dummy" (the default) indicates that
#'  there are dummy values within each dataset for the animals that isGL,
#'  isTelemetry, etc. don't have that data type (FALSE values). If no animals
#'  have a data type, no dummy values are required. If no animals have more than
#'  one type of data, the user can simplify processing their data by choosing
#'  setting "none" here. In this case, there should be no dummy values, and only
#'  the animals with a type of data should be included in that dataset. The
#'  third setting ("named") is not yet implemented, but will eventually allow
#'  another way to allow animals with more than one type of data with named
#'  animals linking records. When there is only one type of data, it is fastest
#'  to leave this on the default. Note that banding data entered through
#'  \code{banded} and \code{reencountered} are assumed to have no
#'  overlap with other data types, so none of this applies to those.
#' @param fixedZero When the user has a priori reasons to believe one or more
#'  transition probabilities are zero, they can indicate those here, and the
#'  model will keep them fixed at zero. This argument should be a matrix with
#'  two columns (for row and column of the transition probability matrix) and
#'  number of transitions being fixed to zero rows. For MCMC modeling,
#'  substantial evidence that a transition fixed to zero isn't zero may
#'  cause an error. For bootstrap modeling, a warning
#'  will come up if any bootstrap runs generate the transition fixed to zero,
#'  and the function will quit with an error if a very large number of runs do
#'  (> 10 * nSamples). Fixing transitions to zero may also slow down the
#'  bootstrap model somewhat.
#' @param targetRelAbund When some/all data have location error at origin sites
#'  (i.e., GL, raster, or probability table data with captured = "target" or
#'  "none"), unless the data were collected in proportion to abundance at target
#'  sites, simulation work indicates substantial bias in transition probability
#'  estimates can result. However, if these data are resampled in proportion to
#'  target site abundance, this bias is removed. This argument allows the user
#'  to provide an estimate of relative abundance at the target sites. Either
#'  a numeric vector of length [number target sites] that sums to 1, or an mcmc
#'  object (such as is produced by \code{\link{modelCountDataJAGS}}) or matrix
#'  with at least \code{nSamples} rows. If there are more than [number target
#'  sites] columns, the relevant columns should be labeled "relN[1]" through
#'  "relN[number target sites]".
#' @param method This important setting lets the user choose the estimation
#'  method used: bootstrap or MCMC (Markov chain Monte Carlo). Bootstrap (the
#'  default) now works with any and all types of data, whereas MCMC currently
#'  only works with banding and telemetry data (enter telemetry data for MCMC
#'  using \code{originAssignment} and \code{targetAssignment}, not
#'  \code{originPoints} and \code{targetPoints}). However, MCMC is
#'  usually faster (and may be a bit more accurate). The third option,
#'  "m-out-of-n-bootstrap", is still under development and should be left alone.
#' @param m We read that the m-out-of-n-bootstrap method may improve the
#'  coverage of confidence intervals for parameters on or near a boundary (0 or
#'  1 in this case). So we're testing that out. This still under development and
#'  not for the end user. In the m-out-of-n-bootstrap, m is the number of
#'  samples taken each time (less than the true sample size, n). If the
#'  "m-out-of-n-bootstrap" is chosen under \code{method} but this is left blank,
#'  currently the default is n/4, rounded up (no idea if that is reasonable).
#' @param psiPrior matrix with same dimensions as psi. Only relevant when
#'  \code{method} is "MCMC". Each row provides a Dirichlet
#'  (https://en.wikipedia.org/wiki/Dirichlet_distribution) prior on the
#'  transition probabilities from that origin site. The default (NULL) supplies
#'  Dirichlet parameters of all 1s, which is a standard uninformative Dirichlet
#'  prior. Setting these to other positive numbers is useful when you think a
#'  priori that certain transitions are unlikely, but don't want to rule them
#'  out altogether using \code{fixedZero}.
#' @param returnAllInput if TRUE (the default) the output includes all of the
#'  inputs. If FALSE, only the inputs currently used by another MigConnectivity
#'  function are included in the output. Switch this if you're worried about
#'  computer memory (and the output will be much slimmer).
#'
#' @return \code{estTransition} returns a list with the elements:
#' \describe{
#'   \item{\code{psi}}{List containing estimates of transition probabilities:
#'   \itemize{
#'    \item{\code{sample}} Array of sampled values for psi. \code{nSamples} x
#'      [number of origin sites] x [number of target sites]. Provided to allow
#'      the user to compute own summary statistics.
#'    \item{\code{mean}} Main estimate of psi matrix. [number of origin sites]
#'      x [number of target sites].
#'    \item{\code{se}} Standard error of psi, estimated from SD of
#'      \code{psi$sample}.
#'    \item{\code{simpleCI}} \code{1 - alpha} confidence interval for psi,
#'      estimated as \code{alpha/2} and \code{1 - alpha/2} quantiles of
#'      \code{psi$sample}.
#'    \item{\code{bcCI}} Bias-corrected \code{1 - alpha} confidence interval
#'      for psi. May be preferable to \code{simpleCI} when \code{mean} is the
#'      best estimate of psi. \code{simpleCI} is preferred when
#'      \code{median} is a better estimator. When the mean and median are equal,
#'      these should be identical.  Estimated as the
#'      \code{pnorm(2 * z0 + qnorm(alpha / 2))} and
#'      \code{pnorm(2 * z0 + qnorm(1 - alpha / 2))} quantiles of \code{sample},
#'      where z0 is the proportion of \code{sample < mean}.
#'    \item{\code{hpdCI}} \code{1 - alpha} credible interval for psi,
#'      estimated using the highest posterior density (HPD) method.
#'    \item{\code{median}} Median estimate of psi matrix.
#'    \item{\code{point}} Simple point estimate of psi matrix, not accounting
#'      for sampling error.
#'   }
#'   }
#'   \item{\code{r}}{List containing estimates of reencounter probabilities at
#'    each target site. NULL except when using direct band/ring reencounter
#'    data.}
#'   \item{\code{input}}{List containing the inputs to \code{estTransition}.}
#'   \item{\code{BUGSoutput}}{List containing \code{R2jags} output. Only present
#'    when using \code{method} of "MCMC".}
#' }
#'
#' @export
#'
#' @seealso \code{\link{estStrength}}, \code{\link{plot.estMigConnectivity}},
#'   \code{\link{estMC}}, \code{\link{estMantel}}
#'
#' @example inst/examples/estTransitionExamples.R
estTransition <- function(originSites = NULL, targetSites = NULL,
                          originPoints = NULL, targetPoints = NULL,
                          originAssignment = NULL, targetAssignment = NULL,
                          originNames = NULL, targetNames = NULL,
                          nSamples = 1000, isGL = FALSE, isTelemetry = FALSE,
                          isRaster = FALSE, isProb = FALSE,
                          captured = "origin", geoBias = NULL, geoVCov = NULL,
                          geoBiasOrigin = geoBias, geoVCovOrigin = geoVCov,
                          targetRaster = NULL, originRaster = NULL,
                          banded = NULL, reencountered = NULL,
                          verbose = 0, alpha = 0.05,
                          resampleProjection = 'ESRI:102010',
                          nSim = ifelse(any(isRaster & isGL) ||
                                          any(isRaster & isProb) ||
                                          any(isGL & isProb), 5000,
                                        ifelse(any(isGL), 1000,
                                               ifelse(any(isRaster), 10, 1))),
                          maxTries = 300,
                          nBurnin = 5000, nChains = 3, nThin = 1,
                          dataOverlapSetting = c("dummy", "none", "named"),
                          fixedZero = NULL,
                          targetRelAbund = NULL,
                          method = c("bootstrap", "MCMC",
                                     "m-out-of-n-bootstrap"),
                          m = NULL,
                          psiPrior = NULL,
                          returnAllInput = TRUE) {
  dataOverlapSetting <- match.arg(dataOverlapSetting)
  method <- match.arg(method)
  if (method != "MCMC") {
    psi <- estTransitionBoot(isGL=isGL, isTelemetry = isTelemetry,
                             isRaster = isRaster, isProb = isProb,
                             geoBias=geoBias, geoVCov=geoVCov,
                             geoBiasOrigin = geoBiasOrigin,
                             geoVCovOrigin=geoVCovOrigin,
                             targetPoints=targetPoints, targetSites=targetSites,
                             targetAssignment=targetAssignment,
                             originPoints=originPoints, originSites=originSites,
                             originAssignment=originAssignment,
                             originNames=originNames, targetNames=targetNames,
                             targetRaster = targetRaster,
                             originRaster = originRaster,
                             captured = captured,
                             nBoot = nSamples, verbose=verbose,
                             nSim = nSim, alpha = alpha,
                             resampleProjection = resampleProjection,
                             maxTries = maxTries,
                             dataOverlapSetting = dataOverlapSetting,
                             fixedZero = fixedZero,
                             targetRelAbund = targetRelAbund,
                             method = method, m = m,
                             banded = banded, reencountered = reencountered,
                             returnAllInput = returnAllInput)
  }
  else {
    psi <- estTransitionJAGS(banded = banded, reencountered = reencountered,
                             originAssignment = originAssignment,
                             targetAssignment = targetAssignment,
                             alpha = alpha,
                             nSamples = nSamples, verbose = verbose,
                             originNames = originNames,
                             targetNames = targetNames,
                             nBurnin = nBurnin, nThin = nThin,
                             nChains = nChains, fixedZero = fixedZero,
                             psiPrior = psiPrior,
                             returnAllInput = returnAllInput)
  }
  class(psi) <- c("estPsi", "estMigConnectivity")
  return(psi)
}

###############################################################################
# Resampling of uncertainty from geolocators and/or GPS data
###############################################################################
estMCGlGps <- function(originDist, targetDist, originRelAbund, isGL,
                       geoBias, geoVCov,
                       targetPoints, targetSites,
                       sampleSize = NULL,
                       targetAssignment=NULL,
                       originPoints=NULL, originSites=NULL,
                       originAssignment=NULL, originNames=NULL,
                       targetNames=NULL, nBoot = 1000, verbose=0,
                       nSim = 1000, calcCorr=TRUE, alpha = 0.05,
                       approxSigTest = FALSE, sigConst = 0,
            resampleProjection = 'ESRI:102010',
                       maxTries = 300,
                       row0=0,
                       maintainLegacyOutput = FALSE) {

  # Input checking and assignment
  if (!(verbose %in% 0:3))
    stop("verbose should be integer 0-3 for level of output during bootstrap: 0 = none, 1 = every 10, 2 = every run, 3 = number of draws")
  if (length(geoBias)!=2 && any(isGL))
    stop("geoBias should be vector of length 2 (expected bias in longitude and latitude of targetPoints, in resampleProjection units, default meters)")
  if (!isTRUE(all.equal(dim(geoVCov), c(2, 2), check.attributes = FALSE)) && any(isGL))
    stop("geoVCov should be 2x2 matrix (expected variance/covariance in longitude and latitude of targetPoints, in resampleProjection units, default meters)")
  if ((is.null(originPoints) || is.null(originSites)) &&
      is.null(originAssignment))
    stop("Need to define either originAssignment or originSites and originPoints")
  if (calcCorr && is.null(originPoints))
    stop('If calcCorr is TRUE, need to define originPoints')
  if ((is.null(targetPoints) || is.null(targetSites)) &&
      is.null(targetAssignment))
    stop("Need to define either targetAssignment or targetSites and targetPoints")
  nAnimals <- max(nrow(targetPoints), length(targetAssignment))
  if (length(isGL)==1)
    isGL <- rep(isGL, nAnimals)
  if(inherits(originSites, "SpatialPolygonsDataFrame") |
     inherits(originSites, "SpatialPolygons")){
    originSites <- sf::st_as_sf(originSites)}
    if(inherits(targetSites, "SpatialPolygonsDataFrame") |
       inherits(targetSites, "SpatialPolygons")){
         targetSites <- sf::st_as_sf(targetSites)}
  if (is.null(originAssignment)) {
    originAssignment <- suppressMessages(unclass(sf::st_intersects(x = originPoints,
                                                                   y = originSites,
                                                                   sparse = TRUE)))
    originAssignment[lengths(originAssignment)==0] <- NA
    if (any(lengths(originAssignment)>1)){
      stop("Overlapping originSites not allowed\n")
    }
    originAssignment <- unlist(originAssignment)
  }
  if (is.null(targetAssignment)) {
    targetAssignment <- suppressMessages(unclass(sf::st_intersects(x = targetPoints,
                                                                   y = targetSites,
                                                                   sparse = TRUE)))
    targetAssignment[lengths(targetAssignment)==0] <- NA
    if (any(lengths(targetAssignment)>1))
      stop("Overlapping targetSites not allowed\n")
    targetAssignment <- unlist(targetAssignment)
  }
  nOriginSites <- ifelse(is.null(originSites), nrow(originDist),
                         nrow(originSites))
  nTargetSites <- ifelse(is.null(targetSites), nrow(targetDist),
                         nrow(targetSites))
  if (nrow(targetPoints)!=nAnimals && length(targetAssignment)!=nAnimals ||
      length(originAssignment)!=nAnimals)
    stop("isGL should be the same length as originAssignment/originPoints and targetPoints/targetAssignment (number of animals)")
  # if (any(is.na(originAssignment)))
  #   stop("NAs in origin sites (make sure all points fall within polygons)")
  if (length(originRelAbund)!=nOriginSites || sum(originRelAbund)!=1)
    stop('originRelAbund should be vector with length number of origin sites that sums to 1')
  if(!is.null(originPoints))
    if(is.na(sf::st_crs(originPoints)) || !is.null(originSites) & is.na(sf::st_crs(originSites))) {
      stop('Coordinate system definition needed for originSites & originPoints')
    }
  if(!is.null(targetSites) & is.na(sf::st_crs(targetSites)) || !is.null(targetPoints) & is.na(sf::st_crs(targetPoints))){
    stop('Coordinate system definition needed for targetSites & targetPoints')
  }
  if(!is.null(targetPoints)){targetPoints <- sf::st_transform(targetPoints, resampleProjection)}
  if(!is.null(targetSites)){targetSites <- sf::st_transform(targetSites, resampleProjection)}
  if (is.null(targetNames))
    targetNames <- as.character(1:nTargetSites)
  if (is.null(originNames)) {
    originNames <- LETTERS[1:nOriginSites]
  }
  if (is.null(sampleSize))
    sampleSize <- nAnimals
  if (!identical(dim(originDist),rep(nOriginSites,2)) ||
      !identical(dim(targetDist),rep(nTargetSites,2)))
    stop('Distance matrices should be square with same number of sites of each type as assignments/points (with distances in meters)')
  sites.array <- psi.array <- array(0, c(nBoot, nOriginSites, nTargetSites),
                                    dimnames = list(1:nBoot, originNames,
                                                    targetNames))


  MC <- corr <- rep(NA, nBoot)

  # determine the number of animals from the input data
  nAnimals <- length(originAssignment)

  # Point estimate of MC
  pointSites <- array(0, c(nOriginSites, nTargetSites),
                       dimnames = list(originNames, targetNames))

  for(i in 1:nAnimals)
      pointSites[originAssignment[i], targetAssignment[i]] <-
      pointSites[originAssignment[i], targetAssignment[i]] + 1

  pointPsi <- prop.table(pointSites, 1)

  if (any(is.na(pointPsi)))
    pointMC <- NA
  else
    pointMC <- calcMC(originDist, targetDist, originRelAbund, pointPsi,
                    sampleSize = sampleSize)

  if (calcCorr) {
    pointMantel <- calcMantel(targetPoints = targetPoints,
                              originPoints = originPoints)
    targetDist1 <- pointMantel$targetDist
    originDistStart <- pointMantel$originDist
    pointCorr <- pointMantel$pointCorr
  }

  boot <- 1
  while (boot <= nBoot) {
    if (verbose > 1 || verbose == 1 && boot %% 100 == 0)
      cat("Bootstrap Run", boot, "of", nBoot, "at", date(), "\n")
    # Make sure have birds from every origin site
    origin.sample <- 'Filler' # Start with one origin site
    while (length(unique(origin.sample)) < nOriginSites) { #2
      # Sample individual animals with replacement
      animal.sample <- sample.int(nAnimals, replace=TRUE)
      # Get origin population for each animal sampled
      origin.sample <- originAssignment[animal.sample]
    }
    tSamp <- targetSample(isGL = isGL, geoBias = geoBias, geoVCov = geoVCov,
                          targetPoints = targetPoints, animal.sample = animal.sample,
                          targetSites = targetSites, targetAssignment = targetAssignment,
                          resampleProjection = resampleProjection, nSim = nSim,
                          maxTries = maxTries)
    target.sample <- tSamp$target.sample
    target.point.sample <- tSamp$target.point.sample
    if (verbose > 2)
      cat(' ', tSamp$draws, 'draw(s) (of length', nSim, 'and of', maxTries, 'possible).\n')
    # Now that we have breeding and non-breeding site for point, add to transition count matrix
    sites <- table(origin.sample, target.sample)
    sites.array[boot, as.integer(rownames(sites)), as.integer(colnames(sites))] <- sites
    # Create psi matrix as proportion of those from each breeding site that went to each NB site
    psi.array[boot, , ] <- prop.table(sites.array[boot, , ], 1)
    # Calculate MC from that psi matrix
    MC[boot] <- calcMC(originDist, targetDist, originRelAbund,
                       psi.array[boot, , ], sampleSize)
    if (verbose > 1 || verbose == 1 && boot %% 10 == 0)
      cat(" MC mean:", mean(MC, na.rm=TRUE), "SD:", sd(MC, na.rm=TRUE),
          "low quantile:", quantile(MC, alpha/2, na.rm=TRUE),
          "high quantile:", quantile(MC, 1-alpha/2, na.rm=TRUE), "\n")
    if (calcCorr) {
      colnames(target.point.sample) <- c("x", "y")
      target.point.sample <- sf::st_as_sf(data.frame(target.point.sample),
                                          coords = c("x","y"),
                                          crs = resampleProjection)
      originDist1 <- originDistStart[animal.sample, animal.sample]
      corr[boot] <- calcMantel(originDist = originDist1,
                               targetPoints = target.point.sample)$pointCorr
      if (verbose > 1 || verbose == 1 && boot %% 10 == 0)
        cat(" Correlation mean:", mean(corr, na.rm=TRUE), "SD:", sd(corr, na.rm=TRUE),
            "low quantile:", quantile(corr, alpha/2, na.rm=TRUE),
            "high quantile:", quantile(corr, 1-alpha/2, na.rm=TRUE), "\n")
    }
    if (!is.na(MC[boot]))
      boot <- boot + 1
  }
  MC.z0 <- qnorm(sum(MC<mean(MC, na.rm = TRUE), na.rm = TRUE)/length(which(!is.na(MC))))
  bcCI <- quantile(MC, pnorm(2*MC.z0+qnorm(c(alpha/2, 1-alpha/2))),
                       na.rm=TRUE, type = 8, names = FALSE)
  MC.mcmc <- coda::as.mcmc(MC) # Ha!
  hpdCI <- as.vector(coda::HPDinterval(MC.mcmc, 1-alpha))
  if (!approxSigTest)
    simpleP <- bcP <- NULL
  else {
    if (pointMC > sigConst)
      simpleP <- sum(MC < sigConst) / nBoot
    else
      simpleP <- sum(MC > sigConst) / nBoot
    if (simpleP == 0)
      simpleP <- 0.5 / nBoot
    bcP <- pnorm(qnorm(simpleP) - 2 * MC.z0)
    if (pointMC < sigConst)
      bcP <- 1 - bcP
  }
  if (calcCorr) {
    meanCorr <- mean(corr, na.rm=TRUE)
    medianCorr <- median(corr, na.rm=TRUE)
    seCorr <- sd(corr, na.rm=TRUE)
    simpleCICorr <- quantile(corr, c(alpha/2, 1-alpha/2), na.rm=TRUE,
                             names = FALSE)
    corr.z0 <- qnorm(sum((corr)<meanCorr)/nBoot)
    bcCICorr <- quantile(corr, pnorm(2*corr.z0+qnorm(c(alpha/2, 1-alpha/2))),
                           na.rm=TRUE, names = FALSE)
  } else
    pointCorr <- meanCorr <- medianCorr <- seCorr <- simpleCICorr <- bcCICorr <- NULL
  meanMC <- mean(MC, na.rm=TRUE)
  medianMC <- median(MC, na.rm=TRUE)
  seMC <- sd(MC, na.rm=TRUE)
  simpleCI <- quantile(MC, c(alpha/2, 1-alpha/2), na.rm=TRUE, names = FALSE)
  meanPsi <- apply(psi.array, 2:3, mean)
  medianPsi <- apply(psi.array, 2:3, median)
  sePsi <- apply(psi.array, 2:3, sd)
  simpleCIPsi <- apply(psi.array, 2:3, quantile, probs = c(alpha/2, 1-alpha/2),
                       na.rm=TRUE, names = FALSE)
  bcCIPsi <- array(NA, dim = c(2, nOriginSites, nTargetSites),
                   dimnames = list(NULL, originNames, targetNames))
  for (i in 1:nOriginSites) {
    for (j in 1:nTargetSites) {
      psi.z0 <- qnorm(sum(psi.array[, i, j] < meanPsi[i, j], na.rm = TRUE) /
                        length(which(!is.na(psi.array[, i, j]))))
      bcCIPsi[ , i, j] <- quantile(psi.array[, i, j],
                                   pnorm(2 * psi.z0 + qnorm(c(alpha/2, 1-alpha/2))),
                                   na.rm=TRUE, type = 8, names = FALSE)
    }
  }
  if (maintainLegacyOutput) {
    return(list(sampleMC=MC, samplePsi = psi.array, pointPsi = pointPsi,
                pointMC=pointMC, meanMC=meanMC,
                medianMC=medianMC, seMC=seMC, simpleCI=simpleCI,
                bcCI=bcCI, hpdCI=hpdCI, simpleP = simpleP, bcP = bcP,
                sampleCorr = corr, pointCorr = pointCorr,
                meanCorr = meanCorr, medianCorr = medianCorr, seCorr=seCorr,
                simpleCICorr=simpleCICorr, bcCICorr=bcCICorr,
                inputSampleSize = sampleSize,
                alpha = alpha, sigConst = sigConst,
                psi = list(sample = psi.array, mean = meanPsi, se = sePsi,
                           simpleCI = simpleCIPsi, bcCI = bcCIPsi,
                           median = medianPsi, point = pointPsi),
                MC = list(sample = MC, mean = meanMC, se = seMC,
                          simpleCI = simpleCI, bcCI = bcCI, hpdCI = hpdCI,
                          median = medianMC, point = pointMC,
                          simpleP = simpleP, bcP = bcP),
                corr = list(sample = corr, mean = meanCorr, se = seCorr,
                            simpleCI = simpleCICorr, bcCI = bcCICorr,
                            median = medianCorr, point = pointCorr),
                input = list(originDist = originDist, targetDist = targetDist,
                             originRelAbund = originRelAbund,
                             sampleSize = sampleSize, originSites = originSites,
                             targetSites = targetSites,
                             originPoints = originPoints,
                             targetPoints = targetPoints,
                             originAssignment = originAssignment,
                             targetAssignment = targetAssignment,
                             originNames = originNames,
                             targetNames = targetNames,
                             nSamples = nBoot, nSim = nSim,
                             isGL=isGL, geoBias=geoBias, geoVCov=geoVCov,
                             row0 = row0,
                             verbose = verbose, calcCorr = calcCorr,
                             alpha = alpha,
                             approxSigTest = approxSigTest, sigConst = sigConst,
                             resampleProjection = resampleProjection,
                             maxTries = maxTries, maintainLegacyOutput = TRUE)))
  }
  else {
    return(list(psi = list(sample = psi.array, mean = meanPsi, se = sePsi,
                           simpleCI = simpleCIPsi, bcCI = bcCIPsi,
                           median = medianPsi, point = pointPsi),
                MC = list(sample = MC, mean = meanMC, se = seMC,
                          simpleCI = simpleCI, bcCI = bcCI, hpdCI = hpdCI,
                          median = medianMC, point = pointMC,
                          simpleP = simpleP, bcP = bcP),
                corr = list(sample = corr, mean = meanCorr, se = seCorr,
                            simpleCI = simpleCICorr, bcCI = bcCICorr,
                            median = medianCorr, point = pointCorr),
                input = list(originDist = originDist, targetDist = targetDist,
                             originRelAbund = originRelAbund,
                             sampleSize = sampleSize, originSites = originSites,
                             targetSites = targetSites,
                             originPoints = originPoints,
                             targetPoints = targetPoints,
                             originAssignment = originAssignment,
                             targetAssignment = targetAssignment,
                             originNames = originNames,
                             targetNames = targetNames,
                             nSamples = nBoot, nSim = nSim,
                             isGL=isGL, geoBias=geoBias, geoVCov=geoVCov,
                             row0 = row0,
                             verbose = verbose, calcCorr = calcCorr,
                             alpha = alpha,
                             approxSigTest = approxSigTest, sigConst = sigConst,
                             resampleProjection = resampleProjection,
                             maxTries = maxTries,
                             maintainLegacyOutput = FALSE)))
  }
}

###############################################################################
#
#
# Resampling of uncertainty from intrinsic markers (Stable-Isotopes)
#
#
###############################################################################
estMCisotope <- function(targetDist=NULL,
                         originRelAbund,
                         targetIntrinsic,
                         targetSites = NULL,
                         sampleSize = NULL,
                         originPoints=NULL,
                         originSites=NULL,
                         originDist = NULL,
                         originAssignment=NULL,
                         originNames=NULL,
                         targetNames=NULL,
                         nBoot = 1000,
                         verbose=0,
                         nSim = NULL,
                         calcCorr=TRUE,
                         alpha = 0.05,
                         approxSigTest = FALSE,
                         sigConst = 0,
                         resampleProjection = sf::st_crs(4326),# MigConnectivity::projections$WGS84,
                         maxTries = 300,
                         maintainLegacyOutput = FALSE) {

  # Input checking and assignment
  if (!(verbose %in% 0:3)){
    stop("verbose should be integer 0-3 for level of output during bootstrap: 0 = none, 1 = every 10, 2 = every run, 3 = number of draws")}

  if ((is.null(originPoints) || is.null(originSites)) && is.null(originAssignment)){
    stop("Need to define either originAssignment or originSites and originPoints")}

  if (calcCorr && is.null(originPoints)){
    stop('If calcCorr is TRUE, need to define originPoints')}

  if (is.null(targetIntrinsic)){
    stop("Need to define targetIntrinsic")
  }

  if (is.null(targetSites))
    targetSites <- targetIntrinsic$targetSites
  if(!inherits(targetSites, "sf")){
    targetSites <- sf::st_as_sf(targetSites)
    targetSites <- sf::st_make_valid(targetSites)
    targetSites <- sf::st_union(targetSites)
  }

  if (is.null(targetDist))
    targetDist <- distFromPos(sf::st_coordinates(sf::st_centroid(targetSites)))


  if (!inherits(targetIntrinsic, 'isoAssign'))
    stop("targetIntrinsic should be output of isoAssign when isIntrinsic == TRUE")

  pointsAssigned <- !(is.null(targetIntrinsic$SingleCell) ||
                        all(is.na(targetIntrinsic$SingleCell)))

  if(inherits(originSites, "SpatialPolygonsDataFrame")){
    originSites <- sf::st_as_sf(originSites)
    originSites <- sf::st_union(originSites)
  }

  if (is.null(originAssignment))
    originAssignment <-
      suppressMessages(as.numeric(unclass(sf::st_intersects(x = originPoints,
                                                            y = originSites,
                                                            sparse = TRUE))))

  nAnimals <- ifelse(pointsAssigned, dim(targetIntrinsic$SingleCell)[3],
                     dim(targetIntrinsic$probassign)[3])
  targetSites <- sf::st_transform(targetSites, resampleProjection)

  if (length(originAssignment)!=nAnimals)
    stop("originAssignment/originPoints should be the same length as targetIntrinsic (number of animals)")

  pointsInSites <- FALSE
  if (pointsAssigned && !is.null(targetSites)) {
    nSamples <- dim(targetIntrinsic$SingleCell)[1]
    targCon <- array(NA, c(nSamples, nAnimals))

    for(i in 1:nSamples) {
      temptargCon <- sf::st_as_sf(data.frame(t(targetIntrinsic$SingleCell[i,,])),
                                  coords = c("Longitude","Latitude"),
                                  crs = sf::st_crs(targetSites))

      targCon[i, ] <- suppressMessages(as.numeric(unclass(sf::st_intersects(x = temptargCon, y = targetSites, sparse = TRUE))))
    }

    if (!any(is.na(targCon)))
      pointsInSites <- TRUE
    else if (verbose > 0)
      cat('Single cell assignment points supplied, but some points (proportion',
        sum(is.na(targCon))/length(targCon), ') not in targetSites\n')
  }
  else
    targCon <- NULL
  nOriginSites <- length(unique(originAssignment))
  nTargetSites <- ifelse(is.null(targetSites), nrow(targetDist), nrow(targetSites))
  # cat("nSites figured\n")

  if (any(is.na(originAssignment)))
    stop("NAs in origin sites (make sure all points fall within polygons)")

  if (length(originRelAbund)!=nOriginSites || sum(originRelAbund)!=1)
    stop('originRelAbund should be vector with length number of origin sites that sums to 1')

  if(!is.null(originPoints))
    if(is.na(sf::st_crs(originPoints)) || is.na(sf::st_crs(originSites))) {
      stop('Coordinate system definition needed for originSites & originPoints')
    }



  if (is.null(targetNames))
    # Need to be able to conserve the names from input
    targetNames <- 1:nrow(targetSites)

  if (is.null(originNames))
    # need to conserve names from inputs
    originNames <- 1:nrow(originSites)
  # cat("names figured", targetNames, originNames, "\n")

  if (is.null(sampleSize))
    sampleSize <- nAnimals

  # cat(nOriginSites, "\n")
  # cat(nTargetSites, "\n")
  if (any(dim(originDist)!=rep(nOriginSites,2)) ||
      any(dim(targetDist)!=rep(nTargetSites,2)))
    stop('Distance matrices should be square with same number of sites of each type as assignments/points (with distances in meters)')


  if (calcCorr) {
    originPoints2 <- sf::st_transform(originPoints, 4326)
    originDistStart <- distFromPos(sf::st_coordinates(originPoints2))
  }



  sites.array <- psi.array <- array(0, c(nBoot, nOriginSites, nTargetSites),
                                    dimnames = list(1:nBoot, originNames,
                                                    targetNames))


  MC <- corr <- rep(NA, nBoot)


  boot <- 1
  while (boot <= nBoot) {
    if (verbose > 1 || verbose == 1 && boot %% 100 == 0)
      cat("Bootstrap Run", boot, "of", nBoot, "at", date(), "\n")
    # Make sure have birds from every origin site
    origin.sample <- 'Filler' # Start with one origin site
    while (length(unique(origin.sample)) < nOriginSites) { #2
      # Sample individual animals with replacement
      animal.sample <- sample.int(nAnimals, replace=TRUE)
      # Get origin points for those animals
      if (calcCorr)
        origin.point.sample <- originPoints[animal.sample,]
      # Get origin population for each animal sampled
      origin.sample <- originAssignment[animal.sample]
    }
    if (verbose > 2)
      cat("Origin sample complete, making target sample with",
          ifelse(pointsAssigned, "points assigned,", "no points assigned,"),
          ifelse(pointsInSites, "all points in sites\n", "not all points in sites\n"))
    # Resample from points for each animal
    tSamp <- targetSampleIsotope(targetIntrinsic = targetIntrinsic,
                                 animal.sample = animal.sample,
                                 targetSites = targetSites,
                                 resampleProjection = resampleProjection, nSim = nSim,
                                 maxTries = maxTries, pointsAssigned = pointsAssigned,
                                 targCon = targCon, pointsInSites = pointsInSites)

    target.sample <- tSamp$target.sample
    target.point.sample <- tSamp$target.point.sample
    if (verbose > 2 & !pointsInSites)
      cat(' ', tSamp$draws, 'draw(s) (of length', nSim, 'and of', maxTries, 'possible).\n')
    # Now that we have breeding and non-breeding site for point, add to transition count matrix
    sites <- table(origin.sample, target.sample)
    sites.array[boot, as.integer(rownames(sites)), as.integer(colnames(sites))] <- sites
    # Create psi matrix as proportion of those from each breeding site that went to each NB site
    psi.array[boot, , ] <- prop.table(sites.array[boot, , ], 1)
    # Calculate MC from that psi matrix
    MC[boot] <- calcMC(originDist, targetDist, originRelAbund,
                       psi.array[boot, , ], sampleSize)
    if (verbose > 1 || verbose == 1 && boot %% 10 == 0)
      cat(" MC mean:", mean(MC, na.rm=TRUE), "SD:", sd(MC, na.rm=TRUE),
          "low quantile:", quantile(MC, alpha/2, na.rm=TRUE),
          "high quantile:", quantile(MC, 1-alpha/2, na.rm=TRUE), "\n")
    if (calcCorr) {
      originDist1 <- originDistStart[animal.sample, animal.sample]
      colnames(target.point.sample) <- c("x","y")
      target.point.sample <- sf::st_as_sf(data.frame(target.point.sample),
                                          coords = c("x","y"),
                                          crs = resampleProjection)
      #target.point.sample <- sf::st_sfc(t.p.s.geom, crs = resampleProjection)
      #target.point.sample <- sp::SpatialPoints(target.point.sample,sp::CRS(resampleProjection))
      corr[boot] <- calcMantel(originDist = originDist1, targetPoints = target.point.sample)$pointCorr
      if (verbose > 1 || verbose == 1 && boot %% 10 == 0)
        cat(" Correlation mean:", mean(corr, na.rm=TRUE), "SD:", sd(corr, na.rm=TRUE),
            "low quantile:", quantile(corr, alpha/2, na.rm=TRUE),
            "high quantile:", quantile(corr, 1-alpha/2, na.rm=TRUE), "\n")
    }
    if (!is.na(MC[boot]))
      boot <- boot + 1
  }
  MC.z0 <- qnorm(sum(MC<mean(MC, na.rm = TRUE), na.rm = TRUE) /
                   length(which(!is.na(MC))))
  bcCI <- quantile(MC, pnorm(2*MC.z0+qnorm(c(alpha/2, 1-alpha/2))),
                   na.rm=TRUE, type = 8, names = FALSE)
  MC.mcmc <- coda::as.mcmc(MC) # Ha!
  hpdCI <- as.vector(coda::HPDinterval(MC.mcmc, 1-alpha))
  if (!approxSigTest)
    simpleP <- bcP <- NULL
  else {
    if (mean(MC, na.rm=TRUE) > sigConst)
      simpleP <- sum(MC < sigConst) / nBoot
    else
      simpleP <- sum(MC > sigConst) / nBoot
    if (simpleP == 0)
      simpleP <- 0.5 / nBoot
    bcP <- pnorm(qnorm(simpleP) - 2 * MC.z0)
    if (mean(MC, na.rm=TRUE) < sigConst)
      bcP <- 1 - bcP
  }
  if (calcCorr) {
    meanCorr <- mean(corr, na.rm=TRUE)
    medianCorr <- median(corr, na.rm=TRUE)
    seCorr <- sd(corr, na.rm=TRUE)
    simpleCICorr <- quantile(corr, c(alpha/2, 1-alpha/2), na.rm=TRUE, type = 8,
                             names = FALSE)
    corr.z0 <- qnorm(sum((corr)<meanCorr)/nBoot)
    bcCICorr <- quantile(corr, pnorm(2*corr.z0+qnorm(c(alpha/2, 1-alpha/2))),
                         na.rm=TRUE, type = 8, names = FALSE)
    pointCorr <- NULL
  } else
    pointCorr <- meanCorr <- medianCorr <- seCorr <- simpleCICorr <- bcCICorr <- NULL
  meanMC <- mean(MC, na.rm=TRUE)
  medianMC <- median(MC, na.rm=TRUE)
  seMC <- sd(MC, na.rm=TRUE)
  simpleCI <- quantile(MC, c(alpha/2, 1-alpha/2), na.rm=TRUE, type = 8,
                       names = FALSE)
  meanPsi <- apply(psi.array, 2:3, mean)
  medianPsi <- apply(psi.array, 2:3, median)
  sePsi <- apply(psi.array, 2:3, sd)
  simpleCIPsi <- apply(psi.array, 2:3, quantile, probs = c(alpha/2, 1-alpha/2),
                       na.rm=TRUE, type = 8, names = FALSE)
  bcCIPsi <- array(NA, dim = c(2, nOriginSites, nTargetSites),
                   dimnames = list(NULL, originNames, targetNames))
  for (i in 1:nOriginSites) {
    for (j in 1:nTargetSites) {
      psi.z0 <- qnorm(sum(psi.array[, i, j] < meanPsi[i, j], na.rm = TRUE) /
                        length(which(!is.na(psi.array[, i, j]))))
      bcCIPsi[ , i, j] <- quantile(psi.array[, i, j],
                                   pnorm(2 * psi.z0 + qnorm(c(alpha/2, 1-alpha/2))),
                                   na.rm=TRUE, type = 8, names = FALSE)
    }
  }
  if (maintainLegacyOutput) {
    return(list(sampleMC = MC, samplePsi = psi.array, pointPsi = NULL,
                pointMC = NULL, meanMC = meanMC,
                medianMC = medianMC, seMC = seMC, simpleCI = simpleCI,
                bcCI = bcCI, hpdCI = hpdCI, simpleP = simpleP, bcP = bcP,
                sampleCorr = corr, pointCorr = pointCorr,
                meanCorr = meanCorr, medianCorr = medianCorr, seCorr=seCorr,
                simpleCICorr=simpleCICorr, bcCICorr=bcCICorr,
                inputSampleSize = sampleSize,
                alpha = alpha, sigConst = sigConst,
                psi = list(sample = psi.array, mean = meanPsi, se = sePsi,
                           simpleCI = simpleCIPsi, bcCI = bcCIPsi,
                           median = medianPsi, point = NULL),
                MC = list(sample = MC, mean = meanMC, se = seMC,
                          simpleCI = simpleCI, bcCI = bcCI, hpdCI = hpdCI,
                          median = medianMC, point = NULL,
                          simpleP = simpleP, bcP = bcP),
                corr = list(sample = corr, mean = meanCorr, se = seCorr,
                            simpleCI = simpleCICorr, bcCI = bcCICorr,
                            median = medianCorr, point = pointCorr),
                input = list(originDist = originDist, targetDist = targetDist,
                             originRelAbund = originRelAbund,
                             sampleSize = sampleSize, originSites = originSites,
                             targetSites = targetSites,
                             originPoints = originPoints,
                             originAssignment = originAssignment,
                             originNames = originNames,
                             targetNames = targetNames,
                             nSamples = nBoot, nSim = nSim,
                             verbose = verbose, calcCorr = calcCorr,
                             alpha = alpha,
                             approxSigTest = approxSigTest, sigConst = sigConst,
                             resampleProjection = resampleProjection,
                             maxTries = maxTries,
                             targetIntrinsic = targetIntrinsic,
                             isIntrinsic = TRUE,
                             maintainLegacyOutput = TRUE)))
  }
  else {
    return(list(psi = list(sample = psi.array, mean = meanPsi, se = sePsi,
                           simpleCI = simpleCIPsi, bcCI = bcCIPsi,
                           median = medianPsi, point = NULL),
                MC = list(sample = MC, mean = meanMC, se = seMC,
                          simpleCI = simpleCI, bcCI = bcCI, hpdCI = hpdCI,
                          median = medianMC, point = NULL,
                          simpleP = simpleP, bcP = bcP),
                corr = list(sample = corr, mean = meanCorr, se = seCorr,
                            simpleCI = simpleCICorr, bcCI = bcCICorr,
                            median = medianCorr, point = pointCorr),
                input = list(originDist = originDist, targetDist = targetDist,
                             originRelAbund = originRelAbund,
                             sampleSize = sampleSize, originSites = originSites,
                             targetSites = targetSites,
                             originPoints = originPoints,
                             originAssignment = originAssignment,
                             originNames = originNames,
                             targetNames = targetNames,
                             nSamples = nBoot, nSim = nSim,
                             verbose = verbose, calcCorr = calcCorr,
                             alpha = alpha,
                             approxSigTest = approxSigTest, sigConst = sigConst,
                             resampleProjection = resampleProjection,
                             maxTries = maxTries,
                             targetIntrinsic = targetIntrinsic,
                             isIntrinsic = TRUE,
                             maintainLegacyOutput = FALSE)))
  }
}

###############################################################################
#' Estimate migratory connectivity
#'
#' Resampling of uncertainty for migratory connectivity strength (MC)
#' and transition probabilities (psi) from RMark psi matrix estimates or
#' samples of psi and/or JAGS
#' relative abundance MCMC samples OR SpatialPoints geolocators and/or GPS
#' data OR intrinsic markers such as isotopes. NOTE: active development of this
#' function is ending. We suggest users estimate psi with
#' \code{\link{estTransition}}, MC with \code{\link{estStrength}}, and Mantel
#' correlations (rM) with \code{\link{estMantel}}.
#'
#' @param originDist Distances between the B origin sites.  Symmetric B by B
#'  matrix
#' @param targetDist Distances between the W target sites.  Symmetric W by W
#'  matrix.  Optional for intrinsic data
#' @param originRelAbund Relative abundance estimates at B origin sites. Either
#'  a numeric vector of length B that sums to 1 or an mcmc object with
#'  \code{nSamples} rows  and columns including 'relN[1]' through 'relN[B]'.
#'  Currently, an mcmc object doesn't work with geolocator, GPS, or intrinsic
#'  data
#' @param psi Transition probabilities between B origin and W target sites.
#'  Either a matrix with B rows and W columns where rows sum to 1, an array with
#'  dimensions x, B, and W (with x samples of the transition probability matrix
#'  from another model), or a MARK object with estimates of transition
#'  probabilities.  If you are estimating MC from GPS, geolocator, or intrinsic
#'  data, leave this as NULL
#' @param sampleSize Total sample size of animals that psi will be estimated
#'  from. Should be the number of animals released in one of the origin sites
#'  and observed in one of the target sites.  Optional, but recommended, unless
#'  you are estimating MC from GPS, geolocator, intrinsic, or direct band return
#'  data (in which case the function can calculate it for you)
#' @param originSites If \code{psi} is a MARK object, this must be a numeric
#'  vector indicating which sites are origin.  If using GPS, geolocator, or
#'  intrinsic data, this can be the geographic definition of sites in the
#'  release season
#' @param targetSites If \code{psi} is a MARK object, this must be a numeric
#'  vector indicating which sites are target.  If using GPS, geolocator, or
#'  intrinsic data, this must be the geographic definition of sites in the
#'  non-release season.  Optional for intrinsic data; if left out, the function
#'  will use the \code{targetSites} defined in \code{targetIntrinsic}
#' @param originPoints A \code{POINT} sf object, with length number of
#'    animals tracked.  Each point indicates the release location of an animal
#' @param targetPoints For GL or GPS data, a \code{POINT} sf object, with
#'    length number of animals tracked.  Each point indicates the point estimate
#'    location in the non-release season
#' @param originAssignment Assignment of \code{originPoints} to release season
#'    sites. Integer vector with length number of animals tracked. Optional,
#'    but if using GL or GPS data, either \code{originAssignment} or
#'    \code{originSites} and \code{originPoints} should be defined
#' @param targetAssignment Optional. Point estimate assignment of
#'    \code{targetPoints} to non-release season sites. Integer vector with
#'    length number of animals tracked
#' @param originNames Optional but recommended. Vector of names for the release season sites
#' @param targetNames Optional but recommended. Vector of names for the non-release season
#'    sites
#' @param nSamples Number of times to resample \code{psi} and/or
#'  \code{originRelAbund} OR number of post-burn-in MCMC samples to store (band
#'  data) OR number of times to resample \code{targetPoints}
#'  for intrinsic data OR number of bootstrap runs for GL or GPS data. In
#'  the last two cases, animals are sampled with replacement for each. For all,
#'  the purpose is to estimate sampling uncertainty
#' @param nSim Tuning parameter for GL or intrinsic data. Affects only the
#'    speed; 1000 seems to work well with our GL data and 10 for our intrinsic
#'    data, but your results may vary.  Should be integer > 0
#' @param isGL Indicates whether or which animals were tracked with geolocators.
#'    Should be either single TRUE or FALSE value, or vector with length of
#'    number of animals tracked, with TRUE for animals in
#'    \code{targetPoints} with geolocators and FALSE for animals with GPS
#' @param geoBias For GL data, vector of length 2 indicating expected bias
#'    in longitude and latitude of \code{targetPoints}, in
#'    \code{resampleProjection} units (default meters)
#' @param geoVCov For GL data, 2x2 matrix with expected variance/covariance
#'    in longitude and latitude of \code{targetPoints}, in
#'    \code{resampleProjection} units (default meters)
#' @param row0 If \code{originRelAbund} is an mcmc object, this can be set
#'  to 0 (default) or any greater integer to specify where to stop ignoring
#'  samples ("burn-in")
#' @param verbose 0 (default) to 3. 0 prints no output during run. 1 prints
#'  a line every 100 samples or bootstraps and a summary every 10.  2 prints a
#'  line and summary every sample or bootstrap. 3 also prints the number of
#'  draws (for tuning nSim for GL/intrinsic data only)
#' @param calcCorr In addition to MC, should function also estimate Mantel
#'    correlation between release and non-release locations (GPS or GL data
#'    only)?  Default is FALSE
#' @param alpha Level for confidence/credible intervals provided
#' @param approxSigTest Should function compute approximate one-sided
#'    significance tests (p-values) for MC from the bootstrap?  Default is
#'    FALSE
#' @param sigConst Value to compare MC to in significance test.
#'    Default is 0
#' @param resampleProjection Projection when sampling from geolocator
#'    bias/error. This projection needs units = m. Default is Equidistant
#'    Conic. The default setting preserves distances around latitude = 0 and
#'    longitude = 0. Other projections may work well, depending on the location
#'    of \code{targetSites}.  Ignored unless data are geolocator or GPS
#' @param maxTries Maximum number of times to run a single GL/intrinsic
#'    bootstrap before exiting with an error.  Default is 300.  Set to NULL to
#'    never stop.  This parameter was added to prevent GL setups where some
#'    sample points never land on target sites from running indefinitely
#' @param targetIntrinsic For intrinsic tracking data, the results of
#'    \code{isoAssign} or a similar function, of class \code{intrinsicAssign}
#' @param isIntrinsic Logical indicating whether the animals are tracked via
#'  intrinsic marker (e.g. isotopes) or not.  Currently estMC will only estimate
#'  connectivity for all intrinsically marked animals or all extrinsic (e.g.,
#'  bands, GL, or GPS), so isIntrinsic should be a single TRUE or FALSE
#' @param maintainLegacyOutput version 0.4.0 of \code{MigConnectivity}
#'  updated the structure of the estimates. If you have legacy code that refers
#'  to elements within a \code{estMigConnectivity} object, you can set this
#'  to TRUE to also keep the old structure. Defaults to FALSE
#'
#' @return NOTE: Starting with version 0.4.0 of \code{MigConnectivity}, we've
#' updated the structure of \code{MigConnectivityEstimate} objects. Below we
#' describe the updated structure. If parameter \code{maintainLegacyOutput} is
#' set to TRUE, the list will start with the old structure: \code{sampleMC},
#' \code{samplePsi}, \code{pointPsi}, \code{pointMC}, \code{meanMC},
#' \code{medianMC}, \code{seMC}, \code{simpleCI}, \code{bcCI}, \code{hpdCI},
#' \code{simpleP}, \code{bcP}, \code{sampleCorr}, \code{pointCorr},
#' \code{meanCorr, medianCorr, seCorr, simpleCICorr, bcCICorr},
#' \code{inputSampleSize}, \code{alpha}, and \code{sigConst}.
#'
#' \code{estMC} returns a list with the elements:
#' \describe{
#'   \item{\code{psi}}{List containing estimates of transition probabilities:
#'   \itemize{
#'    \item{\code{sample}} Array of sampled values for psi. \code{nSamples} x
#'      [number of origin sites] x [number of target sites]. Provided to allow
#'      the user to compute own summary statistics.
#'    \item{\code{mean}} Main estimate of psi matrix. [number of origin sites]
#'      x [number of target sites].
#'    \item{\code{se}} Standard error of psi, estimated from SD of
#'      \code{psi$sample}.
#'    \item{\code{simpleCI}} \code{1 - alpha} confidence interval for psi,
#'      estimated as \code{alpha/2} and \code{1 - alpha/2} quantiles of
#'      \code{psi$sample}.
#'    \item{\code{bcCI}} Bias-corrected \code{1 - alpha} confidence interval
#'      for psi.  Preferable to \code{simpleCI} when \code{mean} is the
#'      best estimate of psi. \code{simpleCI} is preferred when
#'      \code{median} is a better estimator. When \code{meanMC==medianMC},
#'      these should be identical.  Estimated as the
#'      \code{pnorm(2 * z0 + qnorm(alpha / 2))} and
#'      \code{pnorm(2 * z0 + qnorm(1 - alpha / 2))} quantiles of \code{sample},
#'      where z0 is the proportion of \code{sample < mean}.
#'    \item{\code{median}} Median estimate of psi matrix.
#'    \item{\code{point}} Simple point estimate of psi matrix, not accounting
#'      for sampling error. NULL when \code{isIntrinsic == TRUE}.
#'   }
#'   }
#'   \item{\code{MC}}{List containing estimates of migratory connectivity
#'    strength:
#'    \itemize{
#'      \item{\code{sample}} \code{nSamples} sampled values for
#'       MC. Provided to allow the user to compute own summary statistics.
#'      \item{\code{mean}} Mean of \code{MC$sample}. Main estimate of MC,
#'       incorporating parametric uncertainty.
#'      \item{\code{se}} Standard error of MC, estimated from SD of
#'       \code{MC$sample}.
#'      \item{\code{simpleCI}} Default\code{1 - alpha} confidence interval for
#'       MC, estimated as \code{alpha/2} and \code{1 - alpha/2} quantiles of
#'       \code{MC$sample}.
#'      \item{\code{bcCI}} Bias-corrected \code{1 - alpha} confidence interval
#'       for MC.  Preferable to \code{MC$simpleCI} when \code{MC$mean} is the
#'       best estimate of MC. \code{MC$simpleCI} is preferred when
#'       \code{MC$median} is a better estimator. When \code{MC$mean==MC$median},
#'       these should be identical.  Estimated as the
#'       \code{pnorm(2 * z0 + qnorm(alpha / 2))} and
#'       \code{pnorm(2 * z0 + qnorm(1 - alpha / 2))} quantiles of \code{MC$sample},
#'       where z0 is the proportion of \code{MC$sample < MC$mean}.
#'      \item{\code{hpdCI}} \code{1 - alpha} credible interval for MC,
#'       estimated using the highest posterior density (HPD) method.
#'      \item{\code{median}} Median of MC, alternate estimator also including
#'       parametric uncertainty.
#'      \item{\code{point}} Simple point estimate of MC, using the point
#'      estimates of \code{psi} and \code{originRelAbund}, not accounting
#'      for sampling error. NULL when \code{isIntrinsic == TRUE}.
#'      \item{\code{simpleP}} Approximate p-value for MC, estimated as the
#'      proportion of bootstrap iterations where MC < \code{sigConst} (or MC >
#'      \code{sigConst} if \code{pointMC < sigConst}).  Note that if the
#'      proportion is 0, a default value of 0.5 / \code{nSamples} is provided,
#'      but this is best interpreted as p < 1 / \code{nSamples}.  NULL when
#'      \code{approxSigTest==FALSE}.
#'      \item{\code{bcP}} Approximate bias-corrected p-value for MC, estimated as
#'      \code{pnorm(qnorm(simpleP) - 2 * z0)}, where z0 is the proportion of
#'      \code{sampleMC < meanMC}.  May be a better approximation of the p-value
#'      than \code{simpleP}, but many of the same limitations apply.  NULL when
#'      \code{approxSigTest==FALSE}.
#'    }
#'   }
#'   \item{\code{corr}}{List containing estimates of rM, an alternate measure of
#'    migratory connectivity strength. NULL when \code{calcCorr==FALSE} or
#'    \code{!is.null(psi)}:
#'    \itemize{
#'     \item{\code{sample}} \code{nBoot} sampled values for continuous
#'      correlation. Provided to allow the user to compute own summary
#'      statistics.
#'     \item{\code{mean, se, simpleCI, bcCI, median, point}} Summary
#'      statistics for continuous correlation bootstraps.
#'    }
#'   }
#'   \item{\code{input}}{List containing the inputs to \code{estMC}, or at least
#'    the relevant ones, such as sampleSize.}
#' }
#' @example inst/examples/estMCExamples.R
#' @seealso \code{\link{estStrength}}, \code{\link{estTransition}},
#'   \code{\link{estMantel}}, \code{\link{calcMC}}, \code{\link{projections}},
#'   \code{\link{isoAssign}}, \code{\link{plot.estMigConnectivity}}
#' @export
#' @references
#' Cohen, E. B., J. A. Hostetler, M. T. Hallworth, C. S. Rushing, T. S. Sillett,
#' and P. P. Marra. 2018. Quantifying the strength of migratory connectivity.
#' Methods in Ecology and Evolution 9: 513 - 524.
#' \href{http://doi.org/10.1111/2041-210X.12916}{doi:10.1111/2041-210X.12916}
#'
#' Cohen, E. B., C. S. Rushing, F. R. Moore, M. T. Hallworth, J. A. Hostetler,
#' M. Gutierrez Ramirez, and P. P. Marra. 2019. The strength of
#' migratory connectivity for birds en route to breeding through the Gulf of
#' Mexico. Ecography 42: 658–669.
#' \href{https://doi.org/10.1111/ecog.03974}{doi:10.1111/ecog.03974}

estMC <- function(originDist, targetDist = NULL, originRelAbund, psi = NULL,
                  sampleSize = NULL,
                  originSites = NULL, targetSites = NULL,
                  originPoints = NULL, targetPoints = NULL,
                  originAssignment = NULL, targetAssignment = NULL,
                  originNames = NULL, targetNames = NULL,
                  nSamples = 1000, nSim = ifelse(isTRUE(isIntrinsic), 10, 1000),
                  isGL = FALSE, geoBias = NULL, geoVCov = NULL, row0 = 0,
                  verbose = 0, calcCorr = FALSE, alpha = 0.05,
                  approxSigTest = FALSE, sigConst = 0,
                  resampleProjection = 'ESRI:102010',
                  maxTries = 300, targetIntrinsic = NULL,
                  isIntrinsic = FALSE,
                  maintainLegacyOutput = FALSE) {
  if (is.null(psi)) {
    if(isIntrinsic) {
      mc <- estMCisotope(targetDist = targetDist,
                         originRelAbund = originRelAbund,
                         targetIntrinsic = targetIntrinsic,
                         targetSites = targetSites, sampleSize = sampleSize,
                         # targetAssignment = targetAssignment,
                         originPoints=originPoints, originSites=originSites,
                         originDist = originDist,
                         originAssignment = originAssignment,
                         originNames = originNames, targetNames = targetNames,
                         nBoot = nSamples, verbose = verbose, nSim = nSim,
                         calcCorr=calcCorr, alpha = alpha,
                         approxSigTest = approxSigTest, sigConst = sigConst,
                         maxTries = maxTries,
                         maintainLegacyOutput = maintainLegacyOutput)
    }
    else #if (is.null(nReleased)) {
      mc <- estMCGlGps(isGL=isGL, geoBias=geoBias, geoVCov=geoVCov,
                       originRelAbund=originRelAbund, sampleSize = sampleSize,
                       originDist=originDist,
                       targetDist=targetDist,
                       targetPoints=targetPoints, targetSites=targetSites,
                       targetAssignment=targetAssignment,
                       originPoints=originPoints, originSites=originSites,
                       originAssignment=originAssignment,
                       originNames=originNames, targetNames=targetNames,
                       nBoot = nSamples, verbose=verbose,
                       nSim = nSim, calcCorr=calcCorr, alpha = alpha,
                       approxSigTest = approxSigTest, sigConst = sigConst,
                       resampleProjection = resampleProjection,
                       maxTries = maxTries,
                       maintainLegacyOutput = maintainLegacyOutput)
  }
  else {
    mc <- estStrength(originRelAbund = originRelAbund,
                        sampleSize = sampleSize, psi = psi,
                        originDist = originDist, targetDist = targetDist,
                        originSites=originSites, targetSites=targetSites,
                        nSamples = nSamples, row0 = row0, verbose=verbose,
                        alpha = alpha, approxSigTest = approxSigTest,
                        sigConst = sigConst, originNames = originNames,
                        targetNames = targetNames,
                        maintainLegacyOutput = maintainLegacyOutput)
  }
  if (calcCorr && is.null(psi))
    class(mc) <- c("estMC", "estPsi", "estMantel", "estMigConnectivity")
  else if (is.null(psi))
    class(mc) <- c("estMC", "estPsi", "estMigConnectivity")
  else
    class(mc) <- c("estMC", "estMigConnectivity")
  return(mc)
}

#' Estimate Mantel correlation (rM) from geolocator, GPS, and/or raster data.
#'
#' Resampling of uncertainty for migratory connectivity strength, as quantified
#' by Mantel correlation (rM), from geolocators, GPS, and/or raster (e.g.,
#' genoscape or isotope) data.
#'
#' @param targetPoints A \code{POINTS} from sf
#'  object, with length number of animals tracked.  Each point indicates the
#'  point estimate location in the non-release season.
#' @param originPoints A \code{POINTS} from sf
#'  object, with length number of animals tracked.  Each point indicates the
#'  release location of an animal.
#' @param isGL Indicates whether or which animals were tracked with geolocators
#'  Should be either single TRUE or FALSE value, or vector with length of
#'  number of animals tracked, with TRUE for animals in  \code{targetPoints}
#'  with geolocators and FALSE for animals without.
#' @param geoBias For GL data, vector of length 2 indicating expected bias
#'  in longitude and latitude of \code{targetPoints}, in
#'  \code{resampleProjection} units (default meters).
#' @param geoVCov For GL data, 2x2 matrix with expected variance/covariance
#'  in longitude and latitude of \code{targetPoints}, in
#'  \code{resampleProjection} units (default meters).
#' @param targetSites A \code{SpatialPolygons}, \code{SpatialPolygonsDataFrame},
#'  or \code{POLYGONS} sf object indicating valid target location(s). Not
#'  needed unless you want to mask out certain areas (e.g. water) and
#'  \code{captured} is "origin" or you want to use a weighted bootstrap based on
#'  \code{targetRelAbund} for animals captured on the target side.
#' @param nBoot Number of bootstrap runs. Animals are sampled with replacement
#'  for each, to estimate sampling uncertainty.
#' @param nSim Tuning parameter for GL or raster data. Affects only the speed;
#'  1000 seems to work well with our GL data.  Should be integer > 0.
#' @param verbose 0 (default) to 3. 0 prints no output during run. 1 prints
#'  a line every 100 bootstraps.  2 prints a line every bootstrap.
#'  3 also prints the number of draws (for tuning nSim only).
#' @param alpha Level for confidence/credible intervals provided.
#' @param resampleProjection Projection when sampling from geolocator
#'  bias/error. This projection needs units = m. Default is Equidistant
#'  Conic. The default setting preserves distances around latitude = 0 and
#'  longitude = 0. Other projections may work well, depending on the location
#'  of \code{targetPoints}.
#' @param maxTries Maximum number of times to run a single GL bootstrap before
#'  exiting with an error.  Default is 300.  Set to NULL to never stop.  This
#'  parameter was added to prevent GL setups where some sample points never
#'  land on target sites from running indefinitely.
#' @param maintainLegacyOutput version 0.4.0 of \code{MigConnectivity}
#'  updated the structure of the estimates. If you have legacy code that refers
#'  to elements within a \code{estMigConnectivity} object, you can set this
#'  to TRUE to also keep the old structure. Defaults to FALSE.
#' @param originSites A \code{SpatialPolygons}, \code{SpatialPolygonsDataFrame},
#'  or \code{POLYGONS} sf object indicating valid origin location(s). Not
#'  needed unless you want to mask out certain areas (e.g. water) and
#'  \code{captured} is "target" or you want to use a weighted bootstrap based on
#'  \code{originRelAbund} for animals captured on the origin side.
#' @param isTelemetry Indicates whether or which animals were tracked with
#'  telemetry/GPS (no location uncertainty on either end).
#'  Should be either single TRUE or FALSE value, or vector with length of
#'  number of animals tracked, with TRUE or FALSE for each animal in data.
#' @param isRaster Indicates whether or which animals were tracked with
#'  intrinsic markers (e.g., genetics or isotopes), with location uncertainty
#'  expressed as a raster of probabilities by grid cells, either in
#'  \code{targetRaster} or \code{originRaster}. Should be either single TRUE or
#'  FALSE value, or vector with length of number of animals tracked, with TRUE
#'  or FALSE for each animal in data.
#' @param captured Indicates whether or which animals were captured in the
#'  origin sites, the target sites, or neither (another phase of the annual
#'  cycle). Location uncertainty will only be applied where the animal was not
#'  captured. So this doesn't matter for telemetry data. Should be either single
#'  "origin" (default), "target", or "neither" value, or a character vector with
#'  length of number of animals tracked, with "origin", "target", or "neither"
#'  for each animal.
#' @param geoBiasOrigin For GL data where \code{captured}!="origin", vector of
#'  length 2 indicating expected bias in longitude and latitude of
#'  \code{originPoints}, in \code{resampleProjection} units (default meters).
#' @param geoVCovOrigin For GL data where \code{captured}!="origin", 2x2 matrix
#'  with expected variance/covariance in longitude and latitude of
#'  \code{targetPoints}, in \code{resampleProjection} units (default meters).
#' @param targetRaster For intrinsic tracking data, the results of
#'  \code{isoAssign} or a similar function of class \code{intrinsicAssign} or
#'  class \code{RasterBrick}/\code{RasterStack}, for example from the package
#'  \code{assignR}. In any case, it expresses location uncertainty on target
#'  range, through a raster of probabilities by grid cells.
#' @param originRaster For intrinsic tracking data, the results of
#'  \code{isoAssign} or a similar function of class \code{intrinsicAssign} or
#'  class \code{RasterBrick}/\code{RasterStack}, for example from the package
#'  \code{assignR}. In any case, it expresses location uncertainty on origin
#'  range, through a raster of probabilities by grid cells.
#' @param dataOverlapSetting When there is more than one type of data, this
#'  setting allows the user some flexibility for clarifying which type(s) of
#'  data apply to which animals. Setting "dummy" (the default) indicates that
#'  there are dummy values within each dataset for the animals that isGL,
#'  isTelemetry, etc. don't have that data type (FALSE values). If no animals
#'  have a data type, no dummy values are required. If no animals have more than
#'  one type of data, the user can simplify processing their data by choosing
#'  setting "none" here. In this case, there should be no dummy values, and only
#'  the animals with a type of data should be included in that dataset. The
#'  third setting ("named") is not yet implemented, but will eventually allow
#'  another way to allow animals with more than one type of data with named
#'  animals linking records. When there is only one type of data, it is fastest
#'  to leave this on the default.
#' @param originRelAbund the proportion of the total abundance in each of B
#'  \code{originSites}. Used to set up the bootstrap to be weighted by relative
#'  abundance (for animals captured on the origin side). Either a numeric vector
#'  of length B that sums to 1, or an mcmc object (such as is produced by
#'  \code{\link{modelCountDataJAGS}}) or matrix with at least B columns.
#'  If there are more than B columns, the relevant columns should be
#'  labeled "relN[1]" through "relN[B]". Optional, but if you don't set it and
#'  at least some animals are captured on the origin side, there's potential for
#'  rM to be biased (if sampling isn't proportional to abundance).
#' @param targetRelAbund the proportion of the total abundance in each of W
#'  \code{targetSites}. Used to set up the bootstrap to be weighted by relative
#'  abundance (for animals captured on the target side). Either a numeric vector
#'  of length W that sums to 1, or an mcmc object (such as is produced by
#'  \code{\link{modelCountDataJAGS}}) or matrix with at least W columns.
#'  If there are more than W columns, the relevant columns should be
#'  labeled "relN[1]" through "relN[W]". Optional, but if you don't set it and
#'  at least some animals are captured on the target side, there's potential for
#'  rM to be biased (if sampling isn't proportional to abundance).
#'
#' @return \code{estMantel} returns a list with elements:
#' \describe{
#'   \item{\code{corr}}{List containing estimates of rM:
#'    \itemize{
#'     \item{\code{sample}} \code{nBoot} sampled values for Mantel
#'      correlation. Provided to allow the user to compute own summary
#'      statistics.
#'     \item{\code{mean, se, simpleCI, bcCI, median, point}} Summary
#'      statistics for Mantel correlation bootstraps.
#'    }
#'   }
#'   \item{\code{input}}{List containing the inputs to \code{estMantel}}
#' }
#' @export
#'
#' @examples
#' data('OVENdata')
#' rM1 <- estMantel(isGL=OVENdata$isGL,#Logical vector: light-level GL(T)/GPS(F)
#'                  geoBias = OVENdata$geo.bias, # Geolocator location bias
#'                  geoVCov = OVENdata$geo.vcov, # Location covariance matrix
#'                  targetSites = OVENdata$targetSites, # Non-breeding target sites
#'                  originPoints = OVENdata$originPoints, # Capture Locations
#'                  targetPoints = OVENdata$targetPoints, # Target locations
#'                  verbose = 1,   # output options
#'                  nBoot = 100, # This is set low for example
#'                  resampleProjection = sf::st_crs(OVENdata$targetSites))
#' rM1
#' str(rM1, max.level = 2)
#' @seealso \code{\link{estMC}}
#'
#' @references
#' Cohen, E. B., J. A. Hostetler, M. T. Hallworth, C. S. Rushing, T. S. Sillett,
#' and P. P. Marra. 2018. Quantifying the strength of migratory connectivity.
#' Methods in Ecology and Evolution 9: 513 - 524.
#' \href{http://doi.org/10.1111/2041-210X.12916}{doi:10.1111/2041-210X.12916}

estMantel <- function(targetPoints = NULL, originPoints = NULL, isGL,
                      geoBias = NULL, geoVCov = NULL, targetSites = NULL,
                      nBoot = 1000,
                      nSim = ifelse(any(isRaster & isGL), 5000,
                                    ifelse(any(isGL), 1000,
                                           ifelse(any(isRaster), 10, 1))),
                      verbose=0, alpha = 0.05,
                      resampleProjection = 'ESRI:102010',
                      maxTries = 300, maintainLegacyOutput = FALSE,
                      originSites = NULL, isTelemetry = !isGL, isRaster = FALSE,
                      captured = "origin", geoBiasOrigin = geoBias,
                      geoVCovOrigin = geoVCov,
                      targetRaster = NULL, originRaster = NULL,
                      dataOverlapSetting = c("dummy", "none", "named"),
                      originRelAbund = NULL, targetRelAbund = NULL) {
  dataOverlapSetting <- match.arg(dataOverlapSetting)
  # double check that spatial data coming in is in sf format #
  if (inherits(targetPoints, "SpatialPoints"))
    targetPoints <- sf::st_as_sf(targetPoints)
  if (inherits(originPoints, "SpatialPoints"))
    originPoints <- sf::st_as_sf(originPoints)

  # Input checking and assignment
  if (!(verbose %in% 0:3))
    stop("verbose should be integer 0-3 for level of output during bootstrap: 0 = none, 1 = every 10, 2 = every run, 3 = every animal")
  if (length(geoBias)!=2 && any(isGL) && any(captured != "target"))
    stop("geoBias should be vector of length 2 (expected bias in longitude and latitude of targetPoints, in meters)")
  if (!isTRUE(all.equal(dim(geoVCov), c(2, 2), check.attributes = FALSE)) &&
      any(isGL) && any(captured != "target"))
    stop("geoVCov should be 2x2 matrix (expected variance/covariance in longitude and latitude of targetPoints, in meters)")

  targetStats <- assignRasterStats(targetRaster)
  targetPointsAssigned <- targetStats$PointsAssigned
  targetSingleCell <- targetStats$SingleCell
  targetRasterXYZ <- targetStats$RasterXYZ
  targetRasterXYZcrs <- targetStats$RasterXYZcrs

  originStats <- assignRasterStats(originRaster)
  originPointsAssigned <- originStats$PointsAssigned
  originSingleCell <- originStats$SingleCell
  originRasterXYZ <- originStats$RasterXYZ
  originRasterXYZcrs <- originStats$RasterXYZcrs

  if (is.null(targetPoints)) {
    if (any(isGL) || any(isTelemetry) || !all(captured=="origin") ||
        is.null(targetRaster))
      stop("Need to define targetPoints unless all data is raster captured on origin")
  }
  else {
    if(is.na(sf::st_crs(targetPoints))){
      stop('Coordinate system definition needed for targetPoints')
    }
    targetPoints <- sf::st_transform(targetPoints, crs = resampleProjection)
  }

  if (is.null(originPoints)) {
    if (any(isGL) || any(isTelemetry) || any(captured=="origin") ||
        is.null(originRaster))
      stop("Need to define originPoints unless all data is raster captured on target")
  }
  else if (is.na(sf::st_crs(originPoints))) {
    stop('Coordinate system definition needed for originPoints')
  }
  if (dataOverlapSetting != "dummy") {
    temp <- reassignInds(dataOverlapSetting = dataOverlapSetting,
                         originPoints = originPoints,
                         targetPoints = targetPoints,
                         isGL = isGL, isTelemetry = isTelemetry,
                         isRaster = isRaster, isProb = FALSE,
                         captured = captured,
                         originRasterXYZ = originRasterXYZ,
                         #originRasterXYZcrs = originRasterXYZcrs,
                         originSingleCell = originSingleCell,
                         targetRasterXYZ = targetRasterXYZ,
                         #targetRasterXYZcrs = targetRasterXYZcrs,
                         targetSingleCell = targetSingleCell)
    originPoints <- temp$originPoints; targetPoints <- temp$targetPoints
    # originAssignment <- temp$originAssignment
    # targetAssignment <- temp$targetAssignment
    isGL <- temp$isGL; isTelemetry <- temp$isTelemetry
    isRaster <- temp$isRaster; isProb <- temp$isProb
    originRasterXYZ <- temp$originRasterXYZ
    originSingleCell <- temp$originSingleCell
    targetRasterXYZ <- temp$targetRasterXYZ
    targetSingleCell <- temp$targetSingleCell
  }

  nAnimals <- max(nrow(targetPoints), nrow(originPoints), length(isGL),
                  length(isTelemetry), length(isRaster),
                  ifelse(is.null(targetRaster), 0,
                         ifelse(targetPointsAssigned, dim(targetSingleCell)[3],
                                dim(targetRasterXYZ)[2] - 2)),
                  ifelse(is.null(originRaster), 0,
                         ifelse(originPointsAssigned, dim(originSingleCell)[3],
                                dim(originRasterXYZ)[2] - 2)),
                  length(captured))

  if (length(isGL)==1)
    isGL <- rep(isGL, nAnimals)
  if (length(isTelemetry)==1)
    isTelemetry <- rep(isTelemetry, nAnimals)
  if (length(isRaster)==1)
    isRaster <- rep(isRaster, nAnimals)
  if (length(captured)==1)
    captured <- rep(captured, nAnimals)

  if(!is.null(targetSites)){
    #if sp object covert to sf #
    if("SpatialPolygonsDataFrame" %in% class(targetSites)){
      targetSites <- sf::st_as_sf(targetSites)
    }

    if(is.na(sf::st_crs(targetSites))){
      stop('Coordinate system definition needed for targetSites')
    }

    targetSites <- sf::st_transform(targetSites, resampleProjection)
  }
  targetPointsInSites <- FALSE

  if (targetPointsAssigned && !is.null(targetSites) && any(isRaster)) {
    if (verbose > 0){
      cat('Checking if single cell target points in targetSites, may take a moment\n')}
    targetPointSample2 <- apply(targetSingleCell,
                                FUN = function(x)
                                  sf::st_as_sf(data.frame(x),
                                               coords = c("Longitude",
                                                          "Latitude"),
                                               crs = 4326),
                                MARGIN = 3)
    if(!sf::st_crs(targetSites)==sf::st_crs(targetPointSample2[[1]])){
      targetPointSample2 <- sapply(targetPointSample2, sf::st_transform,
                                   crs = resampleProjection)
    }

    targetCon <- sapply(targetPointSample2, FUN = function(z){
      suppressMessages(as.numeric(unclass(sf::st_intersects(x = z,
                                                            y = targetSites,
                                                            sparse = TRUE))))})

    if (!any(is.na(targetCon)))
      targetPointsInSites <- TRUE
    else if (verbose > 0)
      cat('Single cell target points supplied, but some points (proportion',
          format(sum(is.na(targetCon))/length(targetCon), digits = 2),
          ') not in targetSites\n')
  }
  else {
    targetCon <- NULL
  }

  originPointsInSites <- FALSE
  if (originPointsAssigned && !is.null(originSites) && any(isRaster)) {
    if (verbose > 0)
      cat('Checking if single cell origin points in originSites, may take a moment\n')
    nSamples <- dim(originSingleCell)[1]
    originPointSample2 <- apply(originSingleCell,
                                FUN = function(x){sf::st_as_sf(data.frame(x),
                                                               coords = c("Longitude", "Latitude"),
                                                               crs = 4326)},
                                MARGIN = 3)
    if(!sf::st_crs(originSites)==sf::st_crs(originPointSample2[[1]])){
      originPointSample2 <- sapply(originPointSample2, sf::st_transform,
                                   crs = resampleProjection)
    }

    originCon <- sapply(originPointSample2, FUN = function(z){
      suppressMessages(as.numeric(unclass(sf::st_intersects(x = z, y = originSites,
                                                            sparse = TRUE))))})

    if (!any(is.na(originCon)))
      originPointsInSites <- TRUE
    else if (verbose > 0){
      cat('Single cell origin points supplied, but some points (proportion',
          sum(is.na(originCon))/length(originCon), ') not in originSites\n')
    }
  }else{
    originCon <- NULL
  }

  weights <- array(0, c(nBoot, nAnimals))
  # if (is.null(originRelAbund) && any(captured!="target") ||
  #     is.null(targetRelAbund) && any(captured!="origin")) {
  #   warning("rM can be biased if you don't weight by abundance")
  # }
  if (is.null(originSites) && !is.null(originRelAbund)) {
    stop("If you want to set origin relative abundances, you need to define the origin sites")
  }
  if (is.null(targetSites) && !is.null(targetRelAbund)) {
    stop("If you want to set target relative abundances, you need to define the target sites")
  }
  if (is.null(originRelAbund) && is.null(targetRelAbund))
    weights[,] <- 1/nAnimals
  else if (!is.null(originRelAbund) && any(captured!="target")) {
    nOriginSites <- nrow(originSites)
    if (coda::is.mcmc(originRelAbund) || coda::is.mcmc.list(originRelAbund)) {
      originRelAbund <- as.matrix(originRelAbund)
    }
    if (is.matrix(originRelAbund) && dim(originRelAbund)[1]>1) {
      abundFixed <- FALSE
      if (dim(originRelAbund)[2]>nOriginSites)
        abundParams <- paste('relN[', 1:nOriginSites, ']', sep='')
      else if (dim(originRelAbund)[2]==nOriginSites)
        abundParams <- 1:nOriginSites
      else
        stop('Number of origin sites must be constant between sites and abundance')
      if (dim(originRelAbund)[1] >= nBoot)
        abundRows <- round(seq(from = 1, to = dim(originRelAbund)[1],
                               length.out = nBoot))
      else
        abundRows <- sample.int(n = dim(originRelAbund)[1], replace = TRUE,
                                size = nBoot)
      originRelAbund <- as.matrix(originRelAbund[abundRows, abundParams])
      abundBase <- colMeans(originRelAbund)
    }
    else {
      abundFixed <- TRUE
      if (length(originRelAbund)!=nOriginSites)
        stop('Number of origin sites must be constant between sites and abundance')
      abundBase <- originRelAbund
      originRelAbund <- matrix(originRelAbund, nBoot, nOriginSites, TRUE)
    }
    if (verbose > 0)
      cat("Creating origin assignment\n")
    # if geolocator, telemetry and captured in origin then simply get the origin site
    if (!is.null(originPoints)){
      if(!identical(sf::st_crs(originPoints),sf::st_crs(originSites))){
        # project if needed
        originSites <- sf::st_transform(originSites, sf::st_crs(originPoints))
      }
      originAssignment <- suppressMessages(unclass(sf::st_intersects(x = originPoints,
                                                                     y = originSites,
                                                                     sparse = TRUE)))
    }
    else
      originAssignment <- NULL
    if (!is.null(originAssignment)) {
      originAssignment[lengths(originAssignment)==0] <- NA
      if (any(lengths(originAssignment)>1)){
        warning("Overlapping originSites may cause issues\n")
        originAssignment <- lapply(originAssignment, function (x) x[1])
      }
      originAssignment <- array(unlist(originAssignment))
      duds <- is.na(originAssignment) & captured[1:nAnimals] != "target"
      if (any(duds)){
        if (verbose > 0)
          cat("Not all origin capture locations are within originSites. Assigning to closest site\n")
        warning("Not all origin capture locations are within originSites. Assigning to closest site.\n",
                "Affects animals: ", paste(which(duds), collapse = ","))
        originAssignment[duds] <-
          sf::st_nearest_feature(x = originPoints[duds,],
                                 y = originSites)

      }
      nOriginAnimals <- rep(NA, nOriginSites)
      for (i in 1:nOriginSites) {
        nOriginAnimals[i] <- sum(originAssignment==i & captured!="target", na.rm = TRUE)
        if (nOriginAnimals[i] > 0)
          weights[ , originAssignment==i] <- originRelAbund[, i]/nOriginAnimals[i]
      }
    }
    else{
      warning("Can't assign animals to origin sites; using unweighted sampling",
              immediate. = TRUE)
      weights[ , captured!="target"] <- 1/sum(captured!="target")
    }
  }
  if (!is.null(targetRelAbund) && any(captured=="target")) {
    nTargetSites <- nrow(targetSites)
    if (coda::is.mcmc(targetRelAbund) || coda::is.mcmc.list(targetRelAbund)) {
      targetRelAbund <- as.matrix(targetRelAbund)
    }
    if (is.matrix(targetRelAbund) && dim(targetRelAbund)[1]>1) {
      abundFixedT <- FALSE
      if (dim(targetRelAbund)[2]>nTargetSites)
        abundParams <- paste('relN[', 1:nTargetSites, ']', sep='')
      else if (dim(targetRelAbund)[2]==nTargetSites)
        abundParams <- 1:nTargetSites
      else
        stop('Number of target sites must be constant between sites and abundance')
      if (dim(targetRelAbund)[1] >= nBoot)
        abundRows <- round(seq(from = 1, to = dim(targetRelAbund)[1],
                               length.out = nBoot))
      else
        abundRows <- sample.int(n = dim(targetRelAbund)[1], replace = TRUE,
                                size = nBoot)
      targetRelAbund <- as.matrix(targetRelAbund[abundRows, abundParams])
    }
    else {
      abundFixedT <- TRUE
      if (length(targetRelAbund)!=nTargetSites)
        stop('Number of target sites must be constant between sites and abundance')
      targetRelAbund <- matrix(targetRelAbund, nBoot, nTargetSites, TRUE)
    }
    if (verbose > 0)
      cat("Creating target assignment\n")
    # if geolocator, telemetry and captured in target then simply get the target site
    if (!is.null(targetPoints)){
      if(!identical(sf::st_crs(targetPoints),sf::st_crs(targetSites))){
        # project if needed
        targetSites <- sf::st_transform(targetSites, sf::st_crs(targetPoints))
      }
      targetAssignment <- suppressMessages(unclass(sf::st_intersects(x = targetPoints,
                                                                     y = targetSites,
                                                                     sparse = TRUE)))
    }
    else
      targetAssignment <- NULL
    if (!is.null(targetAssignment)) {
      targetAssignment[lengths(targetAssignment)==0] <- NA
      if (any(lengths(targetAssignment)>1)){
        warning("Overlapping targetSites may cause issues\n")
        targetAssignment <- lapply(targetAssignment, function (x) x[1])
      }
      targetAssignment <- array(unlist(targetAssignment))
      duds <- is.na(targetAssignment) & captured[1:nAnimals] == "target"
      if (any(duds)){
        if (verbose > 0)
          cat("Not all target capture locations are within targetSites. Assigning to closest site\n")
        warning("Not all target capture locations are within targetSites. Assigning to closest site.\n",
                "Affects animals: ", paste(which(duds), collapse = ","))
        targetAssignment[duds] <-
          sf::st_nearest_feature(x = targetPoints[duds,],
                                 y = targetSites)

      }
      nTargetAnimals <- rep(NA, nTargetSites)
      for (i in 1:nTargetSites) {
        nTargetAnimals[i] <- sum(targetAssignment==i & captured=="target", na.rm = TRUE)
        if (nTargetAnimals[i] > 0)
          weights[ , targetAssignment==i] <- targetRelAbund[, i]/nTargetAnimals[i]
      }
    }
    else{
      warning("Can't assign animals to targetsites; using unweighted sampling",
              immediate. = TRUE)
      weights[ , captured=="target"] <- 1/sum(captured=="target")
    }
  }

  corr <- rep(NA, nBoot)
  if (!is.null(targetPoints) && !is.null(originPoints)){
    pointMantel <- calcMantel(targetPoints = targetPoints,
                              originPoints = originPoints)
    targetDistStart <- pointMantel$targetDist
    originDistStart <- pointMantel$originDist
    pointCorr <- pointMantel$pointCorr
  }
  else {
    pointCorr <- NA
    if (!is.null(targetPoints)){
      targetPoints2 <- sf::st_transform(targetPoints,4326)
      targetDistStart <- distFromPos(sf::st_coordinates(targetPoints2))
    }
    else { # Can't both be null
      originPoints2 <- sf::st_transform(originPoints,4326)

      originDistStart <- distFromPos(sf::st_coordinates(originPoints2))
    }
  }

  for (boot in 1:nBoot) {
    if (verbose > 1 || verbose == 1 && boot %% 100 == 0)
      cat("Bootstrap Run", boot, "of", nBoot, "at", date(), "\n")
    # Sample individual animals with replacement
    animal.sample <- sample.int(nAnimals, replace=TRUE, prob = weights[boot, ])
    if (any(captured[animal.sample]!='origin')) {
      oSamp <- locSample(isGL = (isGL[animal.sample] & captured[animal.sample]!='origin'),
                         isRaster = (isRaster[animal.sample] & captured[animal.sample]!='origin'),
                         isProb = rep(FALSE, nAnimals),
                         isTelemetry = (isTelemetry[animal.sample] | captured[animal.sample]=='origin'),
                         geoBias = geoBiasOrigin,
                         geoVCov = geoVCovOrigin,
                         points = originPoints[animal.sample, ],
                         matvals = originRasterXYZ[, c(1:2, animal.sample + 2)],
                         matvals_crs = originRasterXYZcrs,
                         singleCell = originSingleCell[,,animal.sample],
                         overlap1 = originCon[,animal.sample],
                         pointsInSites = originPointsInSites,
                         assignment = NULL,
                         sites = originSites,
                         resampleProjection = resampleProjection,
                         nSim = nSim,
                         maxTries = maxTries)
      if (!is.null(oSamp$notfind)) {
        oSamp$notfind$Animal <- animal.sample[oSamp$notfind$Animal]
        notfind <- unique(oSamp$notfind)
        stop('maxTries (',maxTries,') reached during origin location sampling, exiting. ',
             'Animal(s) where location sampling failed to fall in sites:\n',
             paste(utils::capture.output(print(notfind, row.names = FALSE)), collapse = "\n"),
             '\nExamine originSites',
             ifelse(any(notfind$isGL),
                    ', geoBiasOrigin, geoVcovOrigin, originPoints', ''),
             ifelse(any(notfind$isRaster), ', originRaster', ''),
             ifelse(any(notfind$isTelemetry), ', originPoints/captured', ''),
             ', and resampleProjection to determine why sampled points fell outside sites.')
      }
      origin.point.sample <- oSamp$point.sample
      origin.point.sample <- sf::st_as_sf(data.frame(origin.point.sample),
                                          coords = c("x","y"),
                                          crs = resampleProjection)
      if (verbose > 2)
        cat(' ', oSamp$draws, 'origin draw(s) (of length', nSim, 'and of', maxTries, 'possible).\n')
    }
    else {
      # Get origin point for each animal sampled
      origin.point.sample <- originPoints[animal.sample, ]
    }

    if (any(captured[animal.sample]!="target")) {
      tSamp <- locSample(isGL = (isGL[animal.sample] &
                                   captured[animal.sample] != "target"),
                         isRaster = (isRaster[animal.sample] &
                                       captured[animal.sample] != "target"),
                         isProb = rep(FALSE, nAnimals),
                         isTelemetry = (isTelemetry[animal.sample] |
                                          captured[animal.sample] == "target"),
                         geoBias = geoBias, geoVCov = geoVCov,
                         points = targetPoints[animal.sample, ],
                         matvals = targetRasterXYZ[, c(1:2, animal.sample + 2)],
                         matvals_crs = targetRasterXYZcrs,
                         singleCell = targetSingleCell[,,animal.sample],
                         pointsInSites = targetPointsInSites,
                         overlap1 = targetCon[, animal.sample],
                         sites = targetSites,
                         resampleProjection = resampleProjection, nSim = nSim,
                         maxTries = maxTries)
      if (!is.null(tSamp$notfind)) {
        tSamp$notfind$Animal <- animal.sample[tSamp$notfind$Animal]
        notfind <- unique(tSamp$notfind)
        stop('maxTries (',maxTries,') reached during target location sampling, exiting. ',
             'Animal(s) where location sampling failed to fall in sites:\n',
             paste(utils::capture.output(print(notfind, row.names = FALSE)), collapse = "\n"),
             '\nExamine targetSites',
             ifelse(any(notfind$isGL),
                    ', geoBiasOrigin, geoVcovOrigin, targetPoints', ''),
             ifelse(any(notfind$isRaster), ', targetRaster', ''),
             ifelse(any(notfind$isTelemetry), ', targetPoints/captured', ''),
             ', and resampleProjection to determine why sampled points fell outside sites.')
      }
      target.point.sample <- tSamp$point.sample
      target.point.sample <- sf::st_as_sf(data.frame(target.point.sample),
                                          coords = c("x","y"),
                                          crs = resampleProjection)
      if (verbose > 2){
        cat(' ', tSamp$draws, 'target draw(s) (of length', nSim, 'and of', maxTries, 'possible).\n')
      }
    }
    else {
      # Get target point for each animal sampled
      target.point.sample <- targetPoints[animal.sample, ]
    }

    if (all(captured[animal.sample]=="target")) {
      targetDist1 <- targetDistStart[animal.sample, animal.sample]
      corr[boot] <- calcMantel(targetDist = targetDist1,
                               originPoints = origin.point.sample)$pointCorr
    }
    else if (all(captured[animal.sample]=="origin")) {
      originDist1 <- originDistStart[animal.sample, animal.sample]
      corr[boot] <- calcMantel(originDist = originDist1,
                               targetPoints = target.point.sample)$pointCorr
    }
    else {
      corr[boot] <- calcMantel(originPoints = origin.point.sample,
                               targetPoints = target.point.sample)$pointCorr
    }

    if (verbose > 1 || verbose == 1 && boot %% 10 == 0)
      cat(" Correlation mean:", mean(corr, na.rm=TRUE), "SD:", sd(corr, na.rm=TRUE),
          "low quantile:", quantile(corr, alpha/2, na.rm=TRUE),
          "high quantile:", quantile(corr, 1-alpha/2, na.rm=TRUE), "\n")
  }
  meanCorr <- mean(corr, na.rm=TRUE)
  medianCorr <- median(corr, na.rm=TRUE)
  seCorr <- sd(corr, na.rm=TRUE)
  simpleCICorr <- quantile(corr, c(alpha/2, 1-alpha/2), na.rm=TRUE, type = 8,
                           names = FALSE)
  corr.z0 <- qnorm(sum((corr)<meanCorr)/nBoot)
  bcCICorr <- quantile(corr, pnorm(2*corr.z0+qnorm(c(alpha/2, 1-alpha/2))),
                       na.rm=TRUE, names = FALSE)
  if (maintainLegacyOutput)
    return(structure(list(sampleCorr = corr, pointCorr = pointCorr,
                          meanCorr = meanCorr, medianCorr = medianCorr,
                          seCorr=seCorr, simpleCICorr=simpleCICorr,
                          bcCICorr=bcCICorr, alpha = alpha,
                          corr = list(sample = corr, mean = meanCorr,
                                      se = seCorr, simpleCI = simpleCICorr,
                                      bcCI = bcCICorr, median = medianCorr,
                                      point = pointCorr),
                          input = list(targetPoints = targetPoints,
                                       originPoints = originPoints, isGL = isGL,
                                       geoBias = geoBias, geoVCov = geoVCov,
                                       targetSites = targetSites,
                                       targetSites = targetSites, nBoot = nBoot,
                                       nSim = nSim, verbose = verbose,
                                       alpha = alpha,
                                       resampleProjection = resampleProjection,
                                       maxTries = maxTries,
                                       maintainLegacyOutput = TRUE)),
                     class = c("estMantel", "estMigConnectivity")))
  else
    return(structure(list(corr = list(sample = corr, mean = meanCorr,
                                      se = seCorr, simpleCI = simpleCICorr,
                                      bcCI = bcCICorr, median = medianCorr,
                                      point = pointCorr),
                          input = list(targetPoints = targetPoints,
                                       originPoints = originPoints, isGL = isGL,
                                       geoBias = geoBias, geoVCov = geoVCov,
                                       targetSites = targetSites, nBoot = nBoot,
                                       nSim = nSim, verbose = verbose,
                                       alpha = alpha,
                                       resampleProjection = resampleProjection,
                                       maxTries = maxTries,
                                       maintainLegacyOutput = FALSE,
                                       originSites = originSites,
                                       isTelemetry = isTelemetry,
                                       isRaster = isRaster, captured = captured,
                                       geoBiasOrigin = geoBiasOrigin,
                                       geoVCovOrigin = geoVCovOrigin,
                                       targetRaster = targetRaster,
                                       originRaster = originRaster,
                                       dataOverlapSetting = dataOverlapSetting,
                                       originRelAbund = originRelAbund,
                                       targetRelAbund = targetRelAbund)),
                     class = c("estMantel", "estMigConnectivity")))
}

###############################################################################
#' Grab (from https://github.com/SMBC-NZP/MigConnectivity) example RMark
#' transition probability estimates obtained from simulated data
#'
#' Get a dataset containing RMark transition probability estimates from
#' simulated mark-recapture-recovery data from Cohen et al. (2014).  These all
#' represent the intermediate scenario for all settings (moderate connectivity,
#' low re-encounter, 100,000 banded in each breeding area).  Each estimate can
#' be used in \code{estMC} function to estimate MC with uncertainty.  Requires
#' internet connection.
#'
#' @param number Integer 1 - 100, which simulation and RMark estimate you want
#'
#' @return RMark object
#' @export
#' @seealso \code{\link{estMC}}
getCMRexample <- function(number = 1) {
  obj.name <- paste0('psiB.enc2.band100.', number)
  file.name <- paste0('out_', obj.name, '.rds')
  url1 <- paste0('https://github.com/SMBC-NZP/MigConnectivity/blob/master/data-raw/', file.name, '?raw=true')
  temp <- paste(tempdir(), file.name, sep = '/')
  utils::download.file(url1, temp, mode = 'wb')
  fm <- readRDS(temp)
  unlink(temp)
  return(fm)
}


#' Pairwise differences between two or more independent MC estimates
#'
#' Estimates mean (and median) differences in MC, and includes measures of
#' uncertainty (SE and CI).  For those measures of uncertainty to be accurate,
#' only apply this function to MC estimates where all data sources are
#' independent (e.g., different species).
#'
#' @param estimates List of at leat two MC estimates, provided by the estMC
#'    function. If this is a named list (recommended), the function will use
#'    these names in labeling the differences.
#' @param nSamples A positive integer, number of samples (with replacement)
#'    to draw from each pair of MC estimates (default 100000).  If set to NULL,
#'    compares all MC samples from each pair.
#' @param alpha Level for confidence/credible intervals provided.
#' @param returnSamples Should the function return all the sampled differences?
#'    Defaults to FALSE to reduce storage requirements. Change to TRUE to
#'    compute your own summary statistics.
#'
#' @return  \code{diffMC} returns a list with elements:
#' \describe{
#'   \item{\code{meanDiff, medianDiff}}{Vectors with mean and medians of sampled
#'      differences for each pairwise comparison. Estimates of difference
#'      between MC values incorporating parametric uncertainty.}
#'   \item{\code{seDiff}}{Vector with standard errors of MC differences for each
#'      pairwise comparison, estimated from SD of sampled differences.}
#'   \item{\code{simpleCI}}{Matrix of \code{1 - alpha} confidence intervals for
#'      MC differences, estimated as \code{alpha/2} and \code{1 - alpha/2}
#'      quantiles of \code{sampleMC}.}
#'   \item{\code{bcCI}}{Matrix of bias-corrected \code{1 - alpha} confidence
#'      intervals for MC differences for each pairwise comparison. Preferable
#'      to \code{simpleCI} when \code{meanDiff} is the best estimate of the MC
#'      difference. \code{simpleCI} is preferred when
#'      \code{medianDiff} is a better estimator. When \code{meanDiff==medianDiff},
#'      these should be identical.  Estimated as the
#'      \code{pnorm(2 * z0 + qnorm(alpha / 2))} and
#'      \code{pnorm(2 * z0 + qnorm(1 - alpha / 2))} quantiles of sampled
#'      differences, where z0 is the proportion of \code{sampleDiff < meanDiff}.}
#   \item{\code{hpdCI}}{Matrix of \code{1 - alpha} credible intervals for MC
#      differences for each pairwise comparison, etimated using the highest
#      posterior density (HPD) method.}
#'   \item{\code{sampleDiff}}{Only provided if \code{returnSamples} is TRUE.
#'      List of sampled values for each pairwise MC difference.}
#' }
#' @export
#'
#' @references
#' Cohen, E. B., C. S. Rushing, F. R. Moore, M. T. Hallworth, J. A. Hostetler,
#' M. Gutierrez Ramirez, and P. P. Marra. 2019. The strength of
#' migratory connectivity for birds en route to breeding through the Gulf of Mexico.
#'
#' @examples
#' \dontrun{
#' data('OVENdata')
#' ovenPsi <- estTransition(isGL = OVENdata$isGL, #Logical vector:light-level GL(T)
#'                  isTelemetry = !OVENdata$isGL,
#'                  geoBias = OVENdata$geo.bias, # Light-level GL location bias
#'                  geoVCov = OVENdata$geo.vcov, # Location covariance matrix
#'                  targetSites = OVENdata$targetSites, # Non-breeding target sites
#'                  originSites = OVENdata$originSites, # Breeding origin sites
#'                  originPoints = OVENdata$originPoints, # Capture Locations
#'                  targetPoints = OVENdata$targetPoints, # Device target locations
#'                  verbose = 1,   # output options
#'                  nSamples = 1000, # This is set low for example
#'                  resampleProjection = sf::st_crs(OVENdata$targetSites))
#' ovenEst <- estStrength(targetDist = OVENdata$targetDist, # targetSites distance matrix
#'                  originDist = OVENdata$originDist, # originSites distance matrix
#'                  originRelAbund = OVENdata$originRelAbund,#Origin relative abund
#'                  psi = ovenPsi,
#'                  verbose = 1,   # output options
#'                  nSamples = 1000)
#' fm <- getCMRexample()
#' originPos13 <- matrix(c(rep(seq(-99, -81, 2), each = 10),
#'                         rep(seq(49, 31, -2), 10)), 100, 2)
#' targetPos13 <- matrix(c(rep(seq(-79, -61, 2), each = 10),
#'                         rep(seq(9, -9, -2), 10)), 100, 2)
#' originPosCMR <- rowsum(originPos13, c(rep(1:2, 5, each = 5),
#'                                       rep(3:4, 5, each = 5))) / 25
#' targetPosCMR <- rowsum(targetPos13, c(rep(1:2, 5, each = 5),
#'                                       rep(3:4, 5, each = 5))) / 25
#' originDist <- distFromPos(originPosCMR, 'ellipsoid')
#' targetDist <- distFromPos(targetPosCMR, 'ellipsoid')
#' originRelAbundTrue <- rep(0.25, 4)

#' theorEst <- estStrength(originRelAbund = originRelAbundTrue, psi = fm,
#'                   originDist = originDist, targetDist = targetDist,
#'                   originSites = 5:8, targetSites = c(3,2,1,4),
#'                   nSamples = 1000, verbose = 0,
#'                   sampleSize = length(grep("[2-5]", fm$data$data$ch)))
#' ovenEst
#' theorEst
#' diff1 <- diffMC(estimates = list(Ovenbird = ovenEst, Theorybird = theorEst),
#'                 nSamples = 10000, returnSamples = TRUE)
#'
#'}
diffMC <- function(estimates, nSamples = 100000, alpha = 0.05,
                   returnSamples = FALSE) {
  nEst <- length(estimates)
  nComparisons <- choose(nEst, 2)
  for (i in 1:nEst)
    if (is.null(estimates[[i]]$MC))
      estimates[[i]]$MC <- list(sample = estimates[[i]]$sampleMC)
  nSamplesEst <- sapply(estimates, function(x) length(x$MC$sample))
  diffSamples <- vector('list', nComparisons)
  if (is.null(names(estimates)))
    names(estimates) <- 1:nEst
  comparisons <- matrix(c(sequence(1:(nEst - 1)), rep(2:nEst, 1:(nEst - 1))),
                        nComparisons, 2)
  for (i in 1:nComparisons) {
    if (is.null(nSamples))
      diffSamples[[i]] <- rep(estimates[[comparisons[i, 1]]]$MC$sample,
                              times = nSamplesEst[comparisons[i, 2]]) -
        rep(estimates[[comparisons[i, 2]]]$MC$sample,
            each = nSamplesEst[comparisons[i, 1]])
    else
      diffSamples[[i]] <- sample(estimates[[comparisons[i, 1]]]$MC$sample,
                                 nSamples, replace = TRUE) -
        sample(estimates[[comparisons[i, 2]]]$MC$sample, nSamples, replace = TRUE)
  }
  meanDiff <- sapply(diffSamples, mean, na.rm=TRUE)
  medianDiff <- sapply(diffSamples, median, na.rm=TRUE)
  seDiff <- sapply(diffSamples, sd, na.rm=TRUE)
  simpleCI <- sapply(diffSamples, quantile, c(alpha/2, 1-alpha/2), na.rm=TRUE,
                     names = FALSE)
  diff.z0 <- sapply(diffSamples,
                    function(MC) qnorm(sum(MC<mean(MC, na.rm = TRUE),
                                           na.rm = TRUE)/length(which(!is.na(MC)))))
  bcCI <- mapply(function(MC, z0)
    quantile(MC, pnorm(2*z0+ qnorm(c(alpha/2, 1-alpha/2))),
                       na.rm=TRUE, names = FALSE), diffSamples, diff.z0)
  names(diffSamples) <- names(meanDiff) <- paste(names(estimates[comparisons[,1]]),
                                                 '-', names(estimates[comparisons[,2]]))
  names(medianDiff) <- names(seDiff) <- names(diffSamples)
  colnames(simpleCI) <- colnames(bcCI) <- names(diffSamples) #colnames(hpdCI) <-
  sampleDiff <- ifelse(returnSamples, diffSamples, NA)
  return(structure(list(meanDiff = meanDiff, medianDiff = medianDiff,
                        seDiff = seDiff, simpleCI = simpleCI, bcCI = bcCI,
                        sampleDiff = sampleDiff, alpha = alpha),
                   class = c('diffMC', 'diffMigConnectivity')))
}

#' Pairwise differences between two or more independent Mantel correlation
#' estimates
#'
#' Estimates mean (and median) differences in Mantel correlations (rM), and
#' includes measures of uncertainty (SE and CI).  For those measures of
#' uncertainty to be accurate, only apply this function to rM estimates where
#' all data sources are independent (e.g., different species).
#'
#' @param estimates List of at leat two Mantel correlation estimates, provided
#'    by either the estMC or the estMantel functions. If this is a named list
#'    (recommended), the function will use these names in labeling the
#'    differences.
#' @param nSamples A positive integer, number of samples (with replacement)
#'    to draw from each pair of MC estimates (default 100000).  If set to NULL,
#'    compares all Mantel correlation samples from each pair.
#' @param alpha Level for confidence/credible intervals provided.
#' @param returnSamples Should the function return all the sampled differences?
#'    Defaults to FALSE to reduce storage requirements. Change to TRUE to
#'    compute your own summary statistics.
#'
#' @return  \code{diffMantel} returns a list with elements:
#' \describe{
#'   \item{\code{meanDiff, medianDiff}}{Vectors with mean and medians of sampled
#'      differences for each pairwise comparison. Estimates of difference
#'      between rM values incorporating parametric uncertainty.}
#'   \item{\code{seDiff}}{Vector with standard errors of rM differences for each
#'      pairwise comparison, estimated from SD of sampled differences.}
#'   \item{\code{simpleCI}}{Matrix of \code{1 - alpha} confidence intervals for
#'      rM differences, estimated as \code{alpha/2} and \code{1 - alpha/2}
#'      quantiles of \code{sampleCorr}.}
#'   \item{\code{bcCI}}{Matrix of bias-corrected \code{1 - alpha} confidence
#'      intervals for rM differences for each pairwise comparison. Preferable
#'      to \code{simpleCI} when \code{meanDiff} is the best estimate of the rM
#'      difference. \code{simpleCI} is preferred when
#'      \code{medianDiff} is a better estimator. When \code{meanDiff==medianDiff},
#'      these should be identical.  Estimated as the
#'      \code{pnorm(2 * z0 + qnorm(alpha / 2))} and
#'      \code{pnorm(2 * z0 + qnorm(1 - alpha / 2))} quantiles of sampled
#'      differences, where z0 is the proportion of \code{sampleDiff < meanDiff}.}
#   \item{\code{hpdCI}}{Matrix of \code{1 - alpha} credible intervals for rM
#      differences for each pairwise comparison, etimated using the highest
#      posterior density (HPD) method.}
#'   \item{\code{sampleDiff}}{Only provided if \code{returnSamples} is TRUE.
#'      List of sampled values for each pairwise rM difference.}
#' }
#' @export
#'
#' @references
#' Cohen, E. B., C. S. Rushing, F. R. Moore, M. T. Hallworth, J. A. Hostetler,
#' M. Gutierrez Ramirez, and P. P. Marra. 2019. The strength of
#' migratory connectivity for birds en route to breeding through the Gulf of Mexico.
#'
# @examples
diffMantel <- function(estimates, nSamples = 100000, alpha = 0.05,
                       returnSamples = FALSE) {
  nEst <- length(estimates)
  nComparisons <- choose(nEst, 2)
  for (i in 1:nEst)
    if (is.null(estimates[[i]]$corr))
      estimates[[i]]$corr <- list(sample = estimates[[i]]$sampleCorr)
  nSamplesEst <- sapply(estimates, function(x) length(x$corr$sample))
  diffSamples <- vector('list', nComparisons)
  if (is.null(names(estimates)))
    names(estimates) <- 1:nEst
  comparisons <- matrix(c(sequence(1:(nEst - 1)), rep(2:nEst, 1:(nEst - 1))),
                        nComparisons, 2)
  for (i in 1:nComparisons) {
    if (is.null(nSamples))
      diffSamples[[i]] <- rep(estimates[[comparisons[i, 1]]]$corr$sample,
                              times = nSamplesEst[comparisons[i, 2]]) -
        rep(estimates[[comparisons[i, 2]]]$corr$sample,
            each = nSamplesEst[comparisons[i, 1]])
    else
      diffSamples[[i]] <- sample(estimates[[comparisons[i, 1]]]$corr$sample,
                                 nSamples, replace = TRUE) -
        sample(estimates[[comparisons[i, 2]]]$corr$sample, nSamples,
               replace = TRUE)
  }
  meanDiff <- sapply(diffSamples, mean, na.rm=TRUE)
  medianDiff <- sapply(diffSamples, median, na.rm=TRUE)
  seDiff <- sapply(diffSamples, sd, na.rm=TRUE)
  simpleCI <- sapply(diffSamples, quantile, c(alpha/2, 1-alpha/2), na.rm=TRUE,
                     names = FALSE)
  diff.z0 <- sapply(diffSamples,
                    function(MC) qnorm(sum(MC<mean(MC, na.rm = TRUE), na.rm = TRUE)/
                                         length(which(!is.na(MC)))))
  bcCI <- mapply(function(MC, z0) quantile(MC,
                                           pnorm(2*z0+qnorm(c(alpha/2, 1-alpha/2))),
                       na.rm=TRUE, names = FALSE), diffSamples, diff.z0)
  names(diffSamples) <- names(meanDiff) <- paste(names(estimates[comparisons[,1]]),
                                                 '-', names(estimates[comparisons[,2]]))
  names(medianDiff) <- names(seDiff) <- names(diffSamples)
  colnames(simpleCI) <- colnames(bcCI) <- names(diffSamples)
  sampleDiff <- ifelse(returnSamples, diffSamples, NA)
  return(structure(list(meanDiff = meanDiff, medianDiff = medianDiff,
                        seDiff = seDiff, simpleCI = simpleCI, bcCI = bcCI,
                        sampleDiff = sampleDiff, alpha = alpha),
                   class = c('diffMantel', 'diffMigConnectivity')))
}

#' @rdname estTransition
#' @export
estPsi <- estTransition

#' @rdname estMantel
#' @export
estCorr <- estMantel

#' @rdname diffMC
#' @export
diffStrength <- diffMC

#' @rdname diffMantel
#' @export
diffCorr <- diffMantel
