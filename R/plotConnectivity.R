# #
# # data(OVENdata) # Ovenbird
# #
# M<-estMC(isGL=OVENdata$isGL, # Logical vector: light-level geolocator(T)/GPS(F)
#          geoBias = OVENdata$geo.bias, # Light-level geolocator location bias
#          geoVCov = OVENdata$geo.vcov, #Light-level geolocator covariance matrix
#          targetDist = OVENdata$targetDist, # Target location distance matrix
#          originDist = OVENdata$originDist, # Origin location distance matrix
#          targetSites = OVENdata$targetSites, # Non-breeding / target sites
#          originSites = OVENdata$originSites, # Breeding / origin sites
#          originPoints = OVENdata$originPoints, # Capture Locations
#          targetPoints = OVENdata$targetPoints, # Target locations from devices
#          originRelAbund = OVENdata$originRelAbund, # Origin relative abundances
#          resampleProjection = raster::projection(OVENdata$targetPoints),
#          verbose = 0,   # output options - see help ??estMC
#          nSamples = 10000) # This is set low for example
#
# nTarget <- length(OVENdata$targetSites)
# nOrigin <- length(OVENdata$originSites)
# meanPsi <- apply(M$samplePsi, 2:3, mean)
# lowPsi <- apply(M$samplePsi, 2:3, quantile, probs = 0.025)
# highPsi <- apply(M$samplePsi, 2:3, quantile, probs = 0.975)
# library(rgeos)
# library(shape)
# library(raster)
# library(maptools)
# library(rgdal)
# data(wrld_simpl)
# wrld_simple<-sp::spTransform(wrld_simpl,raster::crs(OVENdata$targetSites))
# maxWidth <- 100000
# xOffset <- matrix(c(100000, -100000), nOrigin, nTarget)
# yOffset <- matrix(rep(c(250000, 150000, 150000), each = nOrigin), nOrigin, nTarget)
# alpha <- 0.2
#
# png('psi_plot1.png', width = 6, height = 6, units = 'in', res = 1200)
# par(mar=c(0,0,0,0))
# plot(wrld_simple,xlim=c(raster::extent(OVENdata$targetSites)[1],
#                         raster::extent(OVENdata$targetSites)[2]),
#      ylim=c(raster::extent(OVENdata$targetSites)[3],
#             raster::extent(OVENdata$originSites)[4]))
# plot(OVENdata$originSites[1],add=TRUE,lwd=1.75)
# plot(OVENdata$originSites[2],add=TRUE,lwd=1.75)
# plot(OVENdata$targetSites,add=TRUE,lwd=1.5,col=c("gray70","gray35","gray10"))
#
# legend("topleft",legend=paste("MC =",round(M$meanMC,2), "\u00b1", round(M$seMC,2)),bty="n",cex=1.8,bg="white",xjust=0)
# for (i in 1:nrow(meanPsi)) {
#   for (j in 1:ncol(meanPsi)) {
#     if (highPsi[i,j] > 0) {
#       xO <- gCentroid(OVENdata$originSites[i])@coords[,1]
#       yO <- gCentroid(OVENdata$originSites[i])@coords[,2]
#       xT <- gCentroid(OVENdata$targetSites[j])@coords[,1] + xOffset[i, j]
#       yT <- gCentroid(OVENdata$targetSites[j])@coords[,2] + yOffset[i, j]
#       angle <- atan((yT - yO)/(xT - xO))
#       if (is.nan(angle))
#         angle <- 0
#       if (xT < xO)
#         angle <- angle + pi
#       cosa <- cos(angle)
#       sina <- sin(angle)
#       polygon(c(xO + highPsi[i,j] * sina * maxWidth,
#                 xT + highPsi[i,j] * sina * maxWidth,
#                 xT + lowPsi[i,j] * sina * maxWidth,
#                 xO + lowPsi[i,j] * sina * maxWidth),
#               c(yO - highPsi[i,j] * cosa * maxWidth,
#                 yT - highPsi[i,j] * cosa * maxWidth,
#                 yT - lowPsi[i,j] * cosa * maxWidth,
#                 yO - lowPsi[i,j] * cosa * maxWidth),
#               col = rgb(i/nrow(meanPsi), j/ncol(meanPsi), 1 - (i + j) / (ncol(meanPsi) + nrow(meanPsi)), alpha=alpha), border = NA)
#       polygon(c(xO - highPsi[i,j] * sina * maxWidth,
#                 xT - highPsi[i,j] * sina * maxWidth,
#                 xT - lowPsi[i,j] * sina * maxWidth,
#                 xO - lowPsi[i,j] * sina * maxWidth),
#               c(yO + highPsi[i,j] * cosa * maxWidth,
#                 yT + highPsi[i,j] * cosa * maxWidth,
#                 yT + lowPsi[i,j] * cosa * maxWidth,
#                 yO + lowPsi[i,j] * cosa * maxWidth),
#               col = rgb(i/nrow(meanPsi), j/ncol(meanPsi), 1 - (i + j) / (ncol(meanPsi) + nrow(meanPsi)), alpha=alpha), border = NA)
#       lines(c(xO - meanPsi[i,j] * sina * maxWidth,
#               xT - meanPsi[i,j] * sina * maxWidth),
#             c(yO + meanPsi[i,j] * cosa * maxWidth,
#               yT + meanPsi[i,j] * cosa * maxWidth),
#             col = rgb(i/nrow(meanPsi), j/ncol(meanPsi), 1 - (i + j) / (ncol(meanPsi) + nrow(meanPsi)), alpha=1))
#       lines(c(xO + meanPsi[i,j] * sina * maxWidth,
#               xT + meanPsi[i,j] * sina * maxWidth),
#             c(yO - meanPsi[i,j] * cosa * maxWidth,
#               yT - meanPsi[i,j] * cosa * maxWidth),
#             col = rgb(i/nrow(meanPsi), j/ncol(meanPsi), 1 - (i + j) / (ncol(meanPsi) + nrow(meanPsi)), alpha=1))
#       shape::Arrowhead(xT, yT, angle / pi * 180, arr.width = meanPsi[i,j], arr.length = 1/8,
#                        arr.type = 'curved', npoint = 15,
#                        lcol = rgb(i/nrow(meanPsi), j/ncol(meanPsi), 1 - (i + j) / (ncol(meanPsi) + nrow(meanPsi)), alpha=1), arr.adj = 0)
#     }
#   }
# }
# dev.off()
# shape::Arrows(gCentroid(OVENdata$originSites[1])@coords[,1],
#               gCentroid(OVENdata$originSites[1])@coords[,2],
#               gCentroid(OVENdata$targetSites[2])@coords[,1]+80000,
#               extent(OVENdata$targetSites[2])[4]+150000,
#               arr.length = 0.3,
#               arr.adj = 0.5,
#               arr.lwd = 1,
#               arr.width = 0.4,
#               arr.type = "triangle",
#               lwd=(apply(M$samplePsi[,1,],2,mean)[2]*10),
#               lty=1)
#
# shape::Arrows(gCentroid(OVENdata$originSites[1])@coords[,1],
#               gCentroid(OVENdata$originSites[1])@coords[,2],
#               gCentroid(OVENdata$targetSites[3])@coords[,1],
#               extent(OVENdata$targetSites[3])[4]+150000,
#               arr.length = 0.3,
#               arr.adj = 0.5,
#               arr.lwd = 1,
#               arr.width = 0.4,
#               arr.type = "triangle",
#               lwd=(apply(M$samplePsi[,1,],2,mean)[3]*10),
#               lty=1)
#
# shape::Arrows(gCentroid(OVENdata$originSites[2])@coords[,1],
#               gCentroid(OVENdata$originSites[2])@coords[,2],
#               gCentroid(OVENdata$targetSites[1])@coords[,1],
#               extent(OVENdata$targetSites[1])[4]+150000,
#               arr.length = 0.3,
#               arr.adj = 0.5,
#               arr.lwd = 1,
#               arr.width = 0.4,
#               arr.type = "triangle",
#               lwd=(apply(M$samplePsi[,2,],2,mean)[1]*10),
#               lty=1)
#
# shape::Arrows(gCentroid(OVENdata$originSites[2])@coords[,1],
#               gCentroid(OVENdata$originSites[2])@coords[,2],
#               (gCentroid(OVENdata$targetSites[2])@coords[,1]-80000),
#               extent(OVENdata$targetSites[2])[4]+150000,
#               arr.length = 0.3,
#               arr.adj = 0.5,
#               arr.lwd = 1,
#               arr.width = 0.4,
#               arr.type = "triangle",
#               lwd=(apply(M$samplePsi[,2,],2,mean)[2]*10))
#
# box(which="plot")
# #

# library(ggplot2)
# mc.true <- calcMC(originDist = originDist, targetDist, rep(1/7, 7), psiDeal,
#                   sum(nReleased))
#
# psiDeal.df <- data.frame(Method = "True",
#                       Estimate = "psi",
#                       From = rep(LETTERS[1:nOriginSites], nTargetSites),
#                       To = rep(1:nTargetSites, each = nOriginSites),
#                       Mean = c(psiDeal),
#                       SD = 0,
#                       Lower = c(psiDeal),
#                       Upper = c(psiDeal))
# ests.df <- data.frame(Method = rep(c("DC", "Multinomial"),
#                                    each = nOriginSites*nTargetSites + 1),
#                       Estimate = rep(c("MC", "psi"),
#                                      c(1, nOriginSites*nTargetSites)),
#                       From = c(NA, rep(LETTERS[1:nOriginSites], nTargetSites)),
#                       To = c(NA, rep(1:nTargetSites, each = nOriginSites)),
#                       Mean = c(mc.dc$meanMC, mc.dc$meanPsi,
#                                mc.multi$meanMC,
#                                apply(mc.multi$samplePsi, 2:3, mean)),
#                       SD = c(mc.dc$seMC, apply(mc.dc$samplePsi, 2:3, sd),
#                              mc.multi$seMC, apply(mc.multi$samplePsi, 2:3, sd)),
#                       Lower = c(mc.dc$bcCI[1], mc.dc$bcCIPsi[1,,],
#                                 mc.multi$bcCI[1],
#                                 apply(mc.multi$samplePsi, 2:3, quantile, probs = 0.025)),
#                       Upper = c(mc.dc$bcCI[2], mc.dc$bcCIPsi[2,,],
#                                 mc.multi$bcCI[2],
#                                 apply(mc.multi$samplePsi, 2:3, quantile, probs = 0.975)))
# g.mc <- ggplot(subset(ests.df, Estimate=="MC"),
#                aes(x = Method, y = Mean, ymin = Lower, ymax = Upper)) +
#   geom_hline(yintercept = mc.true, color = "darkgreen") +
#   geom_pointrange() + theme_bw() + ylab("MC")
# g.mc
#
# g.psi <- ggplot(subset(ests.df, Estimate=="psi"),
#              aes(x = Method, y = Mean, ymin = Lower, ymax = Upper)) +
#   geom_hline(data = psiDeal.df, aes(yintercept = Mean), color = "darkgreen") +
#   geom_pointrange() + theme_bw() + ylab("psi") +
#   facet_grid(vars(From), vars(To))
# g.psi
#
# ests.df <- data.frame(Method = "Multinomial",
#                       Estimate = rep(c("MC", "psi"),
#                                      c(1, nOriginSites*nTargetSites)),
#                       From = c(NA, rep(LETTERS[1:nOriginSites], nTargetSites)),
#                       To = c(NA, rep(1:nTargetSites, each = nOriginSites)),
#                       Mean = c(mc.dc$meanMC, mc.dc$meanPsi,
#                                mc.multi$meanMC,
#                                apply(mc.multi$samplePsi, 2:3, mean)),
#                       SD = c(mc.dc$seMC, apply(mc.dc$samplePsi, 2:3, sd),
#                              mc.multi$seMC, apply(mc.multi$samplePsi, 2:3, sd)),
#                       Lower = c(mc.dc$bcCI[1], mc.dc$bcCIPsi[1,,],
#                                 mc.multi$bcCI[1],
#                                 apply(mc.multi$samplePsi, 2:3, quantile, probs = 0.025)),
#                       Upper = c(mc.dc$bcCI[2], mc.dc$bcCIPsi[2,,],
#                                 mc.multi$bcCI[2],
#                                 apply(mc.multi$samplePsi, 2:3, quantile, probs = 0.975)))

plot.estMigConnectivity <- function(x, plot.which = c("MC", "psi", "rM"),
                                    point = c("mean", "median", "point"),
                                    bars = c("bcCI", "simpleCI", "se"),
                                    ylab = plot.which,
                                    originNames = NULL, targetNames = NULL,
                                    map = FALSE) {
  if (map) {
    warning("Map plotting not yet available")
  }
  plot.which <- match.arg(plot.which)
  point <- match.arg(point)
  bars <- match.arg(bars)
  if (is.null(originNames)) {
    if (is.null(dimnames(x$samplePsi)[2])) {
      originNames <- LETTERS[1:dim(x$samplePsi)[2]]
    }
    else {
      originNames <- dimnames(x$samplePsi)[2]
    }
  }
  if (is.null(targetNames)) {
    if (is.null(dimnames(x$samplePsi)[3])) {
      targetNames <- 1:dim(x$samplePsi)[3]
    }
    else {
      targetNames <- dimnames(x$samplePsi)[3]
    }
  }
  nTargetSites <- length(targetNames)
  nOriginSites <- length(originNames)
  if (plot.which=="MC") {
    y <- ifelse(point == "mean", x$meanMC,
                ifelse(point == "median", x$medianMC, x$pointMC))
    ests.df <- data.frame(y = y,
                          lower = ifelse(bars == "bcCI", x$bcCI[1],
                                         ifelse(bars == "simpleCI", x$simpleCI[1],
                                                y - x$seMC)),
                          upper = ifelse(bars == "bcCI", x$bcCI[2],
                                         ifelse(bars == "simpleCI", x$simpleCI[2],
                                                y + x$seMC)))
  }
  else if (plot.which == "rM") {
    y <- ifelse(point == "mean", x$meanCorr,
                ifelse(point == "median", x$medianCorr, x$pointCorr))
    ests.df <- data.frame(y = y,
                          lower = ifelse(bars == "bcCI", x$bcCICorr[1],
                                         ifelse(bars == "simpleCI",
                                                x$simpleCICorr[1],
                                                y - x$seCorr)),
                          upper = ifelse(bars == "bcCI", x$bcCICorr[2],
                                         ifelse(bars == "simpleCI",
                                                x$simpleCICorr[2],
                                                y + x$seCorr)))
  }
  else if (plot.which == "psi") {
    y <- switch(point,
                mean = apply(x$samplePsi, 2:3, mean),
                median = apply(x$samplePsi, 2:3, median),
                point = pointPsi)
    if (is.null(x$bcCIPsi)) {
      x$bcCIPsi <- array(NA, dim = c(2, nOriginSites, nTargetSites),
                       dimnames = list(NULL, originNames, targetNames))
      for (i in 1:nOriginSites) {
        for (j in 1:nTargetSites) {
          psi.z0 <- qnorm(sum(x$samplePsi[, i, j] < mean(x$samplePsi[, i, j],
                                                         na.rm = T)) /
                            length(which(!is.na(x$samplePsi[, i, j]))))
          x$bcCIPsi[ , i, j] <- quantile(x$samplePsi[, i, j],
                                       pnorm(2 * psi.z0 + qnorm(c(x$alpha/2, 1-x$alpha/2))),
                                       na.rm=TRUE, type = 8, names = F)
        }
      }

    }
    yrange <- switch(bars,
                     bcCI = x$bcCIPsi,
                     simpleCI = apply(x$samplePsi, 2:3, quantile,
                                         probs = c(x$alpha/2, 1-x$alpha/2)),
                     se = aperm(array(c(y - apply(x$samplePsi, 2:3, sd),
                              y + apply(x$samplePsi, 2:3, sd)),
                              c(dim(y), 2)), c(3, 1, 2)))
    ests.df <- data.frame(y = c(y),
                          From = rep(originNames, length(targetNames)),
                          To = rep(targetNames, each = length(originNames)),
                          lower = c(yrange[1,,]),
                          upper = c(yrange[2,,]))
  }
  if (plot.which=="psi") {
    g <- ggplot2::ggplot(ests.df, aes(y = y, ymin = lower, ymax = upper, x = To)) +
      ggplot2::geom_pointrange() + ggplot2::theme_bw() +
      ggplot2::facet_grid(vars(From)) + ggplot2::ylab(ylab)
  }
  else {
    g <- ggplot2::ggplot(ests.df, aes(y = y, ymin = lower, ymax = upper, x = plot.which)) +
      ggplot2::geom_pointrange() + ggplot2::theme_bw() + ggplot2::ylab(ylab)
  }
  g
}
