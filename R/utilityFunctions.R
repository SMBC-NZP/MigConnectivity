# Function for generating random transition probability matrix (expression of parametric uncertainty)
make.psi.rand <- function(model,origin.set, target.set) {
  vcv1 <- model$results$beta.vcv
  beta1 <- model$results$beta$estimate
  new.beta <- mvrnorm(1, beta1, vcv1)
  full.beta <- rep(NA, length(model$results$beta$estimate))
  new.real <- get.real(model, "Psi", new.beta, design=model$design.matrix, vcv=F, se=F)
  mat <- TransitionMatrix(new.real)[origin.set, target.set]
  return(mat)
}

# Calculates probability matrix based on exponential decline with distance
mlogit.mat <- function(slope, dist) {
  preMat <- exp(-slope/mean(dist)*dist)
  diag(preMat) <- 0
  nr <- nrow(dist)
  nc <- ncol(dist)
  outMat <- matrix(0, nr, nc)
  for (b in 1:nr) {
    outMat[b,] <- preMat[b,]/(1+sum(preMat[b, ]))
    outMat[b,b] <- 1 - sum(outMat[b, ])
  }
  return(outMat)
}

# Calculates row of probability matrix based on exponential decline with distance
mlogit.row <- function(slope, dist, row.num) {
  preMat <- exp(-slope/mean(dist)*dist)
  diag(preMat) <- 0
  out.row <- preMat[row.num,]/(1+sum(preMat[row.num, ]))
  out.row[row.num] <- 1 - sum(out.row)
  return(out.row)
}

# Crude optimizable function for developing MC pattern based on MC strength
MC.mlogit <- function(slope, MC.in, origin.dist, target.dist, origin.rel.abund) {
  nBreeding <- nrow(origin.dist)
  nWintering <- nrow(target.dist)
  psi <- mlogit.mat(slope, origin.dist)
  if (any(psi<0))
    return(5*slope^2)
  MC <- calcMC(origin.dist, target.dist, psi, origin.rel.abund)
  return((MC.in - MC)^2)
}

