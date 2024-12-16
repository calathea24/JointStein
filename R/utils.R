library(parallel)
library(JuliaCall)


pseudoinverse_new = function (m, tol)
{
  julia <- julia_setup()
  julia_assign("sigma", m)
  msvd<- julia_eval("svd(sigma)")

  S <- msvd$S
  U <- msvd$U
  V <- t(msvd$Vt)

  if (length(S) == 0)
  {
    return(
      array(0, dim(m)[2:1])
    )
  }
  else
  {
    return(
      V %*% (1/S * t(U)) ## slow down
    )
  }
}

pcor2cor_new = function(m, tol)
{
  # negate off-diagonal entries, then invert
  m = -m
  diag(m) = -diag(m)
  m = pseudoinverse_new(m, tol=tol)

  # standardize and return
  return(cov2cor(m))
}

myrmvnorm_new = function(n, mean, sigma)
{
  sigma <- (sigma + t(sigma)) / 2
  julia <- julia_setup()
  julia_assign("sigma", sigma)
  ev <- julia_eval("eigen(sigma)")

  tmp = ev$vectors %*% ( t(ev$vectors) * sqrt(ev$values) ) ## slow down
  tmp = matrix(rnorm(n * ncol(sigma)), nrow = n) %*% tmp ## slow down
  tmp = sweep(tmp, 2, mean, "+")

  return(tmp)
}

ggm_simulate_data = function(sample.size, pcor)
{
  mu = rep(0, dim(pcor)[1])

  cor.mat = pcor2cor_new(pcor) ### slow down

  return( myrmvnorm_new(sample.size, mu, cor.mat) )
}


covcal <- function(data) {
  x = as.matrix(data)
  n = nrow(x)
  p = ncol(x)

  # standardization, No zero standard deviation value
  w = rep(1/n, n)
  sw = sqrt(w)
  h1 = 1/(1-sum(w*w)) ## Bias correction factor based on weight
  m = colSums(w*x) ## vector of gene means based on weight
  v = h1*(colSums(w*x^2)-m^2) ## vector of variance Var=E[X^2] - E[X]^2
  sc = sqrt(v)
  x = sweep(x, 2, m, "-")
  x = sweep(x, 2, sc, "/")

  S = crossprod(sweep(x, MARGIN=1, STATS=sw, FUN="*"))

  return(S)
}

mcc.pcor <- function(est, pcor.truth) {
  pcor.est = est[upper.tri(est)]
  pcor.truth = pcor.truth[upper.tri(pcor.truth)]
  TN = length(which(pcor.est[which(pcor.truth==0)]==0))
  TP = length(which(pcor.est[which(pcor.truth!=0)]!=0))
  FP = length(which(pcor.est[which(pcor.truth==0)]!=0))
  FN = length(which(pcor.est[which(pcor.truth!=0)]==0))
  mcc = mcc(TP=TP,TN=TN,FP=FP,FN=FN)
  return(mcc)
}

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}












