#' Covariance linear shrinkage using GeneNet algorithm
#'
#' Shrinking sample covariance matrix towards identity matrix using GeneNet algorithm from GeneNet R package
#'
#' @param data n x p matrix which n observations in rows and p features in column
#'
#' @return estimated covariance matrix
#' @import GeneNet
#' @references Schäfer, J., & Strimmer, K. (2005). A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statistical applications in genetics and molecular biology, 4(1)

GeneNet <- function(data) {
  return(cor.shrink(data, verbose = FALSE))
}


#' Covariance linear shrinkage using Fisher2011 algorithm
#'
#' Shrinking sample covariance matrix towards identity matrix using Fisher2011 algorithm
#'
#' @param data n x p matrix which n observations in rows and p features in column
#'
#' @return estimated covariance matrix
#' @references Fisher, T. J., & Sun, X. (2011). Improved Stein-type shrinkage estimators for the high-dimensional multivariate normal covariance matrix. Computational Statistics & Data Analysis, 55(5), 1909-1918.

Fisher2011 <- function(data) {
  x = as.matrix(data)
  n = nrow(x)
  p = ncol(x)

  S.emp = covcal(x)
  TS = diag(p)

  #Shrinkage intensity estimate (equations 11)
  a1 = sum(diag(S.emp))/p
  a2 = (n^2/((n-1)*(n+2)*p))*(sum(diag(S.emp %*% S.emp)) - 1/n*(sum(diag(S.emp))^2))

  lambda = (1/n*a2 + p/n*(a1^2))/((n+1)/n*a2+p/n*(a1^2)-2*a1+1)
  lambda = max(0, min(1, lambda)) ## truncate lambda to [0,1]

  S.est = (1-lambda)*S.emp + lambda*TS
  return(S.est)
}



#' Covariance linear shrinkage using HY2014 algorithm
#'
#' Shrinking sample covariance matrix towards identity matrix using HY2014 algorithm
#'
#' @param data n x p matrix which n observations in rows and p features in column
#'
#' @return estimated covariance matrix
#' @references Himeno, T., & Yamada, T. (2014). Estimations for some functions of covariance matrix in high dimension under non-normality and its applications. Journal of Multivariate Analysis, 130, 27-44.

HY2014 <- function(data) {
  x = as.matrix(data)
  N = nrow(x)
  n = N - 1
  p = ncol(x)

  #Shrinkage intensity estimate based on theorem 1 from T.Himeno, T.Yamada 2014 and function of Lambda shrinking
  #towards identity matrix (page 254, A.Touloumis, 2015)
  S = covcal(x)
  Q = 1/(N-1)*sum((diag(x%*%t(x)))^2)
  const = (N-1)/(N*(N-2)*(N-3))
  trSigma2 = const*((N-1)*(N-2)*sum(diag(S%*%S)) + (sum(diag(S)))^2 - N*Q)
  tr2Sigma = const*(2*sum(diag(S%*%S)) +(n^2 - 3*n + 1)*(sum(diag(S)))^2 - n*Q)
  Y2N = trSigma2
  Y1N = sum(diag(S))
  Y1N.2 = tr2Sigma
  beta2 = Y2N + Y1N.2
  delta2 = N*Y2N + Y1N.2 - (N-1)*(2*Y1N - p)
  lambda = beta2/delta2
  lambda = max(0, min(1, lambda))

  #Shrinking covariance matrix
  TS = diag(p)
  S.est = (1-lambda)*S + lambda*TS

  return(S.est)
}



#' Covariance linear shrinkage using IKS algorithm
#'
#' Shrinking sample covariance matrix towards identity matrix using IKS algorithm
#'
#' @param data n x p matrix which n observations in rows and p features in column
#'
#' @return estimated covariance matrix
#' @references Ikeda, Y., Kubokawa, T., & Srivastava, M. S. (2016). Comparison of linear shrinkage estimators of a large covariance matrix in normal and non-normal distributions. Computational Statistics & Data Analysis, 95, 95-108.

IKS <- function(data) {
  x = as.matrix(data)
  N = nrow(x)
  n = N - 1
  p = ncol(x)

  #Shrinkage intensity estimate (equations 11)
  S = covcal(x)
  Q = 1/(N-1)*sum((diag(x%*%t(x)))^2)
  const = (N-1)/(N*(N-2)*(N-3))

  a2c = const*((N-1)*(N-2)*sum(diag(S%*%S)) + (sum(diag(S)))^2 - N*Q)/p
  a1.2 = (sum(diag(S))/p)^2

  numerator = (sum(diag(S%*%S)))/p - a2c
  denominator = (sum(diag(S%*%S)))/p - a1.2

  lambda = numerator/denominator
  lambda = max(0, min(1, lambda))

  TS = diag(p)
  S.est = (1-lambda)*S + lambda*TS

  return(S.est)
}



#' Two-target linear covariance shrinkage
#'
#' Shrinking sample covariance matrix towards identity and shared matrices between groups.
#'
#' @param data n x p matrix which includes observations from all groups in rows and features in column
#' @param group.idx list of indexes for each group. For instance, if data has 2 groups then group.idx list should have 2 objects and each object is vector of indexes of observations belonging to each group.
#' @param group.select a number of which group is selected for the estimation
#' @param covshrinkf name of standard covariance linear shrinkage function to pre-shrink sample covariance matrices from each group before calculating shared target matrix in two-target linear shrinkage
#'
#' @return estimated covariance matrix of selected group


ttls <- function(data, group.idx, group.select = 1, covshrinkf) {
  library(dplyr)
  G <- length(group.idx) # number of groups
  n <- length(group.idx[[group.select]])

  #1. Check if any gene has all zero counts or sd = 0
  x <- data[group.idx[[group.select]],]
  sd_col <- sapply(1:ncol(data), function(i){
    sd(x[,i])
  })
  p_select <- which(sd_col != 0)
  p <- length(p_select)

  #2. Make sure store information of selected genes
  y <- data[,p_select]
  total_p <- ncol(data)
  rm(data)

  #3. Calculate target matrices and other parameters
  S.emp <- lapply(group.idx, function(idx) covshrinkf(y[idx,]))
  S.g <- covcal(y[group.idx[[group.select]],])
  S.all <- 1/(G-1)*Reduce("+", S.emp[-group.select])
  I <- diag(p)

  #4. Calculate A matrix
  aSS <- S.all - S.g
  aIS <- I - S.g
  a11 <- sum(diag(aSS%*%t(aSS)))
  a22 <- sum(diag(aIS%*%t(aIS)))
  a12 <- sum(diag(t(aIS)%*%aSS))
  a21 <- sum(diag(t(aSS)%*%aIS))
  A <- matrix(c(a11,a12,a21,a22), nrow = 2, ncol = 2, byrow = TRUE)

  #5. Calculate b vector
  xs <- wt.scale(y[group.idx[[group.select]],])
  const1 <- n/((n-1)^2*(n-2))
  VS <- const1*(sum(diag(xs%*%t(xs))^2) - 1/n*sum(abs(t(xs)%*%xs)^2))

  b1 = b2 = VS
  b <- c(b1, b2)

  gammaop <- function(g){
    gamma <- c(g[1], g[2])
    t(gamma) %*% A %*% gamma - 2*(t(gamma)%*%b)
  }

  op.gamma <- constrOptim(c(0.4, 0.5), gammaop,
                          ui = cbind(c(1,0,-1), c(0,1,-1)),
                          ci = c(0, 0, -1), method = "Nelder-Mead")

  gamma1 <- op.gamma$par[1]
  gamma2 <- op.gamma$par[2]

  S.est <- (1-gamma1-gamma2)*S.g + gamma1*S.all + gamma2*I

  return(S.est)
}


steinShrink <- function(name){ ## rewrite using match.fun
  covshrinkf <- match.fun(name)
  return(covshrinkf)
}


#' Joint graphical modeling
#'
#' Conducting joint graphical modeling using two-target linear shrinkage
#'
#' @param data n x p matrix which includes observations from all groups in rows and features in column
#' @param group.idx list of indexes for each group. For instance, if data has 2 groups then group.idx list should have 2 objects and each object is vector of indexes of observations belonging to each group.
#' @param group.select a number of which group is selected for the estimation
#' @param covshrinkf name of standard covariance linear shrinkage function to pre-shrink sample covariance matrices from each group before calculating shared target matrix in two-target linear shrinkage
#' @param jointest TRUE (applying joint graphical modeling, default) or FALSE (applying standard graphical modeling)
#'
#'
#' @return estimated partial correlation matrix of selected group


pcorshrink_joint <- function(data, group.idx, group.select = 1, covshrinkf, jointest = TRUE) {
  library(corpcor)


  #========== Step 1: Identity shrinkage ==========#

  G = length(group.idx) # number of groups
  if (group.select > G) {
    stop("select group not in the list")
  }

  S.emp = lapply(group.idx, function(idx) covshrinkf(data[idx,]))

  if (jointest) {
    #=========== Step 2: commonNet shrinkage ==========#
    # Parameter for shrinkage constant
    S.est = ttls(data, group.idx, group.select, covshrinkf)

  } else {

    S.est = S.emp[[group.select]]

  }

  powr = pseudoinverse(S.est)
  pc = -powr
  diag(pc) = -diag(pc)
  pc = cov2cor(pc)
  return(pc)
}














