library(parallel)
library(JuliaCall)

# Coarse-grained structure: 1 shared graph + 1 group-specific non shared graph
coarseGrained_simulate <- function(p, num.group, commonEdge.per = 40, degree = 2, edges.info = TRUE)
{
  #20FEB2024: modify precision.lo[idx.commonEdges] to prevent bias in JGL, networks only shared structure not value of precision matrix
  library(corpcor)
  # p: number of nodes in all networks (shared set of nodes)

  ####################
  eps = 0.0001
  commonEdge.per = commonEdge.per/100
  # eps:        additional additive component to diag(precision_matrix)
  ####################

  # select common edges
  total.edges = p*(p-1)/2
  num.edges = p*degree
  num.commonEdges = round(commonEdge.per*num.edges)
  idx.commonEdges = sample(1:total.edges, num.commonEdges)

  precision.lo = rep(0,total.edges)

  # draw number
  precision.lo[idx.commonEdges] = runif(num.commonEdges,-1.0,+1.0) ## bias in simulation, favor fused than group JGL

  # construct symmetric precision matrix for each network
  pcor.net = lapply(1:num.group, function(i){
    num.nonshared = num.edges - num.commonEdges
    if (length(idx.commonEdges)!=0) {
      idx.nonshared = sample(c(1:total.edges)[-idx.commonEdges], num.nonshared)
    } else {
      idx.nonshared = sample(c(1:total.edges), num.nonshared)
    }

    precision.lo[idx.nonshared] = runif(num.nonshared,-1.0,+1.0)
    precision = matrix(0, nrow = p, ncol = p)
    precision[lower.tri(precision)] = precision.lo
    precision = precision + t(precision)

    for(j in 1:p)
    {
      diag(precision)[j] = sum(abs(precision[,j]))+eps
    }

    # convert to partial correlation matrix
    pcor = cov2cor(precision) # standardized precision matrix (positive definite)
    pcor = -pcor # change signs of the off-diagonal entries to obtain pcor matrix
    diag(pcor) = -diag(pcor) # keep positive sign on diagonal
    return(pcor)

  })

  if (edges.info){
    return(list(pcor.net = pcor.net, common.edges = idx.commonEdges))
  } else {
    return(pcor.net)
  }
}

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

# Shrinkage calculation
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

GeneNet <- function(data) {
  return(cor.shrink(data, verbose = FALSE))
}

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

nlshrink <- function (data) {
  xs <- as.matrix(data)
  n <- nrow(xs)
  p <- ncol(xs)

  sam_cov <- covcal(xs)

  # Get the sorted eigenvalues and eigenvectors
  sam_eig <- eigen(sam_cov, symmetric = TRUE)
  lambda <- sam_eig$values
  u <- sam_eig$vectors

  # Analytical Nonlinear Shrinkage Kernel Formula
  # tolerance value to exclude error eigenvalues
  #  eig_nonzero_tol <- max(n,p)*max(lambda)*.Machine$double.eps
  eig_nonzero_tol <- 0.001
  id_select <- which(lambda > eig_nonzero_tol)

  lambda <- lambda[id_select]
  r <- length(lambda)
  L <- matrix(lambda, nrow = r, ncol = r)

  # LW Equation 4.9
  h <- n^(-1 / 3)
  H <- h * t(L)
  x <- (L - t(L)) / H

  # LW Equation 4.7
  s1 <- (3 / 4) / sqrt(5)
  s2 <- -(3 / 10) / pi
  pos_x <- (1 - (x^2) / 5)
  pos_x <- replace(pos_x, pos_x < 0, 0)
  f_tilde <- s1 * matrixStats::rowMeans2(pos_x / H)

  # LW Equation 4.8
  log_term <- log(abs((sqrt(5) - x) / (sqrt(5) + x)))
  Hftemp <- (s2 * x) + (s1 / pi) * (1 - (x^2) / 5) * log_term
  sq5 <- which(abs(x) == sqrt(5))
  Hftemp[sq5] <- s2 * x[sq5]
  H_tilde <- matrixStats::rowMeans2(Hftemp / H)

  # LW Equation 4.3
  s3 <- pi * (p / n)
  s4 <- 1 / (h^2)
  if (length(id_select) == p) { ## changing
    d_tilde <- lambda / ((s3 * lambda * f_tilde)^2 +
                           (1 - (p / n) - s3 * lambda * H_tilde)^2)
  } else {
    ones <- rep(1, p - length(id_select))
    log_term <- log((1 + sqrt(5) * h) / (1 - sqrt(5) * h))
    m <- mean(1 / lambda)

    # LW Equation C.8
    Hf_tilde0 <- (1 / pi) * ((3 / 10) * s4 +
                               (s1 / h) * (1 - (1 / 5) * s4) * log_term) * m

    # LW Equation C.5
    d_tilde0 <- 1 / (pi * (p - length(id_select)) / length(id_select) * Hf_tilde0)
    #    d_tilde0 <- 1 / (pi * (p - n) / n * Hf_tilde0)

    # LW Equation C.4
    d_tilde1 <- lambda / ((pi^2 * lambda^2) * (f_tilde^2 + H_tilde^2))
    d_tilde <- c(d_tilde1, d_tilde0 * ones)
  }

  # LW Equation 4.4
  sigma_tilde <- u %*% tcrossprod(diag(d_tilde), u)

  return(sigma_tilde)
}

# cv function
cvCovEst_general <- function(dat, estimator_params = list(alpha1 = seq(0,1,0.05), alpha2 = seq(0.1,1,0.05)), v_folds = 10L, S.g, S.all, I) {
  estimator <- function(data, alpha1, alpha2) {
    S.g <- covcal(data)
    estimate <- (1 - alpha1 - alpha2)*S.g + alpha1*S.all + alpha2*I
    return(estimate)
  }
  cv_loss <- function(fold, dat, estimator_fun, estimator_params = NULL) {

    # split the data into training and validation
    train_data <- origami::training(dat, fold = fold)
    valid_data <- origami::validation(dat, fold = fold)

    # compute the sample covariance of the validation set
    sample_cov_valid <- coop::covar(valid_data)

    # extract estimator function and name
    est_fun <- rlang::eval_tidy(estimator_fun)

    # Compute the grid of hyperparameters
    hparam_grid <- expand.grid(estimator_params, stringsAsFactors = FALSE)

    # loop through the estimator hyperparameters
    param_out <- lapply(seq_len(nrow(hparam_grid)), function(idx) {
      # fit the covariance matrix estimator on the training set
      estimator_hparam <- hparam_grid[idx, ]
      if (length(estimator_hparam) > 1) {
        estimator_hparam <- paste(estimator_hparam, collapse = ", ")
      }
      est_mat <- rlang::exec(
        est_fun,
        train_data,
        !!!purrr::flatten(dplyr::slice(hparam_grid, idx))
      )

      # return the results from the fold
      out <- tibble::tibble(
        hyperparameters = estimator_hparam,
        loss = matrixStats::sum2((est_mat - sample_cov_valid)^2),
        fold = origami::fold_index(fold = fold)
      )
      return(out)
    })

    # return data frame of estimator for all considered hyperparameters
    est_out <- dplyr::bind_rows(param_out)
    return(est_out)
  }

  # define folds
  folds <- origami::make_folds(dat, fold_fun = origami::folds_vfold, V = v_folds)

  # apply estimator to each fold
  wrapped_fun <- origami::wrap_in_try(rlang::eval_tidy(cv_loss))
  fold_results <- lapply(folds, wrapped_fun, dat = dat, estimator_fun = estimator, estimator_params = estimator_params)

  # convert results to tibble
  fold_results_concat <- dplyr::bind_rows(fold_results)

  # compute cv risk
  library(magrittr)
  cv_results <- fold_results_concat %>%
    dplyr::group_by(.data$hyperparameters) %>%
    dplyr::summarise(cv_risk = mean(.data$loss)) %>%
    dplyr::arrange(.data$cv_risk) %>%
    dplyr::ungroup()
  library(stringr)
  best_est_hparams <- as.numeric(str_split(cv_results[1, ]$hyperparameters, ", ")[[1]])
  return(best_est_hparams)
}

jointGeneNet <- function(data, group.idx, group.select, covshrinkf) {
  library(corpcor)
  G = length(group.idx) # number of groups
  xs = wt.scale(x=data[group.idx[[group.select]],], center=TRUE, scale=TRUE)

  n = nrow(xs)
  w = rep(1/n, n)
  sw = sqrt(w)
  w2 = sum(w*w)
  h1w2 = w2/(1-w2)

  S.emp = lapply(group.idx, function(idx) covcal(data[idx,]))
  S.g = covcal(data[group.idx[[group.select]],])
  S.all = 1/G*Reduce("+", S.emp)

  E2R = (crossprod(sweep(xs, MARGIN=1, STATS=sw, FUN="*")))^2
  ER2 = crossprod(sweep(xs^2, MARGIN=1, STATS=sw, FUN="*"))

  sE2R = sum(E2R) - sum(diag(E2R))
  sER2 = sum(ER2) - sum(diag(ER2))

  de = (S.g - S.all)^2

  denominator = sum(de) - sum(diag(de))
  numerator = sER2 - sE2R

  if(denominator == 0)
    alpha = 1
  else
    alpha = max(0, min(1, numerator/denominator * h1w2))

  S.est = (1-alpha)*S.g + alpha*S.all
  return(S.est)
}

jointCrossvalidation <- function(data, group.idx, group.select, covshrinkf) {
  library(corpcor)
  G = length(group.idx) # number of groups
  S.emp = lapply(group.idx, function(idx) covshrinkf(data[idx,]))
  S.g = covcal(data[group.idx[[group.select]],])
  I = diag(ncol(data))
  S.all = 1/(G-1)*Reduce("+", S.emp[-group.select])

  est_output <- cvCovEst_general(dat = data[group.idx[[group.select]],], S.g = S.g, S.all = S.all, I = I)
  alpha1 <- est_output[1]
  alpha2 <- est_output[2]
  S.est <- (1 - alpha1 - alpha2)*S.g + alpha1*S.all + alpha2*I

  return(S.est)
}

mtse_v1 <- function(data, group.idx, group.select = 1, covshrinkf) {
  library(dplyr)
  G = length(group.idx) # number of groups
  p = ncol(data)
  n = length(group.idx[[group.select]])
  S.emp = lapply(group.idx, function(idx) covshrinkf(data[idx,]))
  S.g = covcal(data[group.idx[[group.select]],])
  S.all = 1/(G-1)*Reduce("+", S.emp[-group.select])
  I = diag(p)

  #1. Calculate A matrix
  aSS = S.all - S.g
  aIS = I - S.g
  a11 = sum(diag(aSS%*%t(aSS)))
  a22 = sum(diag(aIS%*%t(aIS)))
  a12 = sum(diag(t(aIS)%*%aSS))
  a21 = sum(diag(t(aSS)%*%aIS))
  A = matrix(c(a11,a12,a21,a22), nrow = 2, ncol = 2, byrow = TRUE)

  #2. Calculate b vector
  xs = wt.scale(data[group.idx[[group.select]],])
  const1 = n/((n-1)^2*(n-2))
  VS = const1*(sum(diag(xs%*%t(xs))^2) - 1/n*sum(abs(t(xs)%*%xs)^2))

  b1 = b2 = VS
  b = c(b1, b2)

  gamma1 = c((-a22*b1 + a12*b2)/(a12^2 - a11*a22), 0, b1/a11, (a22-a12)/(a11-2*a12+a22), 0, 0, 1)
  gamma2 = c((b1 - a11*gamma1[1])/a12, b2/a22, 0, 1-gamma1[4], 0, 1, 0)

  g.idx = sapply(1:7, function(i) between(gamma1[i],0,1)&&between(gamma2[i],0,1)&&(gamma1[i]+gamma2[i] <= 1))
  loss = sapply((1:7)[g.idx], function(i){
    gamma <- c(gamma1[i], gamma2[i])
    t(gamma) %*% A %*% gamma - 2*(t(gamma)%*%b)
  })
  gamma1 = gamma1[g.idx][which.min(loss)]
  gamma2 = gamma2[g.idx][which.min(loss)]


  S.est = (1-gamma1-gamma2)*S.g + gamma1*S.all + gamma2*I

  return(S.est)
}

mtse <- function(data, group.idx, group.select = 1, covshrinkf) {
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
  S.emp <- lapply(group.idx, function(idx) covcal(y[idx,]))
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

HC_objective <- function(pi, i, N, method = 'DonohoJin2008') {
  if (method == 'DonohoJin2008') {
    hc <- sqrt(N)*(i/N - pi)/sqrt(i/N*(1-i/N))
  } else if (method == 'BerkJones1979') {
    hc <- sqrt(2*N)*sqrt((i/N)*log(i/(N*pi)) + (1-i/N)*log((1-i/N)/(1-pi)))
  } else if (method == 'LiSiegmund2015') {
    if (pi < i/N) {
      hc <- sqrt(2*N)*sqrt((i/N)*log(i/(N*pi)) - (i/N - pi))
    } else {
      hc <- 0
    }
  }
  return(hc)
}


SEst <- function(method, ...){
  method <- match.fun(method)
  S.est <- method(...)
  return(S.est)
}

steinShrink <- function(name){ ## rewrite using match.fun
  covshrinkf <- match.fun(name)
  return(covshrinkf)
}

cov.shrink.commonNet <- function(data, group.idx, group.select = 1, covshrinkf, alpha = NULL, jointest = TRUE, alpha.method = "jointGeneNet") {
  library(corpcor)
  # data: feature in columns and samples in row
  # group.idx: list of index of samples in each group
  # covshrinkf: shrink function return shrunk correlation matrix


  #========== Step 1: Identity shrinkage ==========#

  G = length(group.idx) # number of groups
  if (group.select > G) {
    stop("select group not in the list")
  }

  S.emp = lapply(group.idx, function(idx) covshrinkf(data[idx,]))

  if (jointest) {
    #=========== Step 2: commonNet shrinkage ==========#
    # Parameter for shrinkage constant

    if (is.null(alpha)){
      S.est = SEst(method = alpha.method, data, group.idx, group.select, covshrinkf)
    } else {
      S.g = covcal(data[group.idx[[group.select]],])
      S.all = 1/(G-1)*Reduce("+", S.emp[-group.select])

      alpha = alpha
      S.est = (1-alpha)*S.g + alpha*S.all
    }

  } else {

    S.est = S.emp[[group.select]]

  }

  powr = pseudoinverse(S.est)
  pc = -powr
  diag(pc) = -diag(pc)
  pc = cov2cor(pc)
  return(pc)
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


