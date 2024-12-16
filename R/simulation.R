# Coarse-grained structure: 1 shared graph + 1 group-specific non shared graph
coarseGrained_simulate <- function(p, num.group, commonEdge.per = 40, degree = 2, edges.info = TRUE) {
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








