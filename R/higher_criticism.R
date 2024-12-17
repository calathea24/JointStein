
#' Second-level significance testing using higher criticism
#'
#' Calculating higher criticism statistics based on calculated p-values using higher criticism. Calculated p values are ordered from lowest to highest values.
#' @param pi p value
#' @param i position of input p value in sorted list
#' @param N number of p values in sorted list
#' @param method name of method to calculate higher criticism statistics. Default is DonohoJin2008. Other methods are BerkJones1979 and LiSiegmund2015.
#'
#' @return value of higher criticism statistics


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









