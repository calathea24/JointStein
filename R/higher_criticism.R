
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









