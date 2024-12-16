library(GeneNet)
library(mltools)
P <- as.numeric(snakemake@params[["P"]])
N <- as.numeric(snakemake@params[["N"]])
G <- as.numeric(snakemake@params[["G"]])
covshrinkf <- snakemake@params[["alg"]]
jointest <- snakemake@params[["joint"]]
padjMet <- snakemake@params[["padj"]]
hc <- snakemake@params[["hc"]]
load(snakemake@input[[1]])
load(snakemake@input[[2]])


source("utils.R")
est <- network.test.edges(pc, plot = FALSE)
padj <- p.adjust(est$pval, method = padjMet)

if (hc) {
  num.edges <- P*(P-1)/2
  alpha <- 0.01*num.edges
  hc_stat <- unlist(lapply(1:num.edges, function(i){
    HC_objective(padj[i], i, num.edges)
  }))
  hc_thres <- which.max(hc_stat[1:alpha])
  idx.pselect <- seq(1, hc_thres)
} else {
  idx.pselect <- which(padj <= 0.01)
}

adj_mat <- matrix(0, P, P)
for (i in idx.pselect) {
  row <- est[i,2]
  col <- est[i,3]
  adj_mat[row, col] <- 1
}
adj_mat <- adj_mat + t(adj_mat)
mcc_val <- mcc.pcor(adj_mat, pcor.net[[1]])

save(mcc_val, file = snakemake@output[[1]])













