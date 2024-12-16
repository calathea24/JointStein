library(GeneNet)
P <- as.numeric(snakemake@params[["P"]])
N <- as.numeric(snakemake@params[["N"]])
G <- as.numeric(snakemake@params[["G"]])
alg <- snakemake@params[["alg"]]
jointest <- snakemake@params[["joint"]]
alpha.method <- snakemake@params[["alphaMet"]]
load(snakemake@input[[1]])

source("utils.R")
group.idx <- split(1:ncol(sim.data), ceiling(seq_along(1:ncol(sim.data))/N))
covshrinkf <- steinShrink(alg)
pc <- cov.shrink.commonNet(t(sim.data), 
                           group.idx = group.idx, 
                           group.select = 1, 
                           covshrinkf = covshrinkf, 
                           jointest = jointest, 
                           alpha.method = alpha.method)

save(pc, file = snakemake@output[[1]])





