library(GeneNet)
P <- as.numeric(snakemake@params[["P"]])
N <- as.numeric(snakemake@params[["N"]])
G <- as.numeric(snakemake@params[["G"]])
commonpercent <- as.numeric(snakemake@params[["coper"]])

source("utils.R")

#----------Simulation from snakemake input------------#
sim.net <- coarseGrained_simulate(P, num.group = G, commonEdge.per = commonpercent, degree = 2, edges.info = TRUE)
pcor.net <- sim.net$pcor.net
sim.data <- do.call(cbind,lapply(pcor.net, function(pcor){
  t(ggm_simulate_data(N, pcor = pcor))
}))

save(pcor.net, file = snakemake@output[[1]])
save(sim.data, file = snakemake@output[[2]])







