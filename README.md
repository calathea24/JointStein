# JointStein

<!-- badges: start -->

<!-- badges: end -->

JointStein is an R package to implement the joint sparse inverse covariance matrix estimation framework for single-cell RNA sequencing (scRNAseq) network analysis using two-target linear shrinkage. Reproducible snakemake workflow is in vignettes

## Installation

You can install the development version of JointStein from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("calathea24/JointStein")
```

## Example

``` r
library(ZINBStein)
## Simulating data of 5 groups with 2000 features/genes and 1000 samples/cells with coarse-grain joint network structure
p <- 2000
n <- 1000
G <- 5
sim.net <- coarseGrained_simulate(p, num.group = G, degree = 2, edges.info = TRUE, commonEdge.per = 0.8)
pcor.net <- sim.net$pcor.net

sim.data <- do.call(cbind,lapply(pcor.net, function(pcor){
  t(ggm_simulate_data(n, pcor = pcor))
}))
covshrinkf <- steinShrink('GeneNet)
group.idx <- split(1:ncol(sim.data), ceiling(seq_along(1:ncol(sim.data))/n))

## Using two-target linear shrinkage algorithm to shrink sample covariance matrix and estimate partial correlation matrix
pc <- pcorshrink_joint(t(sim.data), group.idx, group.select = 1, covshrinkf, jointest = TRUE)
est <- network.test.edges(pc, plot = FALSE)
padj <- p.adjust(est$pval, method = 'BH') # use Benjamini-Hochberg method


## Higher criticism when p > 1000
num.edges <- p*(p-1)/2
alpha <- 0.01*num.edges
hc_stat <- unlist(lapply(1:num.edges, function(i){
  HC_objective(padj[i], i, num.edges)
}))
hc_thres <- which.max(hc_stat[1:alpha])
idx.pselect <- seq(1, hc_thres)


## Generating adjacency matrix and calculate Matthew's correlation coefficient
adj_mat <- matrix(0, p, p)
for (i in idx.pselect) {
  row <- est[i,2]
  col <- est[i,3]
  adj_mat[row, col] <- 1
}
adj_mat <- adj_mat + t(adj_mat)
mcc_val <- mcc.pcor(adj_mat, pcor.net[[1]])





