if (endsWith(getwd(), "sBread")) setwd("Rcpp_code")
import::here(barabasi.game, .from = igraph)
import::here(sourceCpp, .from = Rcpp)
import::here(mclapply, .from = parallel)
import::here(sBread.from.igraph, .from = "sBread_from_igraph.R")
import::here(similarity.sequence, .from = "convergence_diagnostics.R")
import::here(discrete.ess, .from = "discrete_ess.R")
sourceCpp("main.cpp")
sourceCpp("analyze_posterior.cpp")

seed <- 1
set.seed(seed)

# Simulation parameters
T <- 50  # Number of time-steps
N <- 50  # Number of locations
K <- 3   # Number of states
mut.rate <- 0.1

# MCMC parameters
num.chains  <- 2
sample.size.by.chain <- 2000
sample.size.total <- sample.size.by.chain * num.chains
sample.interval <- 1000
burn.in.samples <- 200

# Generate a random neighborhood graph
G <- barabasi.game(N, m = 2)

# Generate random data (i.e. start and end states)
known.data <- matrix(-1, nrow = T, ncol = N)
known.data[1,] <- sample(0:1, N, replace = TRUE)
known.data[T,] <- sample(0:1, N, replace = TRUE)

# Run sBread
run.sBread.chain <- function(chain_id, G, known.data){
  return(
    sBread.from.igraph(G, known.data, K, mut.rate, sample.size.by.chain, 1, sample.interval, seed+chain_id)
  )
}

samples_by_chain <- mclapply(1:num.chains, run.sBread.chain, G=G, known.data=known.data)
# samples <- sBread.from.igraph(G, known.data, K, mut.rate, sample.size.by.chain, num.chains, sample.interval, seed)

for (c in 1:num.chains) {
  # Drop burn-in
  samples <- samples_by_chain[[c]][burn.in.samples:sample.size.by.chain]
  
  # Compute and print ESS
  ess.results <- discrete.ess(samples)
  print(sprintf("The MCMC resulted in an effective sample size of ESS %s %s", ess.results$operator, ess.results$ess))
  
  # Plot the similarity between samples over time
  x <- similarity.sequence(K, samples, known.data, 1)
  plot(x, pch=16, cex=0.2, col=scales::alpha("black", 0.25), ylim=c(0, 3.3))
  abline(h=1)
}
