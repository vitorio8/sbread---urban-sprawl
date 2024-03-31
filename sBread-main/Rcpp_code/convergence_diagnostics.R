#library(tidyverse)
import::here(sourceCpp, .from = Rcpp)
sourceCpp("analyze_posterior.cpp")

#coumpute the expected value of similarity.
#posterior.prob : List of matrices whose [[k]][t][i] element shows the posterior probability of state k at time t in cell i.
#num.known : number of pairs (t, i) such that the state at t in cell i is known.
expected.similarity <- function(posterior.prob, num.known){
  K <- length(posterior.prob)
  T <- nrow(posterior.prob[[1]])
  N <- ncol(posterior.prob[[1]])
  
  coincidence.prob.matrix <- matrix(0, nrow = T, ncol = N)
  for(k in 1:K){
    coincidence.prob.matrix <- coincidence.prob.matrix + (posterior.prob[[k]]) ** 2
  }
  expected.num.coincidence <- sum(coincidence.prob.matrix)
  
  S.expected <- (expected.num.coincidence - num.known) / (T * N - num.known)

  return(S.expected)
}


#coumpute the sequence of the ratio of the expected and realized similarities
similarity.sequence<- function(K, sample.histories, known.data, reference.sample.id){
  T <- nrow(sample.histories[[1]])
  N <- ncol(sample.histories[[2]])
  num.sample <- length(sample.histories)
  num.known <- sum(known.data != -1)
  
  sample.similarity <- function(s1, s2) {
    return((sum(s1 == s2) - num.known) / (T * N - num.known))
  }

  posterior.prob <- compute_posterior(K, T, N, sample.histories)
  S.expected <- expected.similarity(posterior.prob, num.known)

  S.rate.seq <- rep(NA, num.sample)
  for(i in 1:num.sample){
    S.realized <- sample.similarity(sample.histories[[i]], sample.histories[[reference.sample.id]])
    S.rate.seq[i] <- S.realized / S.expected
  }

  return(S.rate.seq)
}

# Coumpute similarities between pairs of samples that are "lag" apart from each
# other in the MCMC chain with varying values of "lag"
autocorr.sequence <- function(K, sample.histories, known.data, max.lag){
  T <- nrow(sample.histories[[1]])
  N <- ncol(sample.histories[[2]])
  num.sample <- length(sample.histories)
  num.known <- sum(known.data != -1)
  
  sample.similarity <- function(s1, s2) {
    return((sum(s1 == s2) - num.known) / (T * N - num.known))
  }
  
  posterior.prob <- compute_posterior(K, T, N, sample.histories)
  S.expected <- expected.similarity(posterior.prob, num.known)
  
  lagged.similarities <- data_frame(row.names = c("lag", "i", "similarity"))
  for(lag in 1:max.lag){
    for (i in seq(1, num.sample-lag, 5)) {
      S.realized <- sample.similarity(sample.histories[[i]], sample.histories[[i+lag]])
      S.rate <- S.realized / S.expected
      lagged.similarities <- lagged.similarities %>% 
        add_row(lag= l, i=i, similarity=S.rate)
    }
  }
  return(lagged.similarities)
}
