import::here("as_adj_list", .from = igraph)

#This function calls the sBread function written in Rcpp from the igraph object.
#G: igraph object
#mut.rate: mutation rate
#sample.size.by.chain: number of samples for each chain.
#num.chains: number of chains
#sample.interval: sampling interval in the MH-algorithm
sBread.from.igraph <- function(G, known.data, K, mut.rate, sample.size.by.chain, num.chains, sample.interval, seed = NULL){
  T <- nrow(known.data)
  N <- ncol(known.data)
  
  if(is.null(seed)){
    seed <- sample(c(1:10^7), 1)
  }
  
  tmp.adj.list <- as_adj_list(G, "total")
  adj.list <- list()
  weight.list <- list()
  
  for(i in 1:N){
    adj.list[[i]] <- as.vector(tmp.adj.list[[i]])
    adj.list[[i]] <- adj.list[[i]] - 1
    adj.list[[i]] <- c(adj.list[[i]], i-1)
    degree <- length(adj.list[[i]])
    weight.list[[i]] <- rep(1/degree, degree)
  }
  
  sample.histories <- sBread(seed, K, mut.rate, num.chains, sample.interval, sample.size.by.chain,
                             sample.interval, adj.list, weight.list, known.data)

  return(sample.histories)
}