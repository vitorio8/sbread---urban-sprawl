sBread.from.ASCIIart <- function(input_directory_name, T, K, mut.rate, sample.size.by.chain, num.chains, sample.interval, seed = NULL){
  if(is.null(seed)){
    seed <- sample(c(1:10^7), 1)
  }
  
  input.data <- ASCIIart_reader(input_directory_name, T)
  
  known.data <- input.data[[1]]
  adj.list <- input.data[[2]]
  weight.list <- input.data[[3]]
  
  sample.histories <- sBread(seed, K, mut.rate, num.chains, sample.interval, sample.size.by.chain,
                             sample.interval, adj.list, weight.list, known.data)
  
  return(sample.histories)
}