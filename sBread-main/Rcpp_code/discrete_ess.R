### This module is adapted from the approximate phylogenetic tree ESS
### in the RWTY package (https://github.com/danlwarren/RWTY).

import::here(detectCores, mclapply, .from = parallel)


discrete.ess <- function(samples, max.sampling.interval = 100, use.all.samples = FALSE){
  N <- length(samples)
  if(N < max.sampling.interval){
    warning("Not enough samples to use your chosen max.sampling.interval")
    warning("Setting it to 90% of the number of samples instead")
    max.sampling.interval = floor(N * 0.9)
  }
  
  # set them equal, so we get every interval.
  autocorr.intervals = max.sampling.interval
  
  print(sprintf("Calculating approximate ESS with sampling intervals from 1 to %d", max.sampling.interval))
  
  autocorr.df = sample.autocorr(samples, max.sampling.interval, autocorr.intervals, use.all.samples = use.all.samples)
  autocorr.m = estimate.autocorr.m(autocorr.df)
  approx.ess.df = approx.ess(autocorr.df, autocorr.m, N)
  
  return(approx.ess.df)
  
}
estimate.autocorr.m <- function(data, ac.cutoff = 0.95){
  
  # Fit an exponential variogram
  this.ac.result <- optim(par = c(1, 1), 
                          function(data, par){sum((par[1] * (1 - exp(-(data$sampling.interval/par[2]))) - data$sample.distance)^2)}, 
                          data = data)
    
  # If the cutoff is exceeded within the set of sampling intervals, return the interval
  # Else return "-1"
  if(any(data$sample.distance/this.ac.result$par[1] > ac.cutoff)){
    autocorr.m <- data$sampling.interval[min(which(data$sample.distance/this.ac.result$par[1] > ac.cutoff))]
  } else {
    autocorr.m <- -1
  }
  
  return(autocorr.m)
}


approx.ess <- function(df, autocorr.time, N){
  # If the autocorrelation time is larger than the maximum sample it is recorded as -1
  if(autocorr.time < 0){
    m = nrow(df) + 1
  }else{
    m = autocorr.time
  }
  
  D = max(df$sample.distance)
  S = 0
  
  if(m>1){
    for(k in 1:(m - 1)){
      f = df$sample.distance[k]
      S = S + ((N - k) * f)
    }
  }
  
  S = S + (N - m + 1) * (N - m) * D / 2
  S = S / 2 / N^2
  ESS = 1 / (1 - 4 * S / D)
  
  # sometimes we can only give an upper bound
  if(autocorr.time<0){
    operator = "<"
  }else{
    operator = "="
  }
  
  return(list("ess" = ESS, "operator" = operator))
}


sample.autocorr <- function(sample.list, max.sampling.interval = NA, autocorr.intervals = 100, use.all.samples = FALSE){
  
  if(!is.numeric(autocorr.intervals)) stop("autocorr.intervals must be a positive integer")
  if(autocorr.intervals<1 | autocorr.intervals%%1!=0) stop("autocorr.intervals must be a positive integer")
  
  # this ensures that we can tell you if your ESS is < some threshold
  # the max(,2) bit is a fallback for extremely short sample lists
  max.thinning <- max.sampling.interval
  
  if(max.thinning > (length(sample.list) - 100)) {
    max.thinning = length(sample.list) - 100
  }
  
  if(max.thinning < 20){
    stop(paste("Too few samples to calculate autocorrelation. Try including more samples or changing the max.sampling.interval. Current value is", 
               max.sampling.interval, ", minimum value is 20."))
  }
  
  # we analyze up to autocorr.intervals thinnings spread evenly, less if there are non-unique numbers
  thinnings <- unique(as.integer(seq(from = 1, to = max.thinning, length.out=autocorr.intervals)))
  
  
  r <- lapply(as.list(thinnings), get.sequential.distances, sample.list, use.all.samples = use.all.samples) 
  r <- data.frame(matrix(unlist(r), ncol=2, byrow=T))
  names(r) = c("sample.distance", "sampling.interval")
  
  return(r)
}

get.sequential.distances <- function(thinning, sample.list, N=500, use.all.samples = FALSE){
  
  starts = 1:(length(sample.list) - thinning)
  ends = starts + thinning
  keep = c(rbind(starts, ends))
  
  pairs = split(keep, ceiling(seq_along(keep)/2))
  
  if(use.all.samples == FALSE){
    if(length(pairs)>N){
      # only subsample if we have enough...
      pairs = sample(pairs, N)
    }
  }
  
  distances <- mclapply(pairs, sample.dist, samples = sample.list, mc.cores =  get.processors(NULL))
  
  distances <- as.numeric(unlist(distances))
  distances <- as.data.frame(distances)
  result <- apply(distances, 2, mean)
  result <- data.frame('distance' = t(t(result)))
  result$sampling.interval <- thinning
  return(result) 
}

sample.dist <- function(pair, samples){
  s1 = samples[[pair[1]]]
  s2 = samples[[pair[2]]]
  return(sum(s1 != s2))
}

get.processors <- function(processors){
  
  if(Sys.info()["sysname"] == 'Windows'){
    # mclapply is not supported on windows
    # so we give a single processor,
    # in which case mclapply calls fall back
    # on lapply
    return(1)
  }
  
  # check for global user-defined variable
  if(exists('rwty.processors')){
    # should be an integer
    if(!is.numeric(rwty.processors)){
      stop("the global rwty.processors variable must be an integer")
    }
    if(rwty.processors%%1==0){
      available_processors = detectCores(all.tests = FALSE, logical = FALSE)
      if(rwty.processors > available_processors){
        rwty.processors = available_processors - 1
      }
      if(rwty.processors < 1){
        rwty.processors = 1
      }
      return(rwty.processors)
    }else{
      stop("the global rwty.processors variable must be an integer")
    }
  }
  
  if(is.null(processors)){ 
    available_processors = detectCores(all.tests = FALSE, logical = FALSE)
    processors = max(c(1, c(available_processors - 1)))
  }
  
  return(processors)
  
}