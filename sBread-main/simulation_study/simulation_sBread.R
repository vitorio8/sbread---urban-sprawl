library(Rcpp)
sourceCpp("main.cpp")
sourceCpp("analyze_posterior.cpp")

seed <- 2
set.seed(seed)

#global variables.
L <- 10
T <- 50
K <- 2
mut.rate <- 0.01
trans.rate <- 0.1

#to create masked history
known.rate <- 0.2

#sBread input
num.sample <- 100
interval <- 10 ^ 4
burn.in <- 10 ^ 6


create.history <- function(adj.list, weight.list, K, num.timestep, mut.rate){
  N <- length(adj.list)
  history <- matrix(nrow = num.timestep, ncol = N)
  
  history[1,] <- sample(0:(K-1), N, replace = TRUE)
  for(t in 2:num.timestep){
    for(i in 1:N){
      model <- sample(adj.list[[i]], 1, prob = weight.list[[i]])
      if(runif(1) < mut.rate){
        history[t,i] <- sample(0:(K-1), 1)
      }else{
        history[t,i] <- history[t-1, model]
      }
    }
  }
  
  return(history)
}

square.grid.adj.list <- function(L, weight){
  N <- L ^ 2
  if(weight > 0.2){
    stop("Weight cannot exceed 0.2")
  }
  
  adj.list <- list()
  for(i in 1:N){
    row <- ceiling(i/L)
    col <- ifelse(i %% L != 0, i %% L, L)
    
    current.neighbors <- c(i)
    if(col > 1){
      current.neighbors <- append(current.neighbors, i-1)
    }
    if(col < L){
      current.neighbors <- append(current.neighbors, i+1)
    }
    if(row > 1){
      current.neighbors <- append(current.neighbors, i-L)
    }
    if(row < L){
      current.neighbors <- append(current.neighbors, i+L)
    }
    
    adj.list[[i]] <- current.neighbors
  }
  
  weight.list <- list()
  for(i in 1:N){
    current.degree <- length(adj.list[[i]])
    weight.list[[i]] <- 
      c(1 - (current.degree - 1) * weight, rep(weight, current.degree - 1))
  }
  
  return(list(adj.list, weight.list))
}

#convert the index in a adj list into the Rcpp format.
index.convert <- function(adj.list){
  len <- length(adj.list)
  for(i in 1:len){
    adj.list[[i]] <- adj.list[[i]] - 1
  }
  return(adj.list)
}

#Censor part of the true history and create the masked history.
#The first and last timestep remain unmasked (two known maps)
#State values in other timesteps are known with prob.known.
mask.history <- function(true.history, prob.known){
  T <- nrow(true.history)
  N <- ncol(true.history)
  result <- matrix(-1, nrow = T, ncol = N)
  is.known <- matrix(runif(T * N) < prob.known, nrow = T, ncol = N)
  result <- ifelse(is.known, true.history, -1)
  result[1,] <- true.history[1,]
  result[T,] <- true.history[T,]
  return(result)
}



##main routine below
adj.weight.list <- square.grid.adj.list(L, trans.rate)
adj.list <- adj.weight.list[[1]]
weight.list <- adj.weight.list[[2]]

#create the ground truth
true.history <- create.history(adj.list, weight.list,
                          K, T * 10, mut.rate)[(9*T+1):(T*10),]
#masked.history <- matrix(-1, nrow = T, ncol = L ^ 2)
#masked.history[1,] <- true.history[1,]
#masked.history[T,] <- true.history[T,]
masked.history <- mask.history(true.history, known.rate)

#run sBread
sampled.histories <- sBread(seed, K, mut.rate, 1, burn.in, num.sample,
                            interval, index.convert(adj.list), 
                            weight.list, masked.history)

num.unknown.timeline <- rowSums(masked.history == -1)
performance.timeline.by.sample <- matrix(NA, nrow = num.sample, ncol = T)
for(i in 1:num.sample){
  performance.timeline.by.sample[i,] <- 
    (rowSums(sampled.histories[[i]] == true.history) 
     - (L^2 - num.unknown.timeline)) / num.unknown.timeline
}

performance.timeline <- colSums(performance.timeline.by.sample) / num.sample
plot(performance.timeline)