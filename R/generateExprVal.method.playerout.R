generateExprVal.method.playerout <- function(probes, weights=F){
  matos <- as.data.frame(lapply(probes,function(x){x$pm}))
  ##names(matos) <- mynames
  matos <- t(as.matrix(matos))
  n <- length(probes[[1]]$pm)
  
  ## skip if only one probe
  if (n == 1) return(unlist(lapply(probes, function(x){x$pm})))
  
  ## I do not know to which extend the use of optim
  ## is really equivalent to the use of nlminb in S-plus
  S1 <- optim(runif(n),
              playerout.costfunction,
              method="L-BFGS-B",
              control=list(maxit=500),
              y=matos)
  ##S1 <- nlm(playerout,runif(20),iterlim=500,y=t(y))
  r <- c(matos %*% S1$par / sum(S1$par))
  if (weights)
    attr(r,"weights") <- S1$par
  return(r)
}


## Code kindly provided by E. Lazaridris
## (I could not help butchering it a bit before putting it into the package)

## The loss function:

playerout.costfunction <- function(w, y) {
  N <- length(w)        # Number of players
  J <- length(y)/N      # Number of games (the number of games is the number of chips used)
  sumw <- sum(w)
  
  tx <- y %*% w    # Full weighted score at each game  
  pl <- matrix(0,J,N)    # Loss at each game due to each player
  
  for(j in 1:J)
    pl[j,] <- w * y[j,] - (tx[j] - w * y[j,]) / (sumw - w)
  
  sum(pl^2)         # Loss
}


