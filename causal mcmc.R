library(nlme)
library(ape)
library(phylolm)

`%notin%` <- Negate(`%in%`)

amat2model <- function(dag) {
  active.pairs <- which(dag == 1, arr.ind=T)
  print(length(active.pairs))
  if (length(active.pairs) == 0) {return("")}
  res <- c()
  for (i in 1:dim(active.pairs)[1]) {
    node1 <- colnames(dag)[active.pairs[i, 1]]
    node2 <- colnames(dag)[active.pairs[i, 2]]
    
    model <- paste(node2, '~', node1)
    res <- c(res, model)
  }
  res <- paste(res, collapse = ',')
  
  return(res)
}

CSTAT <- function(dag, data, tree) {
  q <- dim(dag)[1] + sum(dag)
  n <- length(tree$tip.label)
  p.values <- dsep(dag, data, tree)

  if (is.null(p.values)) {return(Inf)}
  
  C <- -2 * sum(log(p.values))
  CICc <- C + (2*q) * (n / (n - 1 - q))
  
  return(CICc)
}

get.parents <- function(dag, node) {
  parents <- colnames(dag)[dag[,node] == 1]
  return(parents)
}

get.children <- function(dag, node) {
  children <- get.parents(t(dag), node)
  return(children)
}

get.all.descendents <- function(dag, node) {
  descendents <- get.all.ancestors(t(dag), node)
  return(descendents)
}

get.all.ancestors <- function(dag, node) {
  unprocessed <- c(get.parents(dag, node))
  processed <- c()
  ancestors <- unprocessed
  
  repeat {
    if (length(unprocessed) == 0) {break}
    node <- unprocessed[1]
    if (node %notin% processed & !is.na(node)) {
      processed <- c(processed, node)
      
      parents <- get.parents(dag, node)
      ancestors <- unique(c(ancestors, parents))
      for (parent in parents) {
        if (parent %notin% processed) {
          unprocessed <- c(unprocessed, parent)
          }
      }
    }
    unprocessed <- unprocessed[-1]
  }
  
  return(ancestors)
}

mutate <- function(dag, forbidden) {
  legal.adds <- dag == 0
  legal.adds[t(dag) == 1] <- F
  diag(legal.adds) <- F
  legal.adds[forbidden == 1] <- F
  for (node in colnames(dag)) {
    legal.adds[node, get.all.ancestors(dag, node)] <- F
  }
  legal.subs <- dag == 1
  
  legal.inversions <- legal.subs | legal.adds
  
  swappable.pairs <- list()
  active.pairs <- which(dag == 1, arr.ind = T)

  if (length(active.pairs)) {
    for (i in 1:dim(active.pairs)[1]) {
      node1 <- colnames(dag)[active.pairs[i, 1]]
      node2 <- colnames(dag)[active.pairs[i, 2]]
  
      if (legal.inversions[node1, node2] & legal.inversions[node2, node1]) {
      
        tmp.dag <- dag
        tmp.dag[node1, node2] <- 0
        
        all.descendents <- c()
        for (ancestor in get.all.ancestors(tmp.dag, node2)) {
          all.descendents <- c(all.descendents, get.all.descendents(tmp.dag, ancestor))
        }
        if (node2 %notin% all.descendents) {
          swappable.pairs[length(swappable.pairs) + 1] <- list(c(node1, node2))
        }
      }
    }
  }

  can.invert <- sum(legal.inversions) > 0
  can.swap <- length(swappable.pairs) > 0
  
  r <- runif(n=1, min=0, max=1)
  
  if (can.swap & (!can.invert | r > .9)) {
    dag <- swap(dag, swappable.pairs)
  } else if (can.invert & (!can.swap | r < .9)){
    dag <- invert(dag, legal.inversions)
  } else {
    stop('cannot change dag')
  }
 return(dag) 
}

invert <- function(dag, legal.inversions) {
  invertable.coords <- which(legal.inversions)
  to.invert <- sample(invertable.coords, 1)

  debug.dag <- dag
    
  dag[to.invert] <- (dag[to.invert] + 1) %% 2
  
  return(dag)
}

swap <- function(dag, swappable.pairs) {
  to.swap <- sample(swappable.pairs, 1)[[1]]
  node1 <- to.swap[1]
  node2 <- to.swap[2]
  
  dag[node1, node2] <- 0
  dag[node2, node1] <- 1
  
  return(dag)
}

get.basis.set <- function(dag) {
  non.adj <- which(dag == 0, arr.ind = T)
  basis <- list()
  for (i in 1:dim(non.adj)[1]) {
    row <- non.adj[i,]
    node1.i <- row[1]
    node2.i <- row[2]
  
    if (node1.i < node2.i & dag[node2.i, node1.i] == 0) {
      node1 <- colnames(dag)[node1.i]
      node2 <- colnames(dag)[node2.i]
      
      res <- c(node1, node2)
      for (node in c(node1, node2)) {
        parents <- which(dag[,node]==1)
        res <- c(res, colnames(dag)[parents])
      }
      basis[length(basis) + 1] <- list(unique(res))  
      }
    }
  return(basis)
  }
  

dsep <- function(dag, data, tree) {
  basis.set <- get.basis.set(dag)
  p.values <- c()
  for (set in basis.set) {
    if (length(set) > 2 & any(get.parents(dag, set[1]) %in% set[3:length(set)])) {
      y <- set[1]
      interest <- set[2]
      set <- set[-1]
    } else {
      y <- set[2]
      interest <- set[1]
      set <- set[-2]
    }
    x <- paste0(set, collapse='+')
    model <- as.formula(paste(y, '~', x))
    res <- summary(phylolm(model, data=data, phy=tree, model='lambda'))
    p.values <- c(p.values, res$coefficients[interest, 'p.value'])
  }
  return(p.values)
}

mcmc <- function(data, tree, initial.dag=NULL, forbidden.connections=NULL,
                 n.iterations=1e6, sample.every=1, log="") {
  if (is.null(initial.dag)) {
    initial.dag <- matrix(nrow=length(colnames(data)), ncol=length(colnames(data)), dimnames=list(colnames(data), colnames(data)), data=0)
  }
  
  if (is.null(forbidden.connections)) {
    forbidden.connections <- matrix(nrow=length(colnames(data)), ncol=length(colnames(data)), dimnames=list(colnames(initial.dag), rownames(initial.dag)), data=0)
  }
  
  best.dag <- NULL #Only reported for exploratory purposes; not useful in analysis
  best.score <- Inf
  
  next.dag <- NULL
  next.score <- Inf
  
  current.dag <- initial.dag
  current.score <- CSTAT(current.dag, data, tree)
  
  for (i in 1:n.iterations) { 
    if (i%%sample.every == 0) { #Is it time to sample?
      #This is the only line that outputs anything. Th
      write(paste('Generation: ', i, ' | Best: ', best.score, ' | Score: ', current.score, ' | ', amat2model(current.dag)), file=log, append=T)
    }
    
    next.dag <- mutate(current.dag, forbidden.connections)
    next.score <- CSTAT(next.dag, data, tree)
    
    r <- current.score / next.score #Calculate acceptance probablity
    if (is.nan(r) | r > runif(n=1, min=0, max=1)) { #Accept with probability r
      current.score <- next.score
      current.dag <- next.dag
      
      if (current.score < best.score) { #Again, only kept track of for exploratory purposes
        best.score <- current.score
        best.dag <- current.dag
      }
    } 
  } 
}