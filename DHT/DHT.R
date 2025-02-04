library(vegan)
library(transport)
library(pcaPP)
library(ade4)
library(parallel)
library(dplyr)
library(tidyr)
library(Rfast)
library(data.table)
library(psych)
ordimedian.pcaPP <- function(y, groups){
  
  inds <- names(table(groups));
  n.inds <- length(inds); 
  medians <- matrix(0, nrow = n.inds, ncol = ncol(y))
  rownames(medians) <- inds
  for (i in 1:n.inds) {
    X <- y[groups == inds[i], , drop = FALSE]
    medians[i, ] <- l1median_NLM(X)$par;
  }
  medians;
}
cailliez <- function (distmat, print = FALSE, tol = 1e-07, cor.zero = TRUE) {
  if (is.euclid(distmat)) {
    warning("Euclidean distance found : no correction needed")
    return(list(D = distmat, constant = 0))
  }
  distmat <- as.matrix(distmat)
  size <- ncol(distmat)
  m1 <- matrix(0, size, size)
  m1 <- rbind(m1, -diag(size))
  m2 <- -bicenter.wt(distmat * distmat)
  m2 <- rbind(m2, 2 * bicenter.wt(distmat))
  m1 <- cbind(m1, m2)
  lambda <- eigen(m1, only.values = TRUE)$values
  c <- max(Re(lambda)[Im(lambda) < tol])
  if (print) 
    cat(paste("Cailliez constant =", round(c, digits = 5), "\n"))
  if(cor.zero){
    distmat[distmat > tol] <- distmat[distmat > tol] + c
    distmat <- as.dist(distmat)
  } else {      
    distmat <- as.dist(distmat + c)
  }
  attr(distmat, "call") <- match.call()
  attr(distmat, "method") <- "Cailliez"
  return(list(D = distmat, constant = c))
}
p.cct = function(pvals){
  cct = mean(tan(3.14159*(0.5 - pvals))) 
  return(1 - pcauchy(cct))
}
combine.pvals = function(pvals, type){
  if(type == "min") {
    output = 1 - min(pvals)
  }else if (type == "hmp"){
    output = 1/harmonic.mean(pvals)  
  }else if(type == "Fisher"){
    output = -2*sum(log(pvals))
  }else if(type == "cct"){
    output = mean(tan(3.14159*(0.5 - pvals))) 
  }else if(type == "Stouffer"){
    num.pvals = length(pvals)
    Zvals = qnorm(1-pvals)
    output = sum(Zvals)/num.pvals
  }else{output <- NULL;}
  return(output)
}
perm_final = function(obj, n.perm = 999, centered = TRUE, type = "centroid", distance = NULL, binary = TRUE){ 
  ks = obj$ks
  ws = obj$ws
  rel.abund = obj$rel.abund
  n.g = length(rel.abund)
  
  n.sample = unlist(lapply(rel.abund, nrow))
  group = as.factor(rep(1:n.g, times = n.sample))
  D = obj$D
  if(centered == TRUE){
    cail_obj = cailliez(D)
    D.c = cail_obj$D
    const = cail_obj$constant
    if(!is.euclid(D.c)){stop("Cailliez correction failed, try again")}  
    n.tot = sum(n.sample)
    p.1=diag(n.tot) - 1/n.tot
    D.c = as.matrix(D.c)
    a.c=-0.5*D.c^2
    a.1c=p.1 %*% a.c %*% p.1
    eig.a.1c=eigen(a.1c)
    n.var.c=sum( eig.a.1c$values > 10^-12 )
    v.c=eig.a.1c$vectors[,1:n.var.c] %*% diag( sqrt(eig.a.1c$values[1:n.var.c]) )
    if(type == "centroid" ){
      DATAc <- lm(v.c~group)$res
    }else if (type == "median"){
      spatmed <- ordimedian.pcaPP(v.c, group)
      DATAc <- v.c - spatmed[group,]
    }else{
      stop("Centering type not defined. Use either centroid or median.")
    }
    
    
    process_permutation <- function(per) {
      perm <- sample(nrow(DATAc))
      df.perm <- data.frame(DATAc[perm, ], group)
      
      # Compute phi list for each group
      phi <- list()
      for (i in 1:n.g) {
        temp <- df.perm %>% filter(group == i) %>% dplyr::select(-group)
        phi[[i]] <- as.vector(unlist(dist(temp))) - const
      }
      
      # Initialize vectors for ks and wasserstein distances
      ks.perm <- vector()
      ws.perm <- vector()
      count <- 1
      
      # Calculate ks and wasserstein distances
      for (i in 1:(n.g - 1)) {
        for (j in (i + 1):n.g) {
          obj <- suppressWarnings(ks.test(phi[[i]], phi[[j]]))
          ks.perm[count] <- obj$statistic
          ws.perm[count] <- wasserstein1d(phi[[i]], phi[[j]])
          count <- count + 1
        }
      }
      
      list(ks = ks.perm, ws = ws.perm)
    }
    
    # Apply process_permutation across all permutations in parallel
    results <- mclapply(1:n.perm, process_permutation, mc.cores = detectCores()-2)
    
    # Extract results into ks.p and ws.p vectors
    ks.p <- sapply(results, `[[`, "ks")
    ws.p <- sapply(results, `[[`, "ws")
    if(n.g < 3)
    {
      
      
      pval.ks = (sum(ks < ks.p) + 1)/(n.perm + 1)
      pval.ws = (sum(ws < ws.p) + 1)/(n.perm + 1) 
      pval.ks.p = (frank(-ks.p))/n.perm
      pval.ws.p = (frank(-ws.p))/n.perm
      allpvals = rbind(c(pval.ks, pval.ws), cbind(pval.ks.p, pval.ws.p))
      #combining ks and ws: stat = hmp(pvalue)
      H.p = apply(allpvals, 1, combine.pvals, type = "min")
      final.pval = (sum(H.p[1]<= H.p[-1]) + 1)/(n.perm + 1)
      
    }else{
      pval.ks.orig = (rowSums(ks < ks.p) + 1)/(n.perm + 1)
      pval.ws.orig = (rowSums(ws < ws.p) + 1)/(n.perm + 1)
      ks.perm = apply(ks.p, 1, function(x){(frank(-x))/n.perm})
      ws.perm = apply(ws.p, 1, function(x){(frank(-x))/n.perm})
      all.perm = pmin(ks.perm, ws.perm)
      all.orig = pmin(pval.ks.orig, pval.ws.orig)
      all.A = rbind(all.orig, all.perm)
      #typeS1 = c("min", "hmp", "Fisher", "cct", "Stouffer")
      F.A = apply(all.A,1,combine.pvals, type = "Fisher")
      final.pval = (sum(F.A[1]<= F.A[-1]) + 1)/(n.perm + 1)
    }
    
  }else{
    if(is.null(distance)) stop("Need a valid distance type for permutation")
    DATAr = do.call(rbind, rel.abund)
    ks.p <- ws.p <- c()
    process_permutation <- function(per) {
      perm <- sample(nrow(DATAr))
      df.perm <- data.frame(DATAr[perm, ], group)
      
      # Compute phi list for each group
      phi <- list()
      for (i in 1:n.g) {
        temp <- df.perm %>% filter(group == i) %>% dplyr::select(-group)
        phi[[i]] <- as.vector(unlist(vegdist(temp, method = distance, binary = binary)))
      }
      
      # Initialize vectors for ks and wasserstein distances
      ks.perm <- vector()
      ws.perm <- vector()
      count <- 1
      
      # Calculate ks and wasserstein distances
      for (i in 1:(n.g - 1)) {
        for (j in (i + 1):n.g) {
          obj <- suppressWarnings(ks.test(phi[[i]], phi[[j]]))
          ks.perm[count] <- obj$statistic
          ws.perm[count] <- wasserstein1d(phi[[i]], phi[[j]])
          count <- count + 1
        }
      }
      
      list(ks = ks.perm, ws = ws.perm)
    }
    # Apply process_permutation across all permutations in parallel
    results <- mclapply(1:n.perm, process_permutation, mc.cores = detectCores()-2)
    
    # Extract results into ks.p and ws.p vectors
    ks.p <- sapply(results, `[[`, "ks")
    ws.p <- sapply(results, `[[`, "ws")
    if(n.g < 3)
    {
      
      
      pval.ks = (sum(ks < ks.p) + 1)/(n.perm + 1)
      pval.ws = (sum(ws < ws.p) + 1)/(n.perm + 1) 
      
      #For naming continuity, no cct was used here
      pval.ks.cct = pval.ks
      pval.ws.cct = pval.ws
      
      
      pval.ks.p = (frank(-ks.p))/n.perm
      pval.ws.p = (frank(-ws.p))/n.perm
      allpvals = rbind(c(pval.ks, pval.ws), cbind(pval.ks.p, pval.ws.p))
      #combining ks and ws: stat = hmp(pvalue)
      H.p = apply(allpvals, 1, combine.pvals, type = "min")
      final.pval = (sum(H.p[1]<= H.p[-1]) + 1)/(n.perm + 1)
      
    }else{
      pval.ks.orig = (rowSums(ks < ks.p) + 1)/(n.perm + 1)
      pval.ws.orig = (rowSums(ws < ws.p) + 1)/(n.perm + 1)
      pval.ks.cct = p.cct(pval.ks.orig)
      pval.ws.cct = p.cct(pval.ws.orig)
      
      
      
      ks.perm = apply(ks.p, 1, function(x){(frank(-x))/n.perm})
      ws.perm = apply(ws.p, 1, function(x){(frank(-x))/n.perm})
      all.perm = pmin(ks.perm, ws.perm)
      all.orig = pmin(pval.ks.orig, pval.ws.orig)
      all.A = rbind(all.orig, all.perm)
      #typeS1 = c("min", "hmp", "Fisher", "cct", "Stouffer")
      F.A = apply(all.A,1,combine.pvals, type = "Fisher")
      final.pval = (sum(F.A[1]<= F.A[-1]) + 1)/(n.perm + 1)
    }
    
  }
  
  #output = c(pval.ks, pval.ws, pval.comb)
  return(final.pval)
  
}

DHT.pval = function(rel.abund, n.perm = 500, distance, binary = TRUE)
{
  phi = list()
  n.g = length(rel.abund)
  n.taxa <- n.sample <- vector()
  for(i in 1:n.g)
  {
    n.sample[i] = nrow(rel.abund[[i]])
    n.taxa[i] = ncol(rel.abund[[i]])
  }
  
  for(i in 1:n.g)
  {
    phi[[i]] = as.vector(unlist(vegdist(rel.abund[[i]], method = distance, binary = binary)))
  }
  
  k = 1
  ks = vector()
  ws = vector()
  for(i in 1:(n.g-1))
  {
    for(j in (i+1):n.g)
    {
      obj = suppressWarnings(ks.test(phi[[i]], phi[[j]]))
      ks[k] = obj$statistic
      ws[k] = wasserstein1d(phi[[i]], phi[[j]])
      k = k+1
    }
  }
  #Compute distance matrix
  y = do.call(rbind, rel.abund)
  group = factor(rep(1:n.g, times = n.sample))
  D = vegdist(y, method = distance, binary = binary)
  
  #permanova
  perm_obj = adonis2(D~group)
  p_permanova = perm_obj$`Pr(>F)`[1]
  
  
  obj = list(ks = ks, ws = ws, rel.abund = rel.abund, D = D)
  pval = perm_final(obj, n.perm = n.perm, centered = TRUE, type = "median", distance = distance, binary = binary)
  
  
  
  return(pval)
}
