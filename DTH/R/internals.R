#' @importFrom stats as.dist qnorm ks.test lm
#' @importFrom Rfast Dist
#' @importFrom vegan vegdist
#' @importFrom dplyr select filter
#' @importFrom parallel mclapply detectCores
#' @importFrom ade4 bicenter.wt
#' @importFrom pcaPP l1median_NLM
#' @importFrom transport wasserstein1d
#' @importFrom data.table frank
ordimedian.pcaPP <- function(y, groups){

  inds <- names(table(groups));
  n.inds <- length(inds);
  medians <- matrix(0, nrow = n.inds, ncol = ncol(y))
  rownames(medians) <- inds
  for (i in 1:n.inds) {
    X <- y[groups == inds[i], , drop = FALSE]
    medians[i, ] <- pcaPP::l1median_NLM(X)$par;
  }
  medians;
}
cailliez <- function (distmat, print = FALSE, tol = 1e-07, cor.zero = TRUE) {
  if (ade4::is.euclid(distmat)) {
    warning("Euclidean distance found : no correction needed")
    return(list(D = distmat, constant = 0))
  }
  distmat <- as.matrix(distmat)
  size <- ncol(distmat)
  m1 <- matrix(0, size, size)
  m1 <- rbind(m1, -diag(size))
  m2 <- -ade4::bicenter.wt(distmat * distmat)
  m2 <- rbind(m2, 2 * ade4::bicenter.wt(distmat))
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
combine.pvals = function(pvals, type){
  if(type == "min") {
    output = 1 - min(pvals)
  }else if(type == "Fisher"){
    output = -2*sum(log(pvals))
  }else{output <- NULL;}
  return(output)
}
perm = function(obj, n.perm = 999, type, ncores){
  ks = obj$ks
  ws = obj$ws
  group = obj$ident
  n.sample = as.numeric(table(group))
  n.g = length(n.sample)
  n.total = sum(n.sample)
  D = obj$D
  cail_obj = cailliez(D)
  D.c = cail_obj$D
  const = cail_obj$constant
  if(!ade4::is.euclid(D.c)){stop("Cailliez correction failed, try again")}
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
      temp <- df.perm |> filter(group == i) |> dplyr::select(-group)
      phi[[i]] <- as.vector(unlist(Rfast::Dist(temp, method = "euclidean"))) - const
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
        ws.perm[count] <- transport::wasserstein1d(phi[[i]], phi[[j]])
        count <- count + 1
      }
    }

    list(ks = ks.perm, ws = ws.perm)
  }

  # Apply process_permutation across all permutations in parallel
  results <- mclapply(1:n.perm, process_permutation, mc.cores = ncores)

  # Extract results into ks.p and ws.p vectors
  # ks.p <- sapply(results, `[[`, "ks")
  ks.p <- do.call(c, lapply(results, function(xx) xx[[1]]))
  # ws.p <- sapply(results, `[[`, "ws")
  ws.p <- do.call(c, lapply(results, function(xx) xx[[2]]))
  if(n.g < 3){


    pval.ks = (sum(ks < ks.p) + 1)/(n.perm + 1)
    pval.ws = (sum(ws < ws.p) + 1)/(n.perm + 1)
    pval.ks.p = (data.table::frank(-ks.p))/n.perm
    pval.ws.p = (data.table::frank(-ws.p))/n.perm
    allpvals = rbind(c(pval.ks, pval.ws), cbind(pval.ks.p, pval.ws.p))
    H.p = apply(allpvals, 1, combine.pvals, type = "min")
    final.pval = (sum(H.p[1]<= H.p[-1]) + 1)/(n.perm + 1)

  }else{
    pval.ks.orig = (rowSums(ks < ks.p) + 1)/(n.perm + 1)
    pval.ws.orig = (rowSums(ws < ws.p) + 1)/(n.perm + 1)
    ks.perm = apply(ks.p, 1, function(x){(data.table::frank(-x))/n.perm})
    ws.perm = apply(ws.p, 1, function(x){(data.table::frank(-x))/n.perm})
    all.perm = pmin(ks.perm, ws.perm)
    all.orig = pmin(pval.ks.orig, pval.ws.orig)
    all.A = rbind(all.orig, all.perm)
    F.A = apply(all.A,1,combine.pvals, type = "Fisher")
    final.pval = (sum(F.A[1]<= F.A[-1]) + 1)/(n.perm + 1)
  }



  return(final.pval)

}




