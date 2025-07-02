#' Function to implement DTH
#'
#' @param dist_mat The distance matrix, usually obtained from dist() or vegdist()
#' @param ident numeric group identifiers, input as factors
#' @param center either "centroid" or "median"
#' @param n.perm number of permutations
#' @param ncores number of cores to be used for computation
#'
#' @returns A p-value for the DTH test
#' @export
#'
#' @examples
#' n1 = 30
#' n2 = 20
#' p = 300
#' y1 = matrix(rnbinom(n1*p, mu =  4, size = 1), nrow = n1, ncol = p)
#' y2 = matrix(rnbinom(n2*p, mu = 4, size = 10), nrow = n2, ncol = p)
#' y = rbind(y1, y2)
#' ident = rep(1:2, times = c(n1,n2))
#' D = vegan::vegdist(y, method = "bray")
#' DTH(D, ident)
DTH = function(dist_mat, ident, center= "centroid", n.perm = 999, ncores = 2)
{
  n.sample = as.numeric(table(ident))
  n.total = sum(n.sample)
  n.g = length(n.sample)
  mat <- as.matrix(dist_mat)
  qqq <- outer(ident, ident, "!=");
  matB <- mat;
  matB[qqq] <- 0
  phi = list()
  n.cs = c(0, cumsum(n.sample)[-n.g])
  for(i in 1:n.g)
  {
    lower = n.cs[i]+1
    upper= n.cs[i] + n.sample[i]
    group_mat = matB[lower:upper, lower:upper]
    phi[[i]] = group_mat[lower.tri(group_mat)]
  }
  k = 1
  ks = vector()
  ws = vector()
  for(i in 1:(n.g-1))
  {
    for(j in (i+1):n.g)
    {
      obj = suppressWarnings(stats::ks.test(phi[[i]], phi[[j]]))
      ks[k] = obj$statistic
      ws[k] = transport::wasserstein1d(phi[[i]], phi[[j]])
      k = k+1
    }
  }
  obj = list(ks = ks, ws = ws, D = dist_mat, ident = ident)
  pval = perm(obj, n.perm = n.perm, type = center, ncores = ncores)
  return(pval)

}
