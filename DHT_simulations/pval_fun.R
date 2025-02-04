
#source("~/Downloads/Dispersion/centered_permutation/disp_perm.R")
source("~/Downloads/Dispersion/pvalue_combinations_v2/disp_perm.R")
source("~/Downloads/Dispersion/centered_permutation/GO_funs.R")
library(dplyr)
library(tidyr)

#function to generate pvalues from all competing methods
disp.pval = function(rel.abund, n.perm = 500, distance, binary = TRUE)
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
  pval = perm(obj, n.perm = n.perm, centered = TRUE, type = "median", distance = distance, binary = binary)
  
  
  
  #betadisper
  beta_obj = betadisper(D,  group = group, bias.adjust = TRUE)
  p_beta = permutest(beta_obj)$tab[1,6]
  
  
  
  obj.fd = Fd.test(D, group, input = "distances", centring = "median")
  p_go = obj.fd$pval
  p_go_perm = obj.fd$pval.perm
  
  #kmmd
  #obj_kmmd = kmmd(matrix(var1), matrix(var2))
  
  pval_vec = c(p_go,p_go_perm,p_permanova, p_beta, pval)
  #names(pval_vec) = c("GO", "GO-perm","PERMANOVA", "BETADISPER", "KS","WS","min+cct", "min+hmp", "F+sum", "F+min")
  
  
  return(pval_vec)
}



