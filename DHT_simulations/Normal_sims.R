##Data simulation where mean of variance is same but distribution is not
#source("~/Downloads/Dispersion/centered_permutation/pval_fun.R")

n.g = 2
n = c(125,25)
p = 500

sim.test = function(mu_hyper, v_hyper)
{
  rel.abund = list()
  n.g = length(mu_hyper)
  for(i in 1:n.g)
  {
    v = exp(rnorm(n[i], mean = mu_hyper[i], sqrt(v_hyper[i])))
    temp = matrix(nrow = n[i], ncol = p)
    for(j in 1:n[i])
    {
      temp[j,] = rnorm(p, mean = 0, sd = sqrt(v[j]))
    }
    rel.abund[[i]] = temp
  }
  return(disp.pval(rel.abund,n.perm = 999, distance = "euclidean", binary = FALSE))
}

mu_hyper = c(1,0)
v_hyper = c(1,3)
replicate(100, sim.test(mu_hyper, v_hyper)) -> out; rowMeans(out <= 0.05)



#Power simulations

n.g = 5
#n = c(50, 40, 30, 20,10)
#n = c(100, 30, 20)
#n = c(125,25)
n=rep(150/n.g,n.g)
p = 500
mu0 = 1
v0 = 1
theta = seq(0,2, by = 0.5)
phi = theta/2

pow = matrix(nrow = length(theta), ncol = 7)
for(case in 1: length(theta))
{
  mu_case = mu0 - theta[case]
  v_case = 2*(v0/2 + theta[case])
  #((mu_case + v_case/2) == (mu0 + v0/2))
  mu_case2 = mu0 - phi[case]
  v_case2 = 2*(v0/2 + phi[case])
  v_hyper = c(v0,v0,v_case2, v0,v_case)
  mu_hyper = c(mu0, mu0,mu_case2,mu0,mu_case)
  replicate(500, sim.test(mu_hyper, v_hyper)) -> out
  pow[case,] = rowMeans(out <= 0.05, na.rm = TRUE)
  cat(pow[case,])
  if(case==1){
    size_mat = t(out)
    colnames(size_mat) = c("GO", "GO-perm","PERMANOVA", "BETADISPER", "KS.original","WS.original","DHT")
    write.csv(size_mat, "norm_size_5_euc_bal.csv")
  }
}


colnames(pow) = c("GO", "GO-perm","PERMANOVA", "BETADISPER", "KS.original","WS.original","DHT")

pow = cbind(theta, pow)

Sys.sleep(10)
beepr::beep(4)
write.csv(pow, "norm_pow_5_euc_bal.csv")

##Size simulations
library(beepr)
n.g = 3
n = c(20,10, 17)
p = 50
mu0 = 1
v0 = 1
mu_hyper = rep(mu0, n.g)
v_hyper = rep(v0, n.g)
replicate(500, sim.test(mu_hyper, v_hyper)) -> pval

size = rowMeans(pval <= 0.05)

write.csv(pval, "norm_size_5_euc_unbal.csv")
Sys.sleep(5)
beep(4)
size

