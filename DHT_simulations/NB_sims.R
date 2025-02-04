#NB sims
library(beepr)
library(invgamma)


sim.test = function(alpha, beta)
{
  rel.abund = list()
  
  for(i in 1:n.g)
  {
    #size = exp(rnorm(n[i], mean = mu_hyper[i], sqrt(v_hyper[i])))
    size = rinvgamma(n[i], shape = alpha[i], scale = 1/beta[i])
    temp = matrix(nrow = n[i], ncol = p)
    for(j in 1:n[i])
    {
      temp[j,] = rnbinom(p, size[j], mu = 2)
      #zero_mask <- rbinom(n, size = 1, prob = 1 - pi)
    }
    #data_mat = samp_nb(n[i], p = 50, pi = 0.5, r = r, mu = 20)
    #data_mat[data_mat <= 20] <- 0
    rel.abund[[i]] = temp
  }
  return(disp.pval(rel.abund,n.perm = 999, distance = "bray", binary = FALSE))
  #return(sim.fcn(rel.abund, n.perm = 999, distance= "euclidean", binary = FALSE))
}

n.g = 2
#n = c(50, 40, 30, 20,10)
#n = c(100, 30, 20)
#n = c(125,25)
n=rep(150/n.g,n.g)
p = 500
theta = seq(0,3, by = 0.5); 5*exp(-theta)
phi = theta/2; 5*exp(-phi)
pval = matrix(NA,nrow = length(theta), ncol = 7)
a0 = 5; b0 = 5
for(case in 1:length(theta))
{
  
  b_case = 5*exp(-theta[case])
  b_case2 = 5*exp(-phi[case])
  beta = c(b0,b_case)
  alpha = beta
  
  replicate(500, sim.test(alpha, beta)) -> out;
  pval[case,] = rowMeans(out <= 0.05)
  if(case==1){
    size_mat = t(out)
    colnames(size_mat) = c("GO", "GO-perm","PERMANOVA", "BETADISPER", "KS.original","WS.original","DHT")
    write.csv(size_mat, "NB_size_2_bray_bal.csv")
  }
  gc()
}

#colnames(pval) = c("GO", "GO-perm","PERMANOVA", "BETADISPER", "KS.original","WS.original","min+hmp", "min+F", "min+cct", "hmp+hmp", "hmp+F", "hmp+cct", "F+hmp", "F+F","F+cct","cct+hmp", "cct+f", "cct+cct")
colnames(pval) = c("GO", "GO-perm","PERMANOVA", "BETADISPER", "KS.original","WS.original","DHT")
pow_nb = cbind(theta, pval)
Sys.sleep(10)
beepr::beep(4)
write.csv(pow_nb, "NB_pow_2_bray_bal.csv")
pow_nb


































































































































for(i in 1:n.g)
{
  #size = exp(rnorm(n[i], mean = mu_hyper[i], sqrt(v_hyper[i])))
  size = rinvgamma(n[i], shape = alpha[i], scale = 1/beta[i])
  print(mean(1/size))
  print(var(1/size))
}

Sys.sleep(5)
beep(4)


##Simulations for size
n.g = 2
n = c(20, 20)
p = 50
alpha = rep(5,n.g)
beta = rep(5, n.g)
replicate(500, sim.test(alpha, beta)) -> pval
size = rowMeans(pval <= 0.05)

write.csv(pval, "NB_size_2_bray_bal.csv")
Sys.sleep(5)
beep(4)