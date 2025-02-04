##Simulations for dispersion, all on IBD data template
rm(list = ls())
library(gtools)
library(transport)
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(sFFLHD)
library(MIDASim)
library(GUniFrac)
library(MDMR)
library(vegan)
library(tidyr)
library(HMP)
source("~/Downloads/Dispersion/Dispersion_codes/all_funs.R")
load("IBD.Vaginal.simulated.Rdata", verbose = T)

Z1 <- (count.ibd > 0)*1
Z2 <- (count.ibd.sim.sparsedossa > 0)*1

##Size sims




n.sim = 500
pval = matrix(nrow = n.sim, ncol= 5)

for(i in 1:n.sim)
{
  sim_Z1 = generate.Z(Z1)
  sim_Z2 = generate.Z(Z1)
  
  pval[i,] = sim.fcn(list(sim_Z1, sim_Z2), distance = "jaccard", binary = TRUE)
  print(pval[i,])
}
write.csv(pval, "jaccard.csv")
size = apply(pval <= 0.05, 2, mean)
se_size = sqrt(size*(1 - size)/ n.sim)



##Power Sims


alph = seq(0,1,by = 0.1)

rej = matrix(nrow = length(alph), ncol = 7)
n.sample = nrow(Z1)
for(i in 1:length(alph))
{
  
  dil = alph[i]
  n.sim = 200
  pval = matrix(nrow = n.sim, ncol= 7)
 
  for(sim in 1:n.sim)
  {
    
    sim_Z1 = generate.Z(Z1)
    sim1 = generate.Z(Z1)
    sim2 = generate.Z(Z2)
    a = sample(c(0,1), n.sample, replace = TRUE, prob = c(1-dil, dil))
    sim_Z2 = a*sim1 + (1-a)*sim2
    pval[sim,] = disp.pval(list(sim_Z1, sim_Z2), n.perm = 999, distance = "jaccard", binary = TRUE)
    
  }
  if(i == 11){
    write.csv(pval, "IBDsim_size.csv")
  }
  rej[i,] = apply(pval<= 0.05,2,mean, na.rm = TRUE)
  print(rej[i,])
  gc()
  
}
reject = data.frame(alph,rej)
colnames(reject) <- c("dilution", "GO","GO.perm","PERMANOVA", "Betadisper", "KS","WS","DHT")


write.csv(reject, "IBD_pow.csv")

beepr::beep(4)


