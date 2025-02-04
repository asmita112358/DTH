##simulation for continuous covariate

sim.test = function(theta)
{
  x = runif(n, 0, 5)
  data = matrix(nrow = n, ncol = p)
  for(i in 1:n)
  {
    v = exp(rnorm(1, mean = -theta*x[i]^2, sqrt(1 + 2*theta*x[i]^2)))
    data[i,] = rnorm(p, mean = 0, sd = sqrt(v))
  }
  binned <- cut(x, breaks = seq(floor(min(x)), ceiling(max(x)), by = 1), labels = FALSE)
  group = as.factor(binned)
  n.g = nlevels(group)
  
  rel.abund = list()
  for(g in 1:n.g)
  {
    rel.abund[[g]] = data[group == g,]
  }
  return(disp.pval(rel.abund,n.perm = 999, distance = "euclidean", binary = FALSE))
}

safe_sim_test <- function(theta) {
  tryCatch({
    # Begin original function logic
    x <- runif(n, 0, 5)
    data <- matrix(nrow = n, ncol = p)
    for (i in 1:n) {
      v <- exp(rnorm(1, mean = -theta * x[i]^2, sqrt(1 + 2 * theta * x[i]^2)))
      data[i, ] <- rnorm(p, mean = 0, sd = sqrt(v))
    }
    binned <- cut(x, breaks = seq(floor(min(x)), ceiling(max(x)), by = 1), labels = FALSE)
    group <- as.factor(binned)
    n.g <- nlevels(group)
    
    rel.abund <- list()
    for (g in 1:n.g) {
      rel.abund[[g]] <- data[group == g, ]
    }
    # Assuming disp.pval() is defined elsewhere in your code
    result <- disp.pval(rel.abund, n.perm = 999, distance = "euclidean", binary = FALSE)
    
    # Return the result
    return(result)
  }, error = function(e) {
    # Handle the error by returning NA
    rep(NA, 18)
  })
}


n = 150
p = 500
theta = seq(0, 0.15, by = 0.025)


pow = matrix(nrow = length(theta), ncol = 7)
for(case in 1: length(theta))
{
 
  replicate(500, safe_sim_test(theta[case])) -> out; 
 
  if(typeof(out) == "double") {
    out2 = out}else{
      out2=do.call(cbind, out)}
  
  pow[case,] = rowMeans(out2 <= 0.05, na.rm = TRUE)
  if(case == 1){
    size_mat = t(out2)
    colnames(size_mat) = c("GO", "GO-perm","PERMANOVA", "BETADISPER", "KS.original","WS.original","DHT")
    write.csv(size_mat, "cont_size.csv")
  }
}
colnames(pow) <- c("GO", "GO-perm","PERMANOVA", "BETADISPER", "KS.original","WS.original","DHT")
pow = cbind(theta, pow)
write.csv(pow, "cont_pow.csv")
beepr::beep(4)
