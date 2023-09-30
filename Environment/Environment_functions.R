# Initialization data 

dat.gen <- function(zk) { 
  data = data.frame(
    tumor_mass = runif(1, min = 0, max = 2), # random starting tumor masses
    toxicity = runif(1, min = 0, max = 2),   # random starting toxicity measure
    stress = runif(1, min = 0, max = 2)      # random starting stress measure
  )
  if(zk) {
    data$Zk = abs(rnorm(1, mean = 0, sd = 1)) 
  }   
  return(data)
}

# Update functions

updateW <- function(M, W, D, a = 0.1, b = 1.2, c = 1, d = 0.5) {
  W_next <- a * M + b * (c * D - d) + W
  return(replace(W_next, W_next < 0, 0))
}

updateM <- function(M, W, D, a = 0.15, b = 1.2, c = 1, d = 0.5) {
  M_next <- a * W - b * (c * D - d) + M
  return(replace(M_next, M <= 0 | M_next < 0, 0))
}

updateZk <- function(M, W, D, Zk, M_next, a = 0.1, b = 1.2, c = 1, d = 0.5, e = 0.1){
  Zk_next <- - a * W + b * (c * D - d) - e * (M_next - M) + Zk
  return(replace(Zk_next, Zk_next < 0, 0))
}

lambda <- function(M, W, Z = 0, mu0 = -5.5, mu1 = 1, mu2 = 1.2, mu3 = 0.75) {
  return(exp(mu0 + mu1 * W + mu2 * M + mu3 * W * M + Z))
}

updateStress <- function(M, W, STR, doses_values, W_next, M_next){
  stress_next <- STR + 0.1 * (M_next - M) + 0.2 * (W_next - W) + 
    0.8 * doses_values[3] + 0.5 * doses_values[2] + 0.3 * doses_values[1]
  return(replace(stress_next, stress_next < 0, 0)) 
}


# sim.R functions 

updateStress.sim <- function(M, W, STR, doses_values, W_next, M_next){
  stress_next <- STR + 0.1 * (M_next - M) + 0.2 * (W_next - W) + 
    0.8 * doses_values[,3] + 0.5 * doses_values[,2] + 0.3 * doses_values[,1]
  return(replace(stress_next, stress_next < 0, 0))
}

dat.gen.sim <- function(n, seed, zk) { 
  set.seed(seed)
  data = data.frame(
    ID = seq_len(n),
    tumor_mass = runif(n, min = 0, max = 2),  # random starting tumor masses
    toxicity = runif(n, min = 0, max = 2),    # random starting toxicity measure
    stress = runif(n, min = 0, max = 2),      # random starting stress measure
    cW = rep(1, n),
    cM = rep(1, n)
  )
  if(zk) {
    data$Zk = abs(rnorm(n, mean = 0, sd = 1)) # random starting kidney functionality indicator
  }   
  return(data)
}




