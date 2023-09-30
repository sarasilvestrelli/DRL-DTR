pacman::p_load(dplyr, kableExtra, xtable)
rm(list=ls())


# scenarios parameters ----
R = NA       # Reward
beta = NA    # Survival times
method = ""

# exploration strategies
expl_strategy = c("eps_greedy", "softmax")

# update strategies
updt_strategy = c("hard", "soft")

### MODELS ###

# low  --- no noise or interaction
N1 = c(200, 600, 1000)
T1 = c(3, 8, 16)
K1 = c(3, 6, 15)

p1 = 0
q1 = 0
v1 = 0

int1 = F
noise_pred1 = F
zk1 = F

diff1 ="low"


# medium   --- interaction and noise variables
N2 = c(200, 600, 1000)
T2 = c(3, 8, 16)
K2 = c(3, 6, 15)

p2 = 0
q2 = 0
v2 = 30

int2 = T
noise_pred2 = T
zk2 = F

diff2 ="medium"

# high  ---  (p+q+v+zk) noise, interaction and RF variable
N3 = c(200, 600, 1000)
T3 = c(3, 8, 16)
K3 = c(3, 6, 15)

p3 = 5
q3 = 5
v3 = 90

int3 = T
noise_pred3 = T
zk3 = T

diff3 ="high"


# parameters combination

grid_low<-expand.grid(
  complexity = diff1,
  method = method,
  n_cov = 3 + 2*int1 + p1+q1+v1 + zk1,     
  patients=N1,
  stages=T1,
  treatment_options=K1,
  p=p1,
  q=q1,
  v=v1,
  RF = zk1,
  int=int1,
  noise_pred = noise_pred1,
  expl_strategy = expl_strategy,
  updt_strategy = updt_strategy,
  R = R,
  beta = NA
  
  
)

grid_medium<-expand.grid(
  complexity = diff2,
  method = method,
  n_cov = 3 + 2*int2 + p2+q2+v2 + zk2,    
  patients=N2,
  stages=T2,
  treatment_options=K2,
  p=p2,
  q=q2,
  v=v2,
  RF = zk2,
  int=int2,
  noise_pred=noise_pred2,
  expl_strategy = expl_strategy,
  updt_strategy = updt_strategy,
  R = R,
  beta = NA
  
  
)

grid_high<-expand.grid(
  complexity = diff3, 
  method = method,
  n_cov = 3 + 2*int3 + p3+q3+v3 + zk3,
  patients=N3,
  stages=T3,
  treatment_options=K3,
  p=p3,
  q=q3,
  v=v3,
  RF = zk3,
  int=int3,
  noise_pred=noise_pred3,
  expl_strategy = expl_strategy,
  updt_strategy = updt_strategy,
  R = R,
  beta = NA
  
)

grid <-rbind(grid_low, grid_medium, grid_high)

# clean and save ------------------------------
rm(list=grep("^grid", ls(), value = TRUE, invert = TRUE))
rm(list = ls()[grep("^grid_", ls())])

save.image(".//grid_parameters.RData")