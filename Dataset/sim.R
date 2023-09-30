##### MAIN CODE #####

# zk : presence of kidney functionality indicator
# p : number of Z1 variables
# q : number of Z2 variables
# v : number of V variables
# rps : if T reward is available at each time-step

source("Environment_functions.R")

main <- function(seed = 1, n_subjects = NULL, n_stages = NULL,  
                 int = F, noise_pred = F, zk = F,
                 treatment_options = 2, 
                 p = NULL, q = NULL, v = NULL, w = 0.05, rps = F) {
  
  pacman::p_load(purrr, dplyr, stringr)
  
  if((p == 0 & q > 0) | (q == 0 & p > 0)){ 
    stop("p and q must both be zero or an integer greater than zero.", call. = F)      
  }
  
  if((v < 1) & (noise_pred == T)){ 
    stop("p, q and v, or at least the parameter v alone, must be an integer greater than zero if noise options is TRUE.", call. = F)      
  }
  
  possible_actions = rep(0:(treatment_options - 1))/(treatment_options - 1)
  
  l = 0
  output = list()
  
  dat <- dat.gen.sim(n_subjects, seed, zk)
  
  if (int) {
    dat$X1 <- runif(n_subjects, min = 0, max = 1)
    dat$X2 <- runif(n_subjects, min = 0, max = 1)
    dat$cW <- replace(dat$cW, dat$X1 > 0.5, 1.5)
    dat$cM <- replace(dat$cM, dat$X2 > 0.5, 1.5)
  }
  
  dat$noise_chng <- rep(0, n_subjects)
  
  if (noise_pred == T) {
    
    if(p & q > 0){
    Z1 <- lapply(1:p, function(x) rnorm(n_subjects, mean = 1))
    Z2 <- lapply(1:q, function(x) rnorm(n_subjects, mean = -1))
    names(Z1) <- paste0("Z", 1:p)
    names(Z2) <- paste0("Z", (p+1):(p+q))
    z.start <- ncol(dat) + 1
    z1.pos <- c(z.start, z.start + p - 1)            
    z2.pos <- c(z.start + p, z.start + p + q - 1)     
    dat <- cbind(dat, Z1, Z2)
    }
    
    V <- lapply(1:v, function(x) rnorm(n_subjects))    
    names(V) <- paste0("V", 1:v)
    
    dat <- cbind(dat, V)
    
    z.var <- dat[,grep("\\Z",names(dat), value = T)]
    
    if(length(z.var) > 0){
      sum_tot <- rowSums(z.var)
      dat$noise_chng <- w * (sum_tot) 
    }

  }
  
  ## 
  out <- data.frame()
  doses_values <- matrix(0, n_subjects, 3)
  
  for (i in 0:(n_stages - 1)) {
    
    if(zk == T & i > 0){
      sum_tot <- rowSums(dat[, names(z.var)])
      dat$noise_chng <- w * (sum_tot) 
    }
    
    dat$month <- rep(i, nrow(dat))
    
    dat$dose <- sample(possible_actions, nrow(dat), replace = T) 
    
    doses_values <- cbind(doses_values, dat$dose)[, -1]
    colnames(doses_values) <- c("t-3","t-2","t-1")
    
    dat$W_next <- updateW(dat$tumor_mass, dat$toxicity, dat$dose, c = dat$cW)
    dat$W_next <- replace(dat$W_next, dat$dropout , NA_real_)
    dat$M_next <- updateM(dat$tumor_mass, dat$toxicity, dat$dose, c = dat$cM)
    dat$M_next <- replace(dat$M_next, dat$dropout , NA_real_)    
    dat$stress_next <- updateStress.sim(dat$tumor_mass, dat$toxicity, dat$stress, doses_values, 
                                    dat$W_next, dat$M_next) 
    dat$stress_next <- replace(dat$stress_next, dat$dropout , NA_real_)
    
    if(zk == T){
      dat$Zk_next <- updateZk(dat$tumor_mass, dat$toxicity, dat$dose, dat$Zk, dat$M_next)
      dat$Zk_next <- replace(dat$Zk_next, dat$dropout , NA_real_)
    }
    
    dat$beta <- 1 / lambda(dat$M_next, dat$W_next, Z = dat$noise_chng)  
    dat$dropout  <- dat$beta < 1 | is.na(dat$beta) 
    
    dat$reward <- log(i + dat$beta) +  1/(i + dat$stress +1)
    
    if(rps == F){
      if (i < (n_stages - 1)) {
        dat$reward <- replace(log(i + dat$beta) +  1/(i + dat$stress +1), !dat$dropout, 0)  
      } 
    }  
    
    
    out <- rbind(out, dat)
    
    if (i < (n_stages - 1)) {
      dat$tumor_mass <- replace(dat$M_next, dat$dropout , NA_real_)
      dat$toxicity <- replace(dat$W_next, dat$dropout , NA_real_)
      dat$stress <- replace(dat$stress_next, dat$dropout , NA_real_)
      
      if(zk == T){
        dat$Zk <- replace(dat$Zk_next, dat$dropout , NA_real_)
      }
      
    }
  }
  
  out$Qhat <- out$reward
  
      # save
      output[[l+1]] = out 
      names(output)[l+1] <- "out"
      l = l+1
  
  # Dropout  #####
  dropout   <- data.frame(ID = out$ID[which(out$dropout == T)], month = out$month[which(out$dropout == T)])
  dup <- which(duplicated(dropout $ID))
  dropout   <- dropout[- dup,]
  
      # save
      output[[l+1]] = dropout 
      names(output)[l+1] <- "dropout"
      l = l+1
  
  
  
  ifelse(zk == T,  no_var <- c("W_next", "M_next", "stress_next","Zk_next", "noise_chng", "cW", "cM" ),
         no_var <- c("W_next", "M_next", "stress_next", "noise_chng", "cW", "cM" ))
  
  d <- out[, ! names(out) %in% no_var, drop = F]
  
  
      # save
      output[[l+1]] = d
      names(output)[l+1] <- "data"
      l = l+1
      
  
  
  # Organized Data Sequence #####
  
  stage <- seq(0, n_stages-1)
  
  data.seq = as.data.frame(matrix(data = NA, nrow = n_subjects))
  
  for (i in stage) {
    stage.i <- subset(d, month == i)   
    
    O.i <-stage.i %>%              
      select(map(c("tumor", "toxi", "str", 'X', 'Z', 'V', "Q"), 
                 starts_with, 
                 vars = colnames(.)) %>% 
               unlist())
    names(O.i) <- c(paste0(names(O.i),"_",i))
    
    if(i == (n_stages-1)){
      O.i$dropout <- stage.i$dropout  
    }
    
    tr = data.frame(dose = stage.i$dose, reward = stage.i$reward )    
    names(tr) <- c(paste0(names(tr),"_",i))
    
    data.seq <-cbind(data.seq, O.i, tr)
  }
  data.seq = data.seq[,-1]
  
  
      # save
      output[[l+1]] = data.seq
      names(output)[l+1] <- "data_sequence"
      l = l+1 
  
  
  
  # Organized DTR (list) #####
  
  stage <- seq(0, n_stages-1)
  
  for (i in stage) {
    stage.i <- subset(d, month == i)   
    
    O.i <-stage.i %>%                    
      select(map(c("tumor", "toxi", "str", 'X', 'Z', 'V'), 
                 starts_with, 
                 vars = colnames(.)) %>% 
               unlist())
    
    
    A.i <- stage.i$dose                
    R.i <- stage.i$reward              
    
    if (i == stage[1]){
      data <- list(O.i, A.i, R.i)
      names(data) <- c(paste0("X",i), paste0("A",i), paste0("R",i))
    }
    else{
      obj <- list(O.i, A.i, R.i)
      names(obj) <- c(paste0("X",i), paste0("A",i), paste0("R",i))
      
      data <-c(data, obj)
    }
  }
  
      # save
      output[[l+1]] = data
      names(output)[l+1] <- "DTR_list"
      l = l+1
  
  n_patients_state = ncol(data[[1]])
  
  
  # Organized DTR (data) #####
  
  stage <- seq(0, n_stages-1)
  
  data = as.data.frame(matrix(0, nrow = n_subjects))
  
  for (i in stage) {
    stage.i <- subset(d, month == i)   
    
    O.i <-stage.i %>%                      
      select(map(c("tumor", "toxi", "str", 'X', 'Z', 'V'), 
                 starts_with, 
                 vars = colnames(.)) %>% 
               unlist())
    
    
    A.i <- stage.i$dose            
    R.i <- stage.i$reward          
    
    si = cbind(O.i, A.i, R.i)
    names(si) <- c(paste0(names(O.i),paste0("_",i)), paste0("A",i), paste0("R",i))
    data = cbind(data, si)
    
  }
  
  
      # save
      output[[l+1]] = data[,-1]
      names(output)[l+1] <- "DTR_data"
      l = l+1
  
  
  return(output)
  
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#                                                                                       #
# n_subjects = 20                                                                       #
# w = 0.05                                                                              #
# n_stages = 3                                                                          #
# int = F                                                                               #
# noise_pred = T                                                                        #
# zk = F                                                                                #
# treatment_options = 5                                                                 #
# p = 1                                                                                 #
# q = 1                                                                                 #
# v = 1                                                                                 #
#                                                                                       #
# example <- main(n_subjects = n_subjects, seed = 123,                                  #
#           n_stages = n_stages, int = int, noise_pred = noise_pred, zk = zk,           #
#           treatment_options = treatment_options,                                      #
#           p = p, q = q, v = v)$out                                                    #
#                                                                                       #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
