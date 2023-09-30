# MAIN CODE #

# SARA ####
#### Settings ####

rm(list=ls())
pacman::p_load(tictoc, keras, stringr, xtable, dplyr, R6, qtl2, purrr, data.table, tictoc, mltools, 
               magrittr, readr, tibble, openxlsx, extrafont, showtext, ggplot2, forcats, readxl, gridExtra)


setwd(getwd())                          # set the current directory
source("grid_parameters.R")             # load the grid of possible parameter combinations into the data
source("tr_ev_data.R")                  # Inserts initial simulation data for training and validation for each combination into a folder


loadfonts(device = "win")

dir.create(".\\Performance\\QL", recursive = T) 
dir.create(".\\Performance\\DQN", recursive = T)  
dir.create(".\\Performance\\RQN", recursive = T)  
dir.create(".\\Performance\\LSTMQN", recursive = T) 

n_new = 500  # evaluation sample
current_directory <- getwd()

# APPROACHES #

############################################# QL #############################################

rm(list=ls()[!ls() %in% c("grid", "g", "current_directory", "n_new")])

method = "QL"

source("sim.R")
source("train_ql_ols.R")

g[, "method"] <- method

MEMORY_WINDOW = 4L   

w = 0.05

beta_QL <- data.frame(matrix(ncol = nrow(g), nrow = n_new)) 
reward_QL <- data.frame(matrix(ncol = nrow(g), nrow = n_new)) 
colnames(beta_QL) <- colnames(reward_QL) <- c(paste0("comb_", 1:nrow(g))) 

tic.clearlog()

for (j in 1:nrow(g)) {
  
  tic(j)
  Sys.sleep(1)
  
  n_subjects = g[j, "patients"]         
  n_stages = g[j, "stages"]            
  
  p = g[j, "p"]              
  q = g[j, "q"]                 
  v = g[j, "v"]    
  
  zk = g[j, "RF"]
  
  comp <- g[j, "complexity"]
  
  treatment_options = g[j, "treatment_options"]   
  possible_actions = rep(0:(treatment_options-1))/(treatment_options-1)
  
  int = g[j, "int"]                
  noise_pred = g[j, "noise_pred"] 
  
  n_cov = g[j, "n_cov"]
  
  # TRAIN ON 200, 600, 1000 patients # 
  seed = 123
  
  data <- main(seed = seed, n_subjects = n_subjects, n_stages = n_stages , 
               int = int, noise_pred = noise_pred, zk = zk,
               treatment_options = treatment_options, p = p, q = q, v = v, rps = F)$data_sequence
  
  tr <- train_ql_ols(data, n_stages, p, q, v, treatment_options, possible_actions, seed, int, noise_pred, zk)
  
  mod_list <- tr$mod_list
  predictors_name <- tr$predictors_name
  
  rm(data,tr)
  
  # Q LEARNING simulation # 
  seed = 1
  data_q <- as.data.frame(sapply(read_csv(paste0(current_directory, "/eval_dataset_QL/ev", j, ".csv"),
                                          show_col_types = FALSE), as.numeric))
  z.var <- data_q[, grep("\\Z", names(data_q), value = T)]
  
  d <- data_q
  
  # define the initial stage of the main variables 
  var_chng <- c("tumor_mass", "toxicity", "stress", grep("\\X|\\Z|\\V", names(d), value = T))
  id.var_chng <- which(names(d) %in% var_chng)
  
  names(d)[id.var_chng] <- paste0(var_chng, "_0")
  
  out <- data.frame()
  
  # 
  doses_values <- matrix(0,n_new,3)
  len_cov <- (3+int*2+ p+q+v+ zk + treatment_options)
  history = as.data.frame(matrix(0, nrow = n_new, ncol= (MEMORY_WINDOW * len_cov)))
  
  
  for (i in 0:(n_stages - 1)) {
    
    rbind_var <- c("ID", "tumor_mass", "toxicity", "stress", "cW", "cM", "Zk", "X1", "X2", "noise_chng", paste0("Z", 1:(p+q)), paste0("V", 1:v),
                   "month", paste0("dose_", i), "pred", "W_next", "M_next", "stress_next", "Zk_next", "beta", "dropout", "reward")
    
    if (i > 0) {
      
      doses_values <- doses_values[!(d$dropout), ]
      d <- d[!(d$dropout), ]
      
      d$tumor_mass <- d$M_next
      d$toxicity <- d$W_next
      d$stress <- d$stress_next
      
      if(zk){
        d$Zk <- d$Zk_next
        sum_tot <- rowSums(d %>% select(names(z.var)))
        d$noise_chng <- w * (sum_tot) 
      }
      
      names(d)[id.var_chng] <- paste0(var_chng, "_", i)
    }
    
    n <- nrow(d)
    d$month <- rep(i, n)
    
    d_expand <- d[rep(seq_len(n), each = length(possible_actions)), ]
    d_expand[, paste0("dose_",i)] <- rep(possible_actions, n)
    
    m <- mod_list[[paste0('month', i)]]
    predictors <- predictors_name[[paste0('month', i)]]
    train <- cbind(history[rep(seq_len(n), each = length(possible_actions)), ], d_expand)
    
    # Q function prediction
    d_expand <- data.table::as.data.table(d_expand) 
    d_expand$pred <- m %>% predict(subset(train, select = predictors))
    
    d <- d_expand[d_expand[, .I[which.max(pred)], by = ID]$V1]
    
    # update history
    history <- subset(train, select = predictors)
    
    names(d)[id.var_chng] <- var_chng
    
    old_M <- d$tumor_mass
    old_W <- d$toxicity
    old_stress <- d$stress
    
    dose <- d[, get(paste0("dose_", i))]
    doses_values <- cbind(doses_values, as.numeric(as.character(dose)))[, -1]
    
    d$W_next <- updateW(old_M, old_W, dose, c = d$cW)
    d$W_next <- replace(d$W_next, d$dropout, NA_real_)
    d$M_next <- updateM(old_M, old_W, dose, c = d$cM)
    d$M_next <- replace(d$M_next, d$dropout, NA_real_)
    d$stress_next <- updateStress(old_M, old_W, old_stress, doses_values, d$W_next, d$M_next)
    d$stress_next <- replace(d$stress_next, d$dropout, NA_real_)
    
    if(zk){
      d$Zk_next <- updateZk(d$tumor_mass, d$toxicity, dose, d$Zk, d$M_next)
      d$Zk_next <- replace(d$Zk_next, d$dropout, NA_real_)
    }
    
    d$beta <- 1 / lambda(d$M_next, d$W_next, Z = d$noise_chng)
    d$dropout <- d$beta < 1 
    d$reward <- log(i + d$beta) +  1/(i + d$stress +1)
    
    if (i < (n_stages - 1)) {   
      d$reward <- replace(log(i + d$beta) +  1/(i + d$stress +1), !d$dropout, 0)
    }
    
    d_sel <- d[, which(names(d) %in% rbind_var), with=F]
    names(d_sel)  <- gsub(paste0("dose_", (i)), "dose", names(d_sel))
    out <- rbind(out, d_sel)
    
    # cat(
    #   sprintf(
    #     "[%4d of %d] mean rewards: %4f, mean survival rate: %.3f\n",
    #     i,
    #     n_stages-1,
    #     mean(d$reward),
    #     mean(d$beta)
    #   )
    # )
  }
  
  
  R.id <- aggregate(x = out$reward,                
                    by = list(out$ID),              
                    FUN = sum)  
  
  betas <- rbind(subset(out, month < (n_stages-1) & dropout==T), subset(out, month == (n_stages-1)))
  betas <- betas[order(betas$ID, decreasing = F), ] 
  
  g$R[j]<- mean(R.id$x)
  g$beta[j]<- mean(betas$beta)
  
  beta_QL[,j] <- betas$beta 
  reward_QL[,j] <- R.id$x 
  
  cat(
    sprintf(
      "[%4d / %d - %s ] mean rewards: %4f, mean survival rate: %.3f\n",
      j,
      nrow(g),
      comp,
      mean(R.id$x),
      mean(betas$beta)
    )
  )
  
  toc(log = TRUE, quiet = TRUE)
}


log.lst <- tic.log(format = FALSE)
log.txt <- tic.log(format = TRUE)

timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))

# mean(timings)/60
# sum(timings)/60
# writeLines(unlist(log.txt))

g$time <- as.numeric(timings)/60
g$tot_time <- sum(timings)/60

updt_strategy <- expl_strategy <- rep("none", nrow(g))
g <- g %>%
  add_column(expl_strategy, updt_strategy, .after = "noise_pred")

dir <- paste0(current_directory, "/Performance/", method) 
write.xlsx(g, paste0(dir,'/grid_QL.xlsx')) 
write.xlsx(beta_QL, paste0(dir,'/beta.xlsx')) 
write.xlsx(reward_QL, paste0(dir,'/reward.xlsx')) 

save.image(".//QL.RData")

############################################# DQN #############################################

rm(list=ls()[!ls() %in% c("grid", "current_directory", "n_new")])

method = "DQN"
grid[, "method"]  <- method

MEMORY_WINDOW = 4L  

BATCH_SIZE = 32    

beta_DQN <- data.frame(matrix(ncol = nrow(grid), nrow = n_new)) 
reward_DQN <- data.frame(matrix(ncol = nrow(grid), nrow = n_new)) 
CR_DQN <- data.frame(matrix(0, ncol = max(grid$patients), nrow = nrow(grid))) #?
colnames(beta_DQN) <- colnames(reward_DQN) <- c(paste0("comb_", 1:nrow(grid))) 

source("Agent.R")
source("Memory.R")
source("train.R")
source("Environment.R")

#
tic.clearlog()

for (j in 1:nrow(grid)) { 
  tic(j)
  Sys.sleep(1)
  
  # eps_greed = F <- Boltzmann strategy, epsilon greedy strategy otherwise
  expl = grid[j, "expl_strategy"] 
  if(!(expl %in% c("eps_greedy", "softmax"))){ 
    stop("Exploration strategy can only be 'eps_greedy' or 'softmax'. ", call. = F)      
  }
  
  # target update strategy
  update = grid[j, "updt_strategy"]
  if(!(update %in% c("hard", "soft"))){ 
    stop("Target update strategy can only be 'hard' or 'soft'. ", call. = F)      
  }
  
  n_subjects = grid[j,"patients"]         
  n_stages = grid[j,"stages"]            
  
  p = grid[j,"p"]              
  q = grid[j,"q"]                 
  v = grid[j,"v"]    
  
  zk = grid[j,"RF"]
  int = grid[j,"int"]                
  noise_pred = grid[j,"noise_pred"] 
  
  comp <- grid[j,"complexity"]
  
  treatment_options = grid[j,"treatment_options"]
  
  tr <- as.data.frame(sapply(read_csv(paste0(".//tr_dataset/tr", j, ".csv"), show_col_types = FALSE), as.numeric))
  ev <- as.data.frame(sapply(read_csv(paste0(".//eval_dataset/ev", j, ".csv"), show_col_types = FALSE), as.numeric))
  
  # 
  
  N_EPISODE = n_subjects + n_new
  tr_phase = T  # parameter set if an epsilon greedy or softmax exploration strategy is wanted during training phase
  
  env <- Environment$new(seed = 123,
                         int = int,  
                         noise_pred = noise_pred,
                         zk = zk,
                         treatment_options = treatment_options,  
                         p = p,
                         q = q,
                         v = v,
                         n_stages = n_stages)
  
  obs_size <- length(c("tumor_mass","toxicity","stress",
                       grep("\\X|\\Z|\\V", names(tr), value = T)))
  
  agent <-
    Agent$new(
      input_shape = c(MEMORY_WINDOW*obs_size),
      output_dim = env$treatment_options,
      decay_last_episode = 100,
      C = 20,
      tau_softmax = .5,
      tau_soft_updt = 1e-3, 
      updt = update
    )
  
  memory <- Memory$new(capacity = 20000)
  
  rewards <- c()
  betas <- c()
  
  for (episode_i in 1:N_EPISODE) {
    CR <- 0 #?
    if(episode_i <= n_subjects){
      
      s <- env$reset(tr[episode_i, ])
      
    } else{
      
      if (episode_i == n_subjects + 1){
        val_i <- 1
        tr_phase <- F
      } else{
        val_i <- val_i + 1
      }
      
      s <- env$reset(ev[val_i, ])
    }
    
    done <- FALSE
    
    history <- rep(0., MEMORY_WINDOW * obs_size)
    history <- append(history, s)[-c(1:obs_size)]
    
    doses_values <- rep(0,3) 
    
    while (!done) {
      
      a <- agent$get_action(state_ = array(history, dim = c(1, MEMORY_WINDOW*obs_size)), 
                            step = episode_i, tr_phase , expl_str = expl)
      treatment <- as.numeric(a/(env$treatment_options - 1)) 
      doses_values <- append(doses_values, treatment)[-1] 
      ret <- env$step(treatment, doses_values)   
      s2 <- ret$state2
      history2 = append(history, s2)[-c(1:obs_size)]
      done <- ret$dropout
      
      r <- 0 
      
      if(done){
        
        r <- ret$reward  
        
        beta <- ret$beta
        
        betas <- c(betas, beta)
        rewards <- c(rewards, r)
      }
      
      if(episode_i <= n_subjects){
        CR <- CR + ret$reward #?
        memory$push(history, a, r, done, history2)
        
        if (memory$length > BATCH_SIZE) {
          batch <- memory$sample(BATCH_SIZE)
          train(agent, batch, env$i)
        }
        
      }
      s <- s2
      history <- history2
    }
    
    if(episode_i <= n_subjects){ #?
      CR_DQN[j,episode_i] <- CR #?
    } #?
    
    cat(
      sprintf(
        "[%d of %d - %s] %s:  %4d of %d - \t R: %4f \t S: %4f \t %s: %.3f\n",  
        j,
        nrow(grid),
        comp,
        ifelse(episode_i <= n_subjects, "training", "validation"),
        episode_i,
        N_EPISODE,
        r,
        beta,
        ifelse(expl == "eps_greedy", "Epsilon", "Tau"),
        ifelse(expl == "eps_greedy", agent$epsilon, agent$tau_softmax)
      ))
    
    flush.console()
  }
  
  grid$R[j]<- mean(rewards[(n_subjects+1):(n_new + n_subjects)])
  grid$beta[j]<- mean(betas[(n_subjects+1):(n_new + n_subjects)])
  
  beta_DQN[,j] <- betas[(n_subjects+1):(n_new + n_subjects)] 
  reward_DQN[,j] <- rewards[(n_subjects+1):(n_new + n_subjects)] 
  
  cat(
    sprintf(
      "[%4d of %d - %s ] R_mean: %4f \t S_mean: %.3f\n",
      j,
      nrow(grid),
      comp,
      mean(rewards),
      mean(betas)
    )
  )
  
  toc(log = TRUE, quiet = TRUE)
  
}

log.lst <- tic.log(format = FALSE)
log.txt <- tic.log(format = TRUE)

timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))

# mean(timings)/60 
# sum(timings)/60 
# writeLines(unlist(log.txt))

grid$time <- as.numeric(timings)/60
grid$tot_time <- sum(timings)/60

dir <- paste0(current_directory, "/Performance/", method) 
write.xlsx(grid, paste0(dir,'/grid_DQN.xlsx')) 
write.xlsx(beta_DQN, paste0(dir,'/beta.xlsx')) 
write.xlsx(reward_DQN, paste0(dir,'/reward.xlsx')) 
write.xlsx(CR_DQN, paste0(dir,'/CR.xlsx')) #?

save.image(".//DQN.RData")

############################################# RQN #############################################

rm(list=ls()[!ls() %in% c("grid", "current_directory", "n_new")])

method = "RQN"
grid[, "method"]  <- method

# settings 
MEMORY_WINDOW = 4L

BATCH_SIZE = 32  

beta_RQN <- data.frame(matrix(ncol = nrow(grid), nrow = n_new)) 
reward_RQN <- data.frame(matrix(ncol = nrow(grid), nrow = n_new)) 
CR_RQN <- data.frame(matrix(0, ncol = max(grid$patients), nrow = nrow(grid))) #?
colnames(beta_RQN) <- colnames(reward_RQN) <- c(paste0("comb_", 1:nrow(grid))) 


source("RAgent.R")
source("Memory.R")
source("train.R")
source("Environment.R")

#
tic.clearlog()

for (j in 1:nrow(grid)) {
  tic(j)
  Sys.sleep(1)
  
  # eps_greed = F <- Boltzmann strategy, epsilon greedy strategy otherwise
  expl = grid[j, "expl_strategy"] 
  if(!(expl %in% c("eps_greedy", "softmax"))){ 
    stop("Exploration strategy can only be 'eps_greedy' or 'softmax'. ", call. = F)      
  }
  
  # target update strategy
  update = grid[j, "updt_strategy"]
  if(!(update %in% c("hard", "soft"))){ 
    stop("Target update strategy can only be 'hard' or 'soft'. ", call. = F)      
  }   
  
  n_subjects = grid[j, "patients"]         
  n_stages = grid[j, "stages"]            
  
  p = grid[j, "p"]              
  q = grid[j, "q"]                 
  v = grid[j, "v"]    
  
  zk = grid[j, "RF"]
  int = grid[j, "int"]                
  noise_pred = grid[j, "noise_pred"] 
  
  comp <- grid[j, "complexity"]
  
  treatment_options = grid[j, "treatment_options"]
  
  tr <- as.data.frame(sapply(read_csv(paste0(".//tr_dataset/tr", j, ".csv"), show_col_types = FALSE), as.numeric))
  ev <- as.data.frame(sapply(read_csv(paste0(".//eval_dataset/ev", j, ".csv"), show_col_types = FALSE), as.numeric))
  
  # 
  
  N_EPISODE = n_subjects + n_new  
  tr_phase = T  # parameter set if an epsilon greedy or softmax exploration strategy is wanted during training phase
  
  env <- Environment$new(seed = 123,
                         int = int,  
                         noise_pred = noise_pred,
                         zk = zk,
                         treatment_options = treatment_options,  
                         p = p,
                         q = q,
                         v = v,
                         n_stages = n_stages)
  
  obs_size <- length(c("tumor_mass", "toxicity", "stress",
                       grep("\\X|\\Z|\\V", names(tr), value = T)))
  
  agent <-
    Agent$new(
      input_shape = c(MEMORY_WINDOW, obs_size),
      output_dim = env$treatment_options,
      decay_last_episode = 100,
      C = 20,
      tau_softmax = .5,
      tau_soft_updt = 1e-3, 
      updt = update
    )
  
  memory <- Memory$new(capacity = 20000)
  
  rewards <- c()
  betas <- c()
  
  for (episode_i in 1:N_EPISODE) {
    CR <- 0 #?
    if(episode_i <= n_subjects){
      
      s <- env$reset(tr[episode_i, ])
      
    } else{
      
      if (episode_i == n_subjects + 1){
        val_i <- 1
        tr_phase <- F
      } else{
        val_i <- val_i + 1
      }
      
      s <- env$reset(ev[val_i, ])
    }
    
    done <- FALSE
    
    history = rep(0., MEMORY_WINDOW * obs_size)
    history = append(history, s)[-c(1:obs_size)]
    
    doses_values <- rep(0,3)    
    
    while (!done) {
      
      a <- agent$get_action(state_ = array(history, c(1, MEMORY_WINDOW, obs_size)), 
                            step = episode_i, tr_phase , expl_str = expl)
      treatment <- as.numeric(a/(env$treatment_options - 1))  
      doses_values <- append(doses_values, treatment)[-1] 
      ret <- env$step(treatment, doses_values)  
      s2 <- ret$state
      history2 = append(history, s2)[-c(1:obs_size)]
      done <- ret$dropout
      
      r <- 0 
      
      if(done){
        
        r <- ret$reward  
        
        beta <- ret$beta
        
        betas <- c(betas, beta)
        rewards <- c(rewards, r)
      }
      
      if(episode_i <= n_subjects){
        CR <- CR + ret$reward #?
        memory$push(history, a, r, done, history2)
        
        if (memory$length > BATCH_SIZE) {
          batch <- memory$sample(BATCH_SIZE)
          train(agent, batch, env$i)
        }
        
      }
      
      
      s <- s2
      history <- history2
      
      if(episode_i <= n_subjects){ #?
        CR_RQN[j,episode_i] <- CR #?
      } #?
      
    }
    
    cat(
      sprintf(
        "[%d of %d - %s] %s:  %4d of %d - \t R: %4f \t S: %4f \t %s: %.3f\n",  
        j,
        nrow(grid),
        comp,
        ifelse(episode_i <= n_subjects, "training", "validation"),
        episode_i,
        N_EPISODE,
        r,
        beta,
        ifelse(expl == "eps_greedy", "Epsilon", "Tau"),
        ifelse(expl == "eps_greedy", agent$epsilon, agent$tau_softmax)
      ))
    
    flush.console()
  }
  
  grid$R[j]<- mean(rewards[(n_subjects+1):(n_new + n_subjects)])
  grid$beta[j]<- mean(betas[(n_subjects+1):(n_new + n_subjects)])
  
  beta_RQN[,j] <- betas[(n_subjects+1):(n_new + n_subjects)] 
  reward_RQN[,j] <- rewards[(n_subjects+1):(n_new + n_subjects)] 
  
  cat(
    sprintf(
      "[%4d of %d - %s ] mean rewards: %4f, mean survival rate: %.3f\n",
      j,
      nrow(grid),
      comp,
      mean(rewards),
      mean(betas)
    )
  )
  
  toc(log = TRUE, quiet = TRUE)
  
}

log.lst <- tic.log(format = FALSE)
log.txt <- tic.log(format = TRUE)

timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))

# mean(timings)/60 
# sum(timings)/60 
# writeLines(unlist(log.txt))

grid$time <- as.numeric(timings)/60
grid$tot_time <- sum(timings)/60

# write.csv2(grid,".//grid_RQN.csv", row.names = F)

dir <- paste0(current_directory, "/Performance/", method) 
write.xlsx(grid, paste0(dir,'/grid_RQN.xlsx')) 
write.xlsx(beta_RQN, paste0(dir,'/beta.xlsx')) 
write.xlsx(reward_RQN, paste0(dir,'/reward.xlsx')) 
write.xlsx(CR_RQN, paste0(dir,'/CR.xlsx')) #?

save.image(".//RQN.RData")


############################################# LSTM #############################################

rm(list=ls()[!ls() %in% c("grid", "current_directory", "n_new")])

method = "LSTMQN"
grid[, "method"]  <- method

# settings 
MEMORY_WINDOW = 4L  

BATCH_SIZE = 32  

beta_LSTM <- data.frame(matrix(ncol = nrow(grid), nrow = n_new)) 
reward_LSTM <- data.frame(matrix(ncol = nrow(grid), nrow = n_new)) 
CR_LSTM <- data.frame(matrix(0, ncol = max(grid$patients), nrow = nrow(grid))) #?
colnames(beta_LSTM) <- colnames(reward_LSTM) <- c(paste0("comb_", 1:nrow(grid))) 


source("LSTMAgent.R")
source("Memory.R")
source("train.R")
source("Environment.R")

#
tic.clearlog()

for (j in 1:nrow(grid)) {
  tic(j)
  Sys.sleep(1)
  
  # eps_greed = F <- Boltzmann strategy, epsilon greedy strategy otherwise
  expl = grid[j, "expl_strategy"] 
  if(!(expl %in% c("eps_greedy", "softmax"))){ 
    stop("Exploration strategy can only be 'eps_greedy' or 'softmax'. ", call. = F)      
  }
  
  # target update strategy
  update = grid[j, "updt_strategy"]
  if(!(update %in% c("hard", "soft"))){ 
    stop("Target update strategy can only be 'hard' or 'soft'. ", call. = F)      
  }
  
  n_subjects = grid[j, "patients"]         
  n_stages = grid[j, "stages"]            
  
  p = grid[j, "p"]              
  q = grid[j, "q"]                 
  v = grid[j, "v"]    
  
  zk = grid[j, "RF"]
  int = grid[j, "int"]                
  noise_pred = grid[j, "noise_pred"] 
  
  comp <- grid[j, "complexity"]
  
  treatment_options = grid[j, "treatment_options"]
  
  tr <- as.data.frame(sapply(read_csv(paste0(".//tr_dataset/tr", j, ".csv"), show_col_types = FALSE), as.numeric))
  ev <- as.data.frame(sapply(read_csv(paste0(".//eval_dataset/ev", j, ".csv"), show_col_types = FALSE), as.numeric))
  
  # 
  
  N_EPISODE = n_subjects + n_new  
  tr_phase = T  # parameter set if an epsilon greedy or softmax exploration strategy is wanted during training phase
  
  env <- Environment$new(seed = 123,
                         int = int,  
                         noise_pred = noise_pred,
                         zk = zk,
                         treatment_options = treatment_options,  
                         p = p,
                         q = q,
                         v = v,
                         n_stages = n_stages)
  
  obs_size <- length(c("tumor_mass", "toxicity", "stress",
                       grep("\\X|\\Z|\\V", names(tr), value = T)))
  
  agent <-
    Agent$new(
      input_shape = c(MEMORY_WINDOW, obs_size),
      output_dim = env$treatment_options,
      decay_last_episode = 100,
      C = 20,
      tau_softmax = .5,
      tau_soft_updt = 1e-3, 
      updt = update
    )
  
  memory <- Memory$new(capacity = 20000)
  
  rewards <- c()
  betas <- c()
  
  for (episode_i in 1:N_EPISODE) {
    CR <- 0 #?
    if(episode_i <= n_subjects){
      
      s <- env$reset(tr[episode_i, ])
      
    } else{
      
      if (episode_i == n_subjects + 1){
        val_i <- 1
        tr_phase <- F
      } else{
        val_i <- val_i + 1
      }
      
      s <- env$reset(ev[val_i, ])
    }
    
    done <- FALSE
    
    history = rep(0., MEMORY_WINDOW * obs_size)
    history = append(history, s)[-c(1:obs_size)]
    
    doses_values <- rep(0,3)    
    
    while (!done) {
      
      a <- agent$get_action(state_ = array(history, c(1, MEMORY_WINDOW, obs_size)), 
                            step = episode_i, tr_phase , expl_str = expl)
      treatment <- as.numeric(a/(env$treatment_options - 1))  
      doses_values <- append(doses_values, treatment)[-1] 
      ret <- env$step(treatment, doses_values)  
      s2 <- ret$state
      history2 = append(history, s2)[-c(1:obs_size)]
      done <- ret$dropout
      
      r <- 0 
      
      if(done){
        
        r <- ret$reward  
        
        beta <- ret$beta
        
        betas <- c(betas, beta)
        rewards <- c(rewards, r)
      }
      
      if(episode_i <= n_subjects){
        CR <- CR + ret$reward #?
        memory$push(history, a, r, done, history2)
        
        if (memory$length > BATCH_SIZE) {
          batch <- memory$sample(BATCH_SIZE)
          train(agent, batch, env$i)
        }
        
      }
      
      
      s <- s2
      history <- history2
      
      if(episode_i <= n_subjects){ #?
        CR_LSTM[j,episode_i] <- CR #?
      } #?
      
    }
    
    cat(
      sprintf(
        "[%d of %d - %s] %s:  %4d of %d - \t R: %4f \t S: %4f \t %s: %.3f\n",  
        j,
        nrow(grid),
        comp,
        ifelse(episode_i <= n_subjects, "training", "validation"),
        episode_i,
        N_EPISODE,
        r,
        beta,
        ifelse(expl == "eps_greedy", "Epsilon", "Tau"),
        ifelse(expl == "eps_greedy", agent$epsilon, agent$tau_softmax)
      ))
    
    flush.console()
  }
  
  grid$R[j]<- mean(rewards[(n_subjects+1):(n_new + n_subjects)])
  grid$beta[j]<- mean(betas[(n_subjects+1):(n_new + n_subjects)])
  
  beta_LSTM[,j] <- betas[(n_subjects+1):(n_new + n_subjects)] 
  reward_LSTM[,j] <- rewards[(n_subjects+1):(n_new + n_subjects)] 
  
  cat(
    sprintf(
      "[%4d of %d - %s ] mean rewards: %4f, mean survival rate: %.3f\n",
      j,
      nrow(grid),
      comp,
      mean(rewards),
      mean(betas)
    )
  )
  
  toc(log = TRUE, quiet = TRUE)
  
}

log.lst <- tic.log(format = FALSE)
log.txt <- tic.log(format = TRUE)

timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))

# mean(timings)/60 
# sum(timings)/60 
# writeLines(unlist(log.txt))

grid$time <- as.numeric(timings)/60
grid$tot_time <- sum(timings)/60

# write.csv2(grid,".//grid_LSTMQN.csv", row.names = F)
dir <- paste0(current_directory, "/Performance/", method) 
write.xlsx(grid, paste0(dir,'/grid_LSTMQN.xlsx')) 
write.xlsx(beta_LSTM, paste0(dir,'/beta.xlsx')) 
write.xlsx(reward_LSTM, paste0(dir,'/reward.xlsx')) 
write.xlsx(CR_LSTM, paste0(dir,'/CR.xlsx')) #?

save.image(".//LSTMQN.RData")

#### PLOTS #####

source("Plots.R")

