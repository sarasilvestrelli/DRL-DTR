train_ql_ols <- function(data,       
                         n_stages,
                         p,             
                         q,                 
                         v,   
                         treatment_options,  
                         possible_actions,
                         seed,
                         int,              
                         noise_pred,
                         zk){
      
  set.seed(seed)
  
  mod_list <- list()
  predictors_name <- list()
  
  d <- cbind(ID = rownames(data), data)
  
  for (i in (n_stages - 1):0) {
    
    keep <- c(1, which(str_extract(names(d), "([^_]+$)") %in% (paste0(i:(i-MEMORY_WINDOW+1)))))
    
    if(i == (n_stages - 1)){
      # remove dropout patients at the last stage
      x <- d[!(d$dropout), keep]
    }
    
    if(i < (n_stages - 1)){
      # remove dropout patients at step i
      x <- d[!is.na(d[, paste0("tumor_mass_", i + 1)]), keep] 
    }
    
    x <- na.omit(x) 
    n <- nrow(x)
    
    drop <- paste0("reward_", 0:i)
    
    if(i > 0){
      drop <- c(drop, paste0("Qhat_", 0:(i-1)), paste0("dose_", seq(0, i-1))) 
    }
    
    x = x[, !(names(x) %in% drop)]
    
    rem <- c(paste0("Qhat_", i), "ID")
    predictors <- names(x)[!(names(x) %in% rem)]
    predictors_name[[paste0('month', i)]] <- predictors


    form <- as.formula(paste(paste0("Qhat_",i, " ~ ( "), paste(paste(predictors,
                                                                     collapse = " + ")),") *", paste0("dose_",i)))     
    mod <- lm(form, data = x) 
    
    mod_list[[paste0('month', i)]] <- mod
    
    x_expand <- x[rep(seq_len(n), each = length(possible_actions)), ]
    x_expand[, paste0('dose_', i)] <- rep(possible_actions, n)
    x_expand <- data.table::as.data.table(x_expand)
    x_expand[, pred := predict(mod, x_expand)] 
    max_preds <- x_expand[, max(pred), by = ID]
    
    if (i > 0) {
      ind <- (d$ID %in% max_preds$ID)
      q.id <- which(names(d) == paste0("Qhat_", (i-1)))
      r.id <- which(names(d) == paste0('reward_', (i-1)))
      d[ind, q.id] <- d[ind, r.id] + max_preds$V1
    }
  }
  
  return(list(mod_list = mod_list, predictors_name = predictors_name))
  
}
