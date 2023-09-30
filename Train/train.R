train <- function(agent, memory, episode, discount = .999) {
  
  n_samples <- nrow(memory$states)                   
  
  X_memory <- memory$states
  y_memory <- agent$predict(X_memory)                   
  y_future <- agent$predict_target(memory$states2)     
  
  for (i in 1:n_samples) {
    action_idx <- memory$actions[i] + 1L
    
    if (memory$dones[i] == T) {
      y_memory[i, action_idx] <- memory$rewards[i]
    } else {
      y_memory[i, action_idx] <- memory$rewards[i] + discount*max(y_future[i, ])
    }
  }
  
  agent$train(X_memory, y_memory)
  
  # target network hard update every C steps 
  if(episode %% agent$C == 0  & agent$updt == "hard"){
    agent$update_target()
  }
  
  # target network soft update
  if(agent$updt == "soft"){
    agent$soft_target_update()
  }
  
}