
Agent <- R6Class(
    "Agent",
    public = list(
      input_shape = NULL,
      output_dim = NULL,
      epsilon = 0.8,
      model = NULL,
      target_model = NULL, 
      C = 20, 
      decay_last_episode = 100, # epsilon --> 0 always exploit ; epsilon --> 1 always take random actions 
      tau_softmax = 0.5, # tau --> 0 pure exploitation ; tau --> infinity pure exploration   
      tau_soft_updt = 1e-3,
      tau_decay = F,
      updt = NULL,
      
      initialize = function(input_shape,
                            output_dim,
                            decay_last_episode = NULL,
                            C = NULL,
                            tau_softmax = NULL,
                            tau_soft_updt = NULL,
                            updt = NULL) {
          
        self$C <- C 
        self$tau_softmax <- tau_softmax
        self$tau_soft_updt <- tau_soft_updt
        self$updt <- updt
        
        self$input_shape <- input_shape
        self$output_dim <- output_dim
        
        self$model <- self$init_model()
        self$target_model <- self$init_model()
        
        self$update_target()
        self$soft_target_update()
        
        print(self$model)

        },
      
      init_model = function(){
        model <- keras_model_sequential()
        model %>%
          layer_lstm(units = 32,
                     return_sequences=T,
                     input_shape = self$input_shape
                     ) %>%
          layer_lstm(units = 32
                     ) %>%
          layer_dense(units = self$output_dim)
        
        model %>% compile(
          loss = 'mse', 
          optimizer = optimizer_adam(1e-03),
          metrics = c("mse") 
          )
        return(model)
        },
        
      update_target = function(){
        set_weights(self$target_model, get_weights(self$model))
      },
      
      soft_target_update = function(){
        w <- get_weights(self$model)
        w_target <- get_weights(self$target_model)
        
        for (i in 1:length(w)){
          w_target[[i]] = self$tau_soft_updt*w[[i]] + (1-self$tau_soft_updt)*w_target[[i]]
        }
        
        set_weights(self$target_model, w_target)
      },
      
      get_action = function(state_, step, tr_phase, expl_str) {
        
        if (expl_str == "eps_greedy"){                              # get action from model using deterministic epsilon-greedy policy
          
          if (runif(1) < self$epsilon & tr_phase) {
            action <- sample.int(self$output_dim, size = 1) - 1L
          } else {
            action_probs <- self$model$predict(state_)
            action <- which.max(action_probs) - 1L
          }
          
          self$epsilon <-
            max(0.01, -1 / self$decay_last_episode * step + 1.0)
          
        } 
        
        if(expl_str == "softmax"){                                  # get action from model using boltzmann exploration (softmax exploration)
          pred <- self$model$predict(state_)                        # P(a)=exp(Q(s,a)/tau)/sum_t(exp(Q(s,a_t)/tau))
          
          if (tr_phase) { 
            # num <- exp(Q(s,a)/tau)
            num <- exp(pred / self$tau_softmax)
            # den <- sum_i(exp(Q(s,a_i)/tau))
            den <- sum(exp(pred / self$tau_softmax))
            # num/den
            probs <- round(num / den, 10)
            
            # CHECK 
            # cat(sprintf(
            #   "check: %s - Lunghezza : %1.0%f\n",
            #   toString(as.character(probs)),
            #   length(probs)
            # ))
            
            action <- sample.int(n = self$output_dim, size = 1, prob = probs) - 1L
            
            if (self$tau_decay){
              self$tau_softmax <-
                max(0.1, -1 / self$decay_last_episode * step + 1.0)
            }
            
          } else {
            action <- which.max(pred) - 1L
          }
        }
        action
      },
      
      predict = function(states_) {
        b = dim(states_)[1]
        states_ = aperm(array(states_, dim=c(self$input_shape[2], self$input_shape[1], b)), c(3,2,1))
        self$model$predict(states_)
      },
      
      predict_target = function(states_) {
        b = dim(states_)[1]
        states_ = aperm(array(states_, dim=c(self$input_shape[2], self$input_shape[1], b)), c(3,2,1))
        self$target_model$predict(states_)
      },
      
      train = function(states_, targets_) {
        b = dim(states_)[1]   # 3D array (batch_size, time_steps, units)
        states_ = aperm(array(states_, dim=c(self$input_shape[2], self$input_shape[1], b)), c(3,2,1))
        
        self$model %>% fit(
          states_, targets_,
          batch_size = b,
          epochs=1,
          verbose = 0,
          shuffle = F
        )
      }
    )
)
