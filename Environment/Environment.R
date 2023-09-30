
source(".//Environment_functions.R")

Environment <- R6Class(
    "Environment",
    public = list(
        int = NULL,
        noise_pred = NULL,
        zk = NULL,
        treatment_options = NULL,
        p = NULL,
        q = NULL,
        v = NULL,
        w = NULL,
        state = NULL,
        n_stages = NULL,
        i = NULL,
        cW = NULL,
        cM = NULL,
        noise_chng = NULL,
        


        initialize = function(seed=1, 
                              int=T,  
                              noise_pred=T,  
                              zk = T,
                              treatment_options = 5,  #0,0.25,0.5,0.75,1
                              p=1,
                              q=1,
                              v=1,
                              w=0.05,
                              n_stages=5) {
            set.seed(seed)
            self$int <- int
            self$noise_pred <- noise_pred
            self$zk <- zk
            self$treatment_options=treatment_options
            self$p = p
            self$q = q
            self$v = v
            self$w = w
            self$i = 0             
            self$n_stages = n_stages
            self$cW = 1
            self$cM = 1
            self$noise_chng = 0
            
            if((self$p==0 & self$q>0) | (self$q==0 & self$p>0)){ 
                stop("p and q must be both 0 or an integer greater than zero.", call. = F)      
            }
            
            if((self$v<1) & (self$noise_pred==T)){ 
                stop("p, q and v, or at least just the parameter v, must be an integer greater than zero if noise options is TRUE.", call. = F)      
            }
        },
        
        reset = function(vector) {
            
            self$i <- 0 
            
            self$cW <- vector$cW
            self$cM <- vector$cM 
            
            self$noise_chng <- vector$noise_chng
            
            dat <- vector %>% 
              select(map(c("tumor", "tox", "str", 'X', 'Z', 'V'), 
                         starts_with, 
                         vars = colnames(.)) %>%
                       unlist())
            
            self$state <- dat
            
            s = as.numeric(dat)

            return(s)
        }, 
        
        step = function(treatment, doses_values) {   
          
            
            dat <- self$state
            

            old_M <- dat$tumor_mass
            old_W <- dat$toxicity    
            
            dat$toxicity <- updateW(old_M, old_W, treatment, c = self$cW)
            
            dat$tumor_mass <- updateM(old_M, old_W, treatment, c = self$cM)
            
            dat$stress <- updateStress(old_M, old_W, dat$stress, doses_values, dat$toxicity, dat$tumor_mass) 
            
            if(self$zk == T){
              dat$Zk <- updateZk(old_M, old_W, treatment, dat$Zk, dat$tumor_mass)
              self$noise_chng <- self$w * sum(dat[, grep("\\Z", names(dat), value = T)]) 
            }  
            
            
            beta <- 1 / lambda(dat$tumor_mass, dat$toxicity, Z = self$noise_chng)  
            dropout  <- beta < 1 
           
               
            reward <- log(self$i + beta) + 1/(self$i + dat$stress + 1)   
            
            self$i = self$i + 1  
            
            self$state = dat
            
            if(self$i > self$n_stages){dropout = T}
            
            return(list(
                state2 = as.numeric(dat), 
                reward = reward, 
                dropout = dropout, 
                beta = beta))

        }
    ) 
)



