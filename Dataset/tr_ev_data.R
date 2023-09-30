rm(list=ls())

source("sim.R")
load(".\\grid_parameters.RData")


# Duplicated 
fac <- c("complexity", "patients","n_cov", "stages", "treatment_options", "p", "q","v","RF","int","noise_pred")

match = as.data.table(grid)[, c("GRP", "N") := .(.GRP, .N), by = names(grid[,which(colnames(grid)%in%fac)])][
  N > 1, list(list(.I)), by = GRP]
dup = matrix(as.vector(unlist(match$V1)), 
             ncol = length(unlist(match$V1[1])), 
             byrow = T) 

##### TRAIN RESET ##### 

dir.create(".\\tr_dataset" )
dir <- ".\\tr_dataset"

for (j in 1:nrow(grid)) {
  p = grid[j, "p"]              
  q = grid[j, "q"]                 
  v = grid[j, "v"]    
  
  zk = grid[j, "RF"]
  int = grid[j, "int"]                
  noise_pred = grid[j, "noise_pred"] 
  
  obs_size <- 3 + 2*int + p+q+v + zk
  
  comp <- grid[j, "complexity"]
  
  treatment_options = grid[j, "treatment_options"]
  n_stages = grid[j, "stages"]
  n_subjects = grid[j, "patients"]

  dat <-  main(seed = 123, 
               n_subjects = n_subjects, 
               n_stages = 1,
               int = int,  
               noise_pred = noise_pred,
               zk = zk,
               treatment_options = treatment_options,  
               p = p,
               q = q,
               v = v)$out
  
  dat[ , which(names(dat) == "month"):ncol(dat)] <- NULL

  write.csv(dat, paste0(dir, "\\tr", j, ".csv"), row.names=F)
  
  cat(
    sprintf(
      "[%4d / %d - %s ]\n",
      j,
      nrow(grid),
      comp
    )
  )
  
}

##### EVALUATION RESET ##### 

dir.create(".\\eval_dataset" )
dir1 <- ".\\eval_dataset"

# QL
g <- grid[-which(grid$expl_strategy == "softmax" | grid$updt_strategy =="hard"),
          ! names(grid) %in% c("expl_strategy", "updt_strategy")]
rownames(g) <- seq(1:nrow(g))
save.image(".//g.RData")

dir.create(".\\eval_dataset_QL" )
dir2 <- ".\\eval_dataset_QL"

#

n_new = 500
list <- 1:nrow(grid) 
count = 0 
for (j in 1:nrow(grid)) {
  
  if(j %in% list){ 
    
    p = grid[j, "p"]              
    q = grid[j, "q"]                 
    v = grid[j, "v"]    
    
    zk = grid[j, "RF"]
    int = grid[j, "int"]                
    noise_pred = grid[j, "noise_pred"] 
    
    obs_size <- 3 + 2*int + p+q+v + zk
    
    comp <- grid[j, "complexity"]
    
    treatment_options = grid[j, "treatment_options"]
    n_stages = grid[j, "stages"]
    
    
    dat <-  main(seed = 1, 
                 n_subjects = n_new, 
                 n_stages = 1,
                 int = int,  
                 noise_pred = noise_pred,
                 zk = zk,
                 treatment_options = treatment_options,  
                 p = p,
                 q = q,
                 v = v)$out
    
    dat[, which(names(dat) == "month"):ncol(dat)] <- NULL
    
    r.c = which(dup==j, arr.ind=TRUE) 
    for (k in dup[r.c[1],]) { 
      write.csv(dat, paste0(dir1, "\\ev", k, ".csv"), row.names=F)
      list <- list[list != k]
    }
    
    
    count = count + 1 
    write.csv(dat, paste0(dir2, "\\ev", count, ".csv"), row.names=F) 

  }
  
  

  
  cat(
    sprintf(
      "[%4d / %d - %s ]\n",
      j,
      nrow(grid),
      comp
    )
  )
  
}

rm(list=ls()[!ls() %in% c("grid", "g")])

