makePar <- function(migr.spd, grid, L.arr, calcP=FALSE){
  # PROVIDE FIXED KERNEL PARAMETERS
  par0 <- calc.param2(migr.spd = migr.spd, g = grid)
  D1 <- unlist(par0[1:2]) # parameters for kernel 1. this is migratory behavior mode
  D2 <- unlist(par0[3:4]) # parameters for kernel 2. resident behavior mode
  
  # GENERATE MOVEMENT KERNELS. D VALUES ARE MEAN AND SD PIXELS
  K1 <- gausskern(D1[1], D1[2], muadv = 0)
  K2 <- gausskern(D2[1], D2[2], muadv = 0)
  
  if(calcP){
    # MAKE A GUESS AT STATE SWITCHING PROBABILITY
    p <- c(0.7, 0.8)
    
    # RUN EXPECTATION-MAXIMIZATION ROUTINE FOR MATRIX, P (STATE SWITCH PROBABILITY)
    P.init <- matrix(c(p[1], 1 - p[1], 1 - p[2], p[2]), 2, 2, byrow = TRUE)
    P.final <- expmax(P.init, g = grid, L = L.arr, K1, K2, save = T)
    save.p <- P.final[[2]]; P.final <- P.final[[1]]
  } else{
    P.final = NA
  }
  
  return(list(K1 = K1, K2 = K2, P.final = P.final))
}

