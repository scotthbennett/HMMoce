#' HMM filter functions
#'
#' 
#' @param g grid from \code{\link{setup.grid}}
#' @param L final likelihood (2D)
#' @param K1 first movement (diffusion) kernel see \code{\link{gausskern}}
#' @param K2 second movement (diffusion) kernel see \code{\link{gausskern}}
#' @param P 2x2 probability matrix for tranisitons between states (K1 and K2)
#'
#' @return a list: list(phi = phi, pred = pred, psi = psi) where
#' \itemize{
#'  \item phi. is the probability for each state at each trime step 
#'  \item pred. is ....
#'  \item psi. is.... 
#' }
#' @export
#'

hmm.filter <- function(g, L, K1, K2, P){
  
  ## Filter data to estimate locations and behaviour
  
  T <- dim(L)[1] # dimension of time 
  row <- dim(g$lon)[1] # nrows
  col <- dim(g$lon)[2] # ncols
  m <- 2 # Number of behavioural states
  
  pred <- array(1e-15, dim = c(m, T, col, row)) # empty array for prediction step. ordering in col before row emulates lon before lat
  phi  <- array(1e-15, dim = c(m, T, col, row)) # posterior (final) step array
  
  # Start in resident state at the known initial location
  #phi[1,1,,]  <- L[1,,] # first position is known
  phi[2,1,,]  <- L[1,,] # first position is known
  #pred[1,1,,] <- L[1,,] # first position is known
  pred[2,1,,] <- L[1,,] # first position is known
  psi <- rep(0, T - 1) # sum of the probability of both states at each step
  
  # convert movement kernels from matrix to cimg for convolution
  K1 <- imager::as.cimg(K1 / max(K1))
  K2 <- imager::as.cimg(K2 / max(K2))
  
  # Start filter iterations
  for(t in 2:T){
    
    # convolve previous day's likelihood with movement kernels
    p1 = imager::as.cimg(t(phi[1, t-1,,]))
    p2 = imager::as.cimg(t(phi[2, t-1,,]))
    
    #**problem here in which by setting day1, behav1 to all 1e-15, the convolution with migr kernel doesnt allow fish to move away from tagging location
    # solved with if statement
    if(t == 2){
      q1 = imager::convolve(p2, K1)
    } else{
      q1 = imager::convolve(p1, K1)
    }
    
    q2 = imager::convolve(p2, K2)
    q1 = t(as.matrix(q1))
    q2 = t(as.matrix(q2))
    
    # multiply by transition probability 
    pred[1,t,,] <- P[1,1] * q1 + P[2,1] * q2
    pred[2,t,,] <- P[1,2] * q1 + P[2,2] * q2
    
    max1 <- max(pred[1,t,,])
    max2 <- max(pred[2,t,,])
    
    pr1 <- pred[1,t,,] / max1
    pr2 <- pred[2,t,,] / max2
    pr1[pr1 <= .05] <- 0
    pr2[pr2 <= .05] <- 0
    
    # is there a data-based likelihood observation for this day, t?
    sumL = sum(L[t,,])  
    if(sumL > 1e-6){
      
      #post1 <- pred[1,t,,] * L[t,,]
      #post2 <- pred[2,t,,] * L[t,,]
      post1 <- pr1 * L[t,,] + pr1
      post1 <- post1 / max(post1) * max1
      post2 <- pr2 * L[t,,] + pr2
      post2 <- post2 / max(post2) * max2
      
    }else{
      post1 <- pred[1,t,,]
      post2 <- pred[2,t,,]
    }
    
    psi[t-1] <- sum(as.vector(post1), na.rm=T) + sum(as.vector(post2), na.rm=T)
    
    phi[1,t,,] <- post1 / (psi[t-1] + 1e-15)
    phi[2,t,,] <- post2 / (psi[t-1] + 1e-15)
    #phi[phi <= 1e-15] <- 1e-15
    
  }
  
  # End in resident state at the known final location
  #phi[1,T,,]  <- L[T,,] # last position is known
  #phi[2,T,,]  <- L[T,,] # last position is known
  #pred[1,T,,] <- L[T,,] # last position is known
  #pred[2,T,,] <- L[T,,] # last position is known
  
  list(phi = phi, pred = pred, psi = psi)
  
}

