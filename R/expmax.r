#' Expectation-maximization framework for state-switching
#' 
#' @param P.init 2x2 probability matrix for transitions between states (K1 and 
#'   K2). P.init[1,1] is probability of staying in behavior state 1 if currently
#'   in state 1. P.init[1,2] is probability of switching to state 2 if currently
#'   in state 1. Continue this logic for P.init[2,].
#' @param g grid from \code{\link{setup.grid}}
#' @param L final likelihood (2D)
#' @param K1 first movement (diffusion) kernel see \code{\link{gausskern}}
#' @param K2 second movement (diffusion) kernel see \code{\link{gausskern}}
#' @param niter is integer that determines number of iterations to perform
#' @param threshold is threshold of percent change that we consider satisfactory
#'   for convergence. Default is 1%.
#'   
#' @return a 2x2 matrix of state switching probabilities. See P.init input for 
#'   more information.
#' @export
#' @references Woillez M, Fablet R, Ngo TT, et al. (2016) A HMM-based model to
#'   geolocate pelagic fish from high-resolution individual temperature and
#'   depth histories: European sea bass as a case study. Ecol Modell 321:10-22.
#'   

expmax <- function(p.init, g, L, K1, K2, niter = 1000, threshold = .01){
  
  t1 <- Sys.time()
  options(warn = -1)
  
  print(paste('Starting EM for state switching...'))
  
  if (niter < 25){
    stop('Maximum number of iterations (niter) must be > 25.')
  }
  
  if (class(K1) == 'matrix'){
    # convert movement kernels from matrix to cimg for convolution
    K1 <- imager::as.cimg(K1)
    K2 <- imager::as.cimg(K2)
    
  }

  save.p <- data.frame(matrix(NA, ncol=2, nrow=1000))
  save.p[,1] <- p.init[1,1]; save.p[,2] <- p.init[2,2]

  for (i in 1:niter){

    if (i == 1){
      P <- p.init
    }
    
    # RUN THE FILTER STEP
    f <- hmm.filter(g, L, K1, K2, P)
    
    # RUN THE SMOOTHING STEP
    s <- hmm.smoother(f, K1, K2, P)
    
    #------------------------#
    # UPDATE P
    # get states from s
    sv <- apply(s[1,,,], 1, sum) > apply(s[2,,,], 1, sum)
    sv <- -sv + 2
    
    # difference it to find transitions
    svd <- rbind(sv, c(0, diff(sv)))
    
    # get a new transition probability matrix
    r1 <- rev(table(svd[2, svd[1,] == 1]) / sum(svd[1,] == 1))
    r2 <- rev(table(svd[2, svd[1,] == 2]) / sum(svd[1,] == 2))
    P <- matrix(rbind(r1, r2), 2, 2)
    
    # CALCULATE CHANGE IN PARAMETER VALUES OVER LAST SEVERAL ITERATIONS
    # We use "ratio betwen the average over the last 20 values of D and the new value
    # of D below 1%" (Woillez 2016)
    
    save.p[i,] <- c(P[1,1], P[2,2])
    
    if (i > 20){
      thr1 <- (save.p[i,1] - mean(save.p[i:(i-20),1])) / mean(save.p[i:(i-20),1])
      thr2 <- (save.p[i,2] - mean(save.p[i:(i-20),2])) / mean(save.p[i:(i-20),2])
    }
    
    # CHECK IF THRESHOLD IS MET. IF NOT, RE-ITERATE
    if(exists('thr1')){
      if (thr1 <= threshold & thr2 <= threshold){
        t2 <- Sys.time()
        print(paste('Convergence took ', round((t2 - t1) / 60),' minutes...',sep=''))
        break
      } else if (i == niter){
        stop(paste('Maximum iterations reached (', niter, ') without crossing percent change threshold...', sep=''))
      }
    }
  
  }
  
  options(warn = 0)
  return(P)
  
}

