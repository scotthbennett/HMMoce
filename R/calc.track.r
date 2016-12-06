#' Calculate movement track from state estimates
#' 
#' \code{calc.track} uses HMM output via \code{hmm.smoother} to calculate most
#' probable track and behavior state
#' 
#' @param distr is output array from \code{hmm.smoother}
#' @param g is one of the outputs from \code{resample.grid} which denotes what 
#'   spatial scale and grid you're working on
#' @param dateVec is vector of dates from tag to pop-up in 1 day increments.
#' @param method is character indicating what method to use for track 
#'   calculation. Currently only 'mean' is supported.
#'   
#' @return calculated track
#' @export
#' @examples 
#' \dontrun{
#' 
#' # RUN THE SMOOTHING STEP
#' s = hmm.smoother(f, K1, K2, P.final)
#' 
#' # GET AND PLOT THE MOST PROBABLE TRACK
#' tr <- calc.track(s, g, dateVec)
#' 
#' }
#' 

calc.track <- function(distr, g, dateVec, method = 'mean'){
  ## Calculate track from probability distribution of location
  
  if (method == 'mean'){
    T <- dim(distr)[2]
    # Track calculated from mean
    lat <- apply(apply(distr, c(2, 4), sum) * repmat(t(as.matrix(g$lat[,1])), T, 1), 1, sum)
    lon <- apply(apply(distr, c(2, 3), sum) * repmat(t(as.matrix(g$lon[1,])), T, 1), 1, sum)
    
  } else if (method == 'mode'){
    
    stop('Mode is currently unsupported.')
    
    # Track calculated from mode
    row <- dim(g$lon)[1]
    col <- dim(g$lon)[2]
    modelat <- rep(0, T)
    modelon <- rep(0, T)
    
    for(t in 1:T){
      asd <- apply(distr[,t,,], c(2,3), sum)
      ind <- which.max(asd)
      x <- ceiling(ind / col)
      y <- ind %% row
      modelat[t] <- g$lat[y,x]
      modelon[t] <- g$lon[y,x]
    }
  } else if (method == 'Viterbi'){
    
    stop('Viterbi is currently unsupported.')
    # Track calculated with Viterbi
    # --- not included in this script
    
  }
  
  # calculate the estimated behavior state
  p.resid <- apply(distr, c(1,2), sum)[2,]
  
  track <- data.frame(cbind(date = dateVec, lon = lon, lat = lat, p = p.resid))
  track$date <- dateVec

  return(track)

}