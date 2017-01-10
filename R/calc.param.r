#' Calculates movement parameters used for behavior kernels
#' 
#' \code{calc.param} calculates movement parameters used to generate kernels
#' 
#' Movement parameters for each of two behavioral modes are generated using an input grid
#' and user-input speed (m s) for at least the migratory state.
#' 
#' @param migr.spd is numeric indicating movement speed of animal (in m s) while in migratory behavior state.
#' @param resid.spd is numeric indicating movement speed of animal (in m s) while in resident behavior state. If not specified, default is 10 percent of migratory speed. 
#' @param g is grid from \code{setup.grid}
#' 
#' @export
#' @return a list of movement parameters

calc.param <- function(migr.spd, resid.spd = NULL, g){
  
  migr.spd <- migr.spd / 1000 * 3600 * 24 / 111
  migr <- migr.spd / g$dla
  
  if (is.null(resid.spd)){
    resid.spd <- migr.spd / 10
  }
 
  resid.spd <- resid.spd / 1000 * 3600 * 24 / 111
  resid <- resid.spd / g$dla
  
  return(list(migr = migr, sig1 = migr, resid = resid, sig2 = resid / 4))
    
}