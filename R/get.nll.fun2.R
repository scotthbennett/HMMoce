#' Negative Log likelihood of parameters
#' 
#' @param parvec vector of length 6 contianing * kernel 1 * kernel 2 * diagonal
#'   of 2x2 matrix
#' @param g grid from \code{\link{setup.grid}}
#' @param L final likelihood (2D)
#'   
#' @return parameter values
#' @export
#' @references Pedersen MW, Patterson TA, Thygesen UH, Madsen H (2011)
#'   Estimating animal behavior and residency from movement data. Oikos
#'   120:1281-1290. doi: 10.1111/j.1600-0706.2011.19044.x
#'   
#' @examples
#' \dontrun{
#' par0 <- c(log(10), log(10), log(0.5), log(0.5), log(0.95/0.05), log(0.95/0.05))
#' fit <- nlm(get.nll.fun, par0, g, L, dt)
#' D1 <- exp(fit$estimate[1:2])
#' D2 <- exp(fit$estimate[3:4])
#' p <- 1/(1+exp(-fit$estimate[5:6]))
#' fit <- nlm(neg.log.lik.fun, guess, g.mle, L.mle)
#' }

#nll <- nlm(get.nll.fun2, parvec=c(10), g.mle, L.mle, g, L)

get.nll.fun2 <- function(parvec = c(10), g.mle, L.mle, g, L){
  
  # parvec is migr.spd, maskL, bnd 
  migr.spd <- 2 #parvec[1]
  #maskL <- parvec[2]
  bnd <- parvec[1]
  
  # GET MOVEMENT KERNELS AND SWITCH PROB FOR COARSE GRID
  par0 <- makePar(migr.spd=migr.spd, grid=g.mle, L.arr=L.mle, p.guess=c(.9,.9), calcP=T)
  #K1 <- par0$K1; K2 <- par0$K2; 
  P.final <- par0$P.final
  
  # GET MOVEMENT KERNELS AND SWITCH PROB FOR FINER GRID
  par0 <- makePar(migr.spd=migr.spd, grid=g, L.arr=L, calcP=F)
  K1 <- par0$K1; K2 <- par0$K2; #P.final <- par0$P.final
  
  f <- hmm.filter(g, L, K1, K2, maskL = T, P.final, minBounds = bnd)
  
  nllf <- -sum(log(f$psi[f$psi>0]))
  print(paste0("HMM -log(L): ", nllf))
  #flush.console()
  nllf
  
}

