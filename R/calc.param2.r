#' @example #pars=calc.param2(migr.spd=1, resid.frac = 0.1, g=g);siz=pars$siz1;sigma=pars$sig1;muadv=0
calc.param2=function (migr.spd, resid.frac = 0.1, g){
  avg.lon=mean(range(g$lon))
  avg.lat=mean(range(g$lat))
  # cell width in km (at the middle of the study area)
  # dx=g$dlo*sp::spDistsN1(matrix(data=c(avg.lon,avg.lat),ncol=2),c(avg.lon+1,avg.lat),longlat=TRUE)
  # cell height in km
  dy=g$dla*sp::spDistsN1(matrix(data=c(avg.lon,avg.lat),ncol=2),c(avg.lon,avg.lat+1),longlat=TRUE)
  
  #migr.spd.x <- migr.spd/1000 * 3600 * 24 / dx # m/s-1 to cell/day-1
  migr.spd.y <- migr.spd/1000 * 3600 * 24 / dy 
  sigma=migr.spd.y/ sqrt(.5*pi) # cf woillez
  sigmas=c(sigma,sigma*resid.frac)
  dims=sigmas*4 # +-SD (68%) +-2SD (97.5%) +-4 should encompass more than 99%
  return(list(siz1 = dims[1], sig1 = sigmas[1], siz2 = dims[1], sig2 = sigmas[2]))
}