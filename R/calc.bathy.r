# focalDim=3;ncores=4;sens.err=1;bathy.grid=b3
calc.bathy <- function(mmd, bathy.grid, dateVec, focalDim = NULL, sens.err = 1){
  
  # required functions & packages
  # functions needed in parallel
  #ac=function(x){return(as.character(x))}
  #an=function(x){return(as.numeric(ac(x)))}
  #an.=function(x){return(as.numeric(x))}
  
  
  print(paste("Starting bathymetry likelihood calculation..."))
  t0 <- Sys.time()
  mmd$Date <- as.Date(as.POSIXct(mmd$Date, format = findDateFormat(mmd$Date)))
  
  # compute bathy.grid sd
  # bathy.grid=b3;focalDim=3
  sdx = raster::focal(bathy.grid, w = matrix(1, nrow = focalDim, ncol = focalDim), fun = function(x) stats::sd(x,na.rm = T))
  # sdx = t(raster::as.matrix(raster::flip(sdx, 2)))
  #  dat <- base::t(raster::as.matrix(raster::flip(bathy.grid, 2)))
  
  T <- length(dateVec)
  print(paste("Starting iterations through deployment period ", "..."))
  
  L.bathy <- array(0, dim = c(dim(bathy.grid)[1:2], length(dateVec)))
  lon.agg <- seq(raster::extent(bathy.grid)[1], raster::extent(bathy.grid)[2], 
                 length.out = dim(bathy.grid)[2])
  lat.agg <- seq(raster::extent(bathy.grid)[3], raster::extent(bathy.grid)[4], 
                 length.out = dim(bathy.grid)[1])
  
  for (i in 1:T) {
    print(i)
    # i=1; sens.err=1
    idx <- which(mmd$Date %in% dateVec[i])
    if (length(idx) == 0) next
    
    bathy.i <- c(mmd$MaxDepth[idx] * (1 - sens.err / 100), mmd$MaxDepth[idx] * (1 + sens.err / 100))
    #bathy.i = sort(-abs(bathy.i))
    lik.bathy <- likint3(bathy.grid, sdx, bathy.i[1], bathy.i[2])
    
    #idx <- which(dateVec == as.Date(time))
    lik.bathy = as.matrix(lik.bathy) / max(as.matrix(lik.bathy), na.rm = T)
    lik.bathy[is.na(lik.bathy) | is.infinite(lik.bathy)] <- 0
    L.bathy[,,which(dateVec %in% mmd$Date[idx])]<-lik.bathy
  }
  
  L.bathy <- aperm(L.bathy,c(2,1,3))
  L.bathy <- L.bathy[,dim(L.bathy)[2]:1,]
  
  print(paste("Making final likelihood raster..."))
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.bathy <- list(x = lon.agg, y = lat.agg, z = L.bathy)
  ex <- raster::extent(list.bathy)
  L.bathy <- raster::brick(list.bathy$z, xmn = ex[1], xmx = ex[2], 
                           ymn = ex[3], ymx = ex[4], transpose = T, crs)
  L.bathy <- raster::flip(L.bathy, direction = "y")
  L.bathy[L.bathy < 0] <- 0
  t1 <- Sys.time()
  print(paste("Bathymetric calculations took ", round(as.numeric(difftime(t1, t0, units = "mins")), 2), "minutes..."))
  return(L.bathy)
}