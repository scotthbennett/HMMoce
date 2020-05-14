# focalDim=3;ncores=4;sens.err=1;bathy.grid=b3
calc.bathy.par=function (tag.pdt,bathy.grid, dateVec, focalDim = NULL, sens.err = 1,ncores=4){
  
  # required functions & packages
  require(parallel)
  require(doParallel)
  require(foreach)
  # functions needed in parallel
  ac=function(x){return(as.character(x))}
  an=function(x){return(as.numeric(ac(x)))}
  an.=function(x){return(as.numeric(x))}
  likint3=function (w, wsd, minT, maxT){
    widx <- !is.na(w)
    wdf = data.frame(w = as.vector(w[widx]), wsd = as.vector(wsd[widx]))
    wdf$wsd[is.na(wdf$wsd)] = 0.001
    wint = apply(wdf, 1, function(x) stats::integrate(stats::dnorm, 
                                                      # !!!!!! *2 ????
                                                      # minT, maxT, mean = x[1], sd = x[2]*2 )$value)
                                                      minT, maxT, mean = x[1], sd = x[2] )$value)
    w = w * 0
    w[widx] = wint
    w
  }
  
  print(paste("Starting bathymetry likelihood calculation..."))
  t0 <- Sys.time()
  tag.pdt$Date <- as.Date(as.POSIXct(tag.pdt$Date, format = findDateFormat(tag.pdt$Date)))
  
  # compute bathy.grid sd
  # bathy.grid=b3;focalDim=3
  sdx = raster::focal(bathy.grid, w = matrix(1, nrow = focalDim, ncol = focalDim), fun = function(x) stats::sd(x,na.rm = T))
  # sdx = t(raster::as.matrix(raster::flip(sdx, 2)))
  #  dat <- base::t(raster::as.matrix(raster::flip(bathy.grid, 2)))
  
  T <- length(tag.pdt[, 1])
  print(paste("Starting iterations through deployment period ", "..."))
  cl = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl, cores = ncores)
  ans=foreach::foreach(i = 1:T) %dopar% {
    # i=2; sens.err=1
    time <- tag.pdt$Date[i]
    bathy.i <- c(tag.pdt$Depth[i] * (1 - sens.err/100), tag.pdt$Depth[i] *(1 + sens.err/100))
    bathy.i=sort(-abs(bathy.i))
    
    lik.bathy <- likint3(raster::as.matrix(bathy.grid), raster::as.matrix(sdx), bathy.i[1], bathy.i[2])
    #lik.bathy <- likint3(bathy.grid,sdx, bathy.i[1], bathy.i[2])
    
    lik.bathy=raster::as.matrix(lik.bathy)/max(raster::as.matrix(lik.bathy),na.rm=T)
    lik.bathy[is.na(lik.bathy) | is.infinite(lik.bathy)]<-0
    lik.bathy
    #  if (i == 1) {
    #   L.bathy <- array(0, dim = c(dim(lik.bathy)[1:2], length(dateVec)))
    #   lon.agg <- seq(raster::extent(bathy.grid)[1], raster::extent(bathy.grid)[2], 
    #                   length.out = dim(bathy.grid)[2])
    #    lat.agg <- seq(raster::extent(bathy.grid)[3], raster::extent(bathy.grid)[4], 
    #                   length.out = dim(bathy.grid)[1])
    #  }
    #  idx <- which(dateVec == as.Date(time))
    #L.bathy[, , idx] = as.matrix(lik.bathy)/max(as.matrix(lik.bathy), na.rm = T)
  }
  parallel::stopCluster(cl)
  
  L. <- array(0, dim = c(dim(ans[[1]])[1], dim(ans[[2]])[2], length(dateVec)))
  # switch from list to array
  for (i in 1:T) {
    time <- tag.pdt$Date[i]
    idx <- which(dateVec == as.Date(time))
    L.[, , idx] = ans[[i]]
  }
  
  L.bathy=aperm(L.,c(2,1,3))
  L.bathy=L.bathy[,dim(L.bathy)[2]:1,]
  
  lon.agg <- seq(raster::extent(bathy.grid)[1], raster::extent(bathy.grid)[2], 
                 length.out = dim(bathy.grid)[2])
  lat.agg <- seq(raster::extent(bathy.grid)[3], raster::extent(bathy.grid)[4], 
                 length.out = dim(bathy.grid)[1])
  
  print(paste("Making final likelihood raster..."))
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.bathy <- list(x = lon.agg, y = lat.agg, z = L.bathy)
  ex <- raster::extent(list.bathy)
  L.bathy <- raster::brick(list.bathy$z, xmn = ex[1], xmx = ex[2], 
                           ymn = ex[3], ymx = ex[4], transpose = T, crs)
  L.bathy <- raster::flip(L.bathy, direction = "y")
  L.bathy[L.bathy < 0] <- 0
  t1 <- Sys.time()
  print(paste("Bathymetric calculations took ", round(as.numeric(difftime(t1, 
                                                                          t0, units = "mins")), 2), "minutes..."))
  return(L.bathy)
}
