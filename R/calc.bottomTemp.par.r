
calc.bottomTemp.par=function (tag.pdt,bathy.grid, dateVec, focalDim = NULL, sens.err = 1,filename='hycom',hycom.dir,ncores=4){
  # hycom.dir='C:/Users/pgatti/Documents/Halibut/R/functions/whoi_hmmoce_pkg/EnvData/hycom/'
  # filename='hycom'
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
  
  #dateVec=as.Date(substring(ac(tag.pdt$Date,1,10)))
  #--------------------------
  # hycom format
  #--------------------------
  # extract var names in ncdf env (hycom) files
  nc1 = RNetCDF::open.nc(dir(hycom.dir, full.names = T)[1])
  ncnames = NULL
  nmax <- RNetCDF::file.inq.nc(nc1)$nvars - 1
  for (ii in 0:nmax) ncnames[ii + 1] <- RNetCDF::var.inq.nc(nc1,ii)$name
  
  # vars index
  temp.idx <- grep("temp", ncnames, ignore.case = TRUE) - 1
  lat.idx <- grep("lat", ncnames, ignore.case = TRUE) - 1
  lon.idx <- grep("lon", ncnames, ignore.case = TRUE) - 1
  dep.idx <- grep("dep", ncnames, ignore.case = TRUE) - 1
  
  # scale and offset extraction for temperature
  ncatts <- NULL
  nmax <- RNetCDF::var.inq.nc(nc1, temp.idx)$natts - 1
  for (ii in 0:nmax) ncatts[ii + 1] <- RNetCDF::att.inq.nc(nc1,temp.idx, ii)$name
  scale.idx <- grep("scale", ncatts, ignore.case = TRUE) -  1
  if (length(scale.idx) != 0) {
    scale <- RNetCDF::att.get.nc(nc1, temp.idx, attribute = scale.idx)
  }else{
    scale <- 1
  }
  off.idx <- grep("off", ncatts, ignore.case = TRUE) - 1
  if (length(off.idx) != 0) {
    offset <- RNetCDF::att.get.nc(nc1, temp.idx, attribute = off.idx)
  }else{
    offset <- 1
  }
  
  # cell coordinates Depth, lat, lon
  depth <- RNetCDF::var.get.nc(nc1, dep.idx)
  lon <- RNetCDF::var.get.nc(nc1, lon.idx)
  if (length(dim(lon)) == 2)  lon <- lon[, 1]
  if (!any(lon < 180)) lon <- lon - 360
  lat <- RNetCDF::var.get.nc(nc1, lat.idx)
  if (length(dim(lat)) == 2) lat <- lat[1, ]
  
  # max layer
  i=1
  nc <- try(RNetCDF::open.nc(paste(hycom.dir, filename, "_", dateVec[i], ".nc", sep = "")),T)
  while(!class(nc)=="NetCDF"){
    i=i+1
    nc <- try(RNetCDF::open.nc(paste(hycom.dir, filename, "_", dateVec[i], ".nc", sep = "")),T)
  }
  #rm(i)
  dat <- RNetCDF::var.get.nc(nc, temp.idx) * scale + offset
  matz=array(data=NA,dim=dim(dat))
  for(z. in 1:dim(dat)[3]){
    #z.=300
    tt=dat[,,z.]
    ind.=which(is.finite(tt) & !is.na(tt))
    matz[,,z.][ind.]=z.
  }
  
  mz=apply(matz,1:2,function(x)range(x,na.rm=T)[2])
  mz[!is.finite(mz)]<-NA
  
  # index of bottom cells
  # need lon lat  z (indexes) *(nlon*nlat*(z-1))
  ind.=1:length(mz)+dim(mz)[1]*dim(mz)[2]*as.vector(abs(mz)-1)
  #--------------------------
  # hycom format end
  #--------------------------
  
  print(paste("Starting bottom temp likelihood calculation..."))
  t0 <- Sys.time()
  tag.pdt$Date <- as.Date(as.POSIXct(tag.pdt$Date, format = findDateFormat(tag.pdt$Date)))
  
  T <- length(tag.pdt[, 1])
  print(paste("Starting iterations through deployment period ", "..."))
  #L.btemp=array(data=NA,dim=c(dim(dat)[1:2],length(dateVec)))
  # declare cluster parrallel
  cl = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl, cores = ncores)
  ans=foreach::foreach(i = 1:T) %dopar% {
    # i=1; sens.err=1
    time <- tag.pdt$Date[i]
    btemp.i <- c(tag.pdt$MinTemp[i] * (1 - sens.err/100), tag.pdt$MaxTemp[i] *(1 + sens.err/100))
    
    nc <- try(RNetCDF::open.nc(paste(hycom.dir, filename, "_", as.Date(time), ".nc", sep = "")),T)
    # start if class netcdf
    if(class(nc)=="NetCDF"){ 
      # rescale temp data (matrix lon, lat, depth)
      dat <- RNetCDF::var.get.nc(nc, temp.idx) * scale + offset
      btemp=matrix(data=as.vector(dat)[ind.],ncol=dim(dat)[2],nrow=dim(dat)[1])
      btemp.grid=raster::raster(btemp)
      sdx = raster::focal(btemp.grid, w = matrix(1, nrow = focalDim, ncol = focalDim), fun = function(x) stats::sd(x,na.rm = T))
      lik.btemp <- likint3(btemp.grid, sdx, btemp.i[1], btemp.i[2])
      lik.btemp<-raster::as.matrix(lik.btemp)/max(raster::as.matrix(lik.btemp),na.rm=T)
      lik.btemp[is.na(lik.btemp) | is.infinite(lik.btemp)]<-0
    }else{
      lik.btemp=matrix(data=0,nrow=dim(dat)[1],ncol=dim(dat)[2])
    }
    lik.btemp
  }
  parallel::stopCluster(cl)
  
  L. <- array(0, dim = c(dim(ans[[1]])[1], dim(ans[[2]])[2], length(dateVec)))
  # switch from list to array
  for (i in 1:T) {
    time <- tag.pdt$Date[i]
    idx <- which(dateVec == as.Date(time))
    L.[, , idx] = ans[[i]]
  }
  
  print(paste("Making final likelihood raster..."))
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  #list.btemp <- list(x = lon.agg, y = lat.agg, z = L.btemp)
  list.btemp <- list(x = lon, y = lat, z = L.)
  ex <- raster::extent(list.btemp)
  L.btemp <- raster::brick(list.btemp$z, xmn = ex[1], xmx = ex[2], 
                           ymn = ex[3], ymx = ex[4], transpose = T, crs)
  L.btemp <- raster::flip(L.btemp, direction = "y")
  L.btemp[L.btemp < 0] <- 0
  t1 <- Sys.time()
  print(paste("Bottom temp calculations took ", round(as.numeric(difftime(t1, 
                                                                          t0, units = "mins")), 2), "minutes..."))
  return(L.btemp)
}
