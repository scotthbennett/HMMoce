#' Hycom Profile LIkelihood in Parallel
#'
#' Calculate Hycom profile probability surface in parallel
#'
#' @param pdt input PAT data see \code{\link{extract.pdt}}
#' @param ptt name of tag i.e. 123645
#' @param hycom.dir directory of downloaded hycom (or other)data
#' @param focalDim is integer for dimensions of raster::focal used to calculate
#'   sd() of temperature grid cell. Recommend focalDim = 3 if woa.data = woa.one
#'   and 9 if using woa.quarter.
#' @param dateVec vector of complete dates for data range. This should be in 'Date' format
#' @param use.se is logical indicating whether or not to use SE when using regression to predict temperature at specific depth levels.
#' @param ncores specify number of cores, or leave blank and use whatever you have!
#'
#' @return a raster brick of Hycom profile likelihood
#' @export
#' @importFrom foreach %dopar%
#'
#' @examples
#' # load workspace
#' \dontrun{
#' load('~/DATA/blue259.RData')

#' # define hycom.dir
#' hycom.dir = '~/hycom/'
#' # run in parallel
#' res = calc.hycom.par(pdt, ptt, isotherm = '', hycom.dir = hycom.dir, dateVec = dateVec, bathy = T)
#' }


calc.hycom.par <- function(pdt, ptt, hycom.dir, focalDim = 9, dateVec, use.se = TRUE, ncores = parallel::detectCores()){
  
  options(warn=-1)
  
  t0 <- Sys.time()
  print(paste('Starting Hycom profile likelihood calculation...'))
  
  # calculate midpoint of tag-based min/max temps
  pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp) / 2
  
  # get unique time points
  dateVec = lubridate::parse_date_time(dateVec, '%Y-%m-%d')
  
  udates <- unique(lubridate::parse_date_time(pdt$Date, orders = '%Y-%m-%d %H%:%M:%S'))
  T <- length(udates)
  
  print(paste0('Generating profile likelihood for ', udates[1], ' through ', udates[length(udates)]))
  
  nc1 <- RNetCDF::open.nc(dir(hycom.dir, full.names = T)[1])
  depth <- RNetCDF::var.get.nc(nc1, 'depth')
  lon <- RNetCDF::var.get.nc(nc1, 'lon')
  lat <- RNetCDF::var.get.nc(nc1, 'lat')
  # result will be array of likelihood surfaces
  
  L.hycom <- array(0, dim = c(length(lon), length(lat), length(dateVec)))
  
  # BEGIN PARALLEL STUFF  
  
  print('Processing in parallel... ')
  
  # ncores = detectCores()  # should be an input argument
  cl = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl, cores = ncores)
  
  ans = foreach::foreach(i = 1:T) %dopar%{
    
    time <- as.Date(udates[i])
    pdt.i <- pdt[which(pdt$Date == time),]
    
    # open day's hycom data
    nc <- RNetCDF::open.nc(paste(hycom.dir, ptt, '_', as.Date(time), '.nc', sep=''))
    dat <- RNetCDF::var.get.nc(nc, 'water_temp') * RNetCDF::att.get.nc(nc, 'water_temp', attribute='scale_factor') + 
      RNetCDF::att.get.nc(nc, variable='water_temp', attribute='add_offset')
    
    #extracts depth from tag data for day i
    y <- pdt.i$Depth[!is.na(pdt.i$Depth)] 
    y[y < 0] <- 0
    
    #extract temperature from tag data for day i
    x <- pdt.i$MidTemp[!is.na(pdt.i$Depth)]  
    
    # use the which.min
    depIdx = unique(apply(as.data.frame(pdt.i$Depth), 1, FUN = function(x) which.min((x - depth) ^ 2)))
    hycomDep <- depth[depIdx]
    
    # make predictions based on the regression model earlier for the temperature at standard WOA depth levels for low and high temperature at that depth
    suppressWarnings(
      fit.low <- locfit::locfit(pdt.i$MinTemp ~ pdt.i$Depth)
    )
    suppressWarnings(
      fit.high <- locfit::locfit(pdt.i$MaxTemp ~ pdt.i$Depth)
    )
    n = length(hycomDep)
    
    #suppressWarnings(
    pred.low = stats::predict(fit.low, newdata = hycomDep, se = T, get.data = T)
    #suppressWarnings(
    pred.high = stats::predict(fit.high, newdata = hycomDep, se = T, get.data = T)
    
    if (use.se){
      # data frame for next step
      df = data.frame(low = pred.low$fit - pred.low$se.fit * sqrt(n),
                      high = pred.high$fit + pred.high$se.fit * sqrt(n),
                      depth = hycomDep)
    } else{
      # data frame for next step
      df = data.frame(low = pred.low$fit,# - pred.low$se.fit * sqrt(n),
                      high = pred.high$fit,# + pred.high$se.fit * sqrt(n),
                      depth = hycomDep)
    }
    
    # calculate sd using Le Bris neighbor method and focal()
    sd.i = array(NA, dim = c(dim(dat)[1:2], length(depIdx)))
    
    for(ii in 1:length(depIdx)){
      r = raster::flip(raster::raster(t(dat[,,depIdx[ii]])), 2)
      f1 = raster::focal(r, w = matrix(1, nrow = focalDim, ncol = focalDim), fun = function(x) stats::sd(x, na.rm = T))
      f1 = t(raster::as.matrix(raster::flip(f1, 2)))
      sd.i[,,ii] = f1
    }
    
    # make index of dates for filling in lik.prof
    didx = base::match(udates, dateVec)
    
    # setup the likelihood array for each day. Will have length (dim[3]) = n depths
    lik.pdt = array(NA, dim = c(dim(dat)[1], dim(dat)[2], length(depIdx)))
    
    for (b in 1:length(depIdx)) {
      #calculate the likelihood for each depth level, b
      lik.try <- try(likint3(dat[,,depIdx[b]], sd.i[,,b], df[b, 1], df[b, 2]), TRUE)
      
      if(!any(which(lik.try > 0))) class(lik.try) <- 'try-error'
      
      if(class(lik.try) == 'try-error' & use.se == FALSE){
        df[b,1] <- pred.low$fit[b] - pred.low$se.fit[b] * sqrt(n)
        df[b,2] <- pred.high$fit[b] - pred.high$se.fit[b] * sqrt(n)
        
        lik.try <- try(likint3(dat[,,depIdx[b]], sd.i[,,b], df[b, 1], df[b, 2]), TRUE)
        
        if(!any(which(lik.try > 0))) class(lik.try) <- 'try-error'
        
        if (class(lik.try) == 'try-error'){
          lik.try <- dat[,,depIdx[b]] * 0
          warning(paste('Warning: likint3 failed after trying with and without SE prediction of depth-temp profiles. This is most likely a divergent integral for ', time, '...', sep=''))
        }
        
      } else if (class(lik.try) == 'try-error' & use.se == TRUE){
        lik.try <- dat[,,depIdx[b]] * 0
        warning(paste('Warning: likint3 failed after trying with and without SE prediction of depth-temp profiles. This is most likely a divergent integral for ', time, '...', sep=''))
      }
      
      lik.pdt[,,b] <- lik.try
      
    }
    
    lik.pdt0 <- lik.pdt
    lik.pdt0[is.na(lik.pdt0)] <- 0
    use.idx <- unique(which(lik.pdt0 != 0, arr.ind=T)[,3])
    
    # multiply likelihood across depth levels for each day
    lik.pdt <- apply(lik.pdt[,,use.idx], 1:2, prod, na.rm = F)
    
  }
  
  parallel::stopCluster(cl)
  
  # make index of dates for filling in L.hycom
  didx <- base::match(udates, dateVec)
  
  # lapply to normalize
  lik.pdt <- lapply(ans, function(x) x / max(x, na.rm = T))
  
  # fill in L.hycom from the list output
  ii = 1
  for(i in didx){
    L.hycom[,,i] = lik.pdt[[ii]]
    ii = ii+1  
  }
  
  print(paste('Making final likelihood raster...'))
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  L.hycom <- raster::brick(L.hycom, xmn=min(lon-360), xmx=max(lon-360), ymn=min(lat), ymx=max(lat), transpose=T, crs)
  L.hycom <- raster::flip(L.hycom, direction = 'y')
  
  #L.hycom[L.hycom < 0] <- 0
  
  names(L.hycom) = as.character(dateVec)
  
  t1 <- Sys.time()
  print(paste('Hycom profile calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
  
  options(warn=2)
  
  # return hycom likelihood surfaces
  return(L.hycom)
  
}


