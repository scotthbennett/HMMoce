#' OHC Parallel
#' Calculate Ocean Heat Content (OHC) probability surface in parallel
#'
#' @param pdt input PAT data see \code{\link{extract.pdt}}
#' @param filename is the first part of the filename specified to the download 
#'   function \code{\link{get.env}}. For example, if downloaded files were 
#'   specific to a particular dataset, you may want to identify that with a name
#'   like 'tuna' or 'shark1'. This results in a downloaded filename of, for 
#'   example, 'tuna_date.nc'. This filename is required here so the calc
#'   function knows where to get the env data.
#' @param isotherm if specifying a particular isotherm, otherwise leave blank. default value is 
#' @param ohc.dir directory of downloaded hycom (or other)data
#' @param dateVec vector of complete dates for data range. This should be in 'Date' format
#' @param bathy should the land be flagged out? defaults to TRUE
#' @param use.se is logical indicating whether or not to use SE when using regression to predict temperature at specific depth levels.
#' @param ncores specify number of cores, or leave blank and use whatever you have!
#'
#' @return a raster brick of OHC likelihood
#' @seealso \code{\link{calc.ohc}}
#' @export
#' @importFrom foreach "%dopar%"
#'
#' @examples
#' # load workspace
#' \dontrun{
#' load('~/DATA/blue259.RData')

#' # define ohc.dir
#' ohc.dir = '~/hycom/'
#' # run in parallel
#' res = calc.ohc.par(pdt, filename='tuna', isotherm = '', ohc.dir = ohc.dir, dateVec = dateVec, bathy = T)
#' }


calc.ohc.par <- function(pdt, filename, isotherm = '', ohc.dir, dateVec, bathy = TRUE, use.se = TRUE, ncores = parallel::detectCores()){
  
  #max_ohc_date = max(as.Date(substr(dir(ohc.dir), 8, 17)))
  #pdt_idx = as.Date(pdt$Date)<=max_ohc_date
  #pdt = pdt[pdt_idx, ]
  
  #dvidx = dateVec <= max_ohc_date
  
  #dateVec = dateVec[dvidx]
  
  options(warn=1)
  
  t0 <- Sys.time()
  print(paste('Starting OHC likelihood calculation...'))
  
  # constants for OHC calc
  cp <- 3.993 # kJ/kg*C <- heat capacity of seawater
  rho <- 1025 # kg/m3 <- assumed density of seawater
  
  # calculate midpoint of tag-based min/max temps
  pdt$MidTemp <- (pdt$MaxTemp + pdt$MinTemp) / 2
  
  # get unique time points
  dateVec = lubridate::parse_date_time(dateVec, '%Y-%m-%d')
  
  udates <- unique(lubridate::parse_date_time(pdt$Date, orders = '%Y-%m-%d %H%:%M:%S'))
  T <- length(udates)

  if(isotherm != '') iso.def <- TRUE else iso.def <- FALSE
  
  print(paste0('Generating OHC likelihood for ', udates[1], ' through ', udates[length(udates)]))
  
  nc1 =  RNetCDF::open.nc(dir(ohc.dir, full.names = T)[1])
  depth <- RNetCDF::var.get.nc(nc1, 'depth')
  lon <- RNetCDF::var.get.nc(nc1, 'lon')
  lat <- RNetCDF::var.get.nc(nc1, 'lat')
# result will be array of likelihood surfaces
  
  L.ohc <- array(0, dim = c(length(lon), length(lat), length(dateVec)))
  start.t <- Sys.time()
  
# BEGIN PARALLEL STUFF  
  
  print('Processing in parallel... ')
  
  # ncores = detectCores()  # should be an input argument
  cl = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl, cores = ncores)
  
  ans = foreach::foreach(i = 1:T) %dopar%{
    
    time <- as.Date(udates[i])
    pdt.i <- pdt[which(pdt$Date == time),]
    
    # open day's hycom data
    nc <- RNetCDF::open.nc(paste(ohc.dir, filename, '_', as.Date(time), '.nc', sep=''))
    dat <- RNetCDF::var.get.nc(nc, 'water_temp') * RNetCDF::att.get.nc(nc, 'water_temp', attribute='scale_factor') + 
      RNetCDF::att.get.nc(nc, variable='water_temp', attribute='add_offset')
    
    #extracts depth from tag data for day i
    y <- pdt.i$Depth[!is.na(pdt.i$Depth)] 
    y[y<0] <- 0
    
    #extract temperature from tag data for day i
    x <- pdt.i$MidTemp[!is.na(pdt.i$Depth)]  
    
    # use the which.min
    depIdx = unique(apply(as.data.frame(pdt.i$Depth), 1, FUN=function(x) which.min((x - depth) ^ 2)))
    hycomDep <- depth[depIdx]
    
    if(bathy){
      mask <- dat[,,max(depIdx)]
      mask[is.na(mask)] <- NA
      mask[!is.na(mask)] <- 1
      for(bb in 1:length(depth)){
        dat[,,bb] <- dat[,,bb] * mask
      }
    }
    
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

    # isotherm is minimum temperature recorded for that time point
    if(iso.def == FALSE) isotherm <- min(df$low, na.rm = T)
    
    # perform tag data integration at limits of model fits
    minT.ohc <- cp * rho * sum(df$low - isotherm, na.rm = T) / 10000
    maxT.ohc <- cp * rho * sum(df$high - isotherm, na.rm = T) / 10000
    
    # Perform hycom integration
    dat[dat < isotherm] <- NA
    dat <- dat - isotherm
    ohc <- cp * rho * apply(dat[,,depIdx], 1:2, sum, na.rm = T) / 10000 
    ohc[ohc == 0] <- NA
    
    # calc sd of OHC
    # focal calc on mean temp and write to sd var
    r = raster::flip(raster::raster(t(ohc)), 2)
    sdx = raster::focal(r, w = matrix(1, nrow = 9, ncol = 9),
                        fun = function(x) stats::sd(x, na.rm = T))
    sdx = t(raster::as.matrix(raster::flip(sdx, 2)))
    
    # compare hycom to that day's tag-based ohc
    #lik.ohc <- likint3(ohc, sdx, minT.ohc, maxT.ohc)
    
    lik.try <- try(likint3(ohc, sdx, minT.ohc, maxT.ohc), TRUE)
    
    if(class(lik.try) == 'try-error' & use.se == FALSE){
      
      # try ohc again with use.se = T
      df = data.frame(low = pred.low$fit - pred.low$se.fit * sqrt(n),
                      high = pred.high$fit + pred.high$se.fit * sqrt(n),
                      depth = hycomDep)
      
      minT.ohc <- cp * rho * sum(df$low - isotherm, na.rm = T) / 10000
      maxT.ohc <- cp * rho * sum(df$high - isotherm, na.rm = T) / 10000
      
      lik.try <- try(likint3(ohc, sdx, minT.ohc, maxT.ohc), TRUE)
      
      if (class(lik.try) == 'try-error'){
        lik.try <- ohc * 0
        warning(paste('Warning: likint3 failed after trying with and without SE prediction of depth-temp profiles. This is most likely a divergent integral for ', time, '...', sep=''))
      }
      
    } else if (class(lik.try) == 'try-error' & use.se == TRUE){
      lik.try <- ohc * 0
      warning(paste('Warning: likint3 failed after trying with and without SE prediction of depth-temp profiles. This is most likely a divergent integral for ', time, '...', sep=''))
    }
    
    lik.ohc <- lik.try
    
    # if(i == 1){
    #   # result will be array of likelihood surfaces
    #   L.ohc <- array(0, dim = c(dim(lik.ohc), length(dateVec)))
    # }
    
    # idx <- which(dateVec == as.Date(time))
    # L.ohc[,,idx] = (lik.ohc / max(lik.ohc, na.rm=T)) - 0.2
    # 
  }
  
  parallel::stopCluster(cl)
  
  # make index of dates for filling in L.ohc
  
  didx = base::match(udates, dateVec)
  
  
  # lapply 
  lik.ohc = lapply(ans, function(x) x / max(x, na.rm = T))
  
  ii = 1
  for(i in didx){
    L.ohc[,,i] = lik.ohc[[ii]]
    ii = ii+1  
  }

  print(paste('Making final likelihood raster...'))
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.ohc <- list(x = lon - 360, y = lat, z = L.ohc)
  ex <- raster::extent(list.ohc)
  L.ohc <- raster::brick(list.ohc$z, xmn=ex[1], xmx=ex[2], ymn=ex[3], ymx=ex[4], transpose=T, crs)
  L.ohc <- raster::flip(L.ohc, direction = 'y')
  
  L.ohc[L.ohc < 0] <- 0
  
  names(L.ohc) = as.character(dateVec)
  
  t1 <- Sys.time()
  print(paste('OHC calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
  
  # return ohc likelihood surfaces
  return(L.ohc)
  
}


