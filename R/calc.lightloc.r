#' Calculate Light Likelihood from lightloc Output
#' 
#' \code{calc.lightloc} calculates likelihood estimates for each day of animal tag 
#' data.
#' 
#' Light errors are parameterized using elliptical error values output in 
#' '-Locations.csv' (WC tags).
#' 
#' @param lightloc is data frame of light-based location estimates. If errEll is FALSE, only Date (POSIXct), Longitude, Latitude, and Error.Semi.minor.axis are required. If errEll is TRUE, additional required columns are Offset, Offset.orientation, and Error.Semi.major.axis.
#' @param locs.grid is list output from \code{setup.locs.grid}
#' @param dateVec is vector of POSIXct dates for each time step of the likelihood
#' @param errEll is logical indicating whether error ellipses should be 
#'   generated for light-based likelihoods as given from output of WC-GPE. False
#'   if only longitude should be used. If False, standard deviation on light 
#'   measurements is currently fixed at 0.7 deg longitude following Musyl et al 
#'   2011. Default is FALSE and will use longitude only.
#' @references Musyl MK, Domeier ML, Nasby-Lucas N, Brill RW, McNaughton LM, 
#'   Swimmer JY, Lutcavage MS, Wilson SG, Galuardi B, Liddle JB (2011) 
#'   Performance of pop-up satellite archival tags. Mar Ecol Prog Ser
#' @export
#' @return L is an array of lon x lat likelihood surfaces (matrices) for each 
#'   time point (3rd dimension)
#' @seealso \code{\link{calc.srss}}
#' @examples
#' ## Setup for calculating light likelihood
#' # Read the data
#' locsFile <- system.file("extdata", "141259-Locations-lightloc.csv", package = "HMMoce")
#' locs <- read.table(locsFile, sep = ',', header = TRUE, blank.lines.skip = FALSE)
#' 
#' # Set spatial and temporal limits
#' sp.lim <- list(lonmin = -82, lonmax = -25, latmin = 15, latmax = 50)
#' locs.grid <- setup.locs.grid(sp.lim)
#' iniloc <- data.frame(matrix(c(13, 10, 2015, 41.3, -69.27, 10, 4, 2016, 40.251, -36.061),
#'  nrow = 2, ncol = 5, byrow = TRUE))
#' names(iniloc) <- list('day','month','year','lat','lon')
#' tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), 
#' format = '%d/%m/%Y', tz='UTC')
#' pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), 
#' format = '%d/%m/%Y', tz='UTC')
#' dateVec <- as.Date(seq(tag, pop, by = 'day')) 
#' 
#' # Try a calculation
#' L.light <- calc.lightloc(locs[1,], iniloc, locs.grid, dateVec, errEll=TRUE, gpeOnly=TRUE)
#' 
#' \dontrun{
#' # Full example light calculation
#' L.light <- calc.lightloc(locs, iniloc = iniloc, locs.grid = locs.grid,
#'                      dateVec = dateVec, errEll = TRUE, gpeOnly = TRUE)
#' }

calc.lightloc <- function(lightloc, locs.grid, dateVec, errEll = TRUE){
  
  print(paste('Starting light likelihood calculation...'))
  t0 <- Sys.time()
  
  # check date formats match
  lightloc <- lightloc[which(lightloc$Date <= max(dateVec)),]
  lightloc$dateVec <- findInterval(lightloc$Date, dateVec)
  T <- length(dateVec)
  
  # set up results array
  row <- dim(locs.grid$lon)[1]
  col <- dim(locs.grid$lon)[2]
  lat <- locs.grid$lat[,1]
  lon <- locs.grid$lon[1,]
  L.lightloc <- array(0, dim = c(col, row, length(dateVec)))
  
  print(paste('Starting iterations through deployment period...'))
  
  for(t in 1:T){
    
    # data for this time step T
    light.t <- lightloc[which(lightloc$dateVec == t),]

    if (nrow(light.t) > 1){
      
      L.lightloc.try <- array(0, dim = c(col, row, nrow(light.t)))
      
      ## if multiple matches at this time step
      for (ii in 1:nrow(light.t)){
        
        locs.ii <- light.t[ii,]
        
        if(errEll == FALSE){
          if (locs.ii$Error.Semi.minor.axis[1] < 100000) locs.ii$Error.Semi.minor.axis[1] <- 100000
          # create longitude likelihood based on GPE data
          slon.sd <- locs.ii$Error.Semi.minor.axis[1] / 1000 / 111 #semi minor axis
          # use normally distributed error from position using fixed std dev
          L.lightloc.try[,,ii] <- stats::dnorm(t(locs.grid$lon), locs.ii$Longitude[1], slon.sd)
          
        } else if(errEll == TRUE){
          if (locs.ii$Error.Semi.minor.axis[1] < 100000) locs.ii$Error.Semi.minor.axis[1] <- 100000
          L.lightloc.try[,,ii] <- calc.errEll(locs.ii[1,], locs.grid)
          
        } 
        
      }
      
      L.lightloc[,,t] <- apply(L.lightloc.try, 1:2, sum, na.rm = T)
      
      ## normalize
      L.lightloc[,,t] = L.lightloc[,,t] / max(L.lightloc[,,t], na.rm=T)
      
    } else if (nrow(light.t) == 1){
      
      if(errEll == FALSE){
        if (light.t$Error.Semi.minor.axis[1] < 100000) light.t$Error.Semi.minor.axis[1] <- 100000
        # create longitude likelihood based on GPE data
        slon.sd <- light.t$Error.Semi.minor.axis[1] / 1000 / 111 #semi minor axis
        # use normally distributed error from position using fixed std dev
        L.light <- stats::dnorm(t(locs.grid$lon), light.t$Longitude[1], slon.sd)
        
        L.lightloc[,,t] <- L.light
        
      } else if(errEll == TRUE){
        if (light.t$Error.Semi.minor.axis[1] < 100000) light.t$Error.Semi.minor.axis[1] <- 100000
        L.lightloc[,,t] <- calc.errEll(light.t[1,], locs.grid)
        
      }
      
      ## normalize
      L.lightloc[,,t] = L.lightloc[,,t] / max(L.lightloc[,,t], na.rm=T)
      
    } else{
      ## do nothing
    }
    
   
  }
  
  print(paste('Making final likelihood raster...'))
  
  # this performs some transformations to the likelihood array to convert to useable raster
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.locs <- list(x = locs.grid$lon[1,], y = locs.grid$lat[,1], z = L.lightloc)
  ex <- raster::extent(list.locs)
  L.lightloc <- raster::brick(list.locs$z, xmn = ex[1], xmx = ex[2], ymn = ex[3], ymx = ex[4], transpose = T, crs)
  L.lightloc <- raster::flip(L.lightloc, direction = 'y')
  
  L.lightloc[L.lightloc < 0] <- 0
  
  t1 <- Sys.time()
  print(paste('Light calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
  
  return(L.lightloc)
  
}
