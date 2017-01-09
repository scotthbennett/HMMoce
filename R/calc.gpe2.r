#' Calculate Light Likelihood from GPE2 Output
#' 
#' \code{calc.gpe2} calculates likelihood estimates for each day of animal tag 
#' data.
#' 
#' Light errors are parameterized using elliptical error values output in 
#' '-Locations.csv' (WC tags).
#' 
#' @param locs is -Locations file output from DAP/Tag Portal for WC tags and 
#'   contains GPS, Argos, and GPE locations as applicable.
#' @param locDates is vector of dates from locs dataframe
#' @param iniloc is 2 x 5 dataframe containing day, month, year, lat, lon for 
#'   both tag and pop locations
#' @param locs.grid is list output from \code{setup.locs.grid}
#' @param dateVec is vector of dates from tag to pop-up in 1 day increments.
#' @param errEll is logical indicating whether error ellipses should be 
#'   generated for light-based likelihoods as given from output of WC-GPE. False
#'   if only longitude should be used. If False, standard deviation on light 
#'   measurements is currently fixed at 0.7 deg longitude following Musyl et al 
#'   2001. Default is FALSE and will use longitude only.
#' @param gpeOnly is logical. If TRUE, locs input is trimmed to GPE positions 
#'   only. This is most applicable in scenarios with FastGPS data and you're 
#'   adding this as a gps input (see "gps" argument).
#' @export
#' @return L is an array of lon x lat likelihood surfaces (matrices) for each 
#'   time point (3rd dimension)
#' @seealso \code{\link{calc.srss}}
#' @examples
#' \dontrun{
#' # light data as output from GPE2, different filtering algorithm seems to work
#' # better for light likelihood generation
#' locs <- read.table(paste(ptt, '-Locations-GPE2.csv', sep = ''), sep = ',', 
#'                    header = T, blank.lines.skip = F)
#' L.light <- calc.gpe2(locs, iniloc = iniloc, locs.grid = locs.grid,
#'                      dateVec = dateVec, errEll = TRUE, gpeOnly = TRUE)
#' }

calc.gpe2 <- function(locs, locDates, iniloc, locs.grid, dateVec, errEll = TRUE, gpeOnly = TRUE){
  
  print(paste('Starting likelihood calculation...'))
  
  # get rid of non-GPE locations if gpeOnly == TRUE and check length of resulting locs file
  if(gpeOnly == TRUE){
    locs <- locs[which(locs$Type == 'GPE'),]
    if(length(locs[,1]) < 1){
      stop('No GPE positions available in input locs file. Make sure you have calculated GPE positions and exported the resulting -Locations.csv file from WC DAP/GPE2 software.')
    }
  } else{
    stop('Error: Currently gpeOnly must be set to TRUE otherwise duplicate dates are handled improperly.')
  }
  
  #if(any(duplicated(locDates))){
    # run a simplify function
  #  locList <- simplifyLocs(locs, locDates)
  #  locs <- locList$locs; locDates <- locList$locDates
  #}
  
  # set up results array
  row <- dim(locs.grid$lon)[1]
  col <- dim(locs.grid$lon)[2]
  lat <- locs.grid$lat[,1]
  lon <- locs.grid$lon[1,]
  L.gpe2 <- array(0, dim = c(col, row, length(dateVec)))
  
  # add tag/pop locations as known
  #ilo <- which.min(abs(locs.grid$lon[1,] - iniloc$lon[1]))
  #ila <- which.min(abs(locs.grid$lat[,1] - iniloc$lat[1]))
  #L.gpe2[ilo, ila, 1] <- 1   # Initial location is known
  
  #elo <- which.min(abs(locs.grid$lon[1,] - iniloc$lon[2]))
  #ela <- which.min(abs(locs.grid$lat[,1] - iniloc$lat[2]))
  #L.gpe2[elo, ela, length(dateVec)] <- 1  # End location is known
  
  print(paste('Starting iterations through deployment period...'))
  
  for(t in 2:(length(dateVec)) - 1){
    
    if(!is.null(locs) & dateVec[t] %in% locDates){
      # set index to identify position in locs file
      idx <- which(locDates == dateVec[t])
      
      if (nrow(locs[idx,]) > 1){
        
        L.gpe.try <- array(0, dim = c(col, row, nrow(locs[idx,])))
        
        for (ii in 1:nrow(locs[idx,])){
          idx.ii <- idx[ii]
          locs.ii <- locs[idx.ii,]
          
          if(locs.ii$Type == 'GPE'){ #locs includes GPE
            if(errEll == FALSE){
              if (locs.ii$Error.Semi.minor.axis < 100000) locs.ii$Error.Semi.minor.axis <- 100000
              # create longitude likelihood based on GPE data
              slon.sd <- locs.ii$Error.Semi.minor.axis / 1000 / 111 #semi minor axis
              # use normally distributed error from position using fixed std dev
              L.gpe.try[,,ii] <- stats::dnorm(t(locs.grid$lon), locs.ii$Longitude, slon.sd)
              
            } else if(errEll == TRUE){
              if (locs.ii$Error.Semi.minor.axis < 100000) locs.ii$Error.Semi.minor.axis <- 100000
              L.gpe.try[,,ii] <- calc.errEll(locs.ii, locs.grid)
              
            }
            
          } 
          
        }
        
        L.gpe2[,,t] <- apply(L.gpe.try, 1:2, sum, na.rm = T)
        
      } else{
        if(locs$Type[idx] == 'GPE'){ #locs includes GPE
          
          if(errEll == FALSE){
            if (locs$Error.Semi.minor.axis[idx] < 100000) locs$Error.Semi.minor.axis[idx] <- 100000
            # create longitude likelihood based on GPE data
            slon.sd <- locs$Error.Semi.minor.axis[idx] / 1000 / 111 #semi minor axis
            # use normally distributed error from position using fixed std dev
            L.light <- stats::dnorm(t(locs.grid$lon), locs$Longitude[idx], slon.sd)
            
            L.gpe2[,,t] <- L.light
            
          } else if(errEll == TRUE){
            if (locs$Error.Semi.minor.axis[idx] < 100000) locs$Error.Semi.minor.axis[idx] <- 100000
            L.gpe2[,,t] <- calc.errEll(locs[idx,], locs.grid)
            
          }
          
        }  
      }
      
    } 
    
    L.gpe2[,,t] = L.gpe2[,,t] / max(L.gpe2[,,t], na.rm=T) - 0.2
    
  }
  
  print(paste('Making final likelihood raster...'))
  
  # this performs some transformations to the likelihood array to convert to useable raster
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.locs <- list(x = locs.grid$lon[1,], y = locs.grid$lat[,1], z = L.gpe2)
  ex <- raster::extent(list.locs)
  L.gpe2 <- raster::brick(list.locs$z, xmn = ex[1], xmx = ex[2], ymn = ex[3], ymx = ex[4], transpose = T, crs)
  L.gpe2 <- raster::flip(L.gpe2, direction = 'y')
  
  L.gpe2[L.gpe2 < 0] <- 0
  
  return(L.gpe2)
  
}
