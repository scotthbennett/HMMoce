#' Download and Read Oceanographic Data
#' 
#' \code{get.env} accesses oceanographic data like sea surface temperature from
#' a remote server and downloads the temporal and spatial extent of interest for
#' further use
#' 
#' @param uniqueDates is a POSIXct vector of desired dates
#' @param ptt is unique tag identifier
#' @param type is a character string indicating whether you're after sea surface
#'   temperature 'sst' or hybrid coordinate ocean model 'hycom' data
#' @param spatLim is a list of spatial limits as \code{list(xmin, xmax, ymin,
#'   ymax)}
#' @param resol is character describing the desired resolution in degrees if type = 'woa', otherwise NULL.
#' @param save.dir is the directory to save the downloaded data to
#'   
#' @return nothing, just downloads the data to your local machine
#' @export

get.env <- function(uniqueDates = NULL, ptt = NULL, type = NULL, spatLim = NULL, resol = NULL, save.dir = getwd()){
  
  if(is.null(type)){
    
    stop('Type of environmental data desired not specified.')
    
  } else if(type == 'sst'){
    
    for(i in 1:length(uniqueDates)){
      time <- as.Date(uniqueDates[i])
      repeat{
        get.oi.sst(spatLim, time, filename = paste(ptt, '_', time, '.nc', sep = ''), download.file = TRUE, dir = save.dir) # filenames based on dates from above
        tryCatch({
          err <- try(RNetCDF::open.nc(paste(save.dir,'/', ptt, '_', time, '.nc', sep = '')), silent = T)
        }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
        if(class(err) != 'try-error') break
      }
    }

  } else if(type == 'hycom'){
    
    for(i in 1:length(uniqueDates)){
      time <- as.Date(uniqueDates[i])
      repeat{
        get.hycom(spatLim, time, type = 'a', filename = paste(ptt, '_', time, '.nc', sep = ''),
                  download.file = TRUE, dir = save.dir, vars = 'water_temp') 
        tryCatch({
          err <- try(RNetCDF::open.nc(paste(save.dir,'/', ptt, '_', time, '.nc', sep = '')), silent = T)
        }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
        if(class(err) != 'try-error') break
      }
    }
    
  } else if (type == 'woa'){
    
    if(is.null(resol)){
      stop('Error: If type = woa then resol must be specified. See ?get.env for help.')
    }
    
    filename <- get.woa(save.dir = save.dir, resol = resol)
    print(paste('WOA data downloaded to ', filename,'...', sep=''))
    
  }
  
  
}