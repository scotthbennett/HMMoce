#' Download and Read Oceanographic Data
#' 
#' \code{get.env} accesses oceanographic data like sea surface temperature from 
#' a remote server and downloads the temporal and spatial extent of interest for
#' further use
#' 
#' @param uniqueDates is a POSIXct vector of desired dates
#' @param filename is first part of the filename specified to the download 
#'   function. For example, if downloaded files were specific to a particular 
#'   dataset, you may want to identify that with a name like 'tuna' or 'shark1'.
#'   This results in a downloaded filename of, for example, 'tuna_date.nc'.
#' @param type is a character string indicating whether you're after sea surface
#'   temperature 'sst', hybrid coordinate ocean model 'hycom', or world ocean 
#'   atlas 'woa' data
#' @param spatLim is a list of spatial limits as \code{list(xmin, xmax, ymin, 
#'   ymax)}
#' @param resol is character describing the desired resolution in degrees if 
#'   type = 'woa', otherwise NULL.
#' @param save.dir is the directory to save the downloaded data to
#' @param sst.type is character indicating type of desired SST product. Choices
#'   are currently OI and GHR SST.
#' @param depLevels is an integer describing which depth levels to download from Hycom (e.g. 1=surface). Default is NULL and all levels are downloaded.
#'   
#' @return nothing, just downloads the data to your local machine
#' @export

get.env <- function(uniqueDates = NULL, filename = NULL, type = NULL, spatLim = NULL, resol = NULL, save.dir = getwd(), sst.type=NULL, depLevels=NULL){
  
  if(is.null(type)){
    
    stop('Type of environmental data desired not specified.')
    
  } else if(type == 'sst'){
    if(is.null(sst.type)){
      warning('Warning: if type=sst then sst.type should be specified. Using default GHRsst.')
      sst.type = 'ghr'
    }
    
    if(sst.type == 'oi'){
      for(i in 1:length(uniqueDates)){
        time <- as.Date(uniqueDates[i])
        repeat{
          get.oi.sst(spatLim, time, filename = paste(filename, '_', time, '.nc', sep = ''), download.file = TRUE, dir = save.dir) # filenames based on dates from above
          tryCatch({
            err <- try(RNetCDF::open.nc(paste(save.dir, filename, '_', time, '.nc', sep = '')), silent = T)
          }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
          if(class(err) != 'try-error') break
        }
      }
    } else if(sst.type == 'ghr'){
      for(i in 1:length(uniqueDates)){
        time <- as.Date(uniqueDates[i])
        repeat{
          get.ghr.sst(spatLim, time, filename = paste(filename, '_', time, '.nc', sep = ''), download.file = TRUE, dir = save.dir) # filenames based on dates from above
          tryCatch({
            err <- try(RNetCDF::open.nc(paste(save.dir, filename, '_', time, '.nc', sep = '')), silent = T)
          }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
          if(class(err) != 'try-error') break
        }
      }
    }
    

  } else if(type == 'hycom'){
    
    for(i in 1:length(uniqueDates)){
      time <- as.Date(uniqueDates[i])
      repeat{
        get.hycom4(spatLim, time, filename = paste(filename, '_', time, '.nc', sep = ''),
                  download.file = TRUE, dir = save.dir, depLevels=1) 
        tryCatch({
          err <- try(RNetCDF::open.nc(paste(save.dir,'/', filename, '_', time, '.nc', sep = '')), silent = T)
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