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
#'   type = 'woa', otherwise NULL. Choices are 'one' or 'quarter'.
#' @param save.dir is the directory to save the downloaded data to
#' @param sst.type is character indicating type of desired SST product. Choices 
#'   are currently Optimum Interpolation ('oi')
#'   \url{https://www.ncdc.noaa.gov/oisst}, MUR ('mur')\url{https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.html} or a high-resolution composite
#'   ('ghr') \url{https://www.ghrsst.org/}.
#' @param depLevels is an integer describing which depth levels to download from
#'   Hycom (e.g. 1=surface). Default is NULL and all levels are downloaded.
#'   
#' @return nothing, just downloads the data to your local machine
#' @examples
#' \dontrun{
#' # Not run to prevent actual data download
#' sp.lim <- list(lonmin = -82, lonmax = -25, latmin = 15, latmax = 50)
#' # FOR OI SST DATA
#' get.env(as.Date('2015-10-01'), filename='oisst', type = 'sst', 
#' sst.type='oi', spatLim = sp.lim, save.dir = tempdir())
#' 
#' # FOR HYCOM DATA
#' get.env(as.Date('2015-10-01'), filename='hycom', type = 'hycom', 
#' spatLim = sp.lim, save.dir = tempdir())
#' 
#' # FOR WORLD OCEAN ATLAS DATA
#' get.env(type = 'woa', resol = 'quarter', save.dir = woa.dir)
#' }
#' @export

get.env <- function(uniqueDates = NULL, filename = NULL, type = NULL, spatLim = NULL, resol = NULL, save.dir = getwd(), sst.type=NULL, depLevels=NULL, ...){
  
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
          get.oi.sst(spatLim, time, filename = paste(filename, '_', time, '.nc', sep = ''), download.file = TRUE, dir = save.dir, ...) # filenames based on dates from above
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
          get.ghr.sst(spatLim, time, filename = paste(filename, '_', time, '.nc', sep = ''), download.file = TRUE, dir = save.dir, ...) # filenames based on dates from above
          tryCatch({
            err <- try(RNetCDF::open.nc(paste(save.dir, filename, '_', time, '.nc', sep = '')), silent = T)
          }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
          if(class(err) != 'try-error') break
        }
      }
    } else if(sst.type == 'mur'){
      for(i in 1:length(uniqueDates)){
        time <- as.Date(uniqueDates[i])
        
        if (spatLim$lonmax > 180){ ## 0 to 360
          print('Detected input coordinates > 180, downloading multiple files.')
          ex180 <- raster::extent(-180,180,-90,90)
          ex360 <- raster::extent(180,360,-90,90)
          
          ## get part 1 - whatever is above 180
          ex1 <- raster::intersect(raster::extent(unlist(spatLim)), ex360)
          ex1 <- raster::extent(raster::rotate(raster::raster(ex1)))
          if (ex1@xmin == -180) ex1@xmin <- -179.99
          original_dir <- getwd()
          tdir <- tempdir()
          repeat{
            get.mur.sst(c(ex1@xmin, ex1@xmax, ex1@ymin, ex1@ymax), time, filename = paste(filename, '_', time, '_1.nc', sep = ''), download.file = TRUE, dir = tdir)#, ...) # filenames based on dates from above
            tryCatch({
              err <- try(RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_1.nc', sep = '')), silent = T)
            }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
            if(class(err) != 'try-error') break
          }
          
          ## get part 2 - whatever is below 180, if any
          ex2 <- raster::intersect(raster::extent(unlist(spatLim)), ex180)
          #if (ex1@xmax == 180) ex1@xmax <- 179.99
          
          repeat{
            get.mur.sst(c(ex2@xmin, ex2@xmax, ex2@ymin, ex2@ymax), time, filename = paste(filename, '_', time, '_2.nc', sep = ''), download.file = TRUE, dir = tdir)#, ...) # filenames based on dates from above
            tryCatch({
              err <- try(RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_2.nc', sep = '')), silent = T)
            }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
            if(class(err) != 'try-error') break
          }
          
          print('merging those files to a single output file')
          nc1 <- RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_1.nc', sep = ''))
          lon1 <- RNetCDF::var.get.nc(nc1, 'longitude')
          lat1 <- RNetCDF::var.get.nc(nc1, 'latitude')
          dat1 <- RNetCDF::var.get.nc(nc1, 'analysed_sst') # for OI SST
          r1 <- raster::flip(raster::raster(t(dat1), xmn=min(lon1), xmx=max(lon1),
                                           ymn=min(lat1), ymx=max(lat1)), 2)
          
          nc2 <- RNetCDF::open.nc(paste(tdir, '/', filename, '_', time, '_2.nc', sep = ''))
          lon2 <- RNetCDF::var.get.nc(nc2, 'longitude')
          lat2 <- RNetCDF::var.get.nc(nc2, 'latitude')
          dat2 <- RNetCDF::var.get.nc(nc2, 'analysed_sst') # for OI SST
          r2 <- raster::flip(raster::raster(t(dat2), xmn=min(lon2), xmx=max(lon2),
                                            ymn=min(lat2), ymx=max(lat2)), 2)
          
          ## reset original directory after messing w temp files
          setwd(original_dir)
          
          print(paste(save.dir, '/', filename, '_', time, '.nc', sep = ''))
          r1r <- raster::shift(raster::rotate(raster::shift(r1, 180)), 180)
          r3 <- raster::merge(r1r, r2)
          raster::writeRaster(r3, paste(save.dir, '/', filename, '_', time, '.nc', sep = ''), format='CDF')
          if (file.exists(paste(save.dir, '/', filename, '_', time, '.nc', sep = ''))) print(paste0('File output to ', paste(save.dir, '/', filename, '_', time, '.nc', sep = '')))
          
        } else{ ## -180 to 180
          repeat{
            get.mur.sst(spatLim, time, filename = paste(filename, '_', time, '.nc', sep = ''), download.file = TRUE, dir = save.dir, ...) # filenames based on dates from above
            tryCatch({
              err <- try(RNetCDF::open.nc(paste(save.dir, filename, '_', time, '.nc', sep = '')), silent = T)
            }, error=function(e){print(paste('ERROR: Download of data at ', time, ' failed. Trying call to server again.', sep = ''))})
            if(class(err) != 'try-error') break
          }
        }
        
      }
    }
    
    
  } else if(type == 'hycom'){
    
    for(i in 1:length(uniqueDates)){
      time <- as.Date(uniqueDates[i])
      repeat{
        get.hycom(spatLim, time, filename = paste(filename, '_', time, '.nc', sep = ''),
                  download.file = TRUE, dir = save.dir, depLevels=depLevels, ...) 
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
    
    filename <- get.woa(save.dir = save.dir, resol = resol, ...)
    print(paste('WOA data downloaded to ', filename,'...', sep=''))
    
  }
  
  
}
