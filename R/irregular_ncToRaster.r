#' Coerce irregular grid as nc to raster
#' 
#' Main purpose of this function is for converting bathymetry data that is typically acquired as an irregular grid
#' 
#' @param fname is input filename of the .nc file
#' @param varid is character indicating the name of the variable of interest in the input nc variable
#' @param flip is logical indicating whether raster should be flipped along its y axis. Default is TRUE and function pushes warning to tell you it is flipping.
#' @examples
#' \dontrun{
#' bathy <- irregular_ncToRaster(paste0(bathy.dir, 'bathy.nc'), varid = 'topo')
#' }
#' @export
#' @return a raster version of the input .nc file coerced to regular grid

irregular_ncToRaster <- function(fname, varid, flip = TRUE){
  
  nc <- RNetCDF::open.nc(fname)
  lon <- as.numeric(RNetCDF::var.get.nc(nc, variable = "longitude"))
  lat <- as.numeric(RNetCDF::var.get.nc(nc, variable = "latitude"))
  bdata <- RNetCDF::var.get.nc(nc, variable = varid)
  
  if (flip == FALSE){
    if (lat[2] - lat[1] > 0){
      flip = TRUE
    } else{
      lat = lat[order(lat)]
      flip = FALSE
    }
  }
  
  bathy = list(x = lon, y = lat, data = bdata)
  
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  ex <- raster::extent(bathy)
  bathy <- raster::raster(t(bathy$data), xmn = ex[1], xmx = ex[2], ymn = ex[3], ymx = ex[4], crs)

  if (flip){
    bathy <- raster::flip(bathy, 'y')
    warning('Raster is being flipped along its y axis. If you dont want this behavior set flip = FALSE.')
  }
  
  return(bathy)
  
}