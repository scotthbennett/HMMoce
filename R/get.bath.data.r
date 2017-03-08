#' Download bathymetry
#' res is either 30 second, one or two minute resolution (.5, 1)
#' @param lonlow numeric indicating minimum longitude extent of desired download (-180 to 180).
#' @param lonhigh see lonlow
#' @param latlow see lonlow
#' @param lathigh see lonlow
#' @param folder is destination folder. Default is a temporary directory.
#' @param seaonly is logical indicating whether you want to mask out anything on land.
#' @param res is numeric indicating resolution in minutes. Choices are 0.5 or 1.
#' @param raster is logical indicating whether you want the function to return a raster or not (a list will be returned).
#' @export
#' @note Be patient! The download can take a few minutes!

get.bath.data <- function(lonlow, lonhigh, latlow, lathigh, folder = tempdir(), seaonly = T, res = c(.5,1), raster=TRUE){
  #require(RNetCDF)
  
  #rot90 <- function(A) {
  #  n <- dim(A)[2]
  #  A <- t(A)
  #  A[n:1, ]
  #}
  
  #fliplr <- function(A){
  #  A = (A)[(ncol(A)):1,]
  #  A
  #}
  
  fname = paste(folder, "request.nc", sep = "/")
  if(res==1){
    cat('ERDDAP downloading: Topography, Smith & Sandwell v11.1, 1/60-degree \n UCSD   (Dataset ID: usgsCeSS111)')
    opt = "http://coastwatch.pfeg.noaa.gov/erddap/griddap/usgsCeSS111.nc?topo[(LATHIGH):(LATLOW)][(LONLOW):(LONHIGH)]&.draw=surface&.vars=longitude|latitude|topo&.colorBar=|||||&.land=over"
    bathid = 'topo'
  }
  if(res==.5){
    cat('ERDDAP downloading: Topography, SRTM30+ Version 1.0, 30 arc second, Global \n 	Scripps   (Dataset ID: usgsCeSrtm30v1)')
    opt ="http://coastwatch.pfeg.noaa.gov/erddap/griddap/usgsCeSrtm30v1.nc?topo[(LATHIGH):(LATLOW)][(LONLOW):(LONHIGH)]&.draw=surface&.vars=longitude|latitude|topo&.colorBar=|||||"
    #opt = 'http://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.nc?altitude[(LATLOW):1:(LATHIGH)][(LONLOW):1:(LONHIGH)]'
    bathid = 'topo'
  }
  opt <- sub("LATLOW", latlow, opt)
  opt <- sub("LATHIGH", lathigh, opt)
  opt <- sub("LONLOW", lonlow, opt)
  opt <- sub("LONHIGH", lonhigh, opt)
  
  utils::download.file(opt, fname)
  
  nc <- RNetCDF::open.nc(fname)
  lon <- as.numeric(RNetCDF::var.get.nc(nc, variable = "longitude"))
  lat <- as.numeric(RNetCDF::var.get.nc(nc, variable = "latitude"))
  bdata <- RNetCDF::var.get.nc(nc, variable = bathid)

  #if(res==1) bdata = rot90(bdata)
  #if(res==.5) bdata = rot90(bdata)	
  
  lat = lat[order(lat)]
  
  if(seaonly==T) bdata[bdata>=0] = NA
  
  bathy = list(x = lon, y = lat, data = bdata)
  
  if(raster){
    crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
    ex <- raster::extent(bathy)
    bathy <- raster::raster(t(bathy$data), xmn = ex[1], xmx = ex[2], ymn = ex[3], ymx = ex[4], crs)
  }
  
  bathy
  
}

