#' @param bound.thr is numeric indicating the percent threshold that is added
#'   and substracted from the bounding box of the filter output from the
#'   previous day before masking. Default is .05 (5%).

mask.L <- function(pred.t, L.t, lon, lat, bound.thr = .05){
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"
  list.pred <- list(x = lon, y = lat, z = pred.t)
  ex <- raster::extent(list.pred)
  r <- raster::flip(raster::raster(t(list.pred$z), xmn = ex[1], xmx = ex[2], ymn = ex[3], ymx = ex[4], crs), 'y')
  r.m <- r
  lwr.thr <- bound.thr * cellStats(r.m, 'max')
  r.m[r.m <= lwr.thr] <- NA
  r.m <- trim(r.m)
  #r.t <- trim(r.m)
  #plot(r.t)
  
  ex <- extent(r.m)
  #p1 <- as(ex, 'SpatialPolygons')  
  x.diff <- .1 * (ex[2] - ex[1])
  y.diff <- .1 * (ex[4] - ex[3])
  ex[1] <- ex[1] - x.diff
  ex[2] <- ex[2] + x.diff
  ex[3] <- ex[3] - y.diff
  ex[4] <- ex[4] + y.diff
  
  minBounds <- par0$migr * raster::res(r)[1]
  x.diff <- ex[2] - ex[1]
  y.diff <- ex[4] - ex[3]
  
  if(x.diff < minBounds | y.diff < minBounds){
    ex <- extent(r.m)
    mid.x <- (ex[2] + ex[1]) / 2
    ex[1] <- mid.x - minBounds / 2
    ex[2] <- mid.x + minBounds / 2
    
    mid.y <- (ex[4] + ex[3]) / 2
    ex[3] <- mid.y - minBounds / 2
    ex[4] <- mid.y + minBounds / 2
    
  }
  
  new.ex <- as(ex, 'SpatialPolygons')  
  
  p.mask <- raster::mask(r, new.ex, updatevalue=1, inverse=T)
  p.mask[p.mask < 1] <- NA
  #r.mask <- p.mask * r
  #plot(r.mask)
  #plot(p1, add=T, col='green')
  #plot(p2, add=T, col='red')
  
  p.mask <- t(raster::as.matrix(raster::flip(p.mask, 'y')))
  
  new.L <- p.mask * L.t
  if(max(new.L, na.rm=T) > 1e-15){
    post <- pred.t * new.L
  } else{
    post <- pred.t
  }
  
  return(post)
}