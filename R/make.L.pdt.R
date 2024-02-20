#' Combine individual source PDT likelihoods
#' 
#' \code{make.L.pdt} combines individual PDT-based likelihoods from various data sources (e.g.
#' HYCOM OHC, HYCOM depth-levels, GLORYS OHC, GLORYS depth-levels) to make overall 
#' combined PDT-based likelihoods for each time point
#' 
#' @param ras.list a list of likelihood rasters
#' @param dateVec is vector of POSIXct dates for each time step of the likelihood
#' @return an overall PDT-based likelihood RasterStack
#' 
#' @examples
#' \dontrun{
#' L <- make.L(ras.list,dateVec)
#' }
#' @export
#'   
#' @author Martin C. Arostegui

make.L.pdt <- function(ras.list,dateVec){
  
  ## generate blank results likelihood
  L <- ras.list[[1]] * 0
  L[is.na(L)] <- 0
  
  ## COMBINE THE LIKELIHOOD LAYERS
  ## for each day:
  for (i in 1:length(dateVec)){
    
    #print(i)
    
    ## get relevant likelihoods across the list
    for (bb in 1:length(ras.list)){
      if (bb == 1){
        s <- raster::stack(ras.list[[bb]][[i]])
      } else{
        s <- raster::stack(s, ras.list[[bb]][[i]])
      }
    }
    
    ## check for layers with all NA
    sum_NA <- raster::cellStats(!is.na(s), sum, na.rm=T) == 0
    
    ## check for layers that sum to 0
    sum_zero <- rep(NA, raster::nlayers(s))
    for (bb in 1:raster::nlayers(s)){
      s_bb <- s[[bb]]
      s_bb[s_bb == 0] <- 1
      
      sum_zero[bb] <- ifelse(raster::cellStats(s_bb, 'sum') == raster::ncell(s), TRUE, FALSE)
      
    }
    
    ## if all layers are NA, just fill with a 0 raster
    if (length(which(sum_NA)) == raster::nlayers(s)){
      ## nothing left for this iteration
      next
      
      ## if a layer is all NA, drop it
    } else if (length(which(sum_NA)) < raster::nlayers(s) & length(which(sum_NA)) > 0){
      s <- s[[-which(sum_NA)]]
    } else{
      ## do nothing
    }
    
    ## sum & normalize whatever layers remain
    if (raster::nlayers(s) == 1){
      L[[i]] <- s / raster::cellStats(s, 'max') ## do not remove NA yet
    } else{
      L[[i]] <- sum(s) / raster::cellStats(sum(s), 'max') ## do not remove NA yet  
    }
    
  }

  print('Finishing make.L.pdt...', sep='')
  
  return(L)
}

