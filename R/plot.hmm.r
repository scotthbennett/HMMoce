#' Plot track results of HMMoce
#' 
#' \code{plot.hmm} uses HMM output via \code{calc.track} to make simple plots of
#' calculated track and behavior state
#' 
#' @param distr is output array from \code{hmm.smoother}
#' @param track is output dataframe from \code{calc.track}
#' @param dateVec is vector of dates from tag to pop-up location by day.
#' @param known is 3 column data frame containing date, lat, lon of known 
#'   movement track. This is only useful for comparing HMMoce results to known 
#'   track collected by SPOT or GPS, for example. Default is NULL.
#' @param resid is logical indicating whether you want to include a residual plot. This is not yet functional.
#' @param save.plot is logical indicating whether you want the plot written to disk using \code{pdf}.
#' 
#' @return NULL. A plot is rendered on screen or written to disk.
#' @export

plot.hmm <- function(distr, track, dateVec, known = NULL, resid = FALSE, save.plot = FALSE){
  
  ## Show movement as animation
  #if(show.movie) show.movie(s)
  
  ## Calc behaviour
  # get states from s
  sv <- -(apply(distr[1,,,], 1, sum) > apply(distr[2,,,], 1, sum)) + 1
  sv[sv == 0] <- NA
  graphics.off();
  
  if (class(dateVec) != 'Date'){
    stop('Error: dateVec must be of class "Date". See ?as.Date.')
  }
  
  if(save.plot) pdf(paste(ptt, '_track_results.pdf', sep = ''), width = 7, height = 8)
  
  par(mfrow = c(2,1))
  # behavior state plot
  plot(I(sv) ~ dateVec, col = 'grey', type = 'h', lwd=7, ylab = 'Probability of resident', main = 'Estimated behaviour', xlab = '', ylim=c(0,1))
  lines(track$p ~ dateVec, lwd = 2, col = 'red')
  
  # calculated track
  xl <- c(floor(min(lon)), ceiling(max(lon)))
  yl <- c(floor(min(lat)), ceiling(max(lat)))
  plot(track$lon, track$lat, type = 'n', main = 'Estimated movements', ylab = 'Latitude', xlab = 'Longitude', xlim = xl, ylim = yl)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "steelblue1")
  fields::world(add = TRUE, fill = TRUE, col = 'grey60')
  
  if(!is.null(known)){
    lines(known$lon, known$lat, col = 'white')
  }
  
  lines(track$lon, track$lat, col = 'black')
  points(track$lon[1], track$lat[1], bg = 'green', pch = 21)
  TT <- length(track$lon)
  points(track$lon[TT], track$lat[TT], bg = 'red', pch = 21)
  
  if(save.plot) dev.off()
  
  ## NOT YET FUNCTIONAL
  ## Simple diagnostics plot ###
  ## Resample SST
  #if(save.plot) pdf('../plot/sphmmDiagn.pdf',width=7,height=7)
  #if(!save.plot) 
  if (resid){
    dev.new()
    par(mfrow = c(2,2))
    ssterr <- sst - lsst$sst
    sdsst <- sqrt(var(ssterr))
    ssterr <- ssterr / sdsst
    lonerr <- sphmm$meanlon - lsst$lon
    sdlon <- sqrt(var(lonerr))
    lonerr <- lonerr / sdlon
    plot(track$date[ind], ssterr, xlab = 'Date', ylab = 'SST residual', main = 'Residuals through time', ylim = c(-3,3), pch = ".", cex = 3)
    abline(h = c(-2,0,2), col = 'grey', lwd = 2, lty = 2)
    qqnorm(ssterr, main = 'QQ-plot of residuals', pch = ".", cex = 2)
    abline(a = 0, b = 1, col = 'grey', lwd = 2)
    plot(track$date[ind], lonerr, xlab = 'Date', ylab = 'Longitude residual', ylim = c(-3,3), pch = ".", cex = 3)
    abline(h = c(-2,0,2), col = 'grey', lwd = 2, lty = 2)
    qqnorm(lonerr, main = '', pch = ".", cex = 2)
    abline(a = 0, b = 1, col = 'grey', lwd = 2)
  }
  #if(save.plot) dev.off()
    
  }