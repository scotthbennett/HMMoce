# input spot data, hycom directory, pdt, ohc L, hycom L
load('141254_likelihoods_knock.RData')
L.ohc <- L.3; L.hycom <- L.5
hycom.dir <- '~/Documents/WHOI/RCode/HMMoce_run_data/env_data/141254/hycom/'
spot <- read.table('141254_crawl_track.csv', sep=',', header=T)
str(pdt)

pdt$Date <- as.Date(pdt$Date, format=findDateFormat(pdt$Date))

ohc.breaks = seq(0, 1, length.out=201)
ohc.mid = ohc.breaks[1:(length(ohc.breaks)-1)]
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
ohc.col = jet.colors(length(ohc.breaks)-1) #[as.vector((dataT))]


for (i in 62:length(pdt.udates)){
  time <- pdt.udates[i]
  pdt.i <- pdt[which(pdt$Date == time),]
  
  # open day's hycom data
  nc <- RNetCDF::open.nc(paste(hycom.dir, ptt, '_', as.Date(time), '.nc', sep=''))
  dat <- RNetCDF::var.get.nc(nc, 'water_temp') * RNetCDF::att.get.nc(nc, 'water_temp', attribute='scale_factor') + 
    RNetCDF::att.get.nc(nc, variable='water_temp', attribute='add_offset')
  if(i==1){
    lon <- RNetCDF::var.get.nc(nc, 'lon') - 360
    lat <- RNetCDF::var.get.nc(nc, 'lat')
    depth <- RNetCDF::var.get.nc(nc, 'depth')
  }
  
  depIdx = unique(apply(as.data.frame(pdt.i$Depth), 1, FUN = function(x) which.min((x - depth) ^ 2)))
  hycomDep <- depth[depIdx]
  
  s <- raster::flip(raster::brick(dat[,,depIdx], xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), transpose=T), 2)
  
  # get locfit results
  # make predictions based on the regression model earlier for the temperature at standard WOA depth levels for low and high temperature at that depth
  suppressWarnings(fit.low <- locfit::locfit(pdt.i$MinTemp ~ pdt.i$Depth))
  suppressWarnings(fit.high <- locfit::locfit(pdt.i$MaxTemp ~ pdt.i$Depth))
  n = length(hycomDep)
  pred.low = stats::predict(fit.low, newdata = hycomDep, se = T, get.data = T)
  pred.high = stats::predict(fit.high, newdata = hycomDep, se = T, get.data = T)
  df = data.frame(low = pred.low$fit,# - pred.low$se.fit * sqrt(n),
                  high = pred.high$fit,# + pred.high$se.fit * sqrt(n),
                  depth = hycomDep)
  
  idx <- which(dateVec %in% pdt.udates[i])
  spot.prof <- raster::extract(s, cbind(spot$lon[idx], spot$lat[idx]))
  
  # start plotting
  #800 wide x 550 tall
  outDir <- '~/Documents/WHOI/RCode/HMMoce_run_data/data/141254/check_hycom/'
  png(paste(outDir, ptt, '_prof_', dateVec[idx], '.png', sep=''),
      width=16, height=11, units='in', res=400)
  nf <- layout(matrix(c(1,1,2,5,5,
                        3,3,4,5,5), 2, 5, byrow=T), widths=c(4,4,1.8,4,4), heights=c(5,5))
  #layout.show(nf)
  
  # plot 1 - ohc map
  image(L.ohc[[idx]], maxpixels=ncell(L.ohc[[idx]]), #xlim=xlims, ylim=ylims,
        col=ohc.col, breaks=ohc.breaks, xlab='', ylab='', main='OHC')
  fields::world(add=T, fill=T)
  points(spot$lon[idx], spot$lat[idx], pch=21, col='black',bg='white',cex=2)
  
  # plot 2 - ohc legend
  image(1, ohc.mid, t(as.matrix(ohc.mid)), breaks=ohc.breaks, col=ohc.col, axes=FALSE, xlab="",
        ylab='OHC likelihood')
  axis(2);box();
  
  # plot 3 - hycom map
  image(L.hycom[[idx]], maxpixels=ncell(L.hycom[[idx]]), #xlim=xlims, ylim=ylims,
        col=ohc.col, breaks=ohc.breaks, xlab='', ylab='', main='HYCOM')
  fields::world(add=T, fill=T)
  points(spot$lon[idx], spot$lat[idx], pch=21, col='black',bg='white',cex=2)
  
  # plot 4 - hycom legend
  image(1, ohc.mid, t(as.matrix(ohc.mid)), breaks=ohc.breaks, col=ohc.col, axes=FALSE, xlab="",
        ylab='HYCOM likelihood')
  axis(2);box();
  
  # plot 5 - depth profiles
  maxZ <- max(pdt.i$Depth) + 50
  t.lims <- c(min(pdt.i$MinTemp) - 2, max(pdt.i$MaxTemp) + 2)
  
  ohc.prof <- raster::extract(s, cbind(xyFromCell(L.ohc[[idx]], which.max(L.ohc[[idx]]))))
  hyc.prof <- raster::extract(s, cbind(xyFromCell(L.hycom[[idx]], which.max(L.hycom[[idx]]))))
  if(nrow(ohc.prof) > 1) ohc.prof <- ohc.prof[1,]
  if(nrow(hyc.prof) > 1) hyc.prof <- hyc.prof[1,]
  
  yleg <- .9 * maxZ
  xleg <- .9 * t.lims[2]
  
  plot(pdt.i$MinTemp, pdt.i$Depth, type='l', ylim=c(maxZ, 0), xlim=t.lims,
       xlab='Temp', ylab='Depth (m)')
  lines(pdt.i$MaxTemp, pdt.i$Depth)
  lines(spot.prof, hycomDep, col='red')
  lines(ohc.prof, hycomDep, col='blue')
  lines(hyc.prof, hycomDep, col='green')
  lines(df$low, df$depth, lty=2)
  lines(df$high, df$depth, lty=2)
  legend(xleg, yleg,
         c('pdt','spot','ohcmax','hycmax','locfit'),
         lty=c(1,1,1,1,2),
         col=c('black','red','blue','green','black'))
  dev.off()
  
}

