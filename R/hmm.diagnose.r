hmm.diagnose <- function(res, L.idx, L.res, dateVec, locs.grid, iniloc, bathy, pdt, minT=NULL, maxT=NULL, plot=TRUE){
  #res, L.idx, L.res, dateVec, locs.grid, iniloc, bathy, pdt
  require(raster)
  # collapse L.idx so we can get the tt run.idx of interest
  L.idx.coll <- L.idx
  for (ii in 1:length(L.idx)) L.idx.coll[[ii]] <- paste(L.idx[[ii]], collapse='')
  
  bnd <- res$outVec[[2]]; i <- res$outVec[[3]]; tt <- which(L.idx.coll == res$outVec[[4]])
  ptt <- res$outVec[[1]]
  runName <- paste(ptt,'_idx',tt,'_bnd',bnd,'_par',i,sep='')
  
  #----------------------------------------------------------------------------------#
  # COMBINE LIKELIHOOD MATRICES
  #----------------------------------------------------------------------------------#
  L <- HMMoce::make.L(L1 = L.res[[1]][L.idx[[tt]]],
                      L.mle.res = L.res$L.mle.res, dateVec = dateVec,
                      locs.grid = locs.grid, iniloc = iniloc, bathy = bathy,
                      pdt = pdt)
  L.mle <- L$L.mle
  L <- L$L
  g <- L.res$g
  g.mle <- L.res$g.mle
  lon <- g$lon[1,]
  lat <- g$lat[,1]
  
  # GET MOVEMENT KERNELS AND SWITCH PROB FOR COARSE GRID
  par0 <- HMMoce::makePar(migr.spd=i, grid=g.mle, L.arr=L.mle, p.guess=c(.9,.9), calcP=T)
  #K1 <- par0$K1; K2 <- par0$K2; 
  P.final <- par0$P.final
  
  # GET MOVEMENT KERNELS AND SWITCH PROB FOR FINER GRID
  par0 <- HMMoce::makePar(migr.spd=i, grid=g, L.arr=L, p.guess=c(.9,.9), calcP=F)
  K1 <- par0$K1; K2 <- par0$K2; #P.final <- par0$P.final
  
  # RUN THE FILTER STEP
  if(!is.na(bnd)){
    f <- HMMoce::hmm.filter(g, L, K1, K2, maskL=T, P.final, minBounds = bnd)
    maskL.logical <- TRUE
  } else{
    f <- HMMoce::hmm.filter(g, L, K1, K2, P.final, maskL=F)
    maskL.logical <- FALSE
  }
  nllf <- -sum(log(f$psi[f$psi>0]))
  
  # RUN THE SMOOTHING STEP
  #s <- HMMoce::hmm.smoother(f, K1, K2, L, P.final)
  s <- res$s
  
  # GET THE MOST PROBABLE TRACK
  tr <- HMMoce::calc.track(s, g, dateVec)
  
  ## NEED TO WORK ON THE PLOTTING
  if(plot){
    require(fields)
    fname <- paste(runName, '-DIAG.pdf', sep='')
    
    # page 1 is a data summary
    # PLOT PDT USING GGPLOT
    jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    
    pdt$DateOnly <- as.Date(as.POSIXct(pdt$Date, format = '%Y-%m-%d %H:%M:%S'))
    pdt <- pdt[which(pdt$DateOnly >= as.Date(tag) & pdt$DateOnly <= as.Date(pop)),]
    
    # set series colors
    if(is.null(minT)) minT <- min(pdt$MinTemp, na.rm=T)
    if(is.null(maxT)) maxT <- max(pdt$MaxTemp, na.rm=T)
    
    cbrks <- 33
    tz.breaks = seq(minT, maxT, length.out=cbrks)
    tz.mid = tz.breaks[1:(length(tz.breaks)-1)]
    tz.leg.col = jet.colors(length(tz.breaks)-1) #[as.vector((dataT))]
    pdt$fill <- tz.leg.col[findInterval(pdt$MaxTemp, tz.breaks)]
    tag.sst$fill <- tz.leg.col[findInterval(tag.sst$Temperature, tz.breaks)]
    
    # PLOT SST
    tag.sst$DateOnly <- as.Date(as.POSIXct(tag.sst$Date, format = '%H:%M:%S %d-%b-%Y'))
    tag.sst <- tag.sst[which(tag.sst$DateOnly >= as.Date(tag) & tag.sst$DateOnly <= as.Date(pop)),]
    
    locs$Date <- as.Date(as.POSIXct(locs$Date, format=findDateFormat(locs$Date)))
    locs <- locs[which(locs$Date >= as.Date(tag) & locs$Date <= as.Date(pop)),]
    locs <- locs[which(locs$Type == 'GPE'),]
    
    ymax <- max(pdt$Depth, na.rm=T) + max(pdt$Depth, na.rm=T)*.15
    
    if (length(L.idx[[tt]]) == 1){
      pdf(fname, width=8, height=8)
      
      plot(dateVec, rep(0, length.out=length(dateVec)), type='n',
           ylim=c(ymax, -25), ylab='Depth (m)', xlab='')
      points(pdt$DateOnly, pdt$Depth, pch=21, bg=pdt$fill, col=pdt$fill) 
      points(tag.sst$DateOnly, rep(-20, length.out=nrow(tag.sst)), pch=21,
             #bg='black', col='black')
             bg=tag.sst$fill, col='black')
      points(locs$Date, rep(ymax-.05*ymax, length.out=nrow(locs)), pch=15)
      points(as.Date(tag), 0, pch=24, bg='green')
      points(as.Date(pop), 0, pch=25, bg='red')
      
      par(mfrow=c(2,2))
      
      # get plot rasters
      plot.ras.idx <- which(names(L.res[[1]]) %in% paste('L.', L.idx[[tt]], sep=''))
      
    } else if(length(L.idx[[tt]]) == 2){
      pdf(fname, width=8, height=12)
      
      plot(dateVec, rep(0, length.out=length(dateVec)), type='n',
           ylim=c(ymax, -25), ylab='Depth (m)', xlab='')
      points(pdt$DateOnly, pdt$Depth, pch=21, bg=pdt$fill, col=pdt$fill) 
      points(tag.sst$DateOnly, rep(-20, length.out=nrow(tag.sst)), pch=21,
             #bg='black', col='black')
             bg=tag.sst$fill, col='black')
      points(locs$Date, rep(ymax-.05*ymax, length.out=nrow(locs)), pch=15)
      points(as.Date(tag), 0, pch=24, bg='green')
      points(as.Date(pop), 0, pch=25, bg='red')
      
      par(mfrow=c(3,2))
      
      # get plot rasters
      plot.ras.idx <- which(names(L.res[[1]]) %in% paste('L.', L.idx[[tt]], sep=''))
      
    } else if(length(L.idx[[tt]]) == 3){
      pdf(fname, width=8, height=12)
      
      plot(dateVec, rep(0, length.out=length(dateVec)), type='n',
           ylim=c(ymax, -25), ylab='Depth (m)', xlab='')
      points(pdt$DateOnly, pdt$Depth, pch=21, bg=pdt$fill, col=pdt$fill) 
      points(tag.sst$DateOnly, rep(-20, length.out=nrow(tag.sst)), pch=21,
             #bg='black', col='black')
             bg=tag.sst$fill, col='black')
      points(locs$Date, rep(ymax-.05*ymax, length.out=nrow(locs)), pch=15)
      points(as.Date(tag), 0, pch=24, bg='green')
      points(as.Date(pop), 0, pch=25, bg='red')
      
      par(mfrow=c(3,2))
      
      # get plot rasters
      plot.ras.idx <- which(names(L.res[[1]]) %in% paste('L.', L.idx[[tt]], sep=''))
      
    }
    
    for (t in 1:5){#length(dateVec)){
      
      # plot L.x, first one
      plot(L.res[[1]][plot.ras.idx[1]][[1]][[t]])
      world(add=T, fill=T, col='black')
      title(paste(dateVec[t], ' - ', names(L.res[[1]][plot.ras.idx[1]])))
      
      if(length(L.idx[[tt]]) == 2){
        # plot L.x, second
        plot(L.res[[1]][plot.ras.idx[2]][[1]][[t]])
        world(add=T, fill=T, col='black')
        title(paste(names(L.res[[1]][plot.ras.idx[2]])))
        
        plot(0,0, type='n')
        
      } else if(length(L.idx[[tt]]) == 3){
        # plot L.x, second
        plot(L.res[[1]][plot.ras.idx[2]][[1]][[t]])
        world(add=T, fill=T, col='black')
        title(paste(names(L.res[[1]][plot.ras.idx[2]])))
        
        # plot L.x, third
        plot(L.res[[1]][plot.ras.idx[3]][[1]][[t]])
        world(add=T, fill=T, col='black')
        title(paste(names(L.res[[1]][plot.ras.idx[3]])))
        
      }
      
      # plot L
      image.plot(lon, lat, L[t,,], zlim=c(.01, .81))
      world(add=T, fill=T, col='black')
      title('overall L')
      
      # plot f
      image.plot(lon, lat, f$phi[1,t,,])
      world(add=T, fill=T, col='black')
      title('filter output')
      
      # plot s
      image.plot(lon, lat, s[1,t,,])
      world(add=T, fill=T, col='black')
      title('smoother output')
      
    }
    
    dev.off()
    
    
    
  }
}
