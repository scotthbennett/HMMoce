#' TITLE
#' 
#' @param likVec is vector of length 5 and acts as 'switch' to determine which likelihoods are calculated. Default is likVec=rep(1, length.out=5) but change any of the 1s to something else to turn off that likelihood. The vector at position 1 to 5 corresponds to likehoods for light, sst, ohc, woa, hycom respectively.

runHMM <- function(){
  
  parVec <- c(1, 2, 4)
  inilocList <- list()
  pttList <- list()
  
  for (ii in 1:length(pttList)){
    #--------------------------------
    # SET ALL INITIALS AND LOAD DATA
    #--------------------------------
    # READ IN TAG DATA
    # TAG/POPUP DATES AND LOCATIONS (dd, mm, YYYY, lat, lon)
    ptt <- pttList[ii]
    iniloc <- inilocList[ii]
    tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
    pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')
    
    # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
    dateVec <- as.Date(seq(tag, pop, by = 'day')) 
    
    # READ IN DATA FROM WC FILES
    myDir <- paste('/home/rstudio/HMMoce_run/data/', ptt, '/', sep='')
    
    # sst data
    tag.sst <- read.wc(ptt, wd = myDir, type = 'sst', tag=tag, pop=pop); 
    sst.udates <- tag.sst$udates; tag.sst <- tag.sst$data
    
    # depth-temp profile data
    pdt <- read.wc(ptt, wd = myDir, type = 'pdt', tag=tag, pop=pop); 
    pdt.udates <- pdt$udates; pdt <- pdt$data
    
    # light data
    light <- read.wc(ptt, wd = myDir, type = 'light', tag=tag, pop=pop); 
    light.udates <- light$udates; light <- light$data
    
    #----------------------------------------------------------------------------------#
    # FURTHER PREPARATION
    # Set spatial limits and download env data
    #----------------------------------------------------------------------------------#
    
    # SET SPATIAL LIMITS, IF DESIRED
    sp.lim <- list(lonmin = -95, lonmax = -52, latmin = 10, latmax = 55)
    
    if (exists('sp.lim')){
      locs.grid <- setup.locs.grid(sp.lim)
    } else{
      locs.grid <- setup.locs.grid(locs)
      sp.lim <- list(lonmin = min(locs.grid$lon[1,]), lonmax = max(locs.grid$lon[1,]),
                     latmin = min(locs.grid$lat[,1]), latmax = max(locs.grid$lat[,1]))
    }
    
    # IF YOU NEED TO DOWNLOAD SST DATA
    sst.dir <- paste('~/HMMoce_run/env_data/', ptt, '/sst/', sep='')
    get.sst.dates <- sst.udates[!which(sst.udates %in% as.Date(substr(list.files(sst.dir), 8, 17)))]
    if (length(get.sst.dates) > 0) get.env(get.sst.dates, ptt = ptt, type = 'sst', spatLim = sp.lim, save.dir = sst.dir)
    
    # HYCOM DATA
    hycom.dir <- paste('~/HMMoce_run/env_data/', ptt, '/hycom/', sep='')
    get.pdt.dates <- pdt.udates[!which(pdt.udates %in% as.Date(substr(list.files(hycom.dir), 8, 17)))]
    if (length(get.hycom.dates) > 0) get.env(get.pdt.dates, type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)
    
    # AND/OR WOA DATA
    woa.dir <- '~/HMMoce_run/env_data/'
    if(!any(list.files(woa.dir) == 'woa.quarter.rda')){
      download.file('https://raw.githubusercontent.com/camrinbraun/camrinbraun.github.io/master/woa.quarter.rda', 'woa.quarter.rda')
    }
    load(paste(woa.dir,'woa.quarter.rda',sep=''))
    
    # GET BATHYMETRY
    bathy <- get.bath.data(sp.lim$lonmin, sp.lim$lonmax, sp.lim$latmin, sp.lim$latmax, res = c(.5))
    
    #----------------------------------------------------------------------------------#
    # CALCULATE ALL LIKELIHOODS
    #----------------------------------------------------------------------------------#
    if (likVec[1] == 1){
      t0 <- Sys.time()
      L.light <- calc.srss(light, locs.grid = locs.grid, dateVec = dateVec)
      t1 <- Sys.time()
      print(paste('Light calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
    }
    
    if (likVec[2] == 1){
      t0 <- Sys.time()
      L.sst <- calc.sst.par(tag.sst, ptt, sst.dir = sst.dir, dateVec = dateVec, sens.err = 2.5)
      t1 <- Sys.time()
      print(paste('SST calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
    }
    
    if (likVec[3] == 1){
      t0 <- Sys.time()
      L.ohc.par <- calc.ohc.par(pdt, ptt, ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = FALSE)
      t1 <- Sys.time()
      print(paste('OHC calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
    }
    
    if (likVec[4] == 1){
      t0 <- Sys.time()
      L.woa <- calc.woa.par(pdt, ptt, woa.data = woa.quarter, focalDim = 9, dateVec = dateVec, use.se = F)
      t1 <- Sys.time()
      print(paste('WOA calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
    }
    
    if (likVec[5] == 1){
      t0 <- Sys.time()
      L.hycom <- calc.hycom.par(pdt, ptt, hycom.dir, focalDim = 9, dateVec = dateVec, use.se = F)
      t1 <- Sys.time()
      print(paste('HYCOM calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
    }
    
    #----------------------------------------------------------------------------------#
    # LIST, RESAMPLE, SAVE
    #----------------------------------------------------------------------------------#
    
    L.rasters <- mget(ls(pattern = 'L.'))
    resamp.idx <- which.min(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
    L.res <- resample.grid.par(L.rasters, L.rasters[resamp.idx])
    save.image(paste(ptt, '_likelihoods.RData', sep=''))
    
    for (t in xxx){
      for (b in bndVec){
        for (i in parVec){
          
          # GET MOVEMENT KERNELS AND SWITCH PROB FOR COARSE GRID
          par0 <- makePar(migr.spd=i, grid=g.mle, L.arr=L.mle, calcP=TRUE)
          #K1 <- par0$K1; K2 <- par0$K2; 
          P.final <- par0$P.final
          
          # GET MOVEMENT KERNELS AND SWITCH PROB FOR FINER GRID
          par0 <- makePar(migr.spd=i, grid=g, L.arr=L)
          K1 <- par0$K1; K2 <- par0$K2#; P.final <- par0$P.final
          
          # RUN THE FILTER STEP
          f <- hmm.filter.ext(g, L, K1, K2, maskL=T, P.final, minBounds = bnd)
          
          # RUN THE SMOOTHING STEP
          s = hmm.smoother(f, K1, K2, P.final)
          
          # GET THE MOST PROBABLE TRACK
          tr <- calc.track(s, g, dateVec)
          
          # COMPARE HMM, GPE3, SPOT
          
          
          # WRITE OUT RESULTS
          outVec <- c(ptt, minBounds = minBounds, migr.spd = i, fitHMM, fitGPE, names(L.in))
          write.table('outVec_results.csv', sep=',', col.names=F, append=T)
        }
      }
    }
  }
  

}

