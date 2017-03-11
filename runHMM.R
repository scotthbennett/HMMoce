#' TITLE
#' 
#' @param likVec is vector of length 5 and acts as 'switch' to determine which 
#'   likelihoods are calculated. Default is likVec=rep(1, length.out=5) but 
#'   change any of the 1s to something else to turn off that likelihood. The 
#'   vector at position 1 to 5 corresponds to likehoods for light, sst, ohc, 
#'   woa, hycom respectively.
#' @param parVec is vector of possible movement values as input to 
#'   \code{calc.param}
#' @param bndVec is vector of possible minBounds values as input to \code{}
#' @param inilocList is a list in which elements contain individual iniloc 
#'   variables (see \code{make.L.mod}).
#' @param pttList is a list in which elements are ptt IDs that correspond to the
#'   iniloc values in inilocList.
#' @param envDir is path to directory that contains a folder for each ptt value.
#'   The function will look in envDir/ptt/sst for SST data and envDir/ptt/hycom 
#'   for Hycom data.
#' @param dataDir is same as envDir except for input tag data. Function looks in
#'   dataDir/ptt/ for the tag data files.

dataDir <- '/home/rstudio/HMMoce_run/data/'
#dataDir <- '~/Documents/WHOI/RCode/HMMoce_run_data/data/'
#envDir <- '~/Documents/WHOI/RCode/HMMoce_run_data/env_data/'
envDir <- '/home/rstudio/HMMoce_run/env_data/'
source('/home/rstudio/HMMoce_run/formatTracks.r')
source('/home/rstudio/HMMoce_run/compareTracks.r')
parVec <- c(1, 2, 4)
pttList <- c(141259, 141257, 141256, 141254)
gpeNo <- c(4,1,1,5)
bndVec <- c(5, 10)

# TAG/POPUP DATES AND LOCATIONS (dd, mm, YYYY, lat, lon)
inilocList <- list(data.frame(matrix(c(13, 10, 2015, 41.3, -69.27, 
                                       10, 4, 2016, 40.251, -36.061), nrow = 2, ncol = 5, byrow = T)),
                   data.frame(matrix(c(15, 10, 2015, 41.637, -69.706, 
                                       12, 4, 2016, 37.751, -71.649), nrow = 2, ncol = 5, byrow = T)),
                   data.frame(matrix(c(13, 10, 2015, 41.575, -69.423, 
                                       24, 2, 2016, 26.6798, -69.0147), nrow = 2, ncol = 5, byrow = T)),
                   data.frame(matrix(c(21, 10, 2015, 41.597,	-69.445, 
                                       5, 2, 2016, 38.589, -54.874), nrow = 2, ncol = 5, byrow = T)))
for (i in 1:4){
  names(inilocList[[i]]) <- list('day','month','year','lat','lon')
}

sp.limList <- list(list(lonmin = -82, lonmax = -25, latmin = 15, latmax = 50),
                   list(lonmin = -85, lonmax = -50, latmin = 20, latmax = 50),
                   list(lonmin = -95, lonmax = -52, latmin = 10, latmax = 55),
                   list(lonmin = -85, lonmax = -45, latmin = 30, latmax = 50))


runHMM <- function(likVec=c(1,1,1,1,1), inilocList, pttList, sp.limList, bndVec, parVec){
  
  for (ii in 1:length(pttList)){
    #--------------------------------
    # SET ALL INITIALS AND LOAD DATA
    #--------------------------------
    # READ IN TAG DATA
    # TAG/POPUP DATES AND LOCATIONS (dd, mm, YYYY, lat, lon)
    ptt <- pttList[ii]
    iniloc <- inilocList[[ii]]
    tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
    pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')
    
    # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
    dateVec <- as.Date(seq(tag, pop, by = 'day')) 
    
    # READ IN DATA FROM WC FILES
    myDir <- paste(dataDir, ptt, '/', sep='')
    
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
    sp.lim <- sp.limList[[ii]]
    
    if (exists('sp.lim')){
      locs.grid <- setup.locs.grid(sp.lim)
    } else{
      locs.grid <- setup.locs.grid(locs)
      sp.lim <- list(lonmin = min(locs.grid$lon[1,]), lonmax = max(locs.grid$lon[1,]),
                     latmin = min(locs.grid$lat[,1]), latmax = max(locs.grid$lat[,1]))
    }
    
    # IF YOU NEED TO DOWNLOAD SST DATA
    sst.dir <- paste(envDir, ptt, '/sst/', sep='')
    dir.create(file.path(sst.dir), recursive = TRUE, showWarnings = FALSE)
    get.sst.dates <- sst.udates[!(sst.udates %in% as.Date(substr(list.files(sst.dir), 8, 17)))]
    if (length(get.sst.dates) > 0) get.env(get.sst.dates, ptt = ptt, type = 'sst', spatLim = sp.lim, save.dir = sst.dir)
    
    # HYCOM DATA
    hycom.dir <- paste(envDir, ptt, '/hycom/', sep='')
    dir.create(file.path(hycom.dir), recursive = TRUE, showWarnings = FALSE)
    get.pdt.dates <- pdt.udates[!(pdt.udates %in% as.Date(substr(list.files(hycom.dir), 8, 17)))]
    if (length(get.pdt.dates) > 0) get.env(get.pdt.dates, ptt=ptt, type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)
    
    # AND/OR WOA DATA
    woa.dir <- envDir
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
      L.1 <- calc.srss(light, locs.grid = locs.grid, dateVec = dateVec, res=0.25)
      t1 <- Sys.time()
      print(paste('Light calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
    }
    
    if (likVec[2] == 1){
      t0 <- Sys.time()
      L.2 <- calc.sst.par(tag.sst, ptt, sst.dir = sst.dir, dateVec = dateVec, sens.err = 2.5)
      t1 <- Sys.time()
      print(paste('SST calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
    }
    
    if (likVec[3] == 1){
      t0 <- Sys.time()
      L.3 <- calc.ohc.par(pdt, ptt, ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '', use.se = FALSE)
      t1 <- Sys.time()
      print(paste('OHC calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
    }
    
    if (likVec[4] == 1){
      t0 <- Sys.time()
      L.4 <- calc.woa.par(pdt, ptt, woa.data = woa.quarter, focalDim = 9, dateVec = dateVec, use.se = F)
      t1 <- Sys.time()
      print(paste('WOA calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
    }
    
    if (likVec[5] == 1){
      t0 <- Sys.time()
      L.5 <- calc.hycom.par(pdt, ptt, hycom.dir, focalDim = 9, dateVec = dateVec, use.se = F)
      t1 <- Sys.time()
      print(paste('HYCOM calculations took ', round(as.numeric(difftime(t1, t0, units='mins')), 2), 'minutes...'))
    }
    
    #----------------------------------------------------------------------------------#
    # LIST, RESAMPLE, SAVE
    #----------------------------------------------------------------------------------#
    
    L.rasters <- mget(ls(pattern = 'L\\.'))
    resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
    L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])
    save.image(paste(myDir, ptt, '_likelihoods.RData', sep=''))
    
    # Figure out appropriate L combinations
    likVec <- which(likVec == 1)
    if (length(likVec) > 2){
      L.idx <- c(combn(likVec, 2, simplify=F), combn(likVec, 3, simplify=F))
    } else{
      L.idx <- combn(likVec, 2, simplify=F)
    }
    
    for (tt in 1:length(L.idx)){

      #----------------------------------------------------------------------------------#
      # COMBINE LIKELIHOOD MATRICES
      #----------------------------------------------------------------------------------#
      L <- make.L.mod(L1 = L.res[[1]][L.idx[[tt]]],
                      L.mle.res = L.res$L.mle.res, dateVec = dateVec,
                      locs.grid = locs.grid, iniloc = iniloc, bathy=bathy,
                      pdt = pdt)
      
      L.mle <- L$L.mle
      L <- L$L
      g <- L.res$g
      lon <- g$lon[1,]
      lat <- g$lat[,1]
      g.mle <- L.res$g.mle
      
      for (bnd in bndVec){
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
          setwd(myDir)
          write.table(tr, file=paste(ptt, '_HMM_track.csv', sep=''), sep = ',', col.names = T)
          df <- formatTracks(trackDir = myDir, ptt = ptt, gpeNo = gpeNo[ii])
          res <- compareTracks(df)
          res[[4]] <- apply(res$gcd, 2, FUN=function(x) mean(x, na.rm=T))
          res[[5]] <- apply(res$gcd, 2, FUN=function(x) sd(x, na.rm=T))
          
          # WRITE OUT RESULTS
          outVec <- matrix(c(ptt=ptt, minBounds = bnd, migr.spd = i, rmseLon=res$rmse.lon, rmseLat=res$rmse.lat,
                      gcdMean=res[[4]], gcdSD=res[[5]], paste(L.idx[[tt]],collapse='')), ncol=28)
          write.table(outVec,paste(dataDir, 'outVec_results.csv', sep=''), sep=',', col.names=F, append=T)
          #names(out) <- list('rownum','ptt', 'minBnd','migr.spd','rmselon.ti','rmselon.tib','rmselon.kf','rmselon.kfb','rmselon.gpe','rmselon.hmm',
          #                   'rmselat.ti','rmselat.tib','rmselat.kf','rmselat.kfb','rmselat.gpe','rmselat.hmm',
          #                   'gcdm.ti','gcdm.tib','gcdm.kf','gcdm.kfb','gcdm.gpe','gcdm.hmm',
          #                   'gcdsd.ti','gcdsd.tib','gcdsd.kf','gcdsd.kfb','gcdsd.gpe','gcdsd.hmm', 'L.idx')
          
        } # parVec loop
      } # bndVec loop
    } # L.idx loop
  }
  

}

