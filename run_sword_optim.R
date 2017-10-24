#library(HMMoce)
setwd('~/Documents/WHOI/RCode/HMMoce/'); devtools::load_all()
dataDir <- '~/Documents/WHOI/Data/Swordfish/batch/'
#envDir <- '~/EnvData/'
setwd(dataDir)

meta <- read.table(paste(dataDir, 'swords_meta.csv',sep=''), sep=',', header=T)

likVec=c(1,2,3,5)

for (ii in 1:nrow(meta)){

  ptt <- meta$PTT[ii]
  iniloc <- data.frame(matrix(c(meta$TagDay[ii], meta$TagMonth[ii], meta$TagYear[ii], meta$TagLat[ii], meta$TagLong[ii], 
                                meta$PopDay[ii], meta$PopMonth[ii], meta$PopYear[ii], meta$PopLat[ii], meta$PopLong[ii]), nrow = 2, ncol = 5, byrow = T))
  names(iniloc) <- list('day','month','year','lat','lon')
  
  tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
  pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')
  
  # VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
  dateVec <- as.Date(seq(tag, pop, by = 'day')) 
  
  # READ IN DATA FROM WC FILES
  myDir <- paste(dataDir, ptt, '/', sep='')
  load(paste(myDir, ptt,'_likelihoods2.RData', sep=''))
  dataDir <- '~/Documents/WHOI/Data/Swordfish/batch/'
  myDir <- paste(dataDir, ptt, '/', sep='')
  
  #----------------------------------------------------------------------------------#
  # LIST, RESAMPLE, SAVE
  #----------------------------------------------------------------------------------#
  
  bathy <- raster::raster('~/Documents/WHOI/Data/Swordfish/batch/sword_bathy.grd')
  
  #L.rasters <- mget(ls(pattern = 'L\\.'))
  L.rasters <- list(L.1=L.1, L.2=L.2, L.3=L.3f, L.4 = L.1 * 0, L.5=L.5t)
  resamp.idx <- which.max(lapply(L.rasters, FUN=function(x) raster::res(x)[1]))
  L.res <- resample.grid(L.rasters, L.rasters[[resamp.idx]])
  #save.image(paste(myDir, ptt, '_likelihoods.RData', sep=''))
  
  # Figure out appropriate L combinations
  if (length(likVec) > 2){
    L.idx <- c(utils::combn(likVec, 2, simplify=F), utils::combn(likVec, 3, simplify=F))
  } else{
    L.idx <- utils::combn(likVec, 2, simplify=F)
  }
  L.idx <- L.idx[c(1:3,7,8)]
  
  bndVec <- c(NA, 5, 10)
  parVec <- c(2,4)
  
  for (tt in 1:length(L.idx)){
    tt.new <- TRUE
    for (bnd in bndVec){
      for (i in parVec){
        
        #if(bnd == bndVec[1] & i == parVec[1]) tt.new <- TRUE
        if(tt.new){
          #----------------------------------------------------------------------------------#
          # COMBINE LIKELIHOOD MATRICES
          #----------------------------------------------------------------------------------#
          L <- make.L(L1 = L.res[[1]][L.idx[[tt]]],
                      L.mle.res = L.res$L.mle.res, dateVec = dateVec,
                      locs.grid = locs.grid, iniloc = iniloc, bathy = bathy,
                      pdt = pdt)
          
          L.mle <- L$L.mle
          L <- L$L
          g <- L.res$g
          g.mle <- L.res$g.mle
          lon <- g$lon[1,]
          lat <- g$lat[,1]
        }
          # GET MOVEMENT KERNELS AND SWITCH PROB FOR COARSE GRID
          par0 <- makePar(migr.spd=i, grid=g, L.arr=L, calcP=T)
          #K1 <- par0$K1; K2 <- par0$K2; 
          P.final <- par0$P.final
          
          #L.rev <- L
          #for (bb in 1:dim(L)[1]){
          #  iddx <- (dim(L)[1]+1) - bb
          #  L.rev[iddx,,] <- L[bb,,]
          #  print(paste('Replacing', iddx, 'with ', bb))
          #}
          
        #P.final[1,1] <- .86; P.final[2,2] <- .875
        #P.final[1,2] <- 1-P.final[1,1]; P.final[2,1] <- 1-P.final[2,2]
        # GET MOVEMENT KERNELS AND SWITCH PROB FOR FINER GRID
        par0 <- makePar(migr.spd=i, grid=g, L.arr=L.rev, calcP=F)
        K1 <- par0$K1; K2 <- par0$K2; #P.final <- par0$P.final
        
        # RUN THE FILTER STEP
        if(!is.na(bnd)){
          f <- hmm.filter(g, L, K1, K2, maskL=T, P.final, minBounds = bnd)
        } else{
          f <- hmm.filter(g, L, K1, K2, P.final, bathy, maskL=F)
        }
        nllf <- -sum(log(f$psi[f$psi>0]))
        
        # RUN THE SMOOTHING STEP
        s <- hmm.smoother(f, K1, K2, L, P.final)
        #s <- hmm.smoother.mod(f, K1, K2, L.rev, P.final, g, bathy)
        
        # GET THE MOST PROBABLE TRACK
        tr <- calc.track(s, g, dateVec)
        plotHMM(s, tr, dateVec, ptt, save.plot = F, behav.pts=T)

        # COMPARE HMM, GPE3, SPOT
        #setwd(myDir)
        
        # WRITE OUT RESULTS
        outVec <- matrix(c(ptt=paste(ptt,'_nll',sep=''), minBounds = bnd, migr.spd = i, paste(L.idx[[tt]],collapse=''), P1 = P.final[1,1], P2 = P.final[2,2], NLL = nllf), ncol=7)
        write.table(outVec,paste(dataDir, 'outVec_results.csv', sep=''), sep=',', col.names=F, append=T)
        #base::save.image(file=paste(dataDir, ptt, '/', ptt, '-HMMoce.RData', sep=''))
        #res <- list(outVec = outVec, s = s, tr = tr, dateVec=dateVec, iniloc=iniloc, grid = raster::res(L.rasters[[resamp.idx]])[1])
        #base::save(res, file=paste(dataDir, ptt, '/', ptt, '-HMMoce_res.RData', sep=''))
        
        tt.new <- FALSE
        
        print(outVec)
        
      } # parvec
    } # bnd
  } # L.idx
  
  
}


pdf(paste('check_',ptt,'_125_2_5.pdf', sep=''), height=18, width=12)
par(mfrow=c(3,2))
for(i in 1:length(dateVec)){
  #plot(L.res[[1]]$L.1[[i]]); world(add=T); title(paste(dateVec[i], '-light'))#; points(spot$lon[i], spot$lat[i])
  plot(L.res[[1]]$L.2[[i]]); world(add=T); title(paste(dateVec[i], '-sst'))#; points(spot$lon[i], spot$lat[i])
  plot(L.res[[1]]$L.3[[i]]); world(add=T); title(paste(dateVec[i], '-ohc'))#; points(spot$lon[i], spot$lat[i])
  plot(L.res[[1]]$L.5[[i]]); world(add=T); title(paste(dateVec[i], '-hycom'))#; points(spot$lon[i], spot$lat[i])
  image.plot(lon,lat,L[i,,]); world(add=T); title(paste(dateVec[i], '-L'))#; points(spot$lon[i], spot$lat[i])
  image.plot(lon,lat,f$phi[1,i,,]); world(add=T); title(paste(dateVec[i], '-f$phi[1]'))#; points(spot$lon[i], spot$lat[i])
  image.plot(lon,lat,s[1,i,,]); world(add=T); title(paste(dateVec[i], '-s[1]'))#; points(spot$lon[i], spot$lat[i])
}

dev.off()

