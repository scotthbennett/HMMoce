meta$hmmoce <- NA
for (i in 1:nrow(meta)){
  setwd(paste('~/ebs/Data/BaskingSharks/batch/', meta$PTT[i],sep=''))
  fileList <- list.files()
  resFiles <- fileList[grep('_res.rda', fileList)]
  if(length(resFiles) == 36){
    meta$hmmoce[i] <- 1
  } else{
    meta$hmmoce[i] <- 0
    print(warning(paste(meta$PTT[i], 'is not complete.')))
  }
  
  if (meta$hmmoce[i] == 1){
    for (ii in 1:length(resFiles)){
      load(resFiles[ii])
      write.table(res$outVec, file='~/ebs/Data/BaskingSharks/batch/baskbatch_res.csv', sep=',', append=T, col.names = F)
    }
  }
  
}

all <- read.table('~/ebs/Data/BaskingSharks/batch/baskbatch_res.csv', sep=',', header=F)
names(all) <- c('rownum','ptt','bnd','migr.spd','Lidx','P1','P2','xmin','xmax','ymin','ymax','resolx','resoly','maskL','nll','name')

# once we have the "best" model fit, we can look at diagnostics
load(fileList[grep('idx4_bnd5_par2', fileList)[1]])

hmm.diagnose(res, ){
  
  load('check2.rda')
  
  # collapse L.idx so we can get the tt run.idx of interest
  for (ii in 1:length(L.idx)) L.idx[[ii]] <- paste(L.idx[[ii]], collapse='')
  
  bnd <- res$outVec[[2]]; i <- res$outVec[[3]]; tt <- which(L.idx == res$outVec[[4]])
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
  
}
