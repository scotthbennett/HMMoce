# RUN BLUE 259 VIA HMMoce
library(HMMoce)

# SETWD
setwd('~/Documents/WHOI/Data/Blues/2015/141259/') 
load('blue259_runL.RData')

# READ IN TAG DATA
ptt <- 141259

# TAG/POPUP DATES AND LOCATIONS (dd, mm, YYYY, lat, lon)
iniloc <- data.frame(matrix(c(13, 10, 2015, 41.3, -69.27, 
                              10, 4, 2016, 40.251, -36.061), nrow = 2, ncol = 5, byrow = T))
colnames(iniloc) = list('day','month','year','lat','lon')
tag <- as.POSIXct(paste(iniloc[1,1], '/', iniloc[1,2], '/', iniloc[1,3], sep=''), format = '%d/%m/%Y', tz='UTC')
pop <- as.POSIXct(paste(iniloc[2,1], '/', iniloc[2,2], '/', iniloc[2,3], sep=''), format = '%d/%m/%Y', tz='UTC')

# VECTOR OF DATES FROM DATA. THIS WILL BE THE TIME STEPS, T, IN THE LIKELIHOODS
dateVec <- as.Date(seq(tag, pop, by = 'day')) 

# READ IN DATA FROM WC FILES
#myDir <- '~/Documents/WHOI/RCode/HMMoce/inst/extdata/' # WHERE YOUR DATA LIVES, THIS IS THE EXAMPLE DATA
myDir <- getwd()
# sst data
tag.sst <- read.wc(ptt, wd = myDir, type = 'sst', tag=tag, pop=pop); 
sst.udates <- tag.sst$udates; tag.sst <- tag.sst$data

# depth-temp profile data
pdt <- read.wc(ptt, wd = myDir, type = 'pdt', tag=tag, pop=pop); 
pdt.udates <- pdt$udates; pdt <- pdt$data

# light data
light <- read.wc(ptt, wd = myDir, type = 'light', tag=tag, pop=pop); 
light.udates <- light$udates; light <- light$data

# OPTIONAL: light data as output from GPE2, different filtering algorithm seems to work better for light likelihood generation
locs <- read.table(paste(ptt, '-Locations-GPE2.csv', sep = ''), sep = ',', header = T, blank.lines.skip = F)

#----------------------------------------------------------------------------------#
# FURTHER PREPARATION
# Set spatial limits and download env data
#----------------------------------------------------------------------------------#

# SET SPATIAL LIMITS, IF DESIRED
sp.lim <- list(lonmin = -82, lonmax = -25, latmin = 15, latmax = 50)

if (exists('sp.lim')){
  locs.grid <- setup.locs.grid(sp.lim)
} else{
  locs.grid <- setup.locs.grid(locs)
  sp.lim <- list(lonmin = min(locs.grid$lon[1,]), lonmax = max(locs.grid$lon[1,]),
                 latmin = min(locs.grid$lat[,1]), latmax = max(locs.grid$lat[,1]))
}

# IF YOU NEED TO DOWNLOAD SST DATA
sst.dir <- paste('my_sst_dir')
get.env(sst.udates, type = 'sst', spatLim = sp.lim, save.dir = sst.dir)

# HYCOM DATA
hycom.dir <- paste('my_hycom_dir')
get.env(pdt.udates, type = 'hycom', spatLim = sp.lim, save.dir = hycom.dir)

# AND/OR WOA DATA
woa.dir <- paste('my_woa_dir')
get.env(type = 'woa', resol = 'quarter')

#----------------------------------------------------------------------------------#
# CALC LIKELIHOODS
#----------------------------------------------------------------------------------#

# GENERATE LIGHT LIKELIHOOD
# SRSS METHOD
L.light <- calc.srss(light, locs.grid = locs.grid, dateVec = dateVec)
# OR
# GPE2 METHOD
L.light <- calc.gpe2(locs, iniloc = iniloc, locs.grid = locs.grid, dateVec = dateVec, errEll = TRUE, gpeOnly = TRUE)


# GENERATE DAILY SST LIKELIHOODS
L.sst <- calc.sst(tag.sst, sst.dir = sst.dir, dateVec = dateVec)

#-------
# GENERATE DAILY OCEAN HEAT CONTENT (OHC) LIKELIHOODS
L.ohc <- calc.ohc(pdt, ptt, ohc.dir = hycom.dir, dateVec = dateVec, isotherm = '')

#-------
# GENERATE DAILY PROFILE LIKELIHOODS
L.prof <- calc.profile(pdt, ptt, dateVec = dateVec, envType = 'hycom', hycom.dir = hycom.dir)
#L.prof.woa <- calc.profile(pdt, dat = woa, lat = lat, lon = lon, dateVec = dateVec, envType = 'woa')

#----------------------------------------------------------------------------------#
# SETUP A COMMON GRID
#----------------------------------------------------------------------------------#
# create a list of all the desired input likelihood rasters
L.rasters <- list(L.sst = L.sst, L.light = L.light)
# L.sst is the resolution/extent we're sampling everything TO
L.res <- resample.grid(L.rasters, L.rasters$L.sst)

# pull some other helpful variables from the resample.grid() output for later use
L.mle.res <- L.res$L.mle.res
g <- L.res$g; lon <- g$lon[1,]; lat <- g$lat[,1]
g.mle <- L.res$g.mle

#----------------------------------------------------------------------------------#
# LOAD AND FORMAT DATAFRAME OF KNOWN LOCATIONS, IF ANY
#----------------------------------------------------------------------------------#

#colnames(known.locs) <- list('date','lat','lon')
#   where 'date' is from as.Date(known.locs$date)

#----------------------------------------------------------------------------------#
# COMBINE LIKELIHOOD MATRICES
#----------------------------------------------------------------------------------#
# this example just uses L.sst and L.light. You can list up to three (L1, L2, L3 inputs).
L <- make.L(L1 = L.res[[1]]$L.sst,
            L2 = L.res[[1]]$L.light,
            L.mle.res = L.mle.res, dateVec = dateVec,
            locs.grid = locs.grid, iniloc = iniloc)

L.mle <- L$L.mle; L <- L$L

#----------------------------------------------------------------------------------#
# FIGURE OUT MOVEMENT PARAMETERS
#----------------------------------------------------------------------------------#

# PROVIDE FIXED KERNEL PARAMETERS
par0 <- c(8.908,10.27,1.152,0.0472)
D1 <- par0[1:2] # parameters for kernel 1. this is migratory behavior mode
D2 <- par0[3:4] # parameters for kernel 2. resident behavior mode

# GENERATE MOVEMENT KERNELS. D VALUES ARE MEAN AND SD PIXELS
K1 <- gausskern(D1[1], D1[2], muadv = 0)
K2 <- gausskern(D2[1], D2[2], muadv = 0)

# MAKE A GUESS AT STATE SWITCHING PROBABILITY
p <- c(0.7, 0.8)

# RUN EXPECTATION-MAXIMIZATION ROUTINE FOR MATRIX, P (STATE SWITCH PROBABILITY)
P.init <- matrix(c(p[1],1-p[1],1-p[2],p[2]),2,2,byrow=TRUE)
P.final <- expmax(P.init, g = g, L = L, K1, K2)

#----------------------------------------------------------------------------------#
# RUN THE FILTER STEP
f <- hmm.filter(g, L, K1, K2, P.final)

# plot if you want to see confidence limits
#res = apply(f$phi[1,,,],2:3,sum, na.rm=T)
#fields::image.plot(lon, lat, res/max(res), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# RUN THE SMOOTHING STEP
s = hmm.smoother(f, K1, K2, P.final)

# plot if you want to see confidence limits
#sres = apply(s[1,,,], 2:3, sum, na.rm=T)
#fields::image.plot(lon, lat, sres/max(sres), zlim = c(.05,1))

#----------------------------------------------------------------------------------#
# GET THE MOST PROBABLE TRACK
#----------------------------------------------------------------------------------#
T <- dim(s)[2]
meanlat <- apply(apply(s, c(2, 4), sum) * repmat(t(as.matrix(g$lat[,1])), T, 1), 1, sum)
meanlon <- apply(apply(s, c(2, 3), sum) * repmat(t(as.matrix(g$lon[1,])), T, 1), 1, sum)

#**track <- calc.track(distr, g)**

plot(meanlon,meanlat,type='l')
fields::world(add=T, fill=T, col='grey')

#=======================================================================================#
## END
#=======================================================================================#


