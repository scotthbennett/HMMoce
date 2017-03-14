# test sunrise optimization with MTI SRSS data

library(rgdal)
library(oce)
library(RODBC)
library(maptools)
library(maps)
library(lubridate)

# start point
sp = c(-64.700, 43.540)

# OCE methods
# rise <- as.POSIXct("2014-10-17 10:19:00", tz="UTC") #+ 4*3600
# set <- as.POSIXct("2014-10-17 21:44:00", tz="UTC") #+ 4*3600
spm = t(matrix(sp))
sday = as.POSIXct("2014-6-17", tz = 'UTC')

rise <- sunriset(spm,  sday, direction = "sunrise", POSIXct.out=TRUE)$time
set <- sunriset(spm, sday, direction = 'sunset', POSIXct.out = TRUE)$time

mismatch <- function(lonlat) 
{
  sunAngle(rise, lonlat[1], lonlat[2])$altitude^2 + sunAngle(set, lonlat[1], lonlat[2])$altitude^2
}
result <- optim(c(1,1), mismatch)

dist <- geodDist(result$par[1], result$par[2], sp[1], sp[2])
cat(sprintf("Infer start point latitude %.2f and longitude %.2f; distance mismatch %.0f km", 
            result$par[2], result$par[1], dist))

msp = as.matrix(t(sp), nrow=1, ncol = 2)
ssp = SpatialPoints(msp, proj4string=CRS("+proj=longlat +datum=WGS84"))

sdate <- as.POSIXct("2014-10-17", tz="UTC") #+ 4*3600

# correct sunrise time in UTC
crise = sunriset(msp, sdate, direction = 'sunrise', POSIXct.out=TRUE)[2]

# correct sunset time in UTC
cset = sunriset(msp, sdate, direction = 'sunset', POSIXct.out=TRUE)[2]

# how much off is tre first day at liberty from the actual start position
crise[1,1] - rise
# Time difference of 16.27241 mins

cset[1,1] - set
# Time difference of -11.69904 mins

# maptools methods
# look at a matrix of positions
lonvec = seq(sp[1]-5, sp[1]+5, by=.25)
latvec = seq(sp[2]-5, sp[2]+5, by=.25)

llmat = expand.grid(lonvec, latvec)

llmat = SpatialPoints(llmat, proj4string=CRS("+proj=longlat +datum=WGS84"))

srmat = sunriset(llmat, sdate, direction = 'sunrise', POSIXct.out=TRUE)[2]
ssmat = sunriset(llmat, sdate, direction = 'sunset', POSIXct.out=TRUE)[2]

# time difference for sr and ss by grid cell
srdiff = srmat[,1]-rise
ssdiff = ssmat[,1]-set

llmat$srdiff = abs(as.numeric(srdiff))
llmat$ssdiff = abs(as.numeric(ssdiff))


# plot differences 
par(mfrow=c(1,2))
# plot(srmat, pch = 19, cex.=.03, col = 4)
# plot(ssmat, pch = 19, cex.=.03, col=3)


plot(llmat, cex = llmat$srdiff/max(llmat$srdiff), pch=1)
plot(ssp, add=T, col=2, pch=19, cex=1.4)
degAxis(1)
degAxis(2)
title('sunrise difference')
legend('bottom', , legend = round(quantile(llmat$srdiff/60)[2:5]), horiz = T, pch = 1, pt.cex = quantile(llmat$srdiff/max(llmat$srdiff))[2:5], title = 'minutes')
map('world', add=T)

plot(llmat, cex = llmat$ssdiff/max(llmat$ssdiff), pch=1)
plot(ssp, add=T, col=2, pch=19, cex=1.4)
degAxis(1)
degAxis(2)
title('sunset difference')
legend('bottom', , legend = round(quantile(llmat$srdiff/60)[2:5]), horiz = T, pch = 1, pt.cex = quantile(llmat$srdiff/max(llmat$srdiff))[2:5], title = 'minutes')
map('world', add=T)


# look at one point for a year

dates = seq(ymd('2010-1-1', tz= 'UTC'), ymd('2010-12-31', tz= 'UTC'), by = 'day')
newlon = numeric(length = 365)
for(i in 1:365){
  
  spm = t(matrix(sp))
  rise <- sunriset(spm,  dates[i], direction = "sunrise", POSIXct.out=TRUE)$time
  set <- sunriset(spm, dates[i], direction = 'sunset', POSIXct.out = TRUE)$time
  
  mismatch <- function(lonlat) 
  {
    sunAngle(rise, lonlat[1], lonlat[2])$altitude^2 + sunAngle(set, lonlat[1], lonlat[2])$altitude^2
  }
  
  newlon[i] <- optim(c(1,1), mismatch)$par[1]
  
}

plot(dates, newlon)
abline(h = sp[1], col = 2, lty = 2, lwd = 2)
title('one year of optimized longitudes \n with known sr and ss')


# now, train the longitude based on tag measured sunrise and sunset and the pervious days position... 


