light <- read.table('141256-LightLoc.csv', sep=',', header=T, skip=2)

light.m <- light[,c('Day','Time','Type')]
light.m[,4] <- as.POSIXct(paste(light.m$Day, light.m$Time), format='%d-%b-%y %H:%M:%S')
light.m[,5] <- NA
light.m[c(1:(nrow(light.m)-1)),5] <- paste(light.m[2:nrow(light.m),4])
light.m <- light.m[c(4:(nrow(light.m)-1)),]
light.m <- light.m[,c(4,5,3)]
names(light.m) <- list('tFirst', 'tSecond', 'type')

