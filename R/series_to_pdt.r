

series_to_pdt <- function(ts, dateVec, depLevels=NULL){
  
  if(is.null(depLevels)) depLevels <- c(0,2,4,6,8,10,12,15,20,25,30,35,40,45,50,60,70,80,90,100,
                                        125,150,200,250,300,350,400,500,600,700,800,900,1000,
                                        1250,1500,2000,2500,3000,4000,5000) # hycom depth levels 1:40
  
  ts$Day <- as.Date(ts$Day, format=findDateFormat(ts$Day))
  ts$dt <- parse_date_time(paste(ts$Day, ts$Time), orders='Ymd HMS')
  ts$dVec <- findInterval(ts$dt, as.POSIXct(dateVec))
  
  sm0 <- plyr::ddply(ts, c("dVec"), function(x) c(nrecTempTS = nrow(x[which(!is.na(x$Temperature) & 
                                                                              !is.na(x$Depth)), ])))
  if (any(sm0$nrecTempTS == 0)) 
    warning("no Temperature/DepthTS available for dates:\n", 
            paste(sm0$date[which(sm0$nrecTempTS == 0)], collapse = ", "))
  if (!is.null(ts$Corrected.Depth)) {
    warning("Corrected.Depth column found. This data will be used instead!")
    ts$Depth <- ts$Corrected.Depth
  }
  
  ts <- ts[which(!is.na(ts$Temperature) & !is.na(ts$Depth)),]
  sm <- plyr::ddply(ts, c("dVec"), function(x) c(nrecTempTS = nrow(x)))
  u.dVec <- unique(ts$dVec)
  pdt.rec <- c()
  
  for (d in u.dVec) {
    print(dateVec[d])
    i <- which(ts$dVec == d)
    x <- ts[i, ]
    
    x$bin <- findInterval(x$Depth, depLevels)
    add <- data.frame(x %>% group_by(bin) %>% summarise(nrecs = n(), MeanTemp = mean(round(Temperature, 2)), MinTemp = min(round(Temperature, 2)), MaxTemp = max(round(Temperature, 2))))
    add$depth <- depLevels[add$bin]
    add$dVec <- d
    pdt.rec <- rbind(pdt.rec, add)
    
  }
  
  #pdt.rec$MeanPDT <- rowMeans(cbind(pdt.rec$MaxTemp, pdt.rec$MinTemp))
  pdt.rec$date <- dateVec[pdt.rec$dVec]
  names(pdt.rec) <- tolower(names(pdt.rec))
  return(pdt.rec)
}
