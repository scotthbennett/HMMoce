#' Convert time series to PDT-like data for use in HMMoce.
#' 
#' \code{bin_TempTS} converts time series data (such as from MT tags) to PDT-like data for use in HMMoce
#' 
#' @param ts is data frame of time series data. Required inputs cols are Date (POSIX), Depth and Temperature. Note that its up to you to use the appropriate columns with these names. For example, Microwave X-tags can have a "Corrected Depth" column. You should re-name that to "Depth" and provide it as input here.
#' @param out_dates is POSIX vector of desired dates as output. 
#' @param bin_res is numeric indicating desired depth bin resolution. Default is 8 
#'   to match native PDT format on WC tags; however, theres no reason why this 
#'   couldnt be higher resolution, especially if consider we later match this to
#'   a 3D product such as HYCOM that is higher depth resolution (usually).
#' @return a dataframe of summarized depth-temperature information that can be used for likelihood calculations in HMMoce
#'   
#' @export
#' 
#' @examples
#' \dontrun{
#' }

bin_TempTS <- function (ts, out_dates, bin_res = 8) 
{
  ts$time_idx <- findInterval(ts$Date, out_dates)
  
	sm0 <- plyr::ddply(ts, c("time_idx"), function(x) c(nrecTempTS = nrow(x[which(!is.na(x$Temperature) & !is.na(x$Depth)), ])))

	#if (any(sm0$nrecTempTS == 0)) 
	#	warning("no Temperature/DepthTS available for dates:\n", 
	#					paste(sm0$date[which(sm0$nrecTempTS == 0)], collapse = ", "))
	#if (!is.null(ts$Corrected.Depth)) {
	#	warning("Corrected.Depth column found. This data will be used instead!")
	#	ts$Depth <- ts$Corrected.Depth
	#}
	
	ts <- ts[which(!is.na(ts$Temperature) & !is.na(ts$Depth)), ]
	sm <- plyr::ddply(ts, c("Date"), function(x) c(nrecTempTS = nrow(x)))
	fdates <- c()
	pdt.rec <- c()
	ud <- unique(ts$time_idx)
	
	for (d in ud) {
		idx <- which(ts$time_idx %in% d)
	  print(out_dates[d])
		#i <- which(as.character(ts$date) == d)
		x <- ts[idx, ]
		depth.range <- range(x$Depth)
		depth <- x$Depth / bin_res
		ii <- which((x$Depth %% bin_res) / bin_res >= 0.5)
		depth[ii] <- ceiling(depth[ii]) * bin_res
		ii <- which((x$Depth %% bin_res) / bin_res < 0.5)
		depth[ii] <- trunc(depth[ii]) * bin_res
		x$Depth <- depth
		unique_depths <- unique(depth)
		if (length(unique_depths) > 2) {
			h <- hist(depth, breaks = unique_depths, plot = F)
			#identifiers <- c("DeployID", "Serial", "Ptt")
			#identifiers <- identifiers[which(identifiers %in% names(ts))]
			add <- plyr::ddply(x[, which(names(x) %in% c("Depth", "Temperature"))], c("Depth"), function(x) c(nrecs = nrow(x), MeanTemp = mean(round(x$Temp, 2)), MinTemp = min(round(x$Temp, 2)), MaxTemp = max(round(x$Temp, 2))))
			add$Date <- out_dates[d]
			add$bin <- 1:nrow(add)
			pdt.rec <- rbind(pdt.rec, add)
		}
	}
	
	pdt.rec$MeanPDT <- rowMeans(cbind(pdt.rec$MaxTemp, pdt.rec$MinTemp))
	#pdt.rec$date <- as.Date(pdt.rec$date)
	return(pdt.rec)
}