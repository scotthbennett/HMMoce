#########################
#Viterbi most probable track
##########################

###########################################################
#This code finds the most probable sequence of grid cells through the posteriors of HMMoce
###########################################################
#' Original code (using priors, not posteriors) by Dr. Julie Nielsen, adapted to HMMoce by Dr. Martin Arostegui (aka Dr. Mola)
#' Cite: 
#' Nielsen, J. K., F. Mueter, M. Adkison, S. McDermott, T. Loher, and A. C. Seitz. 2019. Effect of study area bathymetric heterogeneity on parameterization and performance of a depth-based geolocation model for demersal fish. Ecological Modelling 402:18-34.
#' Viterbi, A. 1967. Error bounds for convolutional codes and an asymptotically optimum decoding algorithm. IEEE Transactions on Information Theory, 13: 260-269.
#' Fariselli, P., P. L. Martelli, and R. Casadio. 2005. The posterior-Viterbi: a new decoding algorithm for hidden Markov models. arXiv: q-bio/0501006
#' Lember, J., and A. A. Koloydenko. 2014. Bridging Viterbi and posterior decoding: a generalized risk approach to hidden path inference based on hidden Markov models. Journal of Machie Learning Research 15: 1-58.
#'
#' Parameters
#' @param pars HMM parameters output from opt.params()
#' @param s HMM backward smoothed output of hmm.smoother()
#' @param g original grid on which HMM was run, assigned from L.res$g
#' @param iniloc dataframe of tag and popup locations/dates (2 rows, 5 columns)
#' @param dateVec timestep sequence used in HMM
#' @param lims.lat vector of c(min,max) desired latitude for reduced spatial extent in which to run algorithm
#' @param lims.lon vector of c(min,max) desired longitude for reduced spatial extent in which to run algorithm
#' @param percentile uncertainty percentile of the daily posterior to include in admissible paths (default is 99.9% of daily posterior)

hmm.posteriorviterbi <- function(pars, s, g, iniloc, dateVec, lims.lat, lims.lon, percentile = 99.9){
  
  if (length(pars) == 4){ # if two-state model
    sigmas = pars[1:2]
    sizes = rep(ceiling(sigmas[1]*4),2)
    pb = pars[3:4]
    muadvs = c(0,0)
    K1 <- gausskern.pg(sizes[1], sigmas[1], muadv=muadvs[1])
    if (!is.null(pb)) K2 <- gausskern.pg(sizes[2], sigmas[2], muadv=muadvs[2])
    K <- list(K1,K2)
  } else if (length(pars) == 1){ # if one-state model
    sigmas = pars[1]
    sizes = rep(ceiling(sigmas[1]*4),2)
    pb = NULL
    muadvs = c(0)
    K1 <- gausskern.pg(sizes[1], sigmas[1], muadv=muadvs[1])
    K <- list(K1)
  }
  
  
  # normalize the posteriors from 0-1 for each day
  s_idx = s
  for(i in 1:(dim(s)[2])){
    if (dim(s)[1]==2) {
      s_idx[1, i , ,] = imager::renorm(s[1, i, , ], 0, 1)
      s_idx[2, i , ,] = imager::renorm(s[2, i, , ], 0, 1)
    } else {
      s_idx[1, i , ,] = imager::renorm(s[1, i, , ], 0, 1)
    }
    
  }
  
  # set < X% posterior probability to 0
  s[s_idx < (1-percentile/100)] = 0
  
  #Crop s to outer spatial extent of overall UD of focal deployment (degrees lat/lon) [manually determined with plotRD]
  min.lat.crop<- lims.lat[1]
  max.lat.crop<- lims.lat[2]
  min.lon.crop<- lims.lon[1]
  max.lon.crop<- lims.lon[2]
  
  # Latitudes to remove
  crop.lat.index <- which(g$lat[, 1] < min.lat.crop | g$lat[, 1] > max.lat.crop)
  # Longitudes to remove
  crop.lon.index <- which(g$lon[1, ] < min.lon.crop | g$lon[1, ] > max.lon.crop)
  
  # Crop s (spatially)
  if (length(crop.lat.index)==0 & length(crop.lon.index)==0) {
    L.cropped  <- s[,,,,drop=F]
  } else if (length(crop.lat.index)==0 & !length(crop.lon.index)==0) {
    L.cropped  <- s[,,-crop.lon.index,,drop=F]
  } else if (!length(crop.lat.index)==0 & length(crop.lon.index)==0) {
    L.cropped  <- s[,,,-crop.lat.index,drop=F]
  } else {
    L.cropped  <- s[,,-crop.lon.index,-crop.lat.index,drop=F]
  }
  
  # Match HMMoce L format to Nielsen L format (ready for Viterbi)
  L.vit <- aperm(L.cropped, c(4,3,2,1))
  L.vit <- aperm(L.vit,c(2,1,3,4))
  L.vit <- array(reverse(reverse(reverse(apply(L.vit,2,c)))), dim = c(dim(L.vit)[2],dim(L.vit)[1],dim(L.vit)[3],dim(L.vit)[4]))
  
  # Generate other Viterbi-specific objects
  if (length(crop.lon.index)==0) {
    longs.vit <- g$lon[1, ]
  } else {
    longs.vit <- g$lon[1, -crop.lon.index]
  }
  
  if (length(crop.lat.index)==0) {
    lats.vit <- g$lat[ ,1]
  } else {
    lats.vit <- g$lat[-crop.lat.index, 1]
  }
  rows.vit <- dim(L.vit)[1]
  cols.vit <- dim(L.vit)[2]
  timesteps.vit <- dim(L.vit)[3]
  behavs.vit <- dim(L.vit)[4]
  # Identify starting cell
  ilo.vit <- which.min(abs(iniloc[1,"lon"]-longs.vit))
  ila.vit <- which.min(abs(iniloc[1,"lat"]-rev(lats.vit)))
  
  # Determine behavior state sequence if two-state model used
  distr <- s
  if (dim(distr)[1]==2){
    sv <- -(apply(distr[1,,,], 1, sum) > apply(distr[2,,,], 1, sum)) + 1
    sv <- sv + 1 # two-state model, where 1 is migratory and 2 is resident
  } else {
    sv <- NA # one-state model
  }
  
  # Ensure diffusion kernel has a central cell (requires odd dimensions)
  if (is.na(sv)==TRUE) {
    # One-state model
    if (round(pars[1]) == pars[1]) { # if integer
      pars[1] <- pars[1] + 0.000001 # add negligible speed (yields odd dimensions)
      sigmas = pars[1]
      sizes = rep(ceiling(sigmas[1]*4),2)
      pb = NULL
      muadvs = c(0)
      K1 <- gausskern.pg(sizes[1], sigmas[1], muadv=muadvs[1])
      K <- list(K1)
    } else {}
  } else {
    # Two-state model
    if (round(pars[1]) == pars[1] | round(pars[2]) == pars[2]) { # if integer(s)
      if (round(pars[1]) == pars[1]) {pars[1] <- pars[1] + 0.000001} # add negligible speed (yields odd dimensions)
      if (round(pars[2]) == pars[2]) {pars[2] <- pars[2] + 0.000001} # add negligible speed (yields odd dimensions)
      sigmas = pars[1:2]
      sizes = rep(ceiling(sigmas[1]*4),2)
      pb = pars[3:4]
      muadvs = c(0,0)
      K1 <- gausskern.pg(sizes[1], sigmas[1], muadv=muadvs[1])
      if (!is.null(pb)) K2 <- gausskern.pg(sizes[2], sigmas[2], muadv=muadvs[2])
      K <- list(K1,K2)
    } else {}
  }
  
  #mpt.calc: 4-dimensional array to hold pathway and probability of that pathway [x,y,time,index# of most probable previous cell]
  # and [x,y,time,cumulative probability]
  #Keeps track of the highest probability pathway between each grid cell and grid cells from previous day
  
  mpt.calc<- array(NA, dim=c(rows.vit,cols.vit,timesteps.vit,2)) #X,Y,days, pathway (which index of matrix highest), and probability at node
  mpt.calc[ila.vit,ilo.vit,1,1]<- which(row(mpt.calc[,,1,1])==ila.vit & col(mpt.calc[,,1,1])==ilo.vit) #index # of release location
  mpt.calc[ila.vit,ilo.vit,1,2]<- 1 #Probability at release location
  mpt.calc[,,1,2]<- log(mpt.calc[,,1,2]) 
  
  #Day 1: initial probability
  #image.plot(reverse(mpt.calc[,,1,1]))#index# of release location on first day (this is just a single pixel)
  #image.plot(reverse(mpt.calc[,,1,2]))#log probability of release location (this is just a single pixel)
  
  #i=2 #Recommend running this through a day at a time, starting at i=2, to understand how it works, etc.
  #Iterate through each day
  start.time<- Sys.time()
  for(i in 2:timesteps.vit){ #Start at 2
    # One-state model
    if (is.na(sv)==TRUE) {
      kern <- K[[1]]
      z <- 1 # behavior state used in testing path admissibility for percentile uncertainty
    } else {
      # Two-state model
      if (sv[i]==1){
        kern <- K[[1]] #first movement state (HMMoce: migratory), D1, defined during forward filter/reverse smoothing
        z <- 1 # behavior state used in testing path admissibility for percentile uncertainty
      } else if (sv[i]==2) {
        kern <- K[[2]] #second movement state (HMMoce: resident), D2, defined during forward filter/reverse smoothing
        z <- 2 # behavior state used in testing path admissibility for percentile uncertainty
      }
    }
    
    #Define how different kernels will fit into viterbi study area matrix
    kern.dim<- dim(kern)[1] #symmetrical
    kern.y.lower<- ceiling(kern.dim/2) #lowest row where kernel can be fully placed into larger matrix
    kern.y.upper<- rows.vit - floor(kern.dim/2) #largest row where kernel can be fully placed into larger matrix
    kern.x.lower<- ceiling(kern.dim/2) #lowest col where kernel can be fully placed into larger matrix
    kern.x.upper<- cols.vit - floor(kern.dim/2) #largest col where kernel can be fully placed into larger matrix
    kern.offset<- floor(kern.dim/2) #Value subtracted/added from each x and y when kernel fully fits into larger matrix
    
    #For each grid cell, find which previous day's grid cell has highest probability of transition
    #x<- 1 #Look at specific grid cell
    for(x in 1:cols.vit){ 
      #y<- 1 #Look at specific grid cell
      for(y in 1:rows.vit){ 
        if(!is.na(mpt.calc[y,x,i-1,1])){ #For valid grid cells
          temp.mat<- matrix(0,rows.vit,cols.vit)
          if(y>kern.y.lower & y < kern.y.upper & x>kern.x.lower & x < kern.x.upper){ #if kernel fully within larger matrix
            #place kernel in template so that the kernel is centered on the grid cell:
            temp.mat[(y-kern.offset):(y+kern.offset),(x-kern.offset):(x+kern.offset)]<- kern 
          } else { #Modify kernel so that it fits in study area with kernel still centered on grid cell 
            if(y>kern.y.lower & y< kern.y.upper) {ymin<- y-kern.offset; ymax<- y+kern.offset}
            if(y<=kern.y.lower) {ymin<- 1; y.offset<- kern.y.lower-y; kern.y<- (1:kern.dim)[-(1:y.offset)];ymax<- y+kern.offset}
            if(y>= kern.y.upper){ymax<- rows.vit; y.offset<- y-kern.y.upper; kern.y<- 1:(kern.dim-y.offset);ymin<- y-kern.offset}
            if (x>kern.x.lower & x < kern.x.upper){xmin<- x-kern.offset; xmax<- x+kern.offset}
            if(x<=kern.x.lower) {xmin<- 1; x.offset<- kern.x.lower-x;kern.x<- (1:kern.dim)[-(1:x.offset)];xmax<- x+kern.offset} 
            if(x>= kern.x.upper){xmax<- cols.vit; x.offset<- x-kern.x.upper; kern.x<- 1:(kern.dim-x.offset);xmin<- x-kern.offset} 
            if(length(ymin:ymax)==kern.dim){ y.offset<- 0; kern.y<- 1:kern.dim}
            if(length(xmin:xmax)==kern.dim){ x.offset<- 0; kern.x<- 1:kern.dim}
            temp.mat[ymin:ymax,xmin:xmax]<- kern[kern.y,kern.x]
          }
          #image.plot(reverse(temp.mat)) #for i = 2, only one previous grid cell had probability, shows kernel centered on release loc
          if (behavs.vit==1) {
            temp.update<- L.vit[,,i,1]*temp.mat #Multiply diffusion kernel by likelihood
          } else {
            if (sv[i]==1) {
              temp.update<- L.vit[,,i,1]*temp.mat #Multiply migratory diffusion kernel by likelihood
            } else {
              temp.update<- L.vit[,,i,2]*temp.mat #Multiply resident diffusion kernel by likelihood
            }
          }
          
          #image.plot(reverse(temp.update)) #for i=2, cells on right side have more probability
          if(any(temp.update>0)){
            update.index<- which(row(temp.update)==y & col(temp.update)==x) #index number of grid cell in focus
            update.prob<- log(temp.update) + mpt.calc[y,x,i-1,2]  #cumulative probability of cells in temp.update >0
            update.cells<- which(temp.update>0,arr.ind=T) #Cells in temp.update with prob>0
            
            #This is the step (below) that could potentially be done with parallel computing
            #As you are trying to find the grid cell in temp.update that has the highest cumulative probability
            #So this step does not rely on sequential calculations
            for(j in 1:nrow(update.cells)){ #Loop through all cells in temp.update with prob>0
              if(is.na(mpt.calc[update.cells[j,1],update.cells[j,2],i,1])){#If that cell probability is NA
                mpt.calc[update.cells[j,1],update.cells[j,2],i,1]<- update.index #Add the index number of focus grid cell
                mpt.calc[update.cells[j,1],update.cells[j,2],i,2]<- update.prob[update.cells[j,1],update.cells[j,2]] #update cumulative probability
              } 
              
              #If there is existing probability in that cell, and that probability is less than the probability in temp.update
              if(mpt.calc[update.cells[j,1],update.cells[j,2],i,2] < update.prob[update.cells[j,1],update.cells[j,2]]){
                mpt.calc[update.cells[j,1],update.cells[j,2],i,1]<- update.index #Add the index number of the focus grid cell
                mpt.calc[update.cells[j,1],update.cells[j,2],i,2]<- update.prob[update.cells[j,1],update.cells[j,2]] #update cumulative probability
              }
            }
          }
        }
      }
    }
    print (i)
  }
  Sys.time()-start.time
  
  #Reconstruct MPT from grid cell at last time step highest cumulative probability
  mpt<- matrix(NA,nrow=timesteps.vit,ncol=2) #1 = row, 2=col
  end.loc<- arrayInd(which.max(mpt.calc[,,timesteps.vit,2]),c(rows.vit,cols.vit)) #Grid cell with most cumulative probability at the end
  mpt[timesteps.vit,]<- end.loc
  
  prev.loc<- mpt.calc[end.loc[,1],end.loc[,2],timesteps.vit,1] #Index # of grid cell from previous step that had highest probability moving to end loc 
  for(i in (timesteps.vit-1):1){
    mpt[i,]<- arrayInd(prev.loc,c(rows.vit,cols.vit)) #index number of grid cell with highest probability
    prev.loc<-  mpt.calc[mpt[i,1],mpt[i,2],i,1] #re-set prev loc to next time step
  }
  
  mpt.latlong<- matrix(NA,nrow=timesteps.vit,ncol=2) #Convert row/col to lat/long
  for(i in 1:timesteps.vit){
    mpt.latlong[i,1]<- lats.vit[((rows.vit+1)-mpt[i,1])]
    mpt.latlong[i,2]<- longs.vit[mpt[i,2]]
  }
  
  # Generate dataframe of most probable track by date, lat/lon
  vit.mpt<- data.frame("date"=dateVec,"lat"=mpt.latlong[,1],"lon"=mpt.latlong[,2])
  return(vit.mpt)
}
