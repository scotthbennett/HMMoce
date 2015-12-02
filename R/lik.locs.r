lik.locs <- function(locs,iniloc,g){
  ## Calculate the "data" likelihood, i.e. the likelihood field for each
  ## location observation

  #' @param: locs is -Locations file output from DAP for WC tags and contains
  #'        GPS, Argos, and GPE locations as applicable
  #' @param: iniloc is 2 x 5 dataframe containing day, month, year,
  #'          lat, lon for both tag and pop locations
  #' @param: g is output from setup.grid and indicates extent and resolution
  #'        of grid used to calculate likelihoods
  #' @return: L is array of lon x lat likelihood surfaces (matrices)
  #'          for each time point (3rd dimension)
    
  T <- length(locs$Longitude)
  row <- dim(g$lon)[1]
  col <- dim(g$lon)[2]
  
  L <- array(0,dim=c(col, row, T + 2))
  
  # Initial location is known
  ilo <- which.min(abs(g$lon[1,]-iniloc$lon[1]))
  ila <- which.min(abs(g$lat[,1]-iniloc$lat[1]))
  L[ilo, ila, 1] <- 1
  
  # Calculate data likelihood
  # SD for light-based longitude from Musyl et al. (2001)
  sl.sd <- 35/111 # Converting from kms to degrees
  
  for(t in 1:T){
    if(locs$Type[t] == 'GPS'){
      # if GPS exists then other forms of data for that time point are obsolete
      glo <- which.min(abs(g$lon[1,]-locs$Longitude[t]))
      gla <- which.min(abs(g$lat[,1]-locs$Latitude[t]))
      L[glo, gla, (t+1)] <- 1
      
    } else if(locs$Type[t] == 'Argos'){
      # if Argos exists, GPE positions are obsolete
      alo <- which.min(abs(g$lon[1,]-locs$Longitude[t]))
      ala <- which.min(abs(g$lat[,1]-locs$Latitude[t]))
      L[alo, ala, (t+1)] <- 1
      
    } else if(locs$Type[t] == 'GPE'){
     # create longitude likelihood based on GPE data
     # for now, latitude is ignored
       L.light <- dnorm(t(g$lon), locs$Longitude[t], sl.sd) # Longitude data
      L[,, (t + 1)] <- (L.light / max(L.light, na.rm = T)) - .05
      
    } else{}

      }
  
  # End location is known
  elo <- which.min(abs(g$lon[1,]-iniloc$lon[2]))
  ela <- which.min(abs(g$lat[,1]-iniloc$lat[2]))
  L[elo, ela, T + 2] <- 1
  
  return(L)
  
}