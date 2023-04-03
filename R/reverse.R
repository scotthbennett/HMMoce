#Function for manipulating arrays and viewing matrix as image
#Needed for code in hmm.posteriorviterbi.R and hmm.filterExp.R

reverse <- function(m) t(m)[,nrow(m):1] #Reverse lats for using image command
