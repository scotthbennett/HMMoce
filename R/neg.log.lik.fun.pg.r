neg.log.lik.fun.pg <- function(pars,g,L,maskL=TRUE,bound.thr = 0.1, minBounds = 10){
  # pars=c(3,1,.7,.8)
  sigmas=pars[1:2]
  sizes=ceiling(sigmas*4)
  pb=pars[3:4]
  muadvs=c(0,0)
  
  # if(!is.na(pars[5])) muadvs[1]<-pars[5]
  #  if(!is.na(pars[6])) muadvs[2]<-pars[6]
  
  # behav 1
  # ensure that there is one center cell
  if(sizes[1]%%2==0){sizes[1]=sizes[1]+1}
  ss=sum(gausskern.nostd(sizes[1],sigmas[1],muadv=muadvs[1]))
  while(ss<.999){
    sizes[1]=sizes[1]+2
    ss=sum(gausskern.nostd(sizes[1],sigmas[1],muadv=muadvs[1]))
  }
  K1 = gausskern.pg(sizes[1],sigmas[1],muadv=muadvs[1])
  rm(ss)
  #  K1=mask.K(K1)
  
  # behav 2
  # if(sizes[2]%%2==0){sizes[2]=sizes[2]+1}
  ss=sum(gausskern.nostd(sizes[2],sigmas[2],muadv=muadvs[2]))
  while(ss<.999){
    sizes[2]=sizes[2]+2
    ss=sum(gausskern.nostd(sizes[2],sigmas[2],muadv=muadvs[2]))
  }
  K2 = gausskern.pg(sizes[2],sigmas[2],muadv=muadvs[2])
  rm(ss)
  #  K2=mask.K(K2)
  
  P <- matrix(c(pb[1],1-pb[1],1-pb[2],pb[2]),2,2,byrow=TRUE)
  #f. <- hmm.filter.2B.dealNA(g, L, K1,K2, P)
  f. <- hmm.filter(g, L, K1,K2, P,maskL=maskL,bound.thr = bound.thr, minBounds = minBounds)
  nllf. <- -sum(log(f.$psi))
  return(nllf.)
}