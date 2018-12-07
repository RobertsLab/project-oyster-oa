nmds.montejg <-
  function(dissim,rawdata,k,nmds.distance='bray',daisy.distance='gower', trymax=50,autotransform=FALSE,
           trace=0,zerodist='add',perm=99,col.hist='blue',col.line='red',
           lty=2,las=1,lab=c(5,5,4),...){
    
    library(vegan)
    library(MASS)
    
    z<-metaMDS(comm=dissim,k=k,distance=nmds.distance,trymax=trymax,
               autotransform=autotransform,trace=trace,...) #nmds analysis on dissimilarity matrix
    z.stress<-z$stress #get stress
    y.stress<-rep(0,perm)
    
    for(i in 1:perm){
      y<-apply(rawdata,2,sample) #permute raw data matrix
      y<- daisy(y, metric = daisy.distance) # calculate gower's for permuted data
      y<-metaMDS(comm=y,k=k,distance=nmds.distance,trymax=trymax,
                 autotransform=autotransform,trace=trace,...) #nmds analysis on dissimilarity matrix
      y.stress[i]<-y$stress #get stress
    }
    n<-sum(y.stress<=z.stress) #compute number of random runs with stress < observed
    p.value<-(1+n)/(1+perm) #compute p-value
    
    xmin<-min(z.stress,min(y.stress))
    xmax<-max(z.stress,max(y.stress))
    hist(y.stress,col=col.hist,las=las,lab=lab,
         xaxs='i',yaxs='i',xlim=c(xmin,xmax),xlab='Stress',
         main=paste('Random Permutation Distribution of Stress for',k,'Dimensions',sep=' '),...)
    abline(v=z.stress,col=col.line,lty=lty,lwd=2,...)
    
    cat('Randomization Test of Stress:\n')
    cat('Permutation stress values:\n')
    print(y.stress)
    z<-rbind('Observed stress'=z.stress,'P-value'=p.value,'Times Permute Less Than Observed'=n)
    return(z)
  }