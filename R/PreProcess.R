PreProcess <-
function(mydata,spatial.smooth=TRUE, ...){
  Timenames<-dimnames(mydata)[[3]]
  nrows<-dim(mydata)[1];ncols<-dim(mydata)[2];TimeLen<-dim(mydata)[3];Time<-1:TimeLen
  mydata.norm <- array(0,c(nrows,ncols,TimeLen));preproc <- array(0,c(nrows,ncols,TimeLen))

  ## begin clustering
  ##
  X<-Time-mean(Time)           # centered Time
  coef.rq<-apply(mydata,c(1,2),function(Y) coefficients(rq(Y~X,tau=.5)))   
  ## Take log transformation to make variance of two clusters more comparible
  b0.log <- log(as.vector(coef.rq[1,,] - min(coef.rq[1,,]) +1))
  b1.log <- log(-(as.vector(coef.rq[2,,])) - min(-coef.rq[2,,]) + 0.01)
  ## determines the two clusters. The two centers are the initial points
  ## trained from the training set.

  clust2 <- kmeans(cbind(b0.log, b1.log), centers=2)
  group1 <- matrix(clust2$cluster==1, nrow=nrows)
  r1 <- range(mydata[group1])[2]-range(mydata[group1])[1]    # diff in "max" is faster
  r2 <- range(mydata[!group1])[2]-range(mydata[!group1])[1]
  if (r1<r2) { 
      nonmito.ind.array <- group1      # T for nonmito 
} else { 
      nonmito.ind.array <- !group1     # T for nonmito 
}

   ## normalization
   est.drifts <- outer(coef.rq[2,,],X,"*")    
#  est.drifts <- outer(coef.rq[2,,],X,"*") + outer(coef.rq[1,,],rep(1,length(X)),"*")
   mydata.norm <- mydata-est.drifts


  ##  smoothing
  if (spatial.smooth) {
      for (tt in 1:TimeLen){                       
          preproc[,,tt] <- image.smooth(mydata.norm[,,tt], ...)$z
      }
  } else {
      preproc <- mydata.norm
  }
  dimnames(preproc)[[3]]<-Timenames
  return(list(preproc=preproc,mask=nonmito.ind.array))
}

