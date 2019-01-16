kstep <-
function(y,x,ehatk,k,wts2,scores=wscores){

    eps <- 0.000001
#  egt wts
    if(k == 1){
        mat <- hbrrewted(x,y)
        wts2 <- mat[,1]
        ehatk <- mat[,2]
        robdis2 <- mat[,3]
    }

     mstar <- median(ehatk)
     ehatk <- ehatk - mstar
     
     n <- length(ehatk)
#
#   This is not for big data
#

#   For final wts on kth iteration
     r <- rank(ehatk, ties.method = "first")/(length(ehatk) + 1)
     gs <- getScores(scores,r)

     ind <- c()
     wts <- c()

     for(i in 1:n){
          if(abs(ehatk[i]) < eps){
               wts[i] <- 0
               ind <- c(ind,i)
          } else {
               wts[i] <- gs[i]/ehatk[i]
          } 
      }
      wts <- wts/max(wts)
      wts[ind] <- 1
      pwts <- wts2*wts
      yk <- y - mstar

      fitkp1 <- lm(yk~x,weights=pwts)
      if(k == 1){
            list(fitkp1=fitkp1,wts2=wts2,pwts=pwts,robdis2=robdis2,ehatlts=ehatk)
      } else {
            list(fitkp1=fitkp1,wts2=wts2,pwts=pwts)
      }
}
