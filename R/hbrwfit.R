hbrwfit <-
function(xmat,y,numstp=10,epcchk=0.0001,scores=wscores){
##   have y and xmat
     x <- centerx(xmat)

     pchk <- length(x[1,])
     plim <- (pchk+1)*(numstp+1)
     matbeta <- matrix(rep(0,plim),ncol=(numstp+1))


     matwts <- c()
     ehatk  <- c()
     wts2  <- c()
     chks <- 1
     k <- 1

     while(chks > 0){
    
          if(k == 1){
                firstfit <- kstep(y,x,ehatk,1,wts2,scores)
                robdis2 <- firstfit$robdis2
                ehatlts <- firstfit$ehatlts
                betak <- firstfit$fitkp1$coef
                wts2 <- firstfit$wts2
                matwts <- cbind(matwts,firstfit$pwts)
                matbeta[,k] <-betak
                ehatk <- firstfit$fitkp1$resid
                k <- k + 1
          } else {
                kp1fit <- kstep(y,x,ehatk,k,wts2,scores)
                betakp1 <- kp1fit$fitkp1$coef
                matwts <- cbind(matwts,kp1fit$pwts)
                matbeta[,k] <- betakp1
                delta <- betakp1 - betak
                chkdel <- sum(delta^2)/(epcchk + sum(betak^2))
                if(chkdel < epcchk){
                     chks <- -1
                } else {
                     k <- k+1
                     ehatk <- kp1fit$fitkp1$resid
                     betak <- betakp1
                }
          }
          if(k > numstp){
               chks <- -1
               k <- k - 1
          }
      }
                
      findim <- dim(matwts)[2]
      finwts<-cbind(wts2,matwts[,findim])

      if(k == 1){
            coef <- matbeta[,1]
           } else {
             coef <- matbeta[,k]
      }

      fitted.values <- kp1fit$fitkp1$fitted.values
      residuals <- kp1fit$fitkp1$residuals
      
      list(coef=coef,fitted.values=fitted.values,residuals=residuals,kp1fit=kp1fit,matbeta=matbeta,finwts=finwts,robdis2=robdis2,ehatlts=ehatlts)
#       list(coef=coef,matbeta=matbeta)
}
