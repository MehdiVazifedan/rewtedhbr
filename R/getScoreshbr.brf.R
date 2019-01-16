getScoreshbr.brf <-
function(ehat,breaks,scores) {
#  breaks <- unique(breaks)
  ngb <- hist(ehat,br=breaks,plot=FALSE)
  scores1 <- getScores(scores,c(0,cumsum(ngb$counts)/sum(ngb$counts)))
  scoresvec <- (scores1[2:length(scores1)]+scores1[1:(length(scores1)-1)])/2
#  cuts <- as.factor(cut(ehat,ngb$breaks,labels=FALSE))
  cuts <- cut(ehat,ngb$breaks,labels=FALSE)
  lc <- sort(unique(cuts))
  scoresvec <- scoresvec[lc]
#  cutmat0 <- data.table(keys=lc,scores=scoresvec)
#  cutsmat <- data.table(keys=as.numeric(cuts))
  cutmat0 <- data.table(keys=lc,scores=scoresvec)
  cutsmat <- data.table(keys=cuts)
  scoremat <- merge(cutsmat,cutmat0,by='keys',all.x=TRUE,sort=FALSE)
  scorevec <- scoremat$scores - mean(scoremat$scores)
  list(scorevec=scorevec,mids=ngb$mids,counts=ngb$counts)
}
