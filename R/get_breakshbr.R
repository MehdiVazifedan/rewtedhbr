get_breakshbr <-
function(ehat,B,eps=(.Machine$double.eps)^0.5) {
  breaks <- quantile(ehat,seq(0,1,length=B))
  ind <- c(1,length(breaks))
  breaks[ind] <- breaks[ind] + eps*sd(breaks)*c(-1,1)
  unique(breaks)
}
