fitdiag6 <-
function (x, y, est = c("WIL", "HBR"), fitw=fitw,templs=templs,temphbr=temphbr,templts=templts,delta = 0.8, param = 2, 
    conf = 0.95) 
{
    x = as.matrix(centerx(x))
    n = dim(x)[1]
    p = dim(x)[2]
#   tempw = rfit(y ~ x)
#   residw = tempw$residuals
    vcw = vcov(fitw)
    tempw = fitw$coef
#   temphbr = NULL
#   templs = NULL
    if (any("WIL" == est) & any("HBR" == est)) {
#       temphbr = hbrfit(y ~ x)
        diff = tempw - temphbr$coef
        ydiff =  fitw$fitted.value - temphbr$fitted.values
    }
    if (any("WIL" == est) & any("LS" == est)) {
#       templs = lm(y ~ x)
        diff = tempw - templs$coef
        ydiff = fitw$fitted.value - templs$fitted.values
    }
    if (any("HBR" == est) & any("LS" == est)) {
#       temphbr = hbrfit(y ~ x)
#       templs = lm(y ~ x)
        diff = temphbr - templs
    }
    if (any("WIL" == est) & any("LTS" == est)) {
#       templts = ltsreg(x, y)
        diff = tempw - templts
    }
    if (any("HBR" == est) & any("LTS" == est)) {
#       temphbr = hbrfit(y ~ x)
#       templts = ltsreg(x, y)
        diff = temphbr - templts
    }
    if (any("LS" == est) & any("LTS" == est)) {
#       templs = lm(y ~ x)
#       templts = ltsreg(x, y)
        diff = templs - templts
    }
    vcwchol = chol(vcw)
   tdbeta = t(cbind(diff)) %*% chol2inv(vcwchol) %*% cbind(diff)
    bmtd = (4 * (p + 1)^2)/n
    xmat = cbind(rep(1, n), x)
    cfit = c()
    botf <- function(vec){sqrt(sum((vcwchol%*%vec)^2))}
    botres <- apply(xmat,1,botf)
    cfit <- ydiff/botres
    bmcf = 2 * sqrt((p + 1)/n)
#    se = sqrt(diag(vcw))
    list(tdbeta = tdbeta, bmtd = bmtd, cfit = cfit, bmcf = bmcf, 
        est = est)
}