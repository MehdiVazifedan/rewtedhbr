hbrrewted2 <-
function (x, y, percent = 0.95, ehat0 = ltsreg(x, y)$residuals) 
{
    x = as.matrix(x)
    n = dim(x)[1]
    p = dim(x)[2]
    robdis2 <- robdist.hbrfit(x)
    cut = qchisq(percent, p)
    sigma = mad(ehat0)
    m = psi(cut/robdis2)
    a = ehat0/(sigma * m)
    c = (median(a) + 3 * mad(a))^2
    h = sqrt(c)/a
    hbrrewted2 = psi(abs(h))
    mat <- cbind(hbrrewted2,ehat0,robdis2)
    return(mat)
}
