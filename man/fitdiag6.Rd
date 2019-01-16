\name{fitdiag6}
\alias{fitdiag6}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fucntion which computes difference between wted HBR and Wilcoxon fits.
%%  ~~function to do ... ~~
}
\description{
Fucntion which computes difference between wted HBR and Wilcoxon fits.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fitdiag6(x, y, est = c("WIL", "HBR"), fitw = fitw, templs = templs, temphbr = temphbr, templts = templts, delta = 0.8, param = 2, conf = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{est}{
%%     ~~Describe \code{est} here~~
}
  \item{fitw}{
%%     ~~Describe \code{fitw} here~~
}
  \item{templs}{
%%     ~~Describe \code{templs} here~~
}
  \item{temphbr}{
%%     ~~Describe \code{temphbr} here~~
}
  \item{templts}{
%%     ~~Describe \code{templts} here~~
}
  \item{delta}{
%%     ~~Describe \code{delta} here~~
}
  \item{param}{
%%     ~~Describe \code{param} here~~
}
  \item{conf}{
%%     ~~Describe \code{conf} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, est = c("WIL", "HBR"), fitw = fitw, templs = templs, 
    temphbr = temphbr, templts = templts, delta = 0.8, param = 2, 
    conf = 0.95) 
{
    x = as.matrix(centerx(x))
    n = dim(x)[1]
    p = dim(x)[2]
    vcw = vcov(fitw)
    tempw = fitw$coef
    if (any("WIL" == est) & any("HBR" == est)) {
        diff = tempw - temphbr$coef
        ydiff = fitw$fitted.value - temphbr$fitted.values
    }
    if (any("WIL" == est) & any("LS" == est)) {
        diff = tempw - templs$coef
        ydiff = fitw$fitted.value - templs$fitted.values
    }
    if (any("HBR" == est) & any("LS" == est)) {
        diff = temphbr - templs
    }
    if (any("WIL" == est) & any("LTS" == est)) {
        diff = tempw - templts
    }
    if (any("HBR" == est) & any("LTS" == est)) {
        diff = temphbr - templts
    }
    if (any("LS" == est) & any("LTS" == est)) {
        diff = templs - templts
    }
    vcwchol = chol(vcw)
    tdbeta = t(cbind(diff)) \%*\% chol2inv(vcwchol) \%*\% cbind(diff)
    bmtd = (4 * (p + 1)^2)/n
    xmat = cbind(rep(1, n), x)
    cfit = c()
    botf <- function(vec) {
        sqrt(sum((vcwchol \%*\% vec)^2))
    }
    botres <- apply(xmat, 1, botf)
    cfit <- ydiff/botres
    bmcf = 2 * sqrt((p + 1)/n)
    list(tdbeta = tdbeta, bmtd = bmtd, cfit = cfit, bmcf = bmcf, 
        est = est)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
