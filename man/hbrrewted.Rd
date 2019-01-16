\name{hbrrewted}
\alias{hbrrewted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function which computes the weights for the iteratively reweighted LS estimate of the HBR fit.
%%  ~~function to do ... ~~
}
\description{
Function which computes the weights for the iteratively reweighted LS estimate of the HBR fit.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hbrrewted(x, y, percent = 0.95, ehat0 = ltsreg(x, y)$residuals)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{percent}{
%%     ~~Describe \code{percent} here~~
}
  \item{ehat0}{
%%     ~~Describe \code{ehat0} here~~
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
    hbrrewted = psi(abs(h))
    mat <- cbind(hbrrewted, ehat0, robdis2)
    return(mat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
