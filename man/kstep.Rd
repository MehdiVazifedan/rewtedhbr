\name{kstep}
\alias{kstep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function which computes the k steps for the weighted HBR fit.
%%  ~~function to do ... ~~
}
\description{
Function which computes the k steps for the weighted HBR fit.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
kstep(y, x, ehatk, k, wts2, scores = wscores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{ehatk}{
%%     ~~Describe \code{ehatk} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{wts2}{
%%     ~~Describe \code{wts2} here~~
}
  \item{scores}{
%%     ~~Describe \code{scores} here~~
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
function (y, x, ehatk, k, wts2, scores = wscores) 
{
    eps <- 1e-06
    if (k == 1) {
        mat <- hbrrewted(x, y)
        wts2 <- mat[, 1]
        ehatk <- mat[, 2]
        robdis2 <- mat[, 3]
    }
    mstar <- median(ehatk)
    ehatk <- ehatk - mstar
    n <- length(ehatk)
    r <- rank(ehatk, ties.method = "first")/(length(ehatk) + 
        1)
    gs <- getScores(scores, r)
    ind <- c()
    wts <- c()
    for (i in 1:n) {
        if (abs(ehatk[i]) < eps) {
            wts[i] <- 0
            ind <- c(ind, i)
        }
        else {
            wts[i] <- gs[i]/ehatk[i]
        }
    }
    wts <- wts/max(wts)
    wts[ind] <- 1
    pwts <- wts2 * wts
    yk <- y - mstar
    fitkp1 <- lm(yk ~ x, weights = pwts)
    if (k == 1) {
        list(fitkp1 = fitkp1, wts2 = wts2, pwts = pwts, robdis2 = robdis2, 
            ehatlts = ehatk)
    }
    else {
        list(fitkp1 = fitkp1, wts2 = wts2, pwts = pwts)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
