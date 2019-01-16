\name{get_breakshbr}
\alias{get_breakshbr}
\title{
Function that computes breaks for approximate Big data scores.
%%  ~~function to do ... ~~
}
\description{
Computes breaks for approximate Big data scores.
Written by John Kloke for the package bigRfit.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_breakshbr(ehat, B, eps = (.Machine$double.eps)^0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ehat}{
%%     ~~Describe \code{ehat} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
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
function (ehat, B, eps = (.Machine$double.eps)^0.5) 
{
    breaks <- quantile(ehat, seq(0, 1, length = B))
    ind <- c(1, length(breaks))
    breaks[ind] <- breaks[ind] + eps * sd(breaks) * c(-1, 1)
    unique(breaks)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
