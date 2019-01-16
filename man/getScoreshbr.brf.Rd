\name{getScoreshbr.brf}
\alias{getScoreshbr.brf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that computes approximate Big data scores.
%%  ~~function to do ... ~~
}
\description{
Computes approximate (Big data) scores for Big data.
Written by John Kloke for the package bigRfit.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getScoreshbr.brf(ehat, breaks, scores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ehat}{
%%     ~~Describe \code{ehat} here~~
}
  \item{breaks}{
%%     ~~Describe \code{breaks} here~~
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
function (ehat, breaks, scores) 
{
    ngb <- hist(ehat, br = breaks, plot = FALSE)
    scores1 <- getScores(scores, c(0, cumsum(ngb$counts)/sum(ngb$counts)))
    scoresvec <- (scores1[2:length(scores1)] + scores1[1:(length(scores1) - 
        1)])/2
    cuts <- cut(ehat, ngb$breaks, labels = FALSE)
    lc <- sort(unique(cuts))
    scoresvec <- scoresvec[lc]
    cutmat0 <- data.table(keys = lc, scores = scoresvec)
    cutsmat <- data.table(keys = cuts)
    scoremat <- merge(cutsmat, cutmat0, by = "keys", all.x = TRUE, 
        sort = FALSE)
    scorevec <- scoremat$scores - mean(scoremat$scores)
    list(scorevec = scorevec, mids = ngb$mids, counts = ngb$counts)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
