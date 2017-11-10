\name{mattcc}
\alias{mattcc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mattcc(true_classes, predicted_classes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{true_classes}{
%%     ~~Describe \code{true_classes} here~~
}
  \item{predicted_classes}{
%%     ~~Describe \code{predicted_classes} here~~
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
function (true_classes, predicted_classes) 
{
    pos <- true_classes == 1
    neg <- true_classes == 0
    TP <- as.numeric(sum(predicted_classes[pos] == true_classes[pos]))
    TN <- as.numeric(sum(predicted_classes[neg] == true_classes[neg]))
    FP <- as.numeric(sum(predicted_classes[pos] != true_classes[pos]))
    FN <- as.numeric(sum(predicted_classes[neg] != true_classes[neg]))
    N1 <- TP * TN
    N2 <- FP * FN
    D1 <- TP + FP
    D2 <- TP + FN
    D3 <- TN + FP
    D4 <- TN + FN
    N <- N1 - N2
    D <- D1 * D2 * D3 * D4
    if (D == 0) {
        M <- 0
    }
    else {
        M <- N/sqrt(D)
    }
    return(M)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
