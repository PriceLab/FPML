\name{make.pred.df.from.model}
\alias{make.pred.df.from.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.pred.df.from.model(gbdt.model, X.testdata, y.testdata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gbdt.model}{
%%     ~~Describe \code{gbdt.model} here~~
}
  \item{X.testdata}{
%%     ~~Describe \code{X.testdata} here~~
}
  \item{y.testdata}{
%%     ~~Describe \code{y.testdata} here~~
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
function (gbdt.model, X.testdata, y.testdata) 
{
    preds.test <- predict(gbdt.model, X.testdata, missing = NA)
    pred.df <- data.frame(ChIPseq.bound = y.testdata, Prediction = preds.test)
    pred.df$Model.Name = gbdt.model$Model.Name
    return(pred.df)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
