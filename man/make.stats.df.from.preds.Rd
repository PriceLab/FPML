\name{make.stats.df.from.preds}
\alias{make.stats.df.from.preds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.stats.df.from.preds(pred.df, num.roc.points = 101)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pred.df}{
%%     ~~Describe \code{pred.df} here~~
}
  \item{num.roc.points}{
%%     ~~Describe \code{num.roc.points} here~~
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
function (pred.df, num.roc.points = 101) 
{
    preds <- pred.df$Prediction
    cs.vals <- pred.df$ChIPseq.bound
    roc.tf <- roc(cs.vals, preds, ret = c("threshold", "sens", 
        "spec", "ppv", "npv", "acc"))
    stats.tf <- coords((roc.tf), seq(from = 0, to = 1, length.out = num.roc.points), 
        ret = c("threshold", "sens", "spec", "ppv", "npv", "acc"))
    mattcc.tf <- mattcc.curve(cs.vals, preds, num.points = num.roc.points)
    df.stats <- as.data.frame(t(stats.tf))
    df.stats$MattCC <- mattcc.tf
    df.stats$Model.Name <- pred.df$Model.Name[1]
    return(df.stats)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
