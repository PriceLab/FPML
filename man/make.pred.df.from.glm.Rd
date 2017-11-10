\name{make.pred.df.from.glm}
\alias{make.pred.df.from.glm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make.pred.df.from.glm(glm.model, df.testdata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{glm.model}{
%%     ~~Describe \code{glm.model} here~~
}
  \item{df.testdata}{
%%     ~~Describe \code{df.testdata} here~~
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
function (glm.model, df.testdata) 
{
    preds.test <- predict.glm(object = glm.model, newdata = df.testdata, 
        type = "response")
    pred.df <- data.frame(ChIPseq.bound = df.testdata$ChIPseq.bound, 
        Prediction = preds.test)
    pred.df$Model.Name = glm.model$Model.Name
    return(pred.df)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
