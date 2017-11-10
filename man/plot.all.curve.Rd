\name{plot.all.curve}
\alias{plot.all.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.all.curve(stat.dataframe)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stat.dataframe}{
%%     ~~Describe \code{stat.dataframe} here~~
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
function (stat.dataframe) 
{
    stat.dataframe.tidy <- gather(stat.dataframe, "Statistic", 
        "Value", -threshold, -Model.Name)
    ggplot(stat.dataframe.tidy) + geom_line(aes(x = threshold, 
        y = Value, color = Statistic)) + scale_x_continuous(lim = c(0, 
        1)) + scale_y_continuous(lim = c(0, 1)) + coord_fixed(ratio = 1) + 
        labs(x = "Threshold") + theme_minimal()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
