\name{plot.roc.curve}
\alias{plot.roc.curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.roc.curve(stat.dataframe, printAUC = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stat.dataframe}{
%%     ~~Describe \code{stat.dataframe} here~~
}
  \item{printAUC}{
%%     ~~Describe \code{printAUC} here~~
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
function (stat.dataframe, printAUC = FALSE) 
{
    if (printAUC == TRUE) {
        specif.sense.df <- stat.dataframe[, c("specificity", 
            "sensitivity")]
        ordered.specif.sense.df <- specif.sense.df[with(specif.sense.df, 
            order(specificity)), ]
        ss.x <- ordered.specif.sense.df$specificity
        ss.y <- ordered.specif.sense.df$sensitivity
        AUC <- trapz(ss.x, ss.y)
    }
    p <- ggplot(stat.dataframe) + geom_line(aes(x = specificity, 
        y = sensitivity, color = Model.Name)) + geom_abline(intercept = 1, 
        slope = 1, color = "lightgrey") + scale_x_reverse(lim = c(1, 
        0)) + scale_y_continuous(lim = c(0, 1)) + coord_fixed(ratio = 1) + 
        labs(x = "Specificity", y = "Sensitivity") + theme_minimal()
    if (printAUC == TRUE) {
        p <- p + annotate("text", x = 0.25, y = 0.2, label = paste("AUC = ", 
            format(AUC, digits = 5), sep = ""))
    }
    p
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
