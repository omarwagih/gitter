\name{summary.gitter}
\alias{summary.gitter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
View the summary of a gitter data file
}
\description{
This function will show a brief summary of a data frame produced by \code{\link{gitter}}
}
\usage{
\method{summary}{gitter}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  The data.frame produced by \code{\link{gitter}}.
  }
  \item{...}{
    Additional arguments. Not used.
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Summary including the call made to \code{\link{gitter}}, colony size quantiles and more.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Omar Wagih \email{omar.wagih@utoronto.ca}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gitter}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# dat = gitter("/path/to/image")
# summary(dat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ summary }
