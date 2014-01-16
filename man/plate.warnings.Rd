\name{plate.warnings}
\alias{plate.warnings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Show any plate-level warnings associated with a \code{\link{gitter}} data object
}
\description{
This function will show warnings associated with a \code{\link{gitter}} data object.
}
\usage{
plate.warnings(dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
    The data.frame produced by \code{\link{gitter}}.
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Warnings associated with the gitter data object or NULL if no warnings.
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
# plate.warnings(dat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ warnings }
\keyword{ plate }