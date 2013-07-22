\name{gitter.plot}
\alias{gitter.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a gitter dat file 
}
\description{
This function will plot a heatmap or bubble plot of a data.frame produced by \code{\link{gitter}} or a dat file saved to disk.
}
\usage{
gitter.plot(dat, title='', type='heatmap', low='turquoise', mid='black', high='yellow', 
            show.text=F, text.color='white', norm=T, show.circ=T, circ.cutoff=0.6,
            circ.color='white')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
    The data.frame produced by \code{\link{gitter}} or the path to a dat file saved by \code{\link{gitter}}.
  }
  \item{title}{
    Title of plot. Default is blank.
  }
  \item{type}{
    Type of plot. "heatmap" for a heatmap, "bubble" for a bubble plot. Default is "heatmap".
  }
  \item{low}{
    Color for the lower bound of colony sizes. Default is "turquoise".
  }
  \item{mid}{
    Color of the middle value of colony sizes. Default is "black".
  }
  \item{high}{
    Color for the upper bound of colony sizes. Default is "yellow".
  }
  \item{show.text}{
    Logical indicating if text representation of colony sizes should be overlaid on the plot. Default is "TRUE".
  }
  \item{text.color}{
    Color of text if show.text is "TRUE". Default is "white".
  }
  \item{norm}{
    Logical indicating if colony sizes should be normalized by dividing colony sizes the middle mean of values and capping them between 0-2. Default is "TRUE".
  }
  \item{show.circ}{
    Logical indicating if dots should be overlaid on the plot for colonies with a circularity < X, where X is defined by circ.cutoff. Empty colonies are ignored. Default is "TRUE".
  }
  \item{circ.cutoff}{
    Colonies with a circularity of less than this value are displayed as overlaid dots if show.circ is "TRUE". Default is 0.6.
  }
  \item{circ.color}{
    Color of circularity dot if show.circ is "TRUE". Default is "white".
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a ggplot heatmap or bubble plot.
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
# Get dat file path
f = system.file("extdata", "sample.jpg.dat", package="gitter")
# Show a heatmap
gitter.plot(f, type='heatmap')
# Show a bubble plot
gitter.plot(f, type='bubble', low='black', high='red')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }