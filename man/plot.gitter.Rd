\name{plot.gitter}
\alias{plot.gitter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a gitter dat file 
}
\description{
This function will plot a heatmap or bubble plot of a data.frame produced by \code{\link{gitter}} or a dat file saved to disk.
}
\usage{
\method{plot}{gitter}(x, title='', type='heatmap', low='turquoise', mid='black', high='yellow', 
            show.text=F, text.color='white', norm=T, show.flags=T, 
            flag.color='white', ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
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
  \item{show.flags}{
    Logical indicating if dots should be overlaid on the plot for flagged colonies. Default is "TRUE".
  }
  \item{flag.color}{
    Color of flag dot if show.flags is "TRUE". Default is "white".
  }
  \item{...}{
    Additional arguments. Not used.
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
# Read in path as a gitter data object
g = gitter.read(f)
# Plot a heatmap
plot(g, type="heatmap")
# Show a bubble plot
plot(g, type="bubble", low="black", high="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }