\name{gitter}
\alias{gitter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Process a single plate image 
}
\description{
The following function will grid and quantify a single plate image (for batch processing, see \code{\link{gitter.batch}})
}
\usage{
gitter(image.file = file.choose(), plate.format = c(32, 48),
       remove.noise = F, autorotate = F, inverse = F, verbose = "l",
       contrast = NULL, fast = NULL, plot = F, grid.save = getwd(), 
       dat.save = getwd(), .is.ref = F, .params = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image.file}{
  The path to the JPG/JPEG image. Defaults to a file choosing dialog.
}
  \item{plate.format}{
  The plate format, accepted formats: 1536, 768, 384 and 96. Alternatively, you can provide the number of rows and columns on the plate as an integer vector for example "c(32,48)". Default is "1536".
}
  \item{remove.noise}{
  Logical indicating noise/speckles should be remove from the thresholded image prior to analysis. Default is "FALSE". 
}
  \item{autorotate}{
  Logical indicating if image should be auto-rotated prior to processing. Only select this option if image is extremely rotated. gitter is able to handle small variations in rotations (1-2 degrees) without auto-rotating. Default is "FALSE".
}
  \item{inverse}{
  Logical indicating if input image is inverted, meaning colonies are darker compared to their background. Default is "FALSE".
}
  \item{verbose}{
  Shows details about the results of running job. For detailed logs "l", for a progress bar "p" or for no output "n". Default is "l".
}
  \item{contrast}{
  Integer indicating how much contrast should be applied to the image, prior to processing. A value of "NULL" will not apply any contrast. Default is "NULL".
}
  \item{fast}{
  Pixel value between 1500-4000 indicating if images should be resized to a certain width to speed up processing. If "NULL", no resizing occurs. Default is "NULL".
}
  \item{plot}{
 Logical indicating whether intensity profiles should be plotted. Default is "FALSE".
}
  \item{grid.save}{
 Directory path to save gridded/thresholded images. Set to "NULL" if you do not want gridded images saved to disk. Default is the current working directory.
}
  \item{dat.save}{
  Directory path to save dat files. Set to "NULL" if you do not want dat files saved to disk. Default is the current working directory.
}

\item{.is.ref}{
  Specifies if a reference property list is supplied. Warning: NOT for use by casual users.
}
\item{.params}{
  Reference property list. Warning: NOT for use by casual users.
}

}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
\item{ DAT file }{
  Tab delimited file containing quantified colony sizes. 
  \tabular{ll}{
    row: \tab row number\cr
    col: \tab column number\cr
    size: \tab quantified colony size\cr
    circularity: \tab circularity of the colony\cr
    %x: \tab x-coordinate value of the colony centre\cr
    %y: \tab y-coordinate value of the colony centre\cr
    %xl: \tab coordinate defining left bound of the colony (x-left)\cr
    %xr: \tab coordinate defining right bound of the colony (x-right)\cr
    %yt: \tab coordinate defining top bound of the colony (y-top)\cr
    %yb: \tab coordinate defining bottom bound of the colony (y-bottom)\cr
  }
}

\item{ Gridded image }{Thresholded image showing the grid defined over the image}


%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
\code{\link{gitter.batch}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Read sample image
f = system.file("extdata", "sample.jpg", package="gitter")
# Process it
dat = gitter(f)
# View head of the results
head(dat)

# Remember: output files by default are saved to your working directory
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gitter }
