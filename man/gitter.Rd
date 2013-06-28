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
gitter(image.file, plate.format = c(32, 48), remove.noise = F, autorotate = F, 
       inverse = F, method = "kmeans", logging = T, 
       contrast = NA, fast = F, fast.width = 1500, plate.edges = T, plot = F, 
       gridded.save.dir = getwd(), dat.save.dir = getwd(), 
       .is.ref = F, .params = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image.file}{
  The path to the JPG/JPEG image
}
  \item{plate.format}{
  The plate format, which comprises of the number of rows and columns on the plate as an integer vector. Default is "c(32, 48)".
}
  \item{remove.noise}{
  Logical indicating noise/speckles should be remove from the image prior to quantification. Default is "FALSE". 
}
  \item{autorotate}{
  Logical indicating if image should be autorotated prior to processing. Only select this option if image is extremley rotated. gitter is able to handle small variations in rotations (1-3 degrees) without autorotating. Default is "FALSE".
}
  \item{inverse}{
  Logical indicating if input image is inverted, meaning colonies are darker compared to their background. Default is "FALSE".
}
  \item{method}{
  Method used for thesholding the image. If your image is noisy, has an uneven/cloudy background or contains really small colonies, use "tophat", otherwise use "kmeans". Default is "kmeans".
}
%  \item{smooth.factor}{
%  Integer indicating how much intensity profiles are smoothed. If your images are noisy, consider increasing this value as needed. Default is "5".
%}
  \item{logging}{
  Logical indicating whether gitter should show logs while processing images. Default is "TRUE".
}
  \item{contrast}{
  Integer indicating how much contrast should be applied to the image, prior to processing. Default is "NA".
}
  \item{fast}{
  Logical indicating if images should be resized to a certain width to speed up processing. See fast.width to specify a value. Default is "FALSE".
}
  \item{fast.width}{
  Integer indicating what the width of the resized image should be. This value is not used unless fast is "TRUE". Default is "1500".
}
  \item{plate.edges}{
  Logical indicating if plates contain edges. gitter will consider only the central region of the image when computing intensity profiles, to improve gridding. Default is "TRUE".
}
  \item{plot}{
 Logical indicating whether intensity profiles should be plotted. Default is "FALSE".
}
  \item{gridded.save.dir}{
 Directory path to save gridded/thresholded images. Set to "NA" if you do not want gridded images saved to disk. Default is the current working directory.
}
  \item{dat.save.dir}{
  Directory path to save dat files. Set to "NA" if you do not want dat files saved to disk. Default is the current working directory.
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
