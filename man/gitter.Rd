% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{gitter}
\alias{gitter}
\title{Process a single plate image}
\usage{
gitter(
  image.file = file.choose(),
  plate.format = c(32, 48),
  remove.noise = F,
  autorotate = F,
  inverse = F,
  image.align = T,
  verbose = "l",
  contrast = NULL,
  fast = NULL,
  plot = F,
  grid.save = getwd(),
  dat.save = getwd(),
  start.coords = NULL,
  increment.coords = NULL,
  dilation.factor = 0,
  .fx = 2,
  .is.ref = F,
  .params = NULL
)
}
\arguments{
\item{image.file}{The path to the image. Defaults to a file choosing dialog.}

\item{plate.format}{The plate format, accepted formats: 1536, 768, 384 and 96. Alternatively, you can provide the number of rows and columns on the plate as an integer vector for example c(32,48). Default is 1536.}

\item{remove.noise}{Logical indicating noise/speckles should be remove from the thresholded image prior to analysis. Default is \code{FALSE}.}

\item{autorotate}{Logical indicating if image should be auto-rotated prior to processing. Only select this option if image is extremely rotated. gitter is able to handle small variations in rotations (1-2 degrees) without auto-rotating. Default is \code{FALSE}.}

\item{inverse}{Logical indicating if input image is inverted, meaning colonies are darker compared to their background. Default is \code{FALSE}.}

\item{image.align}{Used when reference image is provided. Logical indicating if images should be aligned to reference. This corrects small camera shifts.}

\item{verbose}{Shows details about the results of running job. For detailed logs "l", for a progress bar "p" or for no output "n". Default is "l".}

\item{contrast}{Integer between 1 and 100 indicating how much contrast should be applied to the image, prior to processing. A value of \code{NULL} will not apply any contrast. Default is \code{NULL}.}

\item{fast}{If set to integer value, the image will be resized to this width in pixels to speed up computation. This is useful for very large images that otherwise take a long time to process. We do not recommend resizing to fewer than 1500 pixels or greater that 4000 pixels in width. Default is \code{NULL}.}

\item{plot}{Logical indicating whether intensity profiles should be plotted. Default is \code{FALSE}.}

\item{grid.save}{Directory path to save gridded/thresholded images. Set to \code{NULL} if you do not want gridded images saved to disk. Default is the current working directory.}

\item{dat.save}{Directory path to save resulting data files. Set to \code{NULL} if you do not want resulting data saved to disk. Default is the current working directory.}

\item{start.coords}{Allows the user to specify a grid manually. This is a vector of length two with x,y coordinates of the top left most spot. Must be set in conjunction with increment.coords. By default grid is auto-detected.}

\item{increment.coords}{Allows the user to specify a grid manually. This is a vector of length two with distances between colonies on the x-axis and y-axis. Must be set in conjunction with start.coords. By default grid is auto-detected.}

\item{dilation.factor}{If specified, colonies will be dilated using a kernel of this size in pixels (value must be an odd number and >= 3). This dilation is only for fitting boundaries and will not be applied for counting colony pixels. This option is useful for spotted colonies, which are difficult to grid. By default, factor is 0 and no dialation is applied. See https://homepages.inf.ed.ac.uk/rbf/HIPR2/dilate.htm for more info.}

\item{.fx}{If the central pixel is a zero, set constant boundary with with size of this value}

\item{.is.ref}{Specifies if a reference property list is supplied. Warning: NOT for use by casual users.}

\item{.params}{Reference property list. Warning: NOT for use by casual users.}
}
\value{
\item{ DAT file }{
  Tab delimited file containing quantified colony sizes. There are two types of flags that can be associated with a data file 
  (1) plate-level flags signify possible misgridding of the plate due to a high number of colonies with small size or low circularity. #' These flags can be viewed using the \code{\link{plate.warnings}} function 
  (2) colony-based flags signify warnings associated with individual colonies. These flags can be viewed in the column named flags of #' the data file. 
  \tabular{ll}{
    row: \tab row number\cr
    col: \tab column number\cr
    size: \tab quantified colony size\cr
    circularity: \tab circularity of the colony\cr
    flags: \tab colony-based flags: S - Colony spill or edge interference, C- Low colony circularity\cr
  }
}

\item{ Gridded image }{Thresholded image showing the grid defined over the image}
}
\description{
The following function will grid and quantify a single plate image (for batch processing, see \code{\link{gitter.batch}})
}
\examples{
# Read sample image
f = system.file("extdata", "sample.jpg", package="gitter")
# Process it
dat = gitter(f)
# View head of the results
head(dat)
}
\keyword{gitter}
\keyword{image}
\keyword{process}
\keyword{sga}
\keyword{single}
