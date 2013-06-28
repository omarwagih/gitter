\name{gitter.batch}
\alias{gitter.batch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Process a batch set of plate images
}
\description{
This function will process a directory or list of images in a batch. You can also use this function to process images with sparse to empty rows/columns using a reference image.
}
\usage{
gitter.batch(image.files, ref.image.file = NA, logging = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image.files}{
  Directory containing JPG or JPEG images or a character vector of image paths.
}
  \item{ref.image.file}{
  Specifies path to a reference image, if you are using a reference image to grid other images.
}
  \item{logging}{
  See parameters in \code{\link{gitter}}.
}
  \item{\dots}{
  Additional parameters passed to \code{\link{gitter}}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
gitter.batch does not return any values. DAT and gridded files are saved to their respective directories.
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
\code{\link{gitter}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Processing image using reference image

# This image would typically fail to process, since its missing several rows
f = system.file("extdata", "sample_dead.jpg", package="gitter")
# We will use this image to successfully process the above image
f.ref = system.file("extdata", "sample.jpg", package="gitter")
# Process
gitter.batch(f, f.ref)

# Remember: output files by default are saved to your working directory
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gitter.batch }
