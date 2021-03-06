% Generated by roxygen2 (4.1.1.9000): do not edit by hand
% Please edit documentation in R/get-features-key-functions.R
\name{get.scale}
\alias{get.scale}
\title{Identify map scale}
\usage{
get.scale(site, l = 17)
}
\arguments{
\item{site}{Feature set list, containing a raster of all features, and a matrix assigning each feature to a particular type.}

\item{l}{Length of horizontal strip of pixels to search for. If the image contains a lot of long horizontal objects, this may need to be increased (default is 17).}
}
\value{
New feature set list containing a rescaled raster of numbered features and a matrix of corresponding feature types.
}
\description{
Identify the scale marker on the image. Will either output a list containing the details of the scale identified, or a new feature set object with the raster rescaled as directed by the user. If automatic detection fails, the legend can be selected automatically by the user. (Can be run more than once to rescale map repeatedly if necessary)
}
\details{
Uses a focal raster to find clusters containing large horizontal strips. Function will also rescale map to use life-size coordinate system, but doesn't have to; can be run again at any time to rescale further).
}
\examples{
get.scale(genlis)
# will output a new feature set object, with coordinates rescaled
}

