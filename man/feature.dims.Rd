% Generated by roxygen2 (4.1.1.9000): do not edit by hand
% Please edit documentation in R/get-features-additional-functions.R
\name{feature.dims}
\alias{feature.dims}
\title{Get feature dimensions}
\usage{
feature.dims(site)
}
\arguments{
\item{site}{Feature set list, containing a raster of all features, and a matrix assigning each feature to a particular type.}
}
\value{
A matrix containing the dimensions of every feature identified in the map.
}
\description{
Creates a matrix containing the dimensions of every feature in the map: feature id, midpoint, max and min x and y coordinates; width, height and size in terms of cells covered; and ratio of longest to shortest edge.
Can be used to investigate plausible thresholds for the other functions.
}
\examples{
feature.dims(genlis)
}

