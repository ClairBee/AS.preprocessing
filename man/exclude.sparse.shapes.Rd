% Generated by roxygen2 (4.1.1.9000): do not edit by hand
% Please edit documentation in R/get-features-key-functions.R
\name{exclude.sparse.shapes}
\alias{exclude.sparse.shapes}
\title{Classify sparse shapes}
\usage{
exclude.sparse.shapes(site, density = 0.55, lower = 3, upper, plot = T)
}
\arguments{
\item{site}{Feature set list, containing a raster of all features, and a matrix assigning each feature to a particular type.}

\item{density}{A number from 0 to 1, setting the density threshold below which features are considered 'sparse'. The default is pi/8, based on the assumption that a post-hole is a filled circle, which will take up roughly pi/4 of the cells in a square drawn around it; the threshold is set at half this limit.}

\item{lower}{Number of cells to be considered as just noise. Any cluster of cells below this size is marked as an annotation. Defaults to 3.}

\item{upper}{(Optional) Number of cells above which a feature is considered to be a linear feature of the map, rather than an annotation. No linear features will be identified if this is omitted.
While we are only looking at post-hole positions, the distinction between annotations and linear features is not very important, and only really exists to speed up processing of functions that look at each annotation in turn.}

\item{plot}{Boolean indicator: should the new set of features be displayed on a plot?}
}
\value{
New feature set list containing the original raster of numbered features and a matrix of newly assigned feature types.
}
\description{
Identify shapes that are not compact by calculating the proportion of cells in a bounding square that are coloured, rather than white. Particularly large shapes are marked as linear features, while smaller ones are marked as annotations. Extremely small shapes are also marked as annotations, considering them to be noise.
}
\details{
Density cutoff is the approximate midpoint between the density of a solid circle (pi/8) and the density of the shortest possible straight line (being 3 pixels long, so having density 1/3)
}
\examples{
exclude.sparse.shapes(genlis)
genlis <- sparse.shapes.classified
}

