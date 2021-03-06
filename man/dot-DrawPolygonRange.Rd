% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModAnnot.R
\name{.DrawPolygonRange}
\alias{.DrawPolygonRange}
\title{DrawPolygonRange Function (ModAnnot)}
\usage{
.DrawPolygonRange(nBreaksV, nTextMaxMin, cAxisName, cColor)
}
\arguments{
\item{nBreaksV}{Numeric vector giving the breaks to be used for the polygon.
Should be equal to seq(from = 1, to= c(number of categories), by=1).}

\item{nTextMaxMin}{Numeric vector giving the minimum and maximum values to be plotted on the polygon.
Should be c(min1, max1, min2, max2).}

\item{cAxisName}{The name of the axis to draw.}

\item{cColor}{The color of the polygon to draw.}
}
\description{
Return a polygon representing the range of Mod (or Base) categories of counts (or counts per kbp).
}
\keyword{internal}
