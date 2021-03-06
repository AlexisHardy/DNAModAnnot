% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModAnnot.R
\name{ImportBamExtendedAnnotationTrack}
\alias{ImportBamExtendedAnnotationTrack}
\title{ImportBamExtendedAnnotationTrack Function (ModAnnot)}
\usage{
ImportBamExtendedAnnotationTrack(file, selection)
}
\description{
Subfunction to import a bam file for displaying using Gviz package
while retrieving the content of the "Xa" optional field (will be retrieved as "tag" field).
(See Gviz package documentation to see how to use it while making tracks)
After making the AnnotationTrack,
in order to allow the names of the genomic features to be displayed, the "mapping" sub-list of
the generated annotation track must be completed with the new "id" and "group" values chosen while making the annotationTrack.
Example: trackAnnotation@mapping <- list(id="tag", group="tag")
}
\keyword{internal}
