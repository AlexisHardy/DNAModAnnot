% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DaLoad.R
\name{ImportPacBioGFF}
\alias{ImportPacBioGFF}
\title{ImportPacBioGFF Function (DaLoad)}
\usage{
ImportPacBioGFF(
  cPacBioGFFPath,
  cNameModToExtract,
  cModNameInOutput,
  cContigToBeAnalyzed = NULL,
  lKeepSequence = TRUE
)
}
\arguments{
\item{cPacBioGFFPath}{Path to a PacBio GFF file containing modification detection data.}

\item{cNameModToExtract}{Name of modification to be extracted.}

\item{cModNameInOutput}{Name for the extracted modification in the output.}

\item{cContigToBeAnalyzed}{Names of contigs for which the data will be kept.
If NULL, data from all contigs available will be imported. Defaults to NULL.}

\item{lKeepSequence}{If TRUE, the sequence of the base will be retained in one column.
Otherwise, it will be discarded to reduce object size. Defaults to TRUE.}
}
\description{
Import PacBio GFF file, extract one modification, rename this modification and convert it
as a GRanges object with new colnames similar to PacBio CSV file containing data from all bases sequenced.
}
\examples{
# Loading genome data
myGenome <- Biostrings::readDNAStringSet(system.file(
  package = "DNAModAnnot", "extdata",
  "ptetraurelia_mac_51_sca171819.fa"
))
names(myGenome)

# Loading PacBio data
myGrangesPacBioGFF <-
  ImportPacBioGFF(
    cPacBioGFFPath = system.file(
      package = "DNAModAnnot", "extdata",
      "ptetraurelia.modifications.sca171819.gff"
    ),
    cNameModToExtract = "m6A",
    cModNameInOutput = "6mA",
    cContigToBeAnalyzed = names(myGenome)
  )
myGrangesPacBioGFF

# Loading PacBio data for 2 scaffolds only
myGrangesPacBioGFF <-
  ImportPacBioGFF(
    cPacBioGFFPath = system.file(
      package = "DNAModAnnot", "extdata",
      "ptetraurelia.modifications.sca171819.gff"
    ),
    cNameModToExtract = "m6A",
    cModNameInOutput = "6mA",
    cContigToBeAnalyzed = c("scaffold51_18", "scaffold51_19")
  )
myGrangesPacBioGFF
}
\keyword{ImportPacBioGFF}
