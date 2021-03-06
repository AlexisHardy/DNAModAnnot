% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FDREst.R
\name{GetFdrEstListByThresh}
\alias{GetFdrEstListByThresh}
\title{GetFdrEstListByThresh Function (FDREst)}
\usage{
GetFdrEstListByThresh(
  grangesDataWithSeq,
  grangesDataWithSeqControl = NULL,
  cNameParamToTest,
  nRoundDigits = 1,
  cModMotifsAsForeground = NULL
)
}
\arguments{
\item{grangesDataWithSeq}{A GRanges-like object containing, in the extra columns, the parameter to be tested 
and the sequence of the associated window.}

\item{grangesDataWithSeqControl}{A GRanges-like object to be used as a control containing, in the extra columns, the parameter to be tested
and the sequence of the associated window. This control is usually a non-methylated sample (example: Whole-Genome Amplified/PCR Amplified).
If not NULL, false discovery rate estimation will be calculated using the associated parameter from grangesDataWithSeq and 
grangesDataWithSeqControl. For example, by defining "b" as any adenine and defining "param" as the parameter to be tested, and for each threshold:
\itemize{
  \item foreground = proportion, in the sample, of modified "b" among total "b" = (number of "b" with "param" >= threshold) / total number of "b" in the grangesDataWithSeq
  \item background = proportion, in the control, of modified "b" among total "b" = (number of "b" with "param" >= threshold) / total number of "b" in the grangesDataWithSeqControl
  \item fdr estimation: background / foreground
}
If NULL, only the parameter from grangesDataWithSeq will be used:
here the background and foreground will be estimated using motifs associated to modifications (provided with cModMotifsAsForeground) versus other motifs.
For example, by defining "m" as one motif associated to modifications, defining "b" as any adenine and defining "param" as the parameter to be tested,
for each "m" (provided with cModMotifsAsForeground) and for each threshold:
\itemize{
  \item foreground = proportion of modified "b",     in "m", among total "b" = (number of "b",     in "m", with "param" >= threshold) / total number of "b"     in "m"
  \item background = proportion of modified "b", not in "m", among total "b" = (number of "b", not in "m", with "param" >= threshold) / total number of "b" not in "m"
  \item fdr estimation: background / foreground
}
Defaults to NULL.}

\item{cNameParamToTest}{The name of the column containing the parameter to be tested.}

\item{nRoundDigits}{The number of digits for the thresholds on the parameter. A value of 0 would mean no value to be rounded:
this could results in errors for a continuous variable. Defaults to 1.}

\item{cModMotifsAsForeground}{A character vector of motifs associated to the modifications.
If grangesDataWithSeqControl is NULL, this vector must be given. This vector is not used if grangesDataWithSeqControl is not NULL. Defaults to NULL.}
}
\value{
A list composed of x dataframes: if grangesDataWithSeqControl is NULL, 1 dataframe; if grangesDataWithSeqControl is not NULL,
1 dataframe by motif provided with cModMotifsAsForeground. In each dataframe:
\itemize{
  \item fdr: The false discovery rate estimated for this threshold.
  \item threshold: The threshold on the parameter.
  \item fdr_cummin: The minimum false discovery rate estimated for this threshold and less stringent thresholds (adjusted false discovery rate).
}
}
\description{
Return a list with the false discovery rate estimated by threshold on the parameter provided with the GRanges object(s).
The thresholds to test are determined by all the possible values of the parameter provided.
}
\examples{
library(Biostrings)
myGenome <- Biostrings::readDNAStringSet(system.file(
  package = "DNAModAnnot", "extdata",
  "ptetraurelia_mac_51_sca171819.fa"
))
myGrangesGenome <- GetGenomeGRanges(myGenome)

# Preparing a gposPacBioCSV dataset with sequences
myGposPacBioCSV <-
  ImportPacBioCSV(
    cPacBioCSVPath = system.file(
      package = "DNAModAnnot", "extdata",
      "ptetraurelia.bases.sca171819.csv"
    ),
    cSelectColumnsToExtract = c(
      "refName", "tpl", "strand", "base", "score",
      "ipdRatio", "coverage"
    ),
    lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
    cContigToBeAnalyzed = names(myGenome)
  )
myGrangesBaseCSV <- as(myGposPacBioCSV[myGposPacBioCSV$base == "A"], "GRanges")
myGrangesBaseCSVWithSeq <- GetGRangesWindowSeqandParam(
  grangesData = myGrangesBaseCSV,
  grangesGenome = myGrangesGenome,
  dnastringsetGenome = myGenome,
  nUpstreamBpToAdd = 0,
  nDownstreamBpToAdd = 1
)

# FDR estimation by motif associated to modifications
myFdr_score_per_motif_list <-
  GetFdrEstListByThresh(
    grangesDataWithSeq = myGrangesBaseCSVWithSeq,
    cNameParamToTest = "score",
    nRoundDigits = 1,
    cModMotifsAsForeground = c("AG", "AT")
  )
myFdr_score_per_motif_list

## NOT RUN!
## FDR estimation versus granges control
# myFdr_score_vsCtrl_list <-
#  GetFdrEstListByThresh(grangesDataWithSeq = myGrangesBaseCSVWithSeq,
#                        grangesDataWithSeqControl = myGrangesBaseCSVWithSeq_control,
#                        cNameParamToTest = "score",
#                        nRoundDigits = 1)
# myFdr_score_vsCtrl_list
}
\keyword{GetFdrEstListByThresh}
