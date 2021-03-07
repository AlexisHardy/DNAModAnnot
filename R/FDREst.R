### =========================================================================
### Functions from FDREst category
### -------------------------------------------------------------------------
###
### DNAModAnnot v.0.0.0.9014 - 2020/09/10
### Licence GPL-3
###
### Contributor:
### Alexis Hardy
### - email: "alexis.hardy1994@outlook.fr"
### - Sandra Duharcourt Laboratory
### - Matthieu Defrance Laboratory
###
### -------------------------------------------------------------------------
###
### Functions:
###
### GetFdrEstListByThresh
### .GetCumsumTableParamByMotif
### .EstimateFdr
### DrawFdrEstList
### .DrawFdrEst
### .DrawFdrSegmentsXYAxis
### GetFdrBasedThreshLimit
###
### -------------------------------------------------------------------------

#' GetFdrEstListByThresh Function (FDREst)
#'
#' Return a list with the false discovery rate estimated by threshold on the parameter provided with the GRanges object(s).
#' The thresholds to test are determined by all the possible values of the parameter provided.
#' @param grangesDataWithSeq A GRanges-like object containing, in the extra columns, the parameter to be tested 
#' and the sequence of the associated window.
#' @param grangesDataWithSeqControl A GRanges-like object to be used as a control containing, in the extra columns, the parameter to be tested  
#' and the sequence of the associated window. This control is usually a non-methylated sample (example: Whole-Genome Amplified/PCR Amplified).
#' If not NULL, false discovery rate estimation will be calculated using the associated parameter from grangesDataWithSeq and 
#' grangesDataWithSeqControl. For example, by defining "b" as any adenine and defining "param" as the parameter to be tested, and for each threshold:
#' \itemize{
#'   \item foreground = proportion, in the sample, of modified "b" among total "b" = (number of "b" with "param" >= threshold) / total number of "b" in the grangesDataWithSeq
#'   \item background = proportion, in the control, of modified "b" among total "b" = (number of "b" with "param" >= threshold) / total number of "b" in the grangesDataWithSeqControl
#'   \item fdr estimation: background / foreground
#' }
#' If NULL, only the parameter from grangesDataWithSeq will be used:
#' here the background and foreground will be estimated using motifs associated to modifications (provided with cModMotifsAsForeground) versus other motifs.
#' For example, by defining "m" as one motif associated to modifications, defining "b" as any adenine and defining "param" as the parameter to be tested,
#' for each "m" (provided with cModMotifsAsForeground) and for each threshold:
#' \itemize{
#'   \item foreground = proportion of modified "b",     in "m", among total "b" = (number of "b",     in "m", with "param" >= threshold) / total number of "b"     in "m"
#'   \item background = proportion of modified "b", not in "m", among total "b" = (number of "b", not in "m", with "param" >= threshold) / total number of "b" not in "m"
#'   \item fdr estimation: background / foreground
#' }
#' Defaults to NULL.
#' @param cNameParamToTest The name of the column containing the parameter to be tested.
#' @param nRoundDigits The number of digits for the thresholds on the parameter. A value of 0 would mean no value to be rounded:
#' this could results in errors for a continuous variable. Defaults to 1.
#' @param cModMotifsAsForeground A character vector of motifs associated to the modifications.
#' If grangesDataWithSeqControl is NULL, this vector must be given. This vector is not used if grangesDataWithSeqControl is not NULL. Defaults to NULL.
#' @return A list composed of x dataframes: if grangesDataWithSeqControl is NULL, 1 dataframe; if grangesDataWithSeqControl is not NULL,
#' 1 dataframe by motif provided with cModMotifsAsForeground. In each dataframe:
#' \itemize{
#'   \item fdr: The false discovery rate estimated for this threshold.
#'   \item threshold: The threshold on the parameter.
#'   \item fdr_cummin: The minimum false discovery rate estimated for this threshold and less stringent thresholds (adjusted false discovery rate).
#' }
#' @keywords GetFdrEstListByThresh
#' @export
#' @examples
#' library(Biostrings)
#' myGenome <- Biostrings::readDNAStringSet(system.file(
#'   package = "DNAModAnnot", "extdata",
#'   "ptetraurelia_mac_51_sca171819.fa"
#' ))
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' # Preparing a gposPacBioCSV dataset with sequences
#' myGposPacBioCSV <-
#'   ImportPacBioCSV(
#'     cPacBioCSVPath = system.file(
#'       package = "DNAModAnnot", "extdata",
#'       "ptetraurelia.bases.sca171819.csv"
#'     ),
#'     cSelectColumnsToExtract = c(
#'       "refName", "tpl", "strand", "base", "score",
#'       "ipdRatio", "coverage"
#'     ),
#'     lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'     cContigToBeAnalyzed = names(myGenome)
#'   )
#' myGrangesBaseCSV <- as(myGposPacBioCSV[myGposPacBioCSV$base == "A"], "GRanges")
#' myGrangesBaseCSVWithSeq <- GetGRangesWindowSeqandParam(
#'   grangesData = myGrangesBaseCSV,
#'   grangesGenome = myGrangesGenome,
#'   dnastringsetGenome = myGenome,
#'   nUpstreamBpToAdd = 0,
#'   nDownstreamBpToAdd = 1
#' )
#'
#' # FDR estimation by motif associated to modifications
#' myFdr_score_per_motif_list <-
#'   GetFdrEstListByThresh(
#'     grangesDataWithSeq = myGrangesBaseCSVWithSeq,
#'     cNameParamToTest = "score",
#'     nRoundDigits = 1,
#'     cModMotifsAsForeground = c("AG", "AT")
#'   )
#' myFdr_score_per_motif_list
#'
#' ## NOT RUN!
#' ## FDR estimation versus granges control
#' # myFdr_score_vsCtrl_list <-
#' #  GetFdrEstListByThresh(grangesDataWithSeq = myGrangesBaseCSVWithSeq,
#' #                        grangesDataWithSeqControl = myGrangesBaseCSVWithSeq_control,
#' #                        cNameParamToTest = "score",
#' #                        nRoundDigits = 1)
#' # myFdr_score_vsCtrl_list
GetFdrEstListByThresh <- function(grangesDataWithSeq,
                                  grangesDataWithSeqControl = NULL,
                                  cNameParamToTest,
                                  nRoundDigits = 1,
                                  cModMotifsAsForeground = NULL) {
  lTestBetweenMotifs <- is.null(grangesDataWithSeqControl)
  ct_seq_param <- .GetCumsumTableParamByMotif(
    grangesDataWithSeq = grangesDataWithSeq,
    cNameParamToTest = cNameParamToTest, nRoundDigits = nRoundDigits,
    lTestBetweenMotifs = lTestBetweenMotifs
  )
  if (lTestBetweenMotifs) {
    if (is.null(cModMotifsAsForeground)) {
      print("Error: if grangesDataWithSeqControl is NULL, cModMotifsAsForeground must be given.")
      listFdrEstByThr <- NULL
    } else {
      listFdrEstByThr <- lapply(cModMotifsAsForeground, function(motif_fg) {
        nBg <- rowSums(ct_seq_param[, which(names(ct_seq_param) != motif_fg)]) / as.numeric(sum(tail(ct_seq_param[, which(names(ct_seq_param) != motif_fg)], 1)))
        nFg <- ct_seq_param[, which(names(ct_seq_param) == motif_fg)] / as.numeric(tail(ct_seq_param[, which(names(ct_seq_param) == motif_fg)], 1))
        return(.EstimateFdr(nFg = nFg, nBg = nBg, nThreshold = as.numeric(rownames(ct_seq_param))))
      })
      names(listFdrEstByThr) <- paste("FDRe_", cModMotifsAsForeground, sep = "")
    }
  } else {
    ct_seq_param_ctrl <- .GetCumsumTableParamByMotif(
      grangesDataWithSeq = grangesDataWithSeqControl,
      cNameParamToTest = cNameParamToTest,
      nRoundDigits = nRoundDigits
    )

    ct_seq_param$threshold <- as.numeric(rownames(ct_seq_param))
    ct_seq_param_ctrl$threshold <- as.numeric(rownames(ct_seq_param_ctrl))

    mCounts <- merge(
      x = ct_seq_param, y = ct_seq_param_ctrl,
      by = "threshold", all = TRUE, sort = TRUE
    )
    mCounts <- mCounts[order(mCounts$threshold, decreasing = TRUE), ]

    mCounts[is.na(mCounts[, 2]), 2] <- 0
    mCounts[is.na(mCounts[, 3]), 3] <- 0

    nFg <- mCounts[, 2] / as.numeric(tail(mCounts[, 2], 1))
    nBg <- mCounts[, 3] / as.numeric(tail(mCounts[, 3], 1))
    listFdrEstByThr <- list(FDRe_vsCtrl = .EstimateFdr(
      nFg = nFg, nBg = nBg,
      nThreshold = mCounts$threshold
    ))
  }
  return(listFdrEstByThr)
}

#' GetCumsumTableParamByMotif Function (FDREst)
#'
#' Return a dataframe with cumulative number of values of the parameter provided in the GRanges object (from highest to lowest values).
#' @param grangesDataWithSeq A GRanges-like object containing, in the extra columns, the parameter to be tested and the sequence.
#' @param cNameParamToTest The name of the column containing the parameter to be tested.
#' @param nRoundDigits The number of digits for the thresholds on the parameter. A value of 0 would mean no value to be rounded:
#' this could results in errors for a continuous variable. Defaults to 1.
#' @param lTestBetweenMotifs If TRUE, add one column for each sequence in the GRanges object provided.
#' Else, calculates the cumulative number of values of the parameter provided independantly of the sequence associated. Defaults to FALSE.
#' @keywords internal
.GetCumsumTableParamByMotif <- function(grangesDataWithSeq,
                                        cNameParamToTest,
                                        nRoundDigits = 1,
                                        lTestBetweenMotifs = FALSE) {
  data_t <- grangesDataWithSeq
  if (lTestBetweenMotifs) {
    data_t$sequence <- as.character(data_t$sequence)
    data_t <- elementMetadata(data_t)[c(cNameParamToTest, "sequence")]
  } else {
    data_t <- elementMetadata(data_t)[cNameParamToTest]
  }
  data_t[cNameParamToTest] <- round(as.data.frame(data_t[cNameParamToTest]), digits = nRoundDigits)
  data_t <- table(data_t)
  if (lTestBetweenMotifs) {
    data_t <- data_t[order(as.numeric(rownames(data_t)), decreasing = TRUE), ]
    data_t <- as.data.frame.matrix(apply(data_t, 2, cumsum))
  } else {
    data_t <- data_t[order(as.numeric(names(data_t)), decreasing = TRUE)]
    data_t <- as.data.frame(cumsum(data_t))
  }
  return(data_t)
}

#' EstimateFdr Function (FDREst)
#'
#' Return a dataframe with the false discovery rate estimated by threshold on the parameter provided with the GRanges object(s).
#' @param nFg A numeric vector considered as the foreground signal.
#' @param nBg A numeric vector considered as the background signal.
#' @param nThreshold A numeric vector describing the thresholds used.
#' @keywords internal
.EstimateFdr <- function(nFg, nBg, nThreshold) {
  fdr <- nBg / nFg
  if (any(fdr > 1)) {
    fdr[fdr > 1] <- 1
  }
  fdr <- rev(fdr)
  fdr <- data.frame(fdr = fdr)
  fdr$threshold <- rev(nThreshold)
  fdr$fdr_cummin <- cummin(fdr$fdr)
  return(fdr)
}

#' DrawFdrEstList Function (FDREst)
#'
#' Return a plot describing the false discovery rate estimations by threshold on the parameter provided for each dataframe in the list provided.
#' @param listFdrEstByThr A list composed of x dataframes. In each dataframe:
#' \itemize{
#'   \item fdr: The false discovery rate estimated for this threshold.
#'   \item threshold: The threshold on the parameter.
#'   \item fdr_cummin: The minimum false discovery rate estimated for this threshold and less stringent thresholds (adjusted false discovery rate).
#' }
#' @param cNameParamToTest The name of the column containing the parameter to be tested.
#' @param nFdrPropForFilt A number indicating the false discovery rate to be used for filtering: this will allow to choose the closest threshold
#' below this number and represent it on the plot. Defaults to 0.05 (so fdr of 5\%).
#' @param lAdjustFdr If TRUE, display fdr_cummin (adjusted false discovery rate)
#' instead of fdr column. Defaults to TRUE.
#' @keywords DrawFdrEstList
#' @export
#' @examples
#' myGenome <- Biostrings::readDNAStringSet(system.file(
#'   package = "DNAModAnnot", "extdata",
#'   "ptetraurelia_mac_51_sca171819.fa"
#' ))
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' # Preparing a gposPacBioCSV dataset with sequences
#' myGposPacBioCSV <-
#'   ImportPacBioCSV(
#'     cPacBioCSVPath = system.file(
#'       package = "DNAModAnnot", "extdata",
#'       "ptetraurelia.bases.sca171819.csv"
#'     ),
#'     cSelectColumnsToExtract = c(
#'       "refName", "tpl", "strand", "base", "score",
#'       "ipdRatio", "coverage"
#'     ),
#'     lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'     cContigToBeAnalyzed = names(myGenome)
#'   )
#' myGrangesBaseCSV <- as(myGposPacBioCSV[myGposPacBioCSV$base == "A"], "GRanges")
#' myGrangesBaseCSVWithSeq <- GetGRangesWindowSeqandParam(
#'   grangesData = myGrangesBaseCSV,
#'   grangesGenome = myGrangesGenome,
#'   dnastringsetGenome = myGenome,
#'   nUpstreamBpToAdd = 0,
#'   nDownstreamBpToAdd = 1
#' )
#'
#' # FDR estimation by motif associated to modifications
#' myFdr_score_per_motif_list <-
#'   GetFdrEstListByThresh(
#'     grangesDataWithSeq = myGrangesBaseCSVWithSeq,
#'     cNameParamToTest = "score",
#'     nRoundDigits = 1,
#'     cModMotifsAsForeground = c("AG", "AT")
#'   )
#'
#' DrawFdrEstList(
#'   listFdrEstByThr = myFdr_score_per_motif_list,
#'   cNameParamToTest = "score",
#'   nFdrPropForFilt = 0.05
#' )
DrawFdrEstList <- function(listFdrEstByThr,
                           cNameParamToTest,
                           nFdrPropForFilt = 0.05,
                           lAdjustFdr = TRUE) {
  i <- length(names(listFdrEstByThr))
  layout(matrix(1:i, nrow = 1, ncol = i, byrow = FALSE))
  for (fdr_to_plot in names(listFdrEstByThr)) {
    if (gsub("FDRe_", "", fdr_to_plot) == "vsCtrl") {
      motif_fg <- NULL
    } else {
      motif_fg <- gsub("FDRe_", "", fdr_to_plot)
    }
    .DrawFdrEst(listFdrEstByThr[[fdr_to_plot]],
      cNameParamToTest = cNameParamToTest,
      cMotifFg = motif_fg,
      lAdjustFdr = lAdjustFdr, nFdrPropForFilt = nFdrPropForFilt
    )
  }
  layout(matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2))
}

#' DrawFdrEst Function (FDREst)
#'
#' Return a plot describing the false discovery rate estimations by threshold on the parameter provided for the dataframe provided.
#' @param dFdr A dataframe with 3 columns:
#' \itemize{
#'   \item fdr: The false discovery rate estimated for this threshold.
#'   \item threshold: The threshold on the parameter.
#'   \item fdr_cummin: The minimum false discovery rate estimated for this threshold and less stringent thresholds (adjusted false discovery rate).
#' }
#' @param cNameParamToTest The name of the column containing the parameter to be tested.
#' @param nFdrPropForFilt A number indicating the false discovery rate to be used for filtering: this will allow to choose the closest threshold
#' below this number and represent it on the plot. Defaults to 0.05 (so fdr of 5\%).
#' @param cMotifFg The name of the motif tested. If NULL, the title will display "(Sample vs Control)" instead. Defaults to NULL.
#' @param lAdjustFdr If TRUE, display fdr_cummin (adjusted false discovery rate)
#' instead of fdr column. Defaults to TRUE.
#' @keywords internal
.DrawFdrEst <- function(dFdr, cNameParamToTest, nFdrPropForFilt,
                        cMotifFg = NULL, lAdjustFdr = TRUE) {
  opar <- par()
  par(mar = c(5.1, 5.1, 4.1, 2.1))
  plot(
    x = dFdr$threshold,
    y = if (lAdjustFdr) {
      dFdr$fdr_cummin * 100
    } else {
      dFdr$fdr * 100
    },
    ylim = c(0, 100),
    xlab = paste(cNameParamToTest, "threshold"),
    ylab = "False Discovery Rate Estimation (%)",
    main = ifelse(is.null(cMotifFg),
      paste("FDR estimation (Sample vs Control) ",
        ifelse(lAdjustFdr, "(adjusted FDR)", ""),
        sep = ""
      ),
      paste("FDR estimation (", cMotifFg, " vs non-", cMotifFg, " motifs) ",
        ifelse(lAdjustFdr, "(adjusted FDR)", ""),
        sep = ""
      )
    ),
    cex.main = 1, cex.lab = 1.25, cex.axis = 1.25
  )

  .DrawFdrSegmentsXYAxis(
    dFdr = dFdr,
    nFdrPropForFilt = nFdrPropForFilt
  )

  par(mar = opar$mar)
}

#' DrawFdrSegmentsXYAxis Function (FDREst)
#'
#' Draw a line describing the threshold associated to a value of false discovery rate estimation (defined by nFdrPropForFilt).
#' @param dFdr A dataframe with 3 columns:
#' \itemize{
#'   \item fdr: The false discovery rate estimated for this threshold.
#'   \item threshold: The threshold on the parameter.
#'   \item fdr_cummin: The minimum false discovery rate estimated for this threshold and less stringent thresholds (adjusted false discovery rate).
#' }
#' @param nFdrPropForFilt A number indicating the false discovery rate to be used for filtering: this will allow to choose the closest threshold
#' below this number and represent it on the plot. Defaults to 0.05 (so fdr of 5\%).
#' @keywords internal
.DrawFdrSegmentsXYAxis <- function(dFdr,
                                   nFdrPropForFilt) {
  fdr_limit <- dFdr[head(which(dFdr$fdr_cummin < nFdrPropForFilt), 1), ]
  # draw full line if no fdr below threshold
  if (nrow(fdr_limit) == 0) {
    abline(h = nFdrPropForFilt * 100, col = "red", lwd = 2)
    text(
      x = 0, y = nFdrPropForFilt * 100,
      paste(nFdrPropForFilt * 100, "%", sep = ""),
      col = "red", cex = 1.25, xpd = TRUE, adj = c(0.5, 0)
    )
  } else {
    segments(
      x0 = fdr_limit$threshold, x1 = fdr_limit$threshold,
      y0 = -max(dFdr$fdr_cummin * 100), y1 = fdr_limit$fdr_cummin * 100,
      col = "red", lwd = 2
    )
    segments(
      x0 = -max(dFdr$threshold), x1 = fdr_limit$threshold,
      y0 = fdr_limit$fdr_cummin * 100, y1 = fdr_limit$fdr_cummin * 100,
      col = "red", lwd = 2
    )

    text(
      x = fdr_limit$threshold, y = 0,
      fdr_limit$threshold,
      col = "red", cex = 1.25, xpd = TRUE, adj = c(0, 0.5)
    )
    text(
      x = 0, y = fdr_limit$fdr_cummin * 100,
      paste(round(fdr_limit$fdr_cummin * 100, digits = 2), "%", sep = ""),
      col = "red", cex = 1.25, xpd = TRUE, adj = c(0.5, 0)
    )
  }
}

#' GetFdrBasedThreshLimit Function (FDREst)
#'
#' Return a plot describing the false discovery rate (fdr) estimations by threshold on the parameter provided for each dataframe in the list provided.
#' @param listFdrEstByThr A list composed of x dataframes. In each dataframe:
#' \itemize{
#'   \item fdr: The false discovery rate estimated for this threshold.
#'   \item threshold: The threshold on the parameter.
#'   \item fdr_cummin: The minimum false discovery rate estimated for this threshold and less stringent thresholds (adjusted false discovery rate).
#' }
#' @param nFdrPropForFilt A number indicating the false discovery rate to be used for filtering: this will allow to choose the closest threshold
#' below this number. Defaults to 0.05 (so fdr of 5\%).
#' @param lUseBestThrIfNoFdrThr For fdr calculation by motif: if no fdr-associated threshold can be retrieved for one motif,
#' return the strongest threshold identified for any other motif if lUseBestThrIfNoFdrThr is TRUE; if lUseBestThrIfNoFdrThr is FALSE,
#' return the max value for the threshold (so every modification in that motif will be filtered out automatically). Defaults to TRUE.
#' @keywords GetFdrBasedThreshLimit
#' @export
#' @examples
#' library(Biostrings)
#' myGenome <- readDNAStringSet(system.file(
#'   package = "DNAModAnnot", "extdata",
#'   "ptetraurelia_mac_51_sca171819.fa"
#' ))
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' # Preparing a gposPacBioCSV dataset with sequences
#' myGposPacBioCSV <-
#'   ImportPacBioCSV(
#'     cPacBioCSVPath = system.file(
#'       package = "DNAModAnnot", "extdata",
#'       "ptetraurelia.bases.sca171819.csv"
#'     ),
#'     cSelectColumnsToExtract = c(
#'       "refName", "tpl", "strand", "base",
#'       "score", "ipdRatio", "coverage"
#'     ),
#'     lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'     cContigToBeAnalyzed = names(myGenome)
#'   )
#' myGrangesBaseCSV <- as(myGposPacBioCSV[myGposPacBioCSV$base == "A"], "GRanges")
#' myGrangesBaseCSVWithSeq <- GetGRangesWindowSeqandParam(
#'   grangesData = myGrangesBaseCSV,
#'   grangesGenome = myGrangesGenome,
#'   dnastringsetGenome = myGenome,
#'   nUpstreamBpToAdd = 0,
#'   nDownstreamBpToAdd = 1
#' )
#'
#' # FDR estimation by motif associated to modifications
#' myFdr_score_per_motif_list <-
#'   GetFdrEstListByThresh(
#'     grangesDataWithSeq = myGrangesBaseCSVWithSeq,
#'     cNameParamToTest = "score",
#'     nRoundDigits = 1,
#'     cModMotifsAsForeground = c("AG", "AT")
#'   )
#'
#' GetFdrBasedThreshLimit(
#'   listFdrEstByThr = myFdr_score_per_motif_list,
#'   nFdrPropForFilt = 0.05,
#'   lUseBestThrIfNoFdrThr = TRUE
#' )
GetFdrBasedThreshLimit <- function(listFdrEstByThr, nFdrPropForFilt = 0.05,
                                   lUseBestThrIfNoFdrThr = TRUE) {
  fdr_limit <- lapply(1:length(listFdrEstByThr), function(i) {
    limit <- listFdrEstByThr[[i]][head(which(listFdrEstByThr[[i]]$fdr_cummin < nFdrPropForFilt), 1), "threshold"]
    return(limit)
  })
  if (lUseBestThrIfNoFdrThr) {
    fdr_limit[isEmpty(fdr_limit)] <- max(unlist(fdr_limit[!isEmpty(fdr_limit)]))
  } else {
    fdr_limit[isEmpty(fdr_limit)] <- max(max(do.call(rbind.data.frame, listFdrEstByThr)$threshold))
  }
  return(fdr_limit)
}
