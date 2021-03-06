### =========================================================================
### Functions from ModAn and ModAnPlot categories
### -------------------------------------------------------------------------
###
### DNAModAnnot v.0.0.0.9018 - 2021/02/03
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
### GetMeanParamByContig
### DrawBarplotBothStrands
### DrawDistriHistBox
### GetModReportPacBio
### GetModReportDeepSignal
### GetGRangesWindowSeqandParam
### .AddModMotifPctToDf
### .IncludeModPosInMot
### GetModRatioByContig
### DrawModLogo
### ExtractListModPosByModMotif
###
### -------------------------------------------------------------------------

#' GetMeanParamByContig Function (GloModAn)
#'
#' Return a list with the mean by strand of the parameter provided for all scaffolds of genome assembly provided.
#' @param grangesData A GRanges-like object containing, in the extra columns, the parameter to be analysed.
#' @param dnastringsetGenome A DNAStringSet object containing the sequence for each contig.
#' @param cParamName The name of the column containing the parameter to be analysed.
#' @return A list composed of 3 dataframes: 1 dataframe for each strand and 1 dataframe with both strands. In each dataframe:
#' \itemize{
#'   \item refName: The names of each contig.
#'   \item strand: The strand of each contig.
#'   \item width: The width of each contig.
#'   \item nb_sequenced: The number of bases sequenced by strand for each contig.
#'   \item seqPct: The percentage of bases sequenced by strand for each contig (percentage of sequencing).
#' }
#' @keywords GetMeanParamByContig
#' @importFrom Biobase isUnique
#' @export
#' @examples
#' myGenome <- Biostrings::readDNAStringSet(system.file(
#'   package = "DNAModAnnot", "extdata",
#'   "ptetraurelia_mac_51_sca171819.fa"
#' ))
#'
#' # Preparing a gposPacBioCSV dataset
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
#'
#' myMean_cov_list <- GetMeanParamByContig(
#'   grangesData = myGposPacBioCSV,
#'   dnastringsetGenome = myGenome,
#'   cParamName = "coverage"
#' )
#' myMean_cov_list
GetMeanParamByContig <- function(grangesData,
                                 dnastringsetGenome,
                                 cParamName) {
  count_table <- data.frame(
    refName = as.factor(seqnames(grangesData)),
    strand = as.factor(strand(grangesData)),
    parameter = mcols(grangesData)[[cParamName]]
  )
  mean_table <- aggregate(
    count_table$parameter,
    list(count_table$refName, count_table$strand),
    mean
  )
  colnames(mean_table) <- c("refName", "strand", paste("mean_", cParamName, sep = ""))

  mean_table <- merge(mean_table,
    data.frame(
      refName = names(dnastringsetGenome),
      width = width(dnastringsetGenome)
    ),
    by = "refName",
    all.y = TRUE
  )

  # those who have NA for one or both strands : add missing strands with 0 mean_param
  if (length(which(is.na(mean_table$strand))) > 0) {
    mean_table[is.na(mean_table$strand), paste("mean_", cParamName, sep = "")] <- 0
    mean_table[is.na(mean_table$strand), "strand"] <- "+"
  }
  if (length(which(isUnique(mean_table$refName))) > 0) {
    mean_table2 <- mean_table[isUnique(mean_table$refName), ]
    mean_table2$strand <- ifelse(mean_table2$strand == "+", "-", "+")
    mean_table2[, paste("mean_", cParamName, sep = "")] <- 0
    mean_table <- rbind(mean_table, mean_table2)
  }
  mean_table <- mean_table[with(mean_table, order(width, refName, strand, decreasing = TRUE)), ]

  f_strand <- subset(mean_table, strand == "+")
  r_strand <- subset(mean_table, strand == "-")

  return(list(
    mean_table = mean_table,
    f_strand = f_strand,
    r_strand = r_strand
  ))
}

#' DrawBarplotBothStrands Function (GloModAn)
#'
#' Return a barplot describing some parameter values provided by strand for each contig.
#' @param nParamByContigForward A numeric vector containing the parameter values of the forward strand to be plotted.
#' Must have the same order and same length as cContigNames and nParamByContigReverse.
#' @param nParamByContigReverse A numeric vector containing the parameter values of the reverse strand to be plotted.
#' Must have the same order and same length as nParamByContigForward and cContigNames.
#' @param cContigNames A character vector containing the names of the contigs.
#' Must have the same order and same length as nParamByContigForward and nParamByContigReverse.
#' @param cGraphName The graph name to be displayed on top of the plot.
#' @param lIsOrderedLargestToSmallest If contigs are ordered from largest to smallest contig,
#' add TRUE to display "(from largest to smallest)" below the plot. Defaults to FALSE.
#' @keywords DrawBarplotBothStrands
#' @export
#' @examples
#' DrawBarplotBothStrands(
#'   nParamByContigForward = c(100, 86, 75, 56),
#'   nParamByContigReverse = c(96, 88, 80, 83),
#'   cContigNames = c("chrI", "chrII", "chrIII", "chrIV"),
#'   cGraphName = "Mean Coverage per contig",
#'   lIsOrderedLargestToSmallest = TRUE
#' )
DrawBarplotBothStrands <- function(nParamByContigForward,
                                   nParamByContigReverse,
                                   cContigNames,
                                   cGraphName,
                                   lIsOrderedLargestToSmallest = TRUE) {
  layout(matrix(c(
    4, 4,
    3, 1,
    3, 2
  ),
  nrow = 3, ncol = 2, byrow = TRUE
  ), widths = c(1, 29), heights = c(1, 2, 2))

  nMarLeft <- max(nchar(c(pretty(nParamByContigForward), pretty(nParamByContigReverse)))) * 1.1 + 1
  opar <- par()
  par(mar = c(0, nMarLeft, 4, 2))

  ylim_max <- max(max(nParamByContigForward), max(nParamByContigReverse))

  bp <- barplot(nParamByContigForward,
    ylab = "Forward strand", ylim = c(0, ylim_max),
    col = c("grey30", "grey70"),
    las = 2, xaxs = "i", yaxs = "i", cex.axis = 1.25, cex.lab = 1.5, mgp = c(nMarLeft * 3.5 / 5, 1, 0)
  )

  axis(side = 3, pos = NA, at = bp, labels = cContigNames, las = 2, cex.axis = 1.25)

  par(mar = c(4, nMarLeft, 0, 2))

  barplot(nParamByContigReverse,
    ylab = "Reverse strand", ylim = par("usr")[c(4, 3)],
    col = c("grey30", "grey70"),
    las = 1, xaxs = "i", yaxs = "i", cex.axis = 1.25, cex.lab = 1.5, mgp = c(nMarLeft * 3.5 / 5, 1, 0)
  )
  par(xpd = NA)
  if (lIsOrderedLargestToSmallest) {
    mtext("Contigs (from largest to smallest)", cex = 1.25, side = 1, line = 1.5)
  } else {
    mtext("Contigs", cex = 1.25, side = 1, line = 1.5)
  }
  par(xpd = FALSE)

  par(mar = c(0, 0, 0, 0.5))
  barplot(height = 1, col = "white", border = "white", axes = FALSE)

  par(xpd = NA)
  mtext(cGraphName, cex = 1.25, side = 4, srt = 90)
  par(xpd = FALSE)

  layout(matrix(c(1, 1, 1, 1), nrow = 2, ncol = 2))
  par(mar = opar$mar)
}

#' DrawDistriHistBox Function (GloModAn)
#'
#' Return a line-plot describing the cumulative length of contigs (ordered from largest to smallest contig).
#' @param nParam A numeric vector containing the parameter values to be plotted.
#' @param cGraphName The graph name to be displayed on top of the plot.
#' @param cParamName The name of the parameter to be displayed below the x-axis.
#' @param lTrimOutliers If TRUE, remove the outliers from the boxplot and trim the histogram to the borders
#' of the boxplot. Defaults to FALSE.
#' @param nXLimits A numeric vector giving the limits of the plot on the x-axis. If NULL, the limits will be set
#' to the minimum and the maximum of the nParam data (or to the inner fences (1.5*Interquartile Range) of the
#' boxplot if lTrimmingOutliers is TRUE). Defaults to NULL.
#' @keywords DrawDistriHistBox
#' @export
#' @examples
#' # loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(
#'   package = "DNAModAnnot", "extdata",
#'   "ptetraurelia_mac_51_sca171819.fa"
#' ))
#'
#' # Preparing a gposPacBioCSV dataset
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
#'
#' DrawDistriHistBox(head(myGposPacBioCSV$coverage, 100000),
#'   cGraphName = "Coverage distribution of some bases sequenced",
#'   cParamName = "Coverage",
#'   lTrimOutliers = TRUE
#' )
DrawDistriHistBox <- function(nParam,
                              cGraphName,
                              cParamName,
                              lTrimOutliers = FALSE,
                              nXLimits = NULL) {
  d <- density(nParam)
  h <- hist(nParam, plot = FALSE)
  if (length(seq(from = min(h$breaks), to = max(h$breaks), by = 1)) < 10) {
    breaks_v <- seq(from = min(h$breaks), to = max(h$breaks), length.out = 200)
  } else {
    breaks_v <- seq(from = min(h$breaks), to = max(h$breaks), by = 1)
  }
  h <- hist(nParam, breaks = breaks_v, plot = FALSE)

  if (lTrimOutliers) {
    b <- boxplot(nParam, plot = FALSE)
  }

  if (is.null(nXLimits)) {
    nXLimits <- c()
    nXLimits[1] <- ifelse(!lTrimOutliers, min(h$breaks), b$stats[1])
    nXLimits[2] <- ifelse(!lTrimOutliers, max(h$breaks), b$stats[5])
  }

  opar <- par()
  layout(matrix(c(1, 2), 2, 1), heights = c(1, 9))
  par(mar = c(0.1, 5.1, 2.1, 2.1))
  b <- boxplot(nParam,
    horizontal = TRUE, xaxt = "n", frame = FALSE, outline = !lTrimOutliers,
    main = ifelse(!lTrimOutliers,
      cGraphName,
      paste(cGraphName, " (no outliers)", sep = "")
    ),
    outcol = rgb(0, 0, 0, 0.01),
    ylim = nXLimits
  )

  par(mar = c(5, 5.1, 0.1, 2.1))
  hist(nParam,
    freq = FALSE, probability = TRUE,
    xlab = cParamName, ylab = "Density (%)", main = NULL,
    col = "gray",
    xlim = nXLimits,
    ylim = c(0, ifelse(max(h$density) < max(d$y), max(d$y), max(h$density))),
    breaks = breaks_v
  )
  lines(d, col = "red")
  layout(matrix(c(1), 1, 1))
  par(mar = opar$mar)
}

#' GetModReportPacBio Function (GloModAn)
#'
#' Return a report with global characteristics of DNA modifications (Mod) distribution in the genome assembly provided.
#' (adapted to PacBio data)
#' @param dnastringsetGenome A DNAStringSet object containing the sequence for each contig.
#' @param grangesGenome A GRanges object containing the width of each contig.
#' @param grangesPacBioGFF A GRanges object containing PacBio GFF data.
#' @param gposPacBioCSVBase A GPos object containing PacBio CSV data for sites that can be targeted by the modification only.
#' @param cOrgAssemblyName The name of the genome assembly provided.
#' @param cBaseLetterForMod The name of the base letter of the modified base.
#' @param cModNameInOutput Name for the modification in the output.
#' @keywords GetModReportPacBio
#' @export
#' @examples
#' # loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(
#'   package = "DNAModAnnot", "extdata",
#'   "ptetraurelia_mac_51_sca171819.fa"
#' ))
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' # Preparing a grangesPacBioGFF datasets
#' myGrangesPacBioGFF <-
#'   ImportPacBioGFF(
#'     cPacBioGFFPath = system.file(
#'       package = "DNAModAnnot", "extdata",
#'       "ptetraurelia.modifications.sca171819.gff"
#'     ),
#'     cNameModToExtract = "m6A",
#'     cModNameInOutput = "6mA",
#'     cContigToBeAnalyzed = names(myGenome)
#'   )
#'
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
#' myGposPacBioCSV <- myGposPacBioCSV[myGposPacBioCSV$base == "A"]
#'
#' # Mod report
#' myReport_Mod <- GetModReportPacBio(
#'   grangesGenome = myGrangesGenome,
#'   dnastringsetGenome = myGenome,
#'   grangesPacBioGFF = myGrangesPacBioGFF,
#'   gposPacBioCSVBase = myGposPacBioCSV,
#'   cOrgAssemblyName = "ptetraurelia_mac_51",
#'   cBaseLetterForMod = "A",
#'   cModNameInOutput = "6mA"
#' )
#' myReport_Mod
GetModReportPacBio <- function(dnastringsetGenome,
                               grangesGenome,
                               grangesPacBioGFF,
                               gposPacBioCSVBase,
                               cOrgAssemblyName,
                               cBaseLetterForMod,
                               cModNameInOutput) {
  dReportTable <- data.frame(
    row.names = cOrgAssemblyName,
    nbMod = length(grangesPacBioGFF),
    nbBaseSequenced = length(gposPacBioCSVBase),
    nbBaseAssembly = letterFrequency(dnastringsetGenome, letters = cBaseLetterForMod, collapse = TRUE) +
      letterFrequency(reverseComplement(dnastringsetGenome), letters = cBaseLetterForMod, collapse = TRUE)
  )

  dReportTable$ratioMod <- dReportTable$nbMod / dReportTable$nbBaseSequenced
  dReportTable$ratioMod_corrected <- dReportTable$ratioMod * mean(grangesPacBioGFF$frac)
  dReportTable$ratioModAllAssemblyBases <- dReportTable$nbMod / dReportTable$nbBaseAssembly
  dReportTable$ratioModAllAssemblyBases_corrected <- dReportTable$ratioModAllAssemblyBases * mean(grangesPacBioGFF$frac)
  dReportTable$mean_frac_Mod <- mean(grangesPacBioGFF$frac)
  dReportTable$mean_coverage_Mod <- mean(grangesPacBioGFF$coverage)
  dReportTable$mean_score_Mod <- mean(grangesPacBioGFF$score)
  dReportTable$mean_ipdRatio_Mod <- mean(grangesPacBioGFF$ipdRatio)
  dReportTable$mean_identificationQv_Mod <- mean(grangesPacBioGFF$identificationQv)

  colnames(dReportTable) <- gsub(
    pattern = "Base", replacement = cBaseLetterForMod,
    x = colnames(dReportTable)
  )
  colnames(dReportTable) <- gsub(
    pattern = "Mod", replacement = cModNameInOutput,
    x = colnames(dReportTable)
  )

  dReportTable <- .AddModMotifPctToDf(dReportTable,
    grangesModPos = grangesPacBioGFF, grangesGenome = grangesGenome,
    dnastringsetGenome = dnastringsetGenome,
    nModPositionInMotif = 1, nUpstreamBpToAdd = 0, nDownstreamBpToAdd = 1,
    cBaseLetterForMod = cBaseLetterForMod, cModNameInOutput = cModNameInOutput
  )


  dReportTable <- .AddModMotifPctToDf(dReportTable,
    grangesModPos = grangesPacBioGFF, grangesGenome = grangesGenome,
    dnastringsetGenome = dnastringsetGenome,
    nModPositionInMotif = 2, nUpstreamBpToAdd = 1, nDownstreamBpToAdd = 0,
    cBaseLetterForMod = cBaseLetterForMod, cModNameInOutput = cModNameInOutput
  )


  dReportTable <- .AddModMotifPctToDf(dReportTable,
    grangesModPos = grangesPacBioGFF, grangesGenome = grangesGenome,
    dnastringsetGenome = dnastringsetGenome,
    nModPositionInMotif = 2, nUpstreamBpToAdd = 1, nDownstreamBpToAdd = 1,
    cBaseLetterForMod = cBaseLetterForMod, cModNameInOutput = cModNameInOutput
  )

  return(t(dReportTable))
}

#' GetModReportDeepSignal Function (GloModAn)
#'
#' Return a report with global characteristics of DNA modifications (Mod) distribution in the genome assembly provided.
#' (adapted to data from DeepSignal software)
#' @param dnastringsetGenome A DNAStringSet object containing the sequence for each contig.
#' @param grangesGenome A GRanges object containing the width of each contig.
#' @param gposDeepSignalMod An UnStitched GPos object containing DeepSignal modified sites data.
#' @param gposDeepSignalModBase An UnStitched GPos object containing DeepSignal modification target sites data.
#' @param cOrgAssemblyName The name of the genome assembly provided.
#' @param cBaseLetterForMod The name of the base letter of the modified base.
#' @param cModNameInOutput Name for the modification in the output.
#' @keywords GetModReportDeepSignal
#' @export
#' @examples
#' # preparing genome (simulated)
#' myGenome <- Biostrings::DNAStringSet(paste0(rep("ATCG", 100000), collapse = ""))
#' names(myGenome) <- "NC_000001.11"
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' # Loading Nanopore data
#' myDeepSignalModPath <- system.file(
#'   package = "DNAModAnnot", "extdata",
#'   "FAB39088-288418386-Chr1.CpG.call_mods.frequency.tsv"
#' )
#' mygposDeepSignalModBase <- ImportDeepSignalModFrequency(
#'   cDeepSignalModPath = myDeepSignalModPath,
#'   lSortGPos = TRUE,
#'   cContigToBeAnalyzed = "all"
#' )
#'
#' # Filtering
#' mygposDeepSignalMod <- FiltDeepSignal(
#'   gposDeepSignalModBase = mygposDeepSignalModBase,
#'   cParamNameForFilter = "frac",
#'   nFiltParamLoBoundaries = 0,
#'   nFiltParamUpBoundaries = 1,
#'   cFiltParamBoundariesToInclude = "upperOnly"
#' )$Mod
#'
#' # Mod report
#' myReport_Mod <- GetModReportDeepSignal(
#'   dnastringsetGenome = myGenome,
#'   grangesGenome = myGrangesGenome,
#'   gposDeepSignalMod = as(mygposDeepSignalMod, "GRanges"),
#'   gposDeepSignalModBase = as(mygposDeepSignalModBase, "GRanges"),
#'   cOrgAssemblyName = "Test_function",
#'   cBaseLetterForMod = "C", cModNameInOutput = "5mC"
#' )
#' myReport_Mod
GetModReportDeepSignal <- function(dnastringsetGenome,
                                   grangesGenome,
                                   gposDeepSignalMod,
                                   gposDeepSignalModBase,
                                   cOrgAssemblyName,
                                   cBaseLetterForMod,
                                   cModNameInOutput) {
  dReportTable <- data.frame(
    row.names = cOrgAssemblyName,
    nbMod = length(gposDeepSignalMod),
    nbSitesSequenced = length(gposDeepSignalModBase),
    nbBaseAssembly = letterFrequency(dnastringsetGenome, letters = cBaseLetterForMod, collapse = TRUE) +
      letterFrequency(reverseComplement(dnastringsetGenome), letters = cBaseLetterForMod, collapse = TRUE)
  )

  dReportTable$ratioMod <- dReportTable$nbMod / dReportTable$nbSitesSequenced
  dReportTable$ratioMod_corrected <- dReportTable$ratioMod * mean(gposDeepSignalMod$frac)
  dReportTable$ratioModAllAssemblyBases <- dReportTable$nbMod / dReportTable$nbBaseAssembly
  dReportTable$ratioModAllAssemblyBases_corrected <- dReportTable$ratioModAllAssemblyBases * mean(gposDeepSignalMod$frac)

  dReportTable$mean_prob1sum_Mod <- mean(gposDeepSignalMod$prob_1_sum)
  dReportTable$mean_countModified_Mod <- mean(gposDeepSignalMod$count_modified)
  dReportTable$mean_coverage_Mod <- mean(gposDeepSignalMod$coverage)
  dReportTable$mean_frac_Mod <- mean(gposDeepSignalMod$frac)

  colnames(dReportTable) <- gsub(
    pattern = "Base", replacement = cBaseLetterForMod,
    x = colnames(dReportTable)
  )
  colnames(dReportTable) <- gsub(
    pattern = "Mod", replacement = cModNameInOutput,
    x = colnames(dReportTable)
  )

  dReportTable <- .AddModMotifPctToDf(dReportTable,
    grangesModPos = gposDeepSignalMod, grangesGenome = grangesGenome,
    dnastringsetGenome = dnastringsetGenome,
    nModPositionInMotif = 1, nUpstreamBpToAdd = 0, nDownstreamBpToAdd = 1,
    cBaseLetterForMod = cBaseLetterForMod, cModNameInOutput = cModNameInOutput
  )


  dReportTable <- .AddModMotifPctToDf(dReportTable,
    grangesModPos = gposDeepSignalMod, grangesGenome = grangesGenome,
    dnastringsetGenome = dnastringsetGenome,
    nModPositionInMotif = 2, nUpstreamBpToAdd = 1, nDownstreamBpToAdd = 0,
    cBaseLetterForMod = cBaseLetterForMod, cModNameInOutput = cModNameInOutput
  )


  dReportTable <- .AddModMotifPctToDf(dReportTable,
    grangesModPos = gposDeepSignalMod, grangesGenome = grangesGenome,
    dnastringsetGenome = dnastringsetGenome,
    nModPositionInMotif = 2, nUpstreamBpToAdd = 1, nDownstreamBpToAdd = 1,
    cBaseLetterForMod = cBaseLetterForMod, cModNameInOutput = cModNameInOutput
  )

  return(t(dReportTable))
}

#' GetGRangesWindowSeqandParam Function (GloModAn)
#'
#' Return the GRanges object provided with the sequence associated to each position (and can also retrieve the sequence around each position).
#' @param dnastringsetGenome A DNAStringSet object containing the sequence for each contig.
#' @param grangesGenome A GRanges object containing the width of each contig.
#' @param grangesData A GRanges object containing Modifications Positions data to be extracted with the sequence.
#' @param nUpstreamBpToAdd Number of base pairs to add upstream of the range from the GRanges object provided
#' to obtain some sequence upstream of range. If some new ranges do not fit in the ranges of the contigs (provided with grangesGenome),
#' those new ranges will be removed. New windows with gaps are also removed. Defaults to 0.
#' @param nDownstreamBpToAdd Number of base pairs to add downstream of the range from the GRanges object provided
#' to obtain some sequence downstream of range. If some new ranges do not fit in the ranges of the contigs (provided with grangesGenome),
#' those new ranges will be removed. New windows with gaps are also removed. Defaults to 0.
#' @keywords GetGRangesWindowSeqandParam
#' @importFrom GenomicRanges resize width
#' @importFrom IRanges subsetByOverlaps
#' @importFrom BSgenome getSeq
#' @export
#' @examples
#' # loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(
#'   package = "DNAModAnnot", "extdata",
#'   "ptetraurelia_mac_51_sca171819.fa"
#' ))
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' # Preparing a grangesPacBioGFF datasets
#' myGrangesPacBioGFF <-
#'   ImportPacBioGFF(
#'     cPacBioGFFPath = system.file(
#'       package = "DNAModAnnot", "extdata",
#'       "ptetraurelia.modifications.sca171819.gff"
#'     ),
#'     cNameModToExtract = "m6A",
#'     cModNameInOutput = "6mA",
#'     cContigToBeAnalyzed = names(myGenome)
#'   )
#'
#' # Retrieve GRanges with sequence
#' myPositions_Mod_Granges_wSeq <- GetGRangesWindowSeqandParam(
#'   grangesData = myGrangesPacBioGFF,
#'   grangesGenome = myGrangesGenome,
#'   dnastringsetGenome = myGenome,
#'   nUpstreamBpToAdd = 5,
#'   nDownstreamBpToAdd = 5
#' )
#' myPositions_Mod_Granges_wSeq
GetGRangesWindowSeqandParam <- function(grangesData,
                                        grangesGenome,
                                        dnastringsetGenome,
                                        nUpstreamBpToAdd = 0, nDownstreamBpToAdd = 0) {
  ans <- grangesData

  print("Window adjustment...")
  ans <- resize(ans, width = width(ans) + nUpstreamBpToAdd, fix = "end")
  ans <- resize(ans, width = width(ans) + nDownstreamBpToAdd, fix = "start")

  # remove windows not included completely in contig windows
  gr_a <- grangesGenome
  print("Removing Windows not included completely in contig windows...")
  ans <- subsetByOverlaps(ans, gr_a, ignore.strand = FALSE, type = "within")

  print("Getting sequence...")
  ans$sequence <- as.factor(as.character(getSeq(dnastringsetGenome, ans)))

  print("Window adjustment...")
  ans <- resize(ans, width = width(ans) - nUpstreamBpToAdd, fix = "end")
  ans <- resize(ans, width = width(ans) - nDownstreamBpToAdd, fix = "start")

  return(ans)
}

#' AddModMotifPctToDf Function (GloModAn)
#'
#' Return the dataframe provided with the motif percentages associated to modifications in the GRanges object.
#' @param dReportTable The dataframe in which the motif percentages must be included.
#' @param grangesModPos A GRanges object containing Modifications Positions data.
#' @param grangesGenome A GRanges object containing the width of each contig.
#' @param dnastringsetGenome A DNAStringSet object containing the sequence for each contig.
#' @param cBaseLetterForMod The name of the base letter of the modified base.
#' @param cModNameInOutput Name for the modification in the output.
#' @param nModPositionInMotif Number representing the position of the modification in the ranges
#' (after resizing if some base pairs were added). Defaults to 1.
#' @param nUpstreamBpToAdd Number of base pairs to add upstream of the range from the GRanges object provided
#' to obtain some sequence upstream of range. If some new ranges do not fit in the ranges of the contigs (provided with grangesGenome),
#' those new ranges will be removed. New windows with gaps are also removed. Defaults to 0.
#' @param nDownstreamBpToAdd Number of base pairs to add downstream of the range from the GRanges object provided
#' to obtain some sequence downstream of range. If some new ranges do not fit in the ranges of the contigs (provided with grangesGenome),
#' those new ranges will be removed. New windows with gaps are also removed. Defaults to 0.
#' @keywords internal
#' @importFrom BiocGenerics table
.AddModMotifPctToDf <- function(dReportTable,
                                grangesModPos,
                                grangesGenome,
                                dnastringsetGenome,
                                cBaseLetterForMod,
                                cModNameInOutput,
                                nModPositionInMotif = 1,
                                nUpstreamBpToAdd = 0,
                                nDownstreamBpToAdd = 0) {
  positions_Mod_Granges_wSeq <- GetGRangesWindowSeqandParam(
    grangesModPos, grangesGenome, dnastringsetGenome,
    nUpstreamBpToAdd, nDownstreamBpToAdd
  )

  motif_pct <- as.data.frame.matrix(
    t(table(positions_Mod_Granges_wSeq$sequence)) / sum(t(table(positions_Mod_Granges_wSeq$sequence)))
  )

  motif_pct <- .IncludeModPosInMot(motif_pct, cBaseLetterForMod, cModNameInOutput, nModPositionInMotif)

  names(motif_pct) <- paste("pct_", names(motif_pct), sep = "")

  dReportTable <- cbind(dReportTable, motif_pct)

  return(dReportTable)
}

#' IncludeModPosInMot Function (GloModAn)
#'
#' Return the numeric vector provided with new names to include the position of the modification.
#' @param nParamByMotif A numeric vector named with motifs.
#' @param cBaseLetterForMod The name of the base letter of the modified base.
#' @param cModNameInOutput Name for the modification in the output.
#' @param nModPositionInMotif Number representing the position of the modification in the ranges
#' (after resizing if some base pairs were added). Defaults to 1.
#' @param lExtractOnlyMotifWithMod If TRUE, the numeric vector will be returned with only the
#' motifs that include DNA modifications. Defaults to TRUE.
#' @keywords internal
.IncludeModPosInMot <- function(nParamByMotif,
                                cBaseLetterForMod,
                                cModNameInOutput,
                                nModPositionInMotif = 1,
                                lExtractOnlyMotifWithMod = TRUE) {
  pattern <- paste("^(.{", (nModPositionInMotif - 1), "})", cBaseLetterForMod, "(.*)$", sep = "")
  replacement <- paste("\\1", cModNameInOutput, "\\2", sep = "")
  names(nParamByMotif) <- gsub(
    pattern = pattern, replacement = replacement,
    x = names(nParamByMotif)
  )
  if (lExtractOnlyMotifWithMod) {
    nParamByMotif <- nParamByMotif[grep(cModNameInOutput, names(nParamByMotif))]
  }
  return(nParamByMotif)
}

#' GetModRatioByContig Function (GloModAn)
#'
#' Return a list with the Modification ratio (Mod ratio) by strand for all scaffolds of genome assembly provided.
#' For "b" as the base that can be modified, Mod ratio = Number of modified "b" / Total number of "b".
#' @param grangesModPos A GRanges object containing Modifications Positions data.
#' @param gposModTargetBasePos A GPos object containing Base Positions that can be targeted by the modification.
#' @param dnastringsetGenome A DNAStringSet object containing the sequence for each contig.
#' @param cBaseLetterForMod The name of the base letter of the modified base.
#' @keywords GetModRatioByContig
#' @export
#' @examples
#' # loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(
#'   package = "DNAModAnnot", "extdata",
#'   "ptetraurelia_mac_51_sca171819.fa"
#' ))
#'
#' # Preparing a grangesPacBioGFF and gposPacBioCSV datasets
#' myGrangesPacBioGFF <-
#'   ImportPacBioGFF(
#'     cPacBioGFFPath = system.file(
#'       package = "DNAModAnnot", "extdata",
#'       "ptetraurelia.modifications.sca171819.gff"
#'     ),
#'     cNameModToExtract = "m6A",
#'     cModNameInOutput = "6mA",
#'     cContigToBeAnalyzed = names(myGenome)
#'   )
#'
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
#' myGposPacBioCSV <- myGposPacBioCSV[myGposPacBioCSV$base == "A"]
#'
#' # Mod report
#' myMod_ratio_list <- GetModRatioByContig(
#'   grangesModPos = myGrangesPacBioGFF,
#'   gposModTargetBasePos = myGposPacBioCSV,
#'   dnastringsetGenome = myGenome,
#'   cBaseLetterForMod = "A"
#' )
#' myMod_ratio_list
GetModRatioByContig <- function(grangesModPos,
                                gposModTargetBasePos,
                                dnastringsetGenome,
                                cBaseLetterForMod) {
  mod_count_table <- data.frame(
    refName = as.factor(seqnames(grangesModPos)),
    strand = as.factor(strand(grangesModPos))
  )
  mod_count_table <- as.data.frame(table(mod_count_table))
  colnames(mod_count_table) <- c("refName", "strand", "Mod_count")

  baseS_count_table <- data.frame(
    refName = as.factor(seqnames(gposModTargetBasePos)),
    strand = as.factor(strand(gposModTargetBasePos))
  )
  baseS_count_table <- as.data.frame(table(baseS_count_table))
  colnames(baseS_count_table) <- c("refName", "strand", "BaseSequenced_count")

  dnastringsetGenome_order <- dnastringsetGenome[order(names(dnastringsetGenome)), ]
  baseA_count_table <- data.frame(
    names(dnastringsetGenome_order), "+",
    letterFrequency(dnastringsetGenome_order, letters = cBaseLetterForMod),
    width(dnastringsetGenome_order)
  )
  tmp_table <- data.frame(
    names(dnastringsetGenome_order), "-",
    letterFrequency(reverseComplement(dnastringsetGenome_order), letters = cBaseLetterForMod),
    width(dnastringsetGenome_order)
  )
  baseA_count_table <- rbind(baseA_count_table, tmp_table)
  colnames(baseA_count_table) <- c("refName", "strand", "BaseAssembly_count", "width")

  base_count_table <- merge(baseS_count_table,
    baseA_count_table,
    by = c("refName", "strand"),
    all.y = TRUE
  )

  mod_count_table <- merge(mod_count_table,
    base_count_table,
    by = c("refName", "strand"),
    all.y = TRUE
  )

  # those who have NA for one or both strands : replace NA by 0
  mod_count_table[is.na(mod_count_table$Mod_count), "Mod_count"] <- 0
  mod_count_table[is.na(mod_count_table$BaseSequenced_count), "BaseSequenced_count"] <- 0
  mod_count_table[is.na(mod_count_table$BaseAssembly_count), "BaseAssembly_count"] <- 0

  mod_count_table <- mod_count_table[with(mod_count_table, order(width, refName, strand, decreasing = TRUE)), ]

  mod_count_table$Mod_ratio <- mod_count_table$Mod_count / mod_count_table$BaseSequenced_count
  mod_count_table$Mod_ratio_AllAssemblyBases <- mod_count_table$Mod_count / mod_count_table$BaseAssembly_count

  f_strand <- subset(mod_count_table, strand == "+")
  r_strand <- subset(mod_count_table, strand == "-")

  return(list(
    mod_count_table = mod_count_table,
    f_strand = f_strand,
    r_strand = r_strand
  ))
}

#' DrawModLogo Function (GloModAn)
#'
#' Return a plot describing the sequence motif associated to the sequences provided.
#' Sequences that do not have full width and sequences that have some
#' N or some gaps "-" are automatically removed before drawing the sequence plot.
#'
#' This function reduce the background signal using the genomic composition and also computes
#' the signal of depleted bases among the sequence provided.
#' Positions that are fixed and displaying only one base are taken in account and avoid these two corrections.
#' It is also possible to tag some positions (+/- labels) in the logo.
#'
#' @param dnastringsetSeqAroundMod A DNAStringSet object containing the sequence around each DNA modification.
#' All these sequences must have the same size and the modification must have the same position in each sequence (e.g. at the center).
#' @param cYunit Units to be used for the y axis. Can be "ic" (information content),
#' "ic_hide_bg" (information content + hide low signal) or "prob" (probability).
#' "ic_hide_bg" will hide the letters on the upper panel (OR the lower panel) if these letters
#' have probabilities above (OR below) the background signal based on the genomic composition.
#' Defaults to "ic_hide_bg".
#' @param nGenomicBgACGT A numeric vector giving the background to be corrected with the genomic composition
#' in Adenine (A) then Cytosine (C) then Guanine (G) then Thymine (T). Defaults to c(A=0.25, C=0.25, G=0.25, T=0.25).
#' @param cColorsACGT A character vector giving the color for the Adenine (A) then Cytosine (C) then Guanine (G)
#' then Thymine (T) letters on the plot. Defaults to c(A="red2", C="blue2", G="orange2", T="green3").
#' @param lPlotNegYAxis If TRUE, allow plotting the lower panel to show depletion among the sequences provided. Defaults to TRUE.
#' @param nPositionsToAnnotate A numeric vector giving the positions to highlight on the logo using small triangular tags
#' (e.g. DNA modifications on a fixed position). Defaults to NULL.
#' @param cAnnotationText A character vector. If nPositionsToAnnotate is not NULL, this option can provide labels to be associated to
#' each annotation tag. If the vector's length is 1, it will be used for all tags to be displayed. Defaults to NULL.
#' @param nTagTextFontSize Font size for the text associated to the annotation tags. Defaults to 20.
#' @param nXFontSize Font size for the text associated to the x axis. Defaults to 15.
#' @param nYFontSize Font size for the text associated to the y axis. Defaults to 15.
#' @param lXAxis If TRUE, the x-axis is plotted. Defaults to TRUE.
#' @param lYAxis If TRUE, the y-axis is plotted. Defaults to TRUE.
#' @keywords DrawModLogo
#' @importFrom Biostrings consensusMatrix
#' @importFrom seqLogo makePWM
#' @importFrom stats setNames
#' @export
#' @examples
#' # loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(
#'   package = "DNAModAnnot", "extdata",
#'   "ptetraurelia_mac_51_sca171819.fa"
#' ))
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' # Preparing a grangesPacBioGFF datasets
#' myGrangesPacBioGFF <-
#'   ImportPacBioGFF(
#'     cPacBioGFFPath = system.file(
#'       package = "DNAModAnnot", "extdata",
#'       "ptetraurelia.modifications.sca171819.gff"
#'     ),
#'     cNameModToExtract = "m6A",
#'     cModNameInOutput = "6mA",
#'     cContigToBeAnalyzed = names(myGenome)
#'   )
#'
#' # Retrieve GRanges with sequence
#' myPositions_Mod_Granges_wSeq <- GetGRangesWindowSeqandParam(
#'   grangesData = myGrangesPacBioGFF,
#'   grangesGenome = myGrangesGenome,
#'   dnastringsetGenome = myGenome,
#'   nUpstreamBpToAdd = 5,
#'   nDownstreamBpToAdd = 5
#' )
#'
#' DrawModLogo(
#'   dnastringsetSeqAroundMod = as(myPositions_Mod_Granges_wSeq$sequence, "DNAStringSet"),
#'   nGenomicBgACGT = c(0.35, 0.15, 0.15, 0.35),
#'   nPositionsToAnnotate = c(6), cAnnotationText = c("6mA")
#' )
DrawModLogo <- function(dnastringsetSeqAroundMod,
                        nGenomicBgACGT = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                        cYunit = "ic_hide_bg",
                        lPlotNegYAxis = TRUE,
                        cColorsACGT = c(A = "#4daf4a", C = "#377eb8", G = "#ffd92f", T = "#e41a1c"),
                        nPositionsToAnnotate = NULL,
                        cAnnotationText = NULL,
                        nTagTextFontSize = 20,
                        lXAxis = TRUE, lYAxis = TRUE,
                        nXFontSize = 15, nYFontSize = 15) {
  # remove all sequences that do not have full width
  dnastringsetSeqAroundMod <- dnastringsetSeqAroundMod[which(width(dnastringsetSeqAroundMod) == max(width(dnastringsetSeqAroundMod)))]
  # remove all sequences that have some N or some gaps "-"
  dnastringsetSeqAroundMod <- dnastringsetSeqAroundMod[which(letterFrequency(dnastringsetSeqAroundMod, letters = "N-") == 0)]


  mConsensus <- consensusMatrix(dnastringsetSeqAroundMod, as.prob = TRUE)
  mConsensus <- mConsensus[c("A", "C", "G", "T"), ]
  nFixedPositions <- which(apply(mConsensus, 2, function(x) any(x == 1)))

  # Background correction
  if (cYunit %in% c("ic", "ic_hide_bg")) {
    mConsensus <- t(t(mConsensus / nGenomicBgACGT) / colSums(mConsensus / nGenomicBgACGT))
  }

  pConsensus <- makePWM(mConsensus)

  nConsensus <- round(1 / (mConsensus + 0.00001), 4)
  nConsensus <- t(t(nConsensus) / colSums(nConsensus))

  if (cYunit %in% c("ic", "ic_hide_bg")) {
    nConsensus[, nFixedPositions] <- c(0.25, 0.25, 0.25, 0.25)
  } else {
    # if prob_scale: correct values for fixed positions
    nConsensus <- lapply(1:ncol(nConsensus), function(x) {
      if (any(mConsensus[, x] == 1)) {
        n <- setNames(rep(0, 4), nm = rownames(nConsensus))
        switch(names(mConsensus[mConsensus[, x] == 1, x]),
          "A" = n["T"] <- 1,
          "T" = n["A"] <- 1,
          "G" = n["C"] <- 1,
          "C" = n["G"] <- 1
        )
        return(n)
      }
      return(nConsensus[, x])
    })
    nConsensus <- do.call(cbind, nConsensus)
  }
  nConsensus <- makePWM(nConsensus)

  DrawLogoPosNegAxes(
    pwmUp = pConsensus, pwmDown = nConsensus, cColorsACGT = cColorsACGT,
    nPositionsToAnnotate = nPositionsToAnnotate,
    cAnnotationText = cAnnotationText, cYunit = cYunit,
    lPlotNegYAxis = lPlotNegYAxis,
    nGenomicBgACGT = nGenomicBgACGT,
    nTagTextFontSize = nTagTextFontSize,
    lXAxis = lXAxis, lYAxis = lYAxis, nXFontSize = nXFontSize, nYFontSize = nYFontSize
  )
}

#' DrawLogoPosNegAxes Function (GloModAn)
#'
#' Return a plot describing the sequence motif associated to the sequences provided.
#' Two panels are plotted: if cYunit option is used with "ic_hide_bg", data
#' from the lower panel should correspond to the depleted signal by comparison with upper panel
#' that should described the enriched signal.
#'
#' It is also possible to tag some positions (+/- labels) in the logo.
#'
#' @param pwmUp A seqLogo PWM (position weight matrix) object to be used for the upper panel.
#' @param pwmDown A seqLogo PWM (position weight matrix) object to be used for the lower panel. Should correspond to the depleted signal
#' (by comparison to the data provided with the pwmUp option).
#' @param cYunit Units to be used for the y axis. Can be "ic" (information content),
#' "ic_hide_bg" (information content + hide low signal) or "prob" (probability).
#' "ic_hide_bg" will hide the letters on the upper panel (OR the lower panel) if these letters
#' have probabilities above (OR below) the background signal based on the genomic composition.
#' Defaults to "ic_hide_bg".
#' @param nGenomicBgACGT A numeric vector giving the background to be corrected with the genomic composition
#' in Adenine (A) then Cytosine (C) then Guanine (G) then Thymine (T). Defaults to c(A=0.25, C=0.25, G=0.25, T=0.25).
#' @param cColorsACGT A character vector giving the color for the Adenine (A) then Cytosine (C) then Guanine (G)
#' then Thymine (T) letters on the plot. Defaults to c(A="red2", C="blue2", G="orange2", T="green3").
#' @param lPlotNegYAxis If TRUE, allow plotting the lower panel to show depletion among the sequences provided. Defaults to TRUE.
#' @param nPositionsToAnnotate A numeric vector giving the positions to highlight on the logo using small triangular tags
#' (e.g. DNA modifications on a fixed position). Defaults to NULL.
#' @param cAnnotationText A character vector. If nPositionsToAnnotate is not NULL, this option can provide labels to be associated to
#' each annotation tag. If the vector's length is 1, it will be used for all tags to be displayed. Defaults to NULL.
#' @param nTagTextFontSize Font size for the text associated to the annotation tags. Defaults to 20.
#' @param nXFontSize Font size for the text associated to the x axis. Defaults to 15.
#' @param nYFontSize Font size for the text associated to the y axis. Defaults to 15.
#' @param lXAxis If TRUE, the x-axis is plotted. Defaults to TRUE.
#' @param lYAxis If TRUE, the y-axis is plotted. Defaults to TRUE.
#' @keywords DrawLogoPosNegAxes
#' @importFrom seqLogo pwm ic
#' @importFrom grid grid.newpage plotViewport pushViewport dataViewport grid.abline grid.polygon unit grid.xaxis grid.yaxis grid.text gpar popViewport
#' @export
#' @examples
#' pwm1 <- matrix(c(
#'   0.25, 0.25, 0.25, 0.25,
#'   0.8, 0.05, 0.05, 0.1,
#'   0.33, 0.33, 0.33, 0.01
#' ), nrow = 4, byrow = FALSE)
#' row.names(pwm1) <- c("A", "C", "G", "T")
#' pwm2 <- t(t(1 / pwm1) / colSums((1 / pwm1)))
#'
#' pwm1 <- seqLogo::makePWM(pwm1)
#' pwm2 <- seqLogo::makePWM(pwm2)
#'
#' DrawLogoPosNegAxes(
#'   pwmUp = pwm1, pwmDown = pwm2,
#'   nPositionsToAnnotate = c(1, 3), cAnnotationText = c("Text?", "Depletion of T")
#' )
DrawLogoPosNegAxes <- function(pwmUp, pwmDown,
                               nGenomicBgACGT = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                               cYunit = "ic_hide_bg",
                               lPlotNegYAxis = TRUE,
                               cColorsACGT = c(A = "#4daf4a", C = "#377eb8", G = "#ffd92f", T = "#e41a1c"),
                               nPositionsToAnnotate = NULL,
                               cAnnotationText = NULL,
                               nTagTextFontSize = 20,
                               lXAxis = TRUE, lYAxis = TRUE,
                               nXFontSize = 15, nYFontSize = 15) {
  pUp <- pwm(pwmUp)
  pDo <- pwm(pwmDown)

  if (ncol(pUp) != ncol(pDo)) {
    stop("Error: pwmUp and pwmDown must have the same number of columns!")
  }
  npos <- ncol(pUp)

  if (cYunit == "ic_hide_bg") {
    ylim <- 2
    ylab <- "Information content"
    facs_up <- ic(pwmUp)
    facs_down <- ic(pwmDown)
    hideLetterUp <- pUp - nGenomicBgACGT <= 0
    hideLetterDo <- !hideLetterUp
  } else if (cYunit == "ic") {
    ylim <- 2
    ylab <- "Information content"
    facs_up <- ic(pwmUp)
    facs_down <- ic(pwmDown)
    hideLetterUp <- pUp == 0
    hideLetterDo <- pDo == 0
  } else if (cYunit == "prob") {
    ylim <- 1
    ylab <- "Probability"
    facs_up <- rep(1, npos)
    facs_down <- rep(1, npos)
    hideLetterUp <- pUp == 0
    hideLetterDo <- pDo == 0
  }

  charsUp <- rownames(pUp)
  charsDo <- rownames(pDo)
  lettersUp <- list(x = NULL, y = NULL, id = NULL, fill = NULL)
  lettersDo <- list(x = NULL, y = NULL, id = NULL, fill = NULL)

  wt <- 1
  x.pos <- 0

  for (j in 1:npos) {
    column <- pUp[, j]
    hts <- 0.95 * column * facs_up[j]
    letterOrder <- order(hts)

    y.pos <- 0
    for (i in 1:4) {
      letter <- charsUp[letterOrder[i]]
      ht <- hts[letterOrder[i]]
      if (hideLetterUp[letter, j]) {
        ht <- 0
      }

      lettersUp <- addLetter(
        lettersUp, letter, x.pos,
        y.pos, ht, wt, cColorsACGT
      )
      y.pos <- y.pos + ht + 0.01
    }

    column <- pDo[, j]
    hts <- -0.95 * column * facs_down[j]
    letterOrder <- order(-hts)

    y.pos <- 0
    for (i in 1:4) {
      letter <- charsDo[letterOrder[i]]
      ht <- hts[letterOrder[i]]
      if (hideLetterDo[letter, j]) {
        ht <- 0
      }

      oldYLen <- length(lettersDo$y)
      lettersDo <- addLetter(
        lettersDo, letter, x.pos,
        y.pos, ht, wt, cColorsACGT
      )
      newYLen <- length(lettersDo$y)

      # retrieve index of new values
      newValues <- seq(from = oldYLen + 1, to = newYLen, by = 1)

      # invert y-values of new letter
      old.y <- lettersDo$y[newValues]
      new.y <- (max(old.y) - old.y) + min(old.y)
      lettersDo$y[newValues] <- new.y

      y.pos <- y.pos + ht - 0.01
    }
    x.pos <- x.pos + wt
  }


  grid.newpage()

  bottomMargin <- ifelse(lXAxis, 2 + nXFontSize / 3.5, 2)
  leftMargin <- ifelse(lYAxis, 2 + nYFontSize / 3.5, 2)
  pushViewport(plotViewport(c(bottomMargin, leftMargin, 2, 2)))


  if (lPlotNegYAxis) {
    ylim1 <- min(lettersDo$y)
    ylim2 <- max(lettersUp$y) * 1.25
    if (tail(min(lettersDo$y):ylim2, 1) <= ylim) {
      ylim2 <- ylim * 1.25
      if (tail(min(lettersDo$y):ylim2, 1) <= ylim) {
        ylim2 <- ylim * 1.5
      }
    }
  } else {
    ylim1 <- 0
    ylim2 <- ylim
  }

  pushViewport(dataViewport(0:npos, ylim1:ylim2, name = "vp1"))

  grid.polygon(
    x = unit(lettersUp$x, "native"), y = unit(lettersUp$y, "native"),
    id = lettersUp$id, gp = gpar(fill = lettersUp$fill, col = "transparent")
  )
  if (lPlotNegYAxis) {
    grid.polygon(
      x = unit(lettersDo$x, "native"), y = unit(lettersDo$y, "native"),
      id = lettersDo$id, gp = gpar(fill = lettersDo$fill, col = "transparent")
    )
  }

  grid.abline(intercept = 0, slope = 0)

  if (!is.null(nPositionsToAnnotate)) {
    tag_x <- rep(nPositionsToAnnotate, each = 3) - rep(c(1, 0.5, 0), length(nPositionsToAnnotate))
    tag_y <- rep(ylim * 1.1, 3 * length(nPositionsToAnnotate)) - rep(c(0, 0.1, 0), length(nPositionsToAnnotate))

    # tag
    grid.polygon(
      x = unit(tag_x, "native"), y = unit(tag_y, "native"),
      id = rep(nPositionsToAnnotate, each = 3), gp = gpar(fill = c("orange"), col = "black")
    )

    if (!is.null(cAnnotationText)) {
      # text
      grid.text(
        label = cAnnotationText, x = unit(nPositionsToAnnotate - 0.5, "native"),
        y = unit(rep(ylim * 1.2, each = length(nPositionsToAnnotate)), "native"),
        gp = gpar(fontsize = nTagTextFontSize, col = "black")
      )
    }
  }

  if (lXAxis) {
    grid.xaxis(
      at = seq(0.5, npos - 0.5), label = seq_len(npos),
      gp = gpar(fontsize = nXFontSize)
    )
    grid.text("Position",
      y = unit(-3, "lines"),
      gp = gpar(fontsize = nXFontSize)
    )
  }
  if (lYAxis) {
    grid.yaxis(gp = gpar(fontsize = nYFontSize))
    grid.text(ylab,
      x = unit(-3, "lines"), rot = 90,
      gp = gpar(fontsize = nYFontSize)
    )
  }
  popViewport()
  popViewport()
}

#' addLetter Function from seqLogo package (GloModAn)
#'
#' Internal function from seqLogo package. See seqLogo documentation for more details.
#' User should use the addLetter function from seqLogo and should not use directly this version.
#' @keywords internal
addLetter <- function(letters, which, x.pos, y.pos, ht, wt, fill) {
  if (which == "A") {
    letter <- letterA(x.pos, y.pos, ht, wt, fill = fill["A"])
  }
  else if (which == "C") {
    letter <- letterC(x.pos, y.pos, ht, wt, fill = fill["C"])
  }
  else if (which == "G") {
    letter <- letterG(x.pos, y.pos, ht, wt, fill = fill["G"])
  }
  else if (which == "T") {
    letter <- letterT(x.pos, y.pos, ht, wt, fill = fill["T"])
  }
  else if (which == "U") {
    letter <- letterU(x.pos, y.pos, ht, wt, fill = fill["U"])
  }
  else {
    stop(sprintf("\"which\" must be one of %s", paste(names(fill),
      collapse = ", "
    )))
  }
  letters$x <- c(letters$x, letter$x)
  letters$y <- c(letters$y, letter$y)
  lastID <- ifelse(is.null(letters$id), 0, max(letters$id))
  letters$id <- c(letters$id, lastID + letter$id)
  letters$fill <- c(letters$fill, as.character(letter$fill))
  letters
}

#' letterA Function from seqLogo package (GloModAn)
#'
#' Internal function from seqLogo package. See seqLogo documentation for more details.
#' User should use the letterA function from seqLogo and should not use directly this version.
#' @keywords internal
letterA <- function(x.pos, y.pos, ht, wt, fill = "#61D04F", id = NULL) {
  x <- c(
    0, 4, 6, 2, 0, 4, 6, 10, 8, 4, 3.2, 6.8, 6.4, 3.6,
    3.2
  )
  y <- c(0, 10, 10, 0, 0, 10, 10, 0, 0, 10, 3, 3, 4, 4, 3)
  x <- 0.1 * x
  y <- 0.1 * y
  x <- .rescale(x)
  x <- x.pos + wt * x
  y <- y.pos + ht * y
  if (is.null(id)) {
    id <- c(rep(1L, 5), rep(2L, 5), rep(3L, 5))
  }
  else {
    id <- c(rep(id, 5), rep(id + 1L, 5), rep(id + 2L, 5))
  }
  fill <- rep(fill, 3)
  list(x = x, y = y, id = id, fill = fill)
}

#' letterG Function from seqLogo package (GloModAn)
#'
#' Internal function from seqLogo package. See seqLogo documentation for more details.
#' User should use the letterG function from seqLogo and should not use directly this version.
#' @keywords internal
letterG <- function(x.pos, y.pos, ht, wt, fill = "#F5C710", id = NULL) {
  angle1 <- seq(0.3 + pi / 2, pi, length = 100)
  angle2 <- seq(pi, 1.5 * pi, length = 100)
  x.l1 <- 0.5 + 0.5 * sin(angle1)
  y.l1 <- 0.5 + 0.5 * cos(angle1)
  x.l2 <- 0.5 + 0.5 * sin(angle2)
  y.l2 <- 0.5 + 0.5 * cos(angle2)
  x.l <- c(x.l1, x.l2)
  y.l <- c(y.l1, y.l2)
  x <- c(x.l, rev(x.l))
  y <- c(y.l, 1 - rev(y.l))
  x.i1 <- 0.5 + 0.35 * sin(angle1)
  y.i1 <- 0.5 + 0.35 * cos(angle1)
  x.i1 <- x.i1[y.i1 <= max(y.l1)]
  y.i1 <- y.i1[y.i1 <= max(y.l1)]
  y.i1[1] <- max(y.l1)
  x.i2 <- 0.5 + 0.35 * sin(angle2)
  y.i2 <- 0.5 + 0.35 * cos(angle2)
  x.i <- c(x.i1, x.i2)
  y.i <- c(y.i1, y.i2)
  x1 <- c(x.i, rev(x.i))
  y1 <- c(y.i, 1 - rev(y.i))
  x <- c(x, rev(x1))
  y <- c(y, rev(y1))
  h1 <- max(y.l1)
  r1 <- max(x.l1)
  h1 <- 0.4
  x.add <- c(r1, 0.5, 0.5, r1 - 0.2, r1 - 0.2, r1, r1)
  y.add <- c(h1, h1, h1 - 0.1, h1 - 0.1, 0, 0, h1)
  if (is.null(id)) {
    id <- c(rep(1, length(x)), rep(2, length(x.add)))
  }
  else {
    id <- c(rep(id, length(x)), rep(id + 1, length(x.add)))
  }
  x <- c(rev(x), x.add)
  y <- c(rev(y), y.add)
  x <- .rescale(x)
  x <- x.pos + wt * x
  y <- y.pos + ht * y
  fill <- rep(fill, 2)
  list(x = x, y = y, id = id, fill = fill)
}

#' letterC Function from seqLogo package (GloModAn)
#'
#' Internal function from seqLogo package. See seqLogo documentation for more details.
#' User should use the letterC function from seqLogo and should not use directly this version.
#' @keywords internal
letterC <- function(x.pos, y.pos, ht, wt, fill = "#2297E6", id = NULL) {
  angle1 <- seq(0.3 + pi / 2, pi, length = 100)
  angle2 <- seq(pi, 1.5 * pi, length = 100)
  x.l1 <- 0.5 + 0.5 * sin(angle1)
  y.l1 <- 0.5 + 0.5 * cos(angle1)
  x.l2 <- 0.5 + 0.5 * sin(angle2)
  y.l2 <- 0.5 + 0.5 * cos(angle2)
  x.l <- c(x.l1, x.l2)
  y.l <- c(y.l1, y.l2)
  x <- c(x.l, rev(x.l))
  y <- c(y.l, 1 - rev(y.l))
  x.i1 <- 0.5 + 0.35 * sin(angle1)
  y.i1 <- 0.5 + 0.35 * cos(angle1)
  x.i1 <- x.i1[y.i1 <= max(y.l1)]
  y.i1 <- y.i1[y.i1 <= max(y.l1)]
  y.i1[1] <- max(y.l1)
  x.i2 <- 0.5 + 0.35 * sin(angle2)
  y.i2 <- 0.5 + 0.35 * cos(angle2)
  x.i <- c(x.i1, x.i2)
  y.i <- c(y.i1, y.i2)
  x1 <- c(x.i, rev(x.i))
  y1 <- c(y.i, 1 - rev(y.i))
  x <- c(x, rev(x1))
  y <- c(y, rev(y1))
  x <- .rescale(x)
  x <- x.pos + wt * x
  y <- y.pos + ht * y
  if (is.null(id)) {
    id <- rep(1, length(x))
  }
  else {
    id <- rep(id, length(x))
  }
  fill <- rep(fill, 1)
  list(x = x, y = y, id = id, fill = fill)
}

#' letterT Function from seqLogo package (GloModAn)
#'
#' Internal function from seqLogo package. See seqLogo documentation for more details.
#' User should use the letterT function from seqLogo and should not use directly this version.
#' @keywords internal
letterT <- function(x.pos, y.pos, ht, wt, fill = "#DF536B", id = NULL) {
  x <- c(0, 10, 10, 6, 6, 4, 4, 0)
  y <- c(10, 10, 9, 9, 0, 0, 9, 9)
  x <- 0.1 * x
  y <- 0.1 * y
  x <- .rescale(x)
  x <- x.pos + wt * x
  y <- y.pos + ht * y
  if (is.null(id)) {
    id <- rep(1, 8)
  }
  else {
    id <- rep(id, 8)
  }
  fill <- rep(fill, 1)
  list(x = x, y = y, id = id, fill = fill)
}

#' letterU Function from seqLogo package (GloModAn)
#'
#' Internal function from seqLogo package. See seqLogo documentation for more details.
#' User should use the letterU function from seqLogo and should not use directly this version.
#' @keywords internal
letterU <- function(x.pos, y.pos, ht, wt, fill = "#DF536B", id = NULL) {
  angle1 <- seq(pi / 2, pi, length = 100)
  angle2 <- seq(pi, 1.5 * pi, length = 100)
  x.l1 <- 0.5 + 0.5 * sin(angle1)
  y.l1 <- 0.5 + 0.5 * cos(angle1)
  x.l2 <- 0.5 + 0.5 * sin(angle2)
  y.l2 <- 0.5 + 0.5 * cos(angle2)
  x.l <- c(x.l1, x.l2)
  y.l <- c(y.l1, y.l2)
  x.i1 <- 0.5 + 0.3 * sin(angle1)
  y.i1 <- 0.5 + 0.35 * cos(angle1)
  x.i1 <- x.i1[y.i1 <= max(y.l1)]
  y.i1 <- y.i1[y.i1 <= max(y.l1)]
  y.i1[1] <- max(y.l1)
  x.i2 <- 0.5 + 0.3 * sin(angle2)
  y.i2 <- 0.5 + 0.35 * cos(angle2)
  x.i <- c(x.i1, x.i2)
  y.i <- c(y.i1, y.i2)
  x <- c(x.l, 0, 0.2, 0.2, rev(x.i), 0.8, 0.8, 1, 1)
  y <- c(y.l, 1, 1, 0.5, rev(y.i), 0.5, 1, 1, 0.85)
  x <- .rescale(x)
  x <- x.pos + wt * x
  y <- y.pos + ht * y
  if (is.null(id)) {
    id <- rep(1, length(x))
  }
  else {
    id <- rep(id, length(x))
  }
  fill <- rep(fill, 1)
  list(x = x, y = y, id = id, fill = fill)
}

#' .rescale Function from seqLogo package (GloModAn)
#'
#' Internal function from seqLogo package. See seqLogo documentation for more details.
#' User should use the .rescale function from seqLogo and should not use directly this version.
#' @keywords internal
.rescale <- function(x, to = c(0.025, 0.975)) {
  from <- range(x, na.rm = TRUE, finite = TRUE)
  (x - from[1]) / diff(from) * diff(to) + to[1]
}

#' ExtractListModPosByModMotif Function (GloModAn)
#'
#' Return the GRanges object provided with the sequence associated to each position (and can also retrieve the sequence around each position).
#' @param grangesModPos A GRanges object containing Modifications Positions data to be extracted with the sequence.
#' @param grangesGenome A GRanges object containing the width of each contig.
#' @param dnastringsetGenome A DNAStringSet object containing the sequence for each contig.
#' @param nUpstreamBpToAdd Number of base pairs to add upstream of the range from the GRanges object provided
#' to obtain some sequence upstream of range. If some new ranges do not fit in the ranges of the contigs (provided with grangesGenome),
#' those new ranges will be removed. New windows with gaps are also removed. Defaults to 0.
#' @param nDownstreamBpToAdd Number of base pairs to add downstream of the range from the GRanges object provided
#' to obtain some sequence downstream of range. If some new ranges do not fit in the ranges of the contigs (provided with grangesGenome),
#' those new ranges will be removed. New windows with gaps are also removed. Defaults to 0.
#' @param nModMotifMinProp A number indicating the false discovery rate to be used for filtering: this will allow to choose the closest threshold
#' below this number. Defaults to 0.05 (so fdr of 5\%).
#' @param nModPositionInMotif The position of the modification in the window after resizing with nUpstreamBpToAdd and nDownstreamBpToAdd.
#' If GRanges are 1-bp positions, then 1+nUpstreamBpToAdd will return the right position of the modification.
#' Defaults to 1+nUpstreamBpToAdd.
#' @param cBaseLetterForMod The name of the base letter of the modified base.
#' @param cModNameInOutput Name for the modification in the output.
#' @return A list of 4 objects:
#' \describe{
#'   \item{motifs_to_analyse}{A character vector containing the sequence of motifs associated to the modification.}
#'   \item{mod_motif}{A character vector containing the sequence of motifs associated to the modification with the
#'   modification represented inside those motifs.}
#'   \item{motif_pct}{A table containing the percentage of modifications in each motif tested.}
#'   \item{GRangesbyMotif}{A list of GRanges objects with the sequence: one GRanges object by motif associated to the modification.}
#' }
#' @keywords ExtractListModPosByModMotif
#' @import GenomicRanges
#' @importFrom S4Vectors isEmpty
#' @export
#' @examples
#' # loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(
#'   package = "DNAModAnnot", "extdata",
#'   "ptetraurelia_mac_51_sca171819.fa"
#' ))
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' # Preparing a grangesPacBioGFF dataset
#' myGrangesPacBioGFF <-
#'   ImportPacBioGFF(
#'     cPacBioGFFPath = system.file(
#'       package = "DNAModAnnot", "extdata",
#'       "ptetraurelia.modifications.sca171819.gff"
#'     ),
#'     cNameModToExtract = "m6A",
#'     cModNameInOutput = "6mA",
#'     cContigToBeAnalyzed = names(myGenome)
#'   )
#'
#' # Retrieve GRanges with sequence
#' myMotif_pct_and_GRangesList <- ExtractListModPosByModMotif(
#'   grangesModPos = myGrangesPacBioGFF,
#'   grangesGenome = myGrangesGenome,
#'   dnastringsetGenome = myGenome,
#'   nUpstreamBpToAdd = 0,
#'   nDownstreamBpToAdd = 1,
#'   nModMotifMinProp = 0.05,
#'   cBaseLetterForMod = "A",
#'   cModNameInOutput = "6mA"
#' )
#'
#' myMotif_pct_and_GRangesList$motifs_to_analyse
#' myMotif_pct_and_GRangesList$mod_motif
#' myMotif_pct_and_GRangesList$motif_pct
#' myMotif_pct_and_GRangesList$GRangesbyMotif
ExtractListModPosByModMotif <- function(grangesModPos,
                                        grangesGenome,
                                        dnastringsetGenome,
                                        nUpstreamBpToAdd = 0, nDownstreamBpToAdd = 1,
                                        nModMotifMinProp,
                                        nModPositionInMotif = 1 + nUpstreamBpToAdd,
                                        cBaseLetterForMod,
                                        cModNameInOutput) {
  ans <- GetGRangesWindowSeqandParam(
    grangesData = grangesModPos,
    grangesGenome = grangesGenome,
    dnastringsetGenome = dnastringsetGenome,
    nUpstreamBpToAdd = nUpstreamBpToAdd, nDownstreamBpToAdd = nDownstreamBpToAdd
  )

  motif_pct <- table(ans$sequence) / sum(table(ans$sequence))
  if (isEmpty(which(motif_pct >= nModMotifMinProp))) {
    print(paste("No ", unique(width(ans$sequence)), "nt motif is represented at least ", nModMotifMinProp * 100, "% of 6mA detected.", sep = ""))
    print("Motif width or minimum proportion might be too high.")
  } else {
    motif_to_analyze <- names(motif_pct[motif_pct >= nModMotifMinProp])
    ans <- ans[which(as.character(ans$sequence) %in% motif_to_analyze), ]
    ans <- split(ans, as.factor(ans$sequence))
    ans <- ans[lengths(ans) > 0]
  }
  mod_motif <- .IncludeModPosInMot(motif_pct, cBaseLetterForMod, cModNameInOutput, nModPositionInMotif)

  return(list(
    motifs_to_analyse = names(motif_pct[motif_pct >= nModMotifMinProp]),
    mod_motif = names(mod_motif[motif_pct >= nModMotifMinProp]),
    motif_pct = motif_pct,
    GRangesbyMotif = ans
  ))
}
