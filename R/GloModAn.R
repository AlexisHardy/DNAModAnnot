### =========================================================================
### Functions from ModAn and ModAnPlot categories
### -------------------------------------------------------------------------
###
### DNAModAnnot v.0.0.0.9001 - 2020/04/04
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
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#'
#' #Preparing a gposPacBioCSV dataset
#' myGposPacBioCSV <-
#'    ImportPacBioCSV(cPacBioCSVPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.combinationC.all.corrected_sca171819_sample.csv"),
#'                    cSelectColumnsToExtract = c("refName", "tpl", "strand", "base",
#'                                                "score", "ipdRatio", "coverage"),
#'                    lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'                    cContigToBeAnalyzed = names(myGenome))
#'
#' myMean_cov_list <- GetMeanParamByContig(grangesData = myGposPacBioCSV,
#'                                             dnastringsetGenome = myGenome,
#'                                             cParamName = "coverage")
#' myMean_cov_list
GetMeanParamByContig <- function(grangesData,
                                 dnastringsetGenome,
                                 cParamName) {
  count_table <- data.frame(refName=as.factor(seqnames(grangesData)),
                            strand=as.factor(strand(grangesData)),
                            parameter=mcols(grangesData)[[cParamName]])
  mean_table <- aggregate(count_table$parameter,
                          list(count_table$refName, count_table$strand),
                          mean)
  colnames(mean_table) <- c("refName", "strand", paste("mean_", cParamName, sep="") )

  mean_table <- merge(mean_table,
                      data.frame(refName = names(dnastringsetGenome),
                                 width = width(dnastringsetGenome) ),
                      by = "refName",
                      all.y = TRUE)

  #those who have NA for one or both strands : add missing strands with 0 mean_param
  if(length(which(is.na(mean_table$strand)))>0){
    mean_table[is.na(mean_table$strand), paste("mean_", cParamName, sep="")] <- 0
    mean_table[is.na(mean_table$strand), "strand"] <- "+"
  }
  if(length(which(isUnique(mean_table$refName)))>0){
    mean_table2 <- mean_table[isUnique(mean_table$refName),]
    mean_table2$strand <- ifelse(mean_table2$strand == "+", "-", "+")
    mean_table2[,paste("mean_", cParamName, sep="")] <- 0
    mean_table <- rbind(mean_table, mean_table2)
  }
  mean_table <- mean_table[with(mean_table, order(width, refName, strand, decreasing = TRUE)),]

  f_strand <- subset(mean_table, strand == "+")
  r_strand <- subset(mean_table, strand == "-")

  return(list(mean_table=mean_table,
              f_strand=f_strand,
              r_strand=r_strand))
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
#' DrawBarplotBothStrands(nParamByContigForward = c(100,86,75,56),
#'                               nParamByContigReverse = c(96,88,80,83),
#'                               cContigNames = c("chrI","chrII","chrIII","chrIV"),
#'                               cGraphName = "Mean Coverage per contig",
#'                               lIsOrderedLargestToSmallest = TRUE)
DrawBarplotBothStrands <- function(nParamByContigForward,
                                   nParamByContigReverse,
                                   cContigNames,
                                   cGraphName,
                                   lIsOrderedLargestToSmallest = TRUE){
  layout( matrix(c(4, 4,
                   3, 1,
                   3, 2),
                 nrow = 3, ncol = 2, byrow = TRUE), widths = c(1,29), heights = c(1,2,2) )

  nMarLeft <- max(nchar(c(pretty(nParamByContigForward), pretty(nParamByContigReverse))))*1.1 +1
  opar <- par()
  par(mar = c(0, nMarLeft, 4, 2) )

  ylim_max <- max(  max(nParamByContigForward), max(nParamByContigReverse)  )

  bp <- barplot(nParamByContigForward,
                ylab = "Forward strand", ylim = c(0, ylim_max),
                col = c("grey30", "grey70"),
                las = 2, xaxs = "i", yaxs = "i", cex.axis=1.25, cex.lab=1.5, mgp=c(nMarLeft*3.5/5,1,0)
  )

  axis(side=3, pos=NA, at=bp, labels=cContigNames, las=2, cex.axis=1.25)

  par(mar = c(4, nMarLeft, 0, 2)  )

  barplot(nParamByContigReverse,
          ylab = "Reverse strand", ylim = par("usr")[c(4,3)],
          col = c("grey30", "grey70"),
          las = 1, xaxs = "i", yaxs = "i", cex.axis=1.25, cex.lab=1.5, mgp=c(nMarLeft*3.5/5,1,0)
  )
  par(xpd = NA)
  if(lIsOrderedLargestToSmallest){
    mtext("Contigs (from largest to smallest)", cex = 1.25, side = 1, line = 1.5)
  } else {
    mtext("Contigs", cex = 1.25, side = 1, line = 1.5)
  }
  par(xpd = FALSE)

  par(mar = c(0, 0, 0, 0.5)  )
  barplot(height=1, col="white", border="white", axes = FALSE)

  par(xpd = NA)
  mtext(cGraphName, cex = 1.25, side = 4, srt = 90)
  par(xpd = FALSE)

  layout(matrix(c(1,1,1,1), nrow = 2, ncol = 2))
  par(mar=opar$mar)
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
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#'
#' #Preparing a gposPacBioCSV dataset
#' myGposPacBioCSV <-
#'    ImportPacBioCSV(cPacBioCSVPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.combinationC.all.corrected_sca171819_sample.csv"),
#'                    cSelectColumnsToExtract = c("refName", "tpl", "strand", "base",
#'                                                "score", "ipdRatio", "coverage"),
#'                    lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'                    cContigToBeAnalyzed = names(myGenome))
#'
#' DrawDistriHistBox(head(myGposPacBioCSV$coverage, 100000),
#'                          cGraphName = "Coverage distribution of some bases sequenced",
#'                          cParamName = "Coverage",
#'                          lTrimOutliers = TRUE )
DrawDistriHistBox <- function(nParam,
                              cGraphName,
                              cParamName,
                              lTrimOutliers = FALSE,
                              nXLimits = NULL){
  d <- density(nParam)
  h <- hist(nParam, plot = FALSE)
  if( length( seq(from=min(h$breaks), to=max(h$breaks), by=1) ) < 10) {
    breaks_v = seq(from=min(h$breaks), to=max(h$breaks), length.out = 200)
  } else {
    breaks_v = seq(from=min(h$breaks), to=max(h$breaks), by=1)
  }
  h <- hist(nParam, breaks = breaks_v, plot = FALSE)

  if(lTrimOutliers){ b <- boxplot(nParam, plot=FALSE) }

  if(is.null(nXLimits)) {
    nXLimits <- c()
    nXLimits[1] <- ifelse(!lTrimOutliers, min(h$breaks), b$stats[1] )
    nXLimits[2] <- ifelse(!lTrimOutliers, max(h$breaks), b$stats[5] )
  }

  opar <- par()
  layout(matrix(c(1,2),2,1), heights = c(1,9))
  par(mar = c(0.1, 5.1, 2.1, 2.1) )
  b <- boxplot(nParam, horizontal=TRUE, xaxt="n", frame=FALSE, outline=!lTrimOutliers,
               main= ifelse(!lTrimOutliers,
                            cGraphName,
                            paste(cGraphName, " (no outliers)",sep="") ),
               outcol= rgb(0,0,0,0.01),
               ylim=nXLimits
  )

  par(mar = c(5, 5.1, 0.1, 2.1) )
  hist(nParam, freq = FALSE, probability = TRUE,
       xlab = cParamName, ylab = "Density (%)", main = NULL,
       col="gray",
       xlim = nXLimits,
       ylim = c(0, ifelse(max(h$density) < max(d$y), max(d$y), max(h$density)) ),
       breaks = breaks_v
  )
  lines(d, col="red")
  layout(matrix(c(1),1,1))
  par(mar=opar$mar)
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
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' #Preparing a grangesPacBioGFF datasets
#' myGrangesPacBioGFF <-
#'    ImportPacBioGFF(cPacBioGFFPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.combinationC_merged_sca171819_sample.gff"),
#'                    cNameModToExtract = "m6A",
#'                    cModNameInOutput = "6mA",
#'                    cContigToBeAnalyzed = names(myGenome))
#'
#' myGposPacBioCSV <-
#'    ImportPacBioCSV(cPacBioCSVPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.combinationC.all.corrected_sca171819_sample.csv"),
#'                    cSelectColumnsToExtract = c("refName", "tpl", "strand", "base",
#'                                                "score", "ipdRatio", "coverage"),
#'                    lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <- myGposPacBioCSV[myGposPacBioCSV$base == "A"]
#'
#' #Mod report
#' myReport_Mod <- GetModReportPacBio(grangesGenome = myGrangesGenome,
#'                                  dnastringsetGenome = myGenome,
#'                                  grangesPacBioGFF = myGrangesPacBioGFF,
#'                                  gposPacBioCSVBase = myGposPacBioCSV,
#'                                  cOrgAssemblyName="ptetraurelia_mac_51",
#'                                  cBaseLetterForMod = "A",
#'                                  cModNameInOutput = "6mA")
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
    nbBaseAssembly = letterFrequency(dnastringsetGenome, letters=cBaseLetterForMod, collapse=TRUE) +
      letterFrequency(reverseComplement(dnastringsetGenome), letters=cBaseLetterForMod, collapse=TRUE)
  )

  dReportTable$ratioMod = dReportTable$nbMod / dReportTable$nbBaseSequenced
  dReportTable$ratioMod_corrected = dReportTable$ratioMod * mean(grangesPacBioGFF$frac)
  dReportTable$ratioModAllAssemblyBases = dReportTable$nbMod / dReportTable$nbBaseAssembly
  dReportTable$ratioModAllAssemblyBases_corrected = dReportTable$ratioModAllAssemblyBases * mean(grangesPacBioGFF$frac)
  dReportTable$mean_frac_Mod = mean(grangesPacBioGFF$frac)
  dReportTable$mean_coverage_Mod = mean(grangesPacBioGFF$coverage)
  dReportTable$mean_score_Mod = mean(grangesPacBioGFF$score)
  dReportTable$mean_ipdRatio_Mod = mean(grangesPacBioGFF$ipdRatio)
  dReportTable$mean_identificationQv_Mod = mean(grangesPacBioGFF$identificationQv)

  colnames(dReportTable) <- gsub(pattern="Base", replacement = cBaseLetterForMod,
                                 x = colnames(dReportTable))
  colnames(dReportTable) <- gsub(pattern="Mod", replacement = cModNameInOutput,
                                 x = colnames(dReportTable))

  dReportTable <- .AddModMotifPctToDf(dReportTable,
                                      grangesModPos = grangesPacBioGFF, grangesGenome=grangesGenome,
                                      dnastringsetGenome=dnastringsetGenome,
                                      nModPositionInMotif = 1, nUpstreamBpToAdd=0, nDownstreamBpToAdd=1,
                                      cBaseLetterForMod=cBaseLetterForMod, cModNameInOutput=cModNameInOutput)


  dReportTable <- .AddModMotifPctToDf(dReportTable,
                                      grangesModPos = grangesPacBioGFF, grangesGenome=grangesGenome,
                                      dnastringsetGenome=dnastringsetGenome,
                                      nModPositionInMotif = 2, nUpstreamBpToAdd=1, nDownstreamBpToAdd=0,
                                      cBaseLetterForMod=cBaseLetterForMod, cModNameInOutput=cModNameInOutput)


  dReportTable <- .AddModMotifPctToDf(dReportTable,
                                      grangesModPos = grangesPacBioGFF, grangesGenome=grangesGenome,
                                      dnastringsetGenome=dnastringsetGenome,
                                      nModPositionInMotif = 2, nUpstreamBpToAdd=1, nDownstreamBpToAdd=1,
                                      cBaseLetterForMod=cBaseLetterForMod, cModNameInOutput=cModNameInOutput)

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
#' #preparing genome (simulated)
#' myGenome <- Biostrings::DNAStringSet(paste0(rep("ATCG", 100000), collapse = ""))
#' names(myGenome) <- "NC_000001.11"
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' #Loading Nanopore data
#' myDeepSignalModPath <- system.file(
#'   package="DNAModAnnot", "extdata",
#'   "Notts.FAB39088-288418386-NC_000001.11_Sample.CpG.call_mods.frequency.tsv")
#' mygposDeepSignalModBase <- ImportDeepSignalModFrequency(cDeepSignalModPath=myDeepSignalModPath,
#'                                                         lSortGPos=TRUE,
#'                                                         cContigToBeAnalyzed = "all")
#'
#' #Filtering
#' mygposDeepSignalMod <- FiltDeepSignal(gposDeepSignalModBase = mygposDeepSignalModBase,
#'                                       cParamNameForFilter = "frac",
#'                                       lFiltParam = TRUE,
#'                                       nFiltParamLoBoundaries = 0,
#'                                       nFiltParamUpBoundaries = 1,
#'                                       cFiltParamBoundariesToInclude = "upperOnly")$Mod
#'
#' #Mod report
#' myReport_Mod <- GetModReportDeepSignal(
#'   dnastringsetGenome = myGenome,
#'   grangesGenome = myGrangesGenome,
#'   gposDeepSignalMod = as(mygposDeepSignalMod, "GRanges"),
#'   gposDeepSignalModBase = as(mygposDeepSignalModBase, "GRanges"),
#'   cOrgAssemblyName = "Test_function",
#'   cBaseLetterForMod = "C", cModNameInOutput = "5mC")
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
    nbBaseAssembly = letterFrequency(dnastringsetGenome, letters=cBaseLetterForMod, collapse=TRUE) +
      letterFrequency(reverseComplement(dnastringsetGenome), letters=cBaseLetterForMod, collapse=TRUE)
  )

  dReportTable$ratioMod = dReportTable$nbMod / dReportTable$nbSitesSequenced
  dReportTable$ratioMod_corrected = dReportTable$ratioMod * mean(gposDeepSignalMod$frac)
  dReportTable$ratioModAllAssemblyBases = dReportTable$nbMod / dReportTable$nbBaseAssembly
  dReportTable$ratioModAllAssemblyBases_corrected = dReportTable$ratioModAllAssemblyBases * mean(gposDeepSignalMod$frac)

  dReportTable$mean_prob1sum_Mod = mean(gposDeepSignalMod$prob_1_sum)
  dReportTable$mean_countModified_Mod = mean(gposDeepSignalMod$count_modified)
  dReportTable$mean_coverage_Mod = mean(gposDeepSignalMod$coverage)
  dReportTable$mean_frac_Mod = mean(gposDeepSignalMod$frac)

  colnames(dReportTable) <- gsub(pattern="Base", replacement = cBaseLetterForMod,
                                 x = colnames(dReportTable))
  colnames(dReportTable) <- gsub(pattern="Mod", replacement = cModNameInOutput,
                                 x = colnames(dReportTable))

  dReportTable <- .AddModMotifPctToDf(dReportTable,
                                      grangesModPos = gposDeepSignalMod, grangesGenome=grangesGenome,
                                      dnastringsetGenome=dnastringsetGenome,
                                      nModPositionInMotif = 1, nUpstreamBpToAdd=0, nDownstreamBpToAdd=1,
                                      cBaseLetterForMod=cBaseLetterForMod, cModNameInOutput=cModNameInOutput)


  dReportTable <- .AddModMotifPctToDf(dReportTable,
                                      grangesModPos = gposDeepSignalMod, grangesGenome=grangesGenome,
                                      dnastringsetGenome=dnastringsetGenome,
                                      nModPositionInMotif = 2, nUpstreamBpToAdd=1, nDownstreamBpToAdd=0,
                                      cBaseLetterForMod=cBaseLetterForMod, cModNameInOutput=cModNameInOutput)


  dReportTable <- .AddModMotifPctToDf(dReportTable,
                                      grangesModPos = gposDeepSignalMod, grangesGenome=grangesGenome,
                                      dnastringsetGenome=dnastringsetGenome,
                                      nModPositionInMotif = 2, nUpstreamBpToAdd=1, nDownstreamBpToAdd=1,
                                      cBaseLetterForMod=cBaseLetterForMod, cModNameInOutput=cModNameInOutput)

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
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' #Preparing a grangesPacBioGFF datasets
#' myGrangesPacBioGFF <-
#'    ImportPacBioGFF(cPacBioGFFPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.combinationC_merged_sca171819_sample.gff"),
#'                    cNameModToExtract = "m6A",
#'                    cModNameInOutput = "6mA",
#'                    cContigToBeAnalyzed = names(myGenome))
#'
#' #Retrieve GRanges with sequence
#' myPositions_Mod_Granges_wSeq <- GetGRangesWindowSeqandParam(grangesData=myGrangesPacBioGFF,
#'                                                             grangesGenome = myGrangesGenome,
#'                                                             dnastringsetGenome = myGenome,
#'                                                             nUpstreamBpToAdd = 5,
#'                                                             nDownstreamBpToAdd = 5)
#' myPositions_Mod_Granges_wSeq
GetGRangesWindowSeqandParam <- function(grangesData,
                                               grangesGenome,
                                               dnastringsetGenome,
                                               nUpstreamBpToAdd=0, nDownstreamBpToAdd=0){

  ans <- grangesData

  print("Window adjustment...")
  ans <- resize(ans, width = width(ans)+nUpstreamBpToAdd, fix = 'end')
  ans  <- resize(ans, width = width(ans)+nDownstreamBpToAdd, fix = 'start')

  # remove windows not included completely in contig windows
  gr_a <- grangesGenome
  print("Removing Windows not included completely in contig windows...")
  ans <- subsetByOverlaps(ans, gr_a, ignore.strand = FALSE, type="within")

  print("Getting sequence...")
  ans$sequence <- as.factor(as.character(getSeq(dnastringsetGenome, ans)))

  print("Window adjustment...")
  ans <- resize(ans, width = width(ans)-nUpstreamBpToAdd, fix = 'end')
  ans  <- resize(ans, width = width(ans)-nDownstreamBpToAdd, fix = 'start')

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
                                         nUpstreamBpToAdd=0,
                                         nDownstreamBpToAdd=0){
  positions_Mod_Granges_wSeq <- GetGRangesWindowSeqandParam(grangesModPos, grangesGenome, dnastringsetGenome,
                                                                   nUpstreamBpToAdd, nDownstreamBpToAdd)

  motif_pct <- as.data.frame.matrix(
    t(table(positions_Mod_Granges_wSeq$sequence))/sum(t(table(positions_Mod_Granges_wSeq$sequence)))
    )

  motif_pct <- .IncludeModPosInMot(motif_pct, cBaseLetterForMod, cModNameInOutput, nModPositionInMotif)

  names(motif_pct) <- paste("pct_", names(motif_pct), sep="")

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
                                       lExtractOnlyMotifWithMod = TRUE){
  pattern=paste("^(.{",(nModPositionInMotif-1),"})",cBaseLetterForMod,"(.*)$",sep="")
  replacement = paste("\\1",cModNameInOutput,"\\2",sep="")
  names(nParamByMotif) <- gsub(pattern=pattern, replacement = replacement,
                               x = names(nParamByMotif))
  if(lExtractOnlyMotifWithMod) {
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
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#'
#' #Preparing a grangesPacBioGFF and gposPacBioCSV datasets
#' myGrangesPacBioGFF <-
#'    ImportPacBioGFF(cPacBioGFFPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.combinationC_merged_sca171819_sample.gff"),
#'                    cNameModToExtract = "m6A",
#'                    cModNameInOutput = "6mA",
#'                    cContigToBeAnalyzed = names(myGenome))
#'
#' myGposPacBioCSV <-
#'    ImportPacBioCSV(cPacBioCSVPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.combinationC.all.corrected_sca171819_sample.csv"),
#'                    cSelectColumnsToExtract = c("refName", "tpl", "strand", "base",
#'                                                "score", "ipdRatio", "coverage"),
#'                    lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <- myGposPacBioCSV[myGposPacBioCSV$base == "A"]
#'
#' #Mod report
#' myMod_ratio_list <- GetModRatioByContig(grangesModPos = myGrangesPacBioGFF,
#'                                         gposModTargetBasePos = myGposPacBioCSV,
#'                                         dnastringsetGenome = myGenome,
#'                                         cBaseLetterForMod = "A")
#' myMod_ratio_list
GetModRatioByContig <- function(grangesModPos,
                                gposModTargetBasePos,
                                dnastringsetGenome,
                                cBaseLetterForMod) {

  mod_count_table <- data.frame(refName=as.factor(seqnames(grangesModPos)),
                                strand=as.factor(strand(grangesModPos))  )
  mod_count_table <- as.data.frame(table(mod_count_table))
  colnames(mod_count_table) <- c("refName", "strand", "Mod_count")

  baseS_count_table <- data.frame(refName=as.factor(seqnames(gposModTargetBasePos)),
                                  strand=as.factor(strand(gposModTargetBasePos))  )
  baseS_count_table <- as.data.frame(table(baseS_count_table))
  colnames(baseS_count_table) <- c("refName", "strand", "BaseSequenced_count")

  dnastringsetGenome_order <- dnastringsetGenome[order(names(dnastringsetGenome)),]
  baseA_count_table <- data.frame(names(dnastringsetGenome_order), "+",
                                  letterFrequency(dnastringsetGenome_order, letters=cBaseLetterForMod),
                                  width(dnastringsetGenome_order))
  tmp_table <- data.frame(names(dnastringsetGenome_order), "-",
                          letterFrequency(reverseComplement(dnastringsetGenome_order), letters=cBaseLetterForMod),
                          width(dnastringsetGenome_order))
  baseA_count_table <- rbind(baseA_count_table, tmp_table)
  colnames(baseA_count_table) <- c("refName", "strand", "BaseAssembly_count", "width")

  base_count_table <- merge(baseS_count_table,
                            baseA_count_table,
                            by = c("refName", "strand"),
                            all.y = TRUE)

  mod_count_table <- merge(mod_count_table,
                           base_count_table,
                           by = c("refName", "strand"),
                           all.y = TRUE)

  #those who have NA for one or both strands : replace NA by 0
  mod_count_table[is.na(mod_count_table$Mod_count),"Mod_count"] <- 0
  mod_count_table[is.na(mod_count_table$BaseSequenced_count),"BaseSequenced_count"] <- 0
  mod_count_table[is.na(mod_count_table$BaseAssembly_count),"BaseAssembly_count"] <- 0

  mod_count_table <- mod_count_table[with(mod_count_table, order(width, refName, strand, decreasing = TRUE)),]

  mod_count_table$Mod_ratio <- mod_count_table$Mod_count / mod_count_table$BaseSequenced_count
  mod_count_table$Mod_ratio_AllAssemblyBases <- mod_count_table$Mod_count / mod_count_table$BaseAssembly_count

  f_strand <- subset(mod_count_table, strand == "+")
  r_strand <- subset(mod_count_table, strand == "-")

  return(list(mod_count_table=mod_count_table,
              f_strand=f_strand,
              r_strand=r_strand))
}

#' DrawModLogo Function (GloModAn)
#'
#' Return a plot describing the sequence motif associated to the sequences provided.
#' Sequences that do not have full width and sequences that have some
#' N or some gaps "-" are automatically removed before drawing the sequence plot.
#' If a position is fixed with one base only, that position will not be corrected with
#' the background in order to avoid other letters to appear below this base if using cLogoType="EDLogo".
#' @param dnastringsetSeqAroundMod A DNAStringSet object containing the sequence around each DNA modification.
#' All these sequences must have the same size and the modification must have the same position in each sequence (example: at the center).
#' @param cLogoType Can be "Logo" or "EDLogo". Defaults to "EDLogo".
#' \itemize{
#'   \item Logo: A classic logo showing sequence enrichment around the modification.
#'   \item EDLogo: A logo showing sequence enrichment and depletion around the modification.
#' }
#' @param nGenomicBgACGT A numeric vector giving the background to be corrected with the genomic composition
#' in Adenine (A) then Cytosine (C) then Guanine (G) then Thymine (T). Defaults to c(0.25, 0.25, 0.25, 0.25).
#' @param cColorsCTAG A character vector giving the color for the Cytosine (C) then Thymine (T) then
#' Adenine (A) then Guanine (G) letters on the plot. Defaults to c("orange2","green3","red2","blue2").
#' @keywords DrawModLogo
#' @importFrom Logolas logomaker
#' @importFrom Biostrings consensusMatrix
#' @export
#' @examples
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' #Preparing a grangesPacBioGFF datasets
#' myGrangesPacBioGFF <-
#'    ImportPacBioGFF(cPacBioGFFPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.combinationC_merged_sca171819_sample.gff"),
#'                                  cNameModToExtract = "m6A",
#'                                  cModNameInOutput = "6mA",
#'                                  cContigToBeAnalyzed = names(myGenome))
#'
#' #Retrieve GRanges with sequence
#' myPositions_Mod_Granges_wSeq <- GetGRangesWindowSeqandParam(grangesData=myGrangesPacBioGFF,
#'                                                             grangesGenome = myGrangesGenome,
#'                                                             dnastringsetGenome = myGenome,
#'                                                             nUpstreamBpToAdd = 5,
#'                                                             nDownstreamBpToAdd = 5)
#'
#' DrawModLogo(dnastringsetSeqAroundMod = as(myPositions_Mod_Granges_wSeq$sequence,"DNAStringSet"),
#'                    cLogoType = "EDLogo",
#'                    nGenomicBgACGT = c(0.35, 0.15, 0.15, 0.35))
DrawModLogo <- function(dnastringsetSeqAroundMod,
                               cLogoType = "EDLogo",
                               nGenomicBgACGT=c(0.25, 0.25, 0.25, 0.25),
                               cColorsCTAG=c("orange2","green3","red2","blue2")){
  names(nGenomicBgACGT) <- c("A", "C", "G", "T")

  mConsensus <- consensusMatrix(dnastringsetSeqAroundMod, as.prob=TRUE)
  nFixedPositions <- which(apply(mConsensus, 2, function(x) any(x == 1)))

  #remove all sequences that do not have full width
  dnastringsetSeqAroundMod <- dnastringsetSeqAroundMod[which( width(dnastringsetSeqAroundMod) == max(width(dnastringsetSeqAroundMod)) )]
  #remove all sequences that have some N or some gaps "-"
  dnastringsetSeqAroundMod <- dnastringsetSeqAroundMod[which(letterFrequency(dnastringsetSeqAroundMod, letters="N-" ) == 0)]

  bgmatrix <- matrix(rep(nGenomicBgACGT, max(width(dnastringsetSeqAroundMod)) ),
                     ncol = max(width(dnastringsetSeqAroundMod)), byrow = FALSE)
  #remove background correction for fixed positions (example: Position Of The Modification)
  bgmatrix[,nFixedPositions] <- c(0.25,0.25,0.25,0.25)

  logomaker(as.character(dnastringsetSeqAroundMod),
            type = cLogoType,
            color_type = "per_row",
            colors = cColorsCTAG,
            bg=bgmatrix)

  # #alternative method
  # c_m <- consensusMatrix(dnastringsetSeqAroundMod)
  # f_m_corrected <- prop.table(c_m, 2)/nGenomicBgACGT
  # logomaker(f_m_corrected,
  #           type = cLogoType,
  #           color_type = "per_row",
  #           colors = cColorsCTAG)

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
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' #Preparing a grangesPacBioGFF dataset
#' myGrangesPacBioGFF <-
#'    ImportPacBioGFF(cPacBioGFFPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.combinationC_merged_sca171819_sample.gff"),
#'                    cNameModToExtract = "m6A",
#'                    cModNameInOutput = "6mA",
#'                    cContigToBeAnalyzed = names(myGenome))
#'
#' #Retrieve GRanges with sequence
#' myMotif_pct_and_GRangesList <- ExtractListModPosByModMotif(grangesModPos=myGrangesPacBioGFF,
#'                                                         grangesGenome = myGrangesGenome,
#'                                                         dnastringsetGenome = myGenome,
#'                                                         nUpstreamBpToAdd = 0,
#'                                                         nDownstreamBpToAdd = 1,
#'                                                         nModMotifMinProp = 0.05,
#'                                                         cBaseLetterForMod = "A",
#'                                                         cModNameInOutput = "6mA")
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
                                             nModPositionInMotif = 1+nUpstreamBpToAdd,
                                             cBaseLetterForMod,
                                             cModNameInOutput) {
  ans <- GetGRangesWindowSeqandParam(grangesData = grangesModPos,
                                            grangesGenome=grangesGenome,
                                            dnastringsetGenome = dnastringsetGenome,
                                            nUpstreamBpToAdd=nUpstreamBpToAdd, nDownstreamBpToAdd=nDownstreamBpToAdd)

  motif_pct <- table(ans$sequence)/sum(table(ans$sequence))
  if( isEmpty(which(motif_pct >= nModMotifMinProp)) ){
    print(paste("No ", unique(width(ans$sequence)),"nt motif is represented at least ",nModMotifMinProp*100,"% of 6mA detected.",sep=""))
    print("Motif width or minimum proportion might be too high.")
  } else {
    motif_to_analyze <- names(motif_pct[motif_pct >= nModMotifMinProp])
    ans <- ans[which(as.character(ans$sequence) %in% motif_to_analyze), ]
    ans <- split(ans, as.factor(ans$sequence))
    ans <- ans[lengths(ans) > 0]
  }
  mod_motif <- .IncludeModPosInMot(motif_pct, cBaseLetterForMod, cModNameInOutput, nModPositionInMotif)

  return(list(motifs_to_analyse = names(motif_pct[motif_pct >= nModMotifMinProp]),
              mod_motif         = names(mod_motif[motif_pct >= nModMotifMinProp]),
              motif_pct         = motif_pct,
              GRangesbyMotif    = ans))
}
