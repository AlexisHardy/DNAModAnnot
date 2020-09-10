### =========================================================================
### Functions from ModAnnot
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
### GetModBaseCountsByFeature
### DrawModBasePropByFeature
### .GetModBaseCountsCategories
### DrawParamPerModBaseCategories
### .DrawPolygonRange
### GetModBaseCountsWithinFeature
### DrawModBaseCountsWithinFeature
### GetDistFromFeaturePos
### .GetPosDistFromFeaturePos
### GetListCountsByDist
### .GetCountsByDist
### DrawModBasePropDistFromFeature
### AddToModBasePropDistFromFeaturePlot
### AdaptedIdeogramTrackWithoutBandsData
### ExportFilesForGViz
### .ExportFastaForGviz
### .ExportBigwigForGviz
### .ExportBamForGviz
### ImportBamExtendedAnnotationTrack
###
### -------------------------------------------------------------------------

#' GetModBaseCountsByFeature Function (ModAnnot)
#'
#' Return the Annotation provided with the counts (and counts per KiloBase pairs (kbp)) of
#' the base modified (Mod) and the base letter of the modified base (Base).
#' Example: for Mod="6mA", Base="A"; for Mod="5mC", Base="C".
#' @param grangesAnnotations A GRanges object containing the annotation for the genome assembly
#' corresponding to the grangesModPos and gposModTargetBasePos provided. The Genomic features categories must be in a column named "type".
#' @param grangesModPos A GRanges object containing Modifications Positions data to be counted.
#' @param gposModTargetBasePos A GRanges or GPos object containing Base Positions (which can be targeted by the modification) to be counted.
#' @param lIgnoreStrand If TRUE, Mod and Base will be counted independently of the strand of each feature. If FALSE, only Mod and Base
#' on the same strand as the feature will be counted. Defaults to FALSE.
#' @return A GRanges object based on grangesAnnotations with the counts:
#' \itemize{
#'   \item Modcount:           The number of "Mod"  within this feature.
#'   \item Modcount_perkbp:   (The number of "Mod"  within this feature divided by the size of the feature)*1000.
#'   \item Basecount:          The number of "Base" within this feature.
#'   \item Basecount_perkbp:  (The number of "Base" within this feature divided by the size of the feature)*1000.
#' }
#' @keywords GetModBaseCountsByFeature
#' @importFrom S4Vectors countQueryHits
#' @export
#' @examples
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#'
#' #loading annotation
#' library(rtracklayer)
#' myAnnotations <- readGFFAsGRanges(system.file(package="DNAModAnnot", "extdata",
#'                                           "ptetraurelia_mac_51_annotation_v2.0_sca171819.gff3"))
#'
#' #Preparing a grangesPacBioGFF and a grangesPacBioCSV datasets
#' myGrangesPacBioGFF <-
#'    ImportPacBioGFF(cPacBioGFFPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.modifications.sca171819.gff"),
#'                    cNameModToExtract = "m6A",
#'                    cModNameInOutput = "6mA",
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <-
#'    ImportPacBioCSV(cPacBioCSVPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.bases.sca171819.csv"),
#'                    cSelectColumnsToExtract = c("refName", "tpl", "strand", "base",
#'                                                "score", "ipdRatio", "coverage"),
#'                    lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'                    cContigToBeAnalyzed = names(myGenome))
#'
#' #Retrieve annotations with "Mod" and "Base" counts (and counts per kbp)
#' myAnn_ModBase_counts <- GetModBaseCountsByFeature(grangesAnnotations=myAnnotations,
#'                                                       grangesModPos=myGrangesPacBioGFF,
#'                                                       gposModTargetBasePos=myGposPacBioCSV)
#' myAnn_ModBase_counts
GetModBaseCountsByFeature <- function(grangesAnnotations,
                                          grangesModPos,
                                          gposModTargetBasePos,
                                          lIgnoreStrand = FALSE){
  gr_a <- grangesAnnotations
  features_count <- as.data.frame(table(gr_a$type))
  features_count <- features_count[order(features_count$Var1, decreasing = FALSE),]
  gr_b <- grangesModPos

  hits <- findOverlaps(gr_a, gr_b, ignore.strand = lIgnoreStrand)
  ans <- gr_a

  ans$Modcount <- countQueryHits(hits)
  ans$Modcount_perkbp <- 1000*ans$Modcount/width(ans)

  gr_a <- ans
  gr_b <- gposModTargetBasePos

  hits <- findOverlaps(gr_a, gr_b, ignore.strand = lIgnoreStrand)
  ans <- gr_a
  ans$Basecount <- countQueryHits(hits)
  ans$Basecount_perkbp <- 1000*ans$Basecount/width(ans)

  return(ans)
}

#' DrawModBasePropByFeature Function (ModAnnot)
#'
#' Return a plot describing, between the features to be displayed, the proportion
#' (or proportion per KiloBase pairs (kbp)) of
#' the base modified (Mod) and the base letter of the modified base (Base).
#' Example: for Mod="6mA", Base="A"; for Mod="5mC", Base="C".
#' @param grangesAnnotationsWithCounts A GRanges object based on grangesAnnotations with the counts:
#' \itemize{
#'   \item Modcount:           The number of "Mod"  within this feature.
#'   \item Modcount_perkbp:   (The number of "Mod"  within this feature divided by the size of the feature)*1000.
#'   \item Basecount:          The number of "Base" within this feature.
#'   \item Basecount_perkbp:  (The number of "Base" within this feature divided by the size of the feature)*1000.
#' }
#' The Genomic features categories must be in a column named "type".
#' @param cFeaturesToCompare Names of the features to be displayed and compared. Defaults to c("gene", "intergenic").
#' @param lUseCountsPerkbp If TRUE, use counts per kbp (Modcount_perkbp and Basecount_perkbp) to calculate
#' the proportion between features displayed. If FALSE, use counts (Modcount and Basecount) to calculate
#' the proportion between features displayed. Defaults to FALSE.
#' @param cBaseMotif The name of the motif with the base letter of the modified base.
#' @param cModMotif The name of the motif with the modification in the output.
#' @keywords DrawModBasePropByFeature
#' @export
#' @examples
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' #loading annotation
#' library(rtracklayer)
#' myAnnotations <- readGFFAsGRanges(system.file(package="DNAModAnnot", "extdata",
#'                                           "ptetraurelia_mac_51_annotation_v2.0_sca171819.gff3"))
#' myAnnotations <- PredictMissingAnnotation(grangesAnnotations = myAnnotations,
#'                                           grangesGenome = myGrangesGenome,
#'                                           cFeaturesColName = "type",
#'                                           cGeneCategories = c("gene"),
#'                                           lAddIntronRangesUsingExon = TRUE)
#'
#' #Preparing a grangesPacBioGFF and a grangesPacBioCSV datasets
#' myGrangesPacBioGFF <-
#'    ImportPacBioGFF(cPacBioGFFPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.modifications.sca171819.gff"),
#'                    cNameModToExtract = "m6A",
#'                    cModNameInOutput = "6mA",
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <-
#'    ImportPacBioCSV(cPacBioCSVPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.bases.sca171819.csv"),
#'                    cSelectColumnsToExtract = c("refName", "tpl", "strand", "base",
#'                                                "score", "ipdRatio", "coverage"),
#'                    lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <- myGposPacBioCSV[myGposPacBioCSV$base == "A"]
#'
#' #Retrieve annotations with "Mod" and "Base" counts (and counts per kbp)
#' myAnn_ModBase_counts <- GetModBaseCountsByFeature(grangesAnnotations=myAnnotations,
#'                                                   grangesModPos=myGrangesPacBioGFF,
#'                                                   gposModTargetBasePos=myGposPacBioCSV)
#'
#' DrawModBasePropByFeature(grangesAnnotationsWithCounts=myAnn_ModBase_counts,
#'                                 cFeaturesToCompare=c("gene", "intergenic"),
#'                                 lUseCountsPerkbp=TRUE,
#'                                 cBaseMotif="A",
#'                                 cModMotif="6mA")
DrawModBasePropByFeature <- function(grangesAnnotationsWithCounts,
                                            cFeaturesToCompare = c("gene", "intergenic"),
                                            lUseCountsPerkbp = FALSE,
                                            cBaseMotif,
                                            cModMotif) {
  opar <- par()

  dataToPlot <- subset(grangesAnnotationsWithCounts, type %in% cFeaturesToCompare)
  dataToPlot$type <- as.factor(as.character(dataToPlot$type))

  if(lUseCountsPerkbp){
    dataToPlot$Modprop  <- 100 * dataToPlot$Modcount_perkbp  / sum(dataToPlot$Modcount_perkbp)
    dataToPlot$Baseprop <- 100 * dataToPlot$Basecount_perkbp / sum(dataToPlot$Basecount_perkbp)
  } else {
    dataToPlot$Modprop  <- 100 * dataToPlot$Modcount  / sum(dataToPlot$Modcount)
    dataToPlot$Baseprop <- 100 * dataToPlot$Basecount / sum(dataToPlot$Basecount)
  }

  #Barplot
  dataToPlot2 <- aggregate(dataToPlot$Modprop, by=list(dataToPlot$type), sum)
  dataToPlot3 <- aggregate(dataToPlot$Baseprop, by=list(dataToPlot$type), sum)
  dataToPlot4 <- matrix(c(dataToPlot2$x, dataToPlot3$x), nrow=2, byrow = TRUE )
  rownames(dataToPlot4) <- c(cModMotif, cBaseMotif)
  colnames(dataToPlot4) <- levels(dataToPlot$type)

  layout(mat=matrix(1:2, ncol = 2), widths = c(8,2))
  par(mar=c(5.1,5.1,3.1,0))
  barplot(height = dataToPlot4, col=c("red3", "grey"),
          beside = TRUE,
          cex.names = 1.5, cex.axis = 1.5, cex.lab=1.5,
          ylab = "Cumulated proportion (%)",
          # main = paste0(cModMotif," (or ", cBaseMotif,
          #               ") proportion between displayed features (% of counts",
          #               ifelse(lUseCountsPerkbp, " per kbp)",")")),
          ylim=c(0,100))

  par( mar=c(0,0,0,0))
  plot(1:10, col="transparent", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
  legend("left",
         legend = paste0(c(cModMotif, cBaseMotif),"%"),
         col=c("red3","grey"),
         pch=15, cex=1.25,
         bty="n")

  par(mar=opar$mar)
  layout(mat=matrix(1:1, ncol = 2), widths = c(1,1))
  title(main = paste0(cModMotif," (or ", cBaseMotif,
                      ") proportion between displayed features (% of counts",
                      ifelse(lUseCountsPerkbp, " per kbp)",")")),)
}

#' GetModBaseCountsCategories Function (ModAnnot)
#'
#' Subfunction to return categories of counts (or counts per KiloBase pairs (kbp)) of
#' the base modified (Mod) and the base letter of the modified base (Base). These categories
#' should have a comparable amount of Mod (or Base) between them.
#' Example: for Mod="6mA", Base="A"; for Mod="5mC", Base="C".
#' @param grangesAnnotationsWithCounts A GRanges object based on grangesAnnotations with the counts:
#' \itemize{
#'   \item Modcount:           The number of "Mod"  within this feature.
#'   \item Modcount_perkbp:   (The number of "Mod"  within this feature divided by the size of the feature)*1000.
#'   \item Basecount:          The number of "Base" within this feature.
#'   \item Basecount_perkbp:  (The number of "Base" within this feature divided by the size of the feature)*1000.
#' }
#' The Genomic features categories must be in a column named "type".
#' @keywords internal
.GetModBaseCountsCategories <- function(grangesAnnotationsWithCounts){

  breaks_v <- unique(quantile(grangesAnnotationsWithCounts$Modcount, probs = seq(0, 1, length.out = 11)))
  nLoop <- 0
  while(length(breaks_v) < 5 & nLoop < 10){
    nLoop <- nLoop + 1
    breaks_v <- unique(quantile(grangesAnnotationsWithCounts$Modcount,
                                probs = seq(0, 1, length.out = 11+10*nLoop)))
  }
  grangesAnnotationsWithCounts$Modcount_category <- cut(grangesAnnotationsWithCounts$Modcount, breaks= breaks_v,
                                                        include.lowest = TRUE, right = TRUE  )


  breaks_v <- unique(quantile(grangesAnnotationsWithCounts$Basecount, probs = seq(0, 1, length.out = 11)))
  nLoop <- 0
  while(length(breaks_v) < 5 & nLoop < 10){
    nLoop <- nLoop + 1
    breaks_v <- unique(quantile(grangesAnnotationsWithCounts$Basecount,
                                probs = seq(0, 1, length.out = 11+10*nLoop)))
  }
  grangesAnnotationsWithCounts$Basecount_category <- cut(grangesAnnotationsWithCounts$Basecount, breaks= breaks_v,
                                                         include.lowest = TRUE, right = TRUE  )

  breaks_v <- unique(quantile(grangesAnnotationsWithCounts$Modcount_perkbp, probs = seq(0, 1, length.out = 11)))
  nLoop <- 0
  while(length(breaks_v) < 5 & nLoop < 10){
    nLoop <- nLoop + 1
    breaks_v <- unique(quantile(grangesAnnotationsWithCounts$Modcount_perkbp,
                                probs = seq(0, 1, length.out = 11+10*nLoop)))
  }
  grangesAnnotationsWithCounts$Modcount_perkbp_category <- cut(grangesAnnotationsWithCounts$Modcount_perkbp,
                                                               breaks= breaks_v,
                                                               include.lowest = TRUE, right = TRUE  )


  breaks_v <- unique(quantile(grangesAnnotationsWithCounts$Basecount_perkbp, probs = seq(0, 1, length.out = length(breaks_v) )))
  nLoop <- 0
  while(length(breaks_v) < 5 & nLoop < 10){
    nLoop <- nLoop + 1
    breaks_v <- unique(quantile(grangesAnnotationsWithCounts$Basecount_perkbp,
                                probs = seq(0, 1, length.out = 11+10*nLoop)))
  }
  grangesAnnotationsWithCounts$Basecount_perkbp_category <- cut(grangesAnnotationsWithCounts$Basecount_perkbp,
                                                                breaks= breaks_v,
                                                                include.lowest = TRUE, right = TRUE  )

  return(grangesAnnotationsWithCounts)
}

#' DrawParamPerModBaseCategories Function (ModAnnot)
#'
#' Return boxplots describing, for each feature provided, the distribution of a parameter provided by categories of counts
#' (or counts per KiloBase pairs (kbp)) of
#' the base modified (Mod) and the base letter of the modified base (Base). These categories
#' should have a comparable amount of Mod (or Base) between them.
#' Example: for Mod="6mA", Base="A"; for Mod="5mC", Base="C".
#' @param grangesAnnotationsWithCounts A GRanges object based on grangesAnnotations with the counts:
#' \itemize{
#'   \item Modcount:           The number of "Mod"  within this feature.
#'   \item Modcount_perkbp:   (The number of "Mod"  within this feature divided by the size of the feature)*1000.
#'   \item Basecount:          The number of "Base" within this feature.
#'   \item Basecount_perkbp:  (The number of "Base" within this feature divided by the size of the feature)*1000.
#' }
#' The Genomic features categories must be in a column named "type".
#' @param cParamColname Name of the column, in the grangesAnnotationsWithCounts provided, containing
#' the parameter to be compared with Mod|Base categories.
#' @param cParamFullName Name of the parameter to be displayed in the title of the plot. Defaults to cParamColname.
#' @param cParamYLabel Name of the parameter to be displayed on the Y-axis of the plot. Defaults to cParamColname.
#' @param cSelectFeature The name of the feature from the annotation to be analysed along counts and parameter provided.
#' Defaults to NULL (No subsetting: counts and parameter provided from all features in the annotation provided will
#' be used).
#' @param lUseCountsPerkbp If TRUE, use counts per kbp (Modcount_perkbp and Basecount_perkbp) to calculate
#' the Mod|Base categories to be displayed. If FALSE, use counts (Modcount and Basecount) to calculate
#' the Mod|Base categories to be displayed. Defaults to TRUE.
#' @param lKeepOutliers If FALSE, remove outliers before plotting. Defaults to FALSE.
#' @param lUseSameYAxis If TRUE, the 2 plots will use the same range for the
#' y-axis. If FALSE, y-axis of the 2 plots will be independant.
#' Defaults to FALSE.
#' @param cBaseMotif The name of the motif with the base letter of the modified base.
#' @param cModMotif The name of the motif with the modification in the output.
#' @param lBoxPropToCount If TRUE, the width of each boxplot depends
#' on the number of Mod (or Base) in the categories. If FALSE, all boxplots will have the
#' same size. Defaults to TRUE.
#' @keywords DrawParamPerModBaseCategories
#' @export
#' @examples
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#'
#' #loading annotation
#' library(rtracklayer)
#' myAnnotations <- readGFFAsGRanges(system.file(package="DNAModAnnot", "extdata",
#'                                           "ptetraurelia_mac_51_annotation_v2.0_sca171819.gff3"))
#'
#' #Preparing a grangesPacBioGFF and a grangesPacBioCSV datasets
#' myGrangesPacBioGFF <-
#'    ImportPacBioGFF(cPacBioGFFPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.modifications.sca171819.gff"),
#'                    cNameModToExtract = "m6A",
#'                    cModNameInOutput = "6mA",
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <-
#'    ImportPacBioCSV(cPacBioCSVPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.bases.sca171819.csv"),
#'                    cSelectColumnsToExtract = c("refName", "tpl", "strand", "base",
#'                                                "score", "ipdRatio", "coverage"),
#'                    lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <- myGposPacBioCSV[myGposPacBioCSV$base == "A"]
#'
#' #Retrieve annotations with "Mod" and "Base" counts (and counts per kbp)
#' myAnn_ModBase_counts <- GetModBaseCountsByFeature(grangesAnnotations=myAnnotations,
#'                                                   grangesModPos=myGrangesPacBioGFF,
#'                                                   gposModTargetBasePos=myGposPacBioCSV)
#'
#' #Add Parameter by feature to annotation file
#' myAnn_ModBase_counts$ParamToPlot <- width(myAnn_ModBase_counts)
#'
#' DrawParamPerModBaseCategories(myAnn_ModBase_counts,
#'                                      cParamColname="ParamToPlot",
#'                                      cParamFullName = "Gene Width",
#'                                      cParamYLabel = "Gene Width (bp)",
#'                                      cSelectFeature = c("gene"),
#'                                      lUseCountsPerkbp = TRUE,
#'                                      lKeepOutliers = FALSE,
#'                                      lUseSameYAxis = FALSE,
#'                                      cBaseMotif = "A",
#'                                      cModMotif = "6mA",
#'                                      lBoxPropToCount=FALSE)
DrawParamPerModBaseCategories <- function(grangesAnnotationsWithCounts,
                                                 cParamColname,
                                                 cParamFullName = cParamColname,
                                                 cParamYLabel = cParamColname,
                                                 cSelectFeature = NULL,
                                                 lUseCountsPerkbp = TRUE,
                                                 lKeepOutliers = FALSE,
                                                 lUseSameYAxis = FALSE,
                                                 cBaseMotif,
                                                 cModMotif,
                                                 lBoxPropToCount=TRUE){
  if(is.null(cSelectFeature)){
    dataToPlot <- grangesAnnotationsWithCounts
  } else {
    dataToPlot <- subset(grangesAnnotationsWithCounts, type %in% c(cSelectFeature))
    dataToPlot$type <- as.factor(as.character(dataToPlot$type))
  }
  dataToPlot <- .GetModBaseCountsCategories(dataToPlot)
  cParamColname=mcols(dataToPlot)[[cParamColname]]

  if(lUseCountsPerkbp){
    Modcount_vect <- mcols(dataToPlot)[["Modcount_perkbp_category"]]
    Basecount_vect <- mcols(dataToPlot)[["Basecount_perkbp_category"]]
    nTextMaxMin <- round(c(min(dataToPlot$Modcount_perkbp),  max(dataToPlot$Modcount_perkbp),
                           min(dataToPlot$Basecount_perkbp), max(dataToPlot$Basecount_perkbp)), 2)
  } else {
    Modcount_vect <- mcols(dataToPlot)[["Modcount_category"]]
    Basecount_vect <- mcols(dataToPlot)[["Basecount_category"]]
    nTextMaxMin <- c(min(dataToPlot$Modcount),  max(dataToPlot$Modcount),
                     min(dataToPlot$Basecount), max(dataToPlot$Basecount))
  }

  nBreaksVMod <- 1:length(levels(Modcount_vect))
  nBreaksVBase <- 1:length(levels(Basecount_vect))

  if(lBoxPropToCount){
    w1 <- boxplot( cParamColname~Modcount_vect, plot=FALSE)$n
    w1 <- w1/sum(w1)
    w2 <- boxplot( cParamColname~Basecount_vect, plot=FALSE)$n
    w2 <- w2/sum(w2)
  } else {
    w1 <- rep(1, length(levels(Modcount_vect )))
    w2 <- rep(1, length(levels(Basecount_vect)))
  }

  if(lUseSameYAxis){
    if(lKeepOutliers){
      ylimits <- c(min(cParamColname), max(cParamColname))
    } else {
      b_stats1 <- boxplot( cParamColname~Modcount_vect, plot=FALSE)$stats
      b_stats2 <- boxplot( cParamColname~Basecount_vect, plot=FALSE)$stats
      ylimits <- c(min( b_stats1[1,], b_stats2[1,] ),
                   max( b_stats1[5,], b_stats2[5,] ) )
    }
  } else{
    ylimits <- NULL
  }

  nMarBottom <- max(nchar(c(levels(Modcount_vect),levels(Basecount_vect))))*8.1/12

  opar <- par()
  layout(mat = matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE), heights = c(8,2))
  par(mar=c(nMarBottom,5.1,4.1,2.1))
  boxplot(cParamColname~Modcount_vect, col = "red3", outline=lKeepOutliers, width = w1, ylab=cParamYLabel,
          cex.lab=1.5, cex.axis=1.5, xlab="", las=3,
          ylim = ylimits )

  boxplot(cParamColname~Basecount_vect, col = "grey", outline=lKeepOutliers, width = w2, ylab=cParamYLabel,
          cex.lab=1.5, cex.axis=1.5, xlab="", las=3,
          ylim =  ylimits)

  par(mar=c(1.1,4.1,0.1,2.1))
  .DrawPolygonRange(nBreaksV=nBreaksVMod, nTextMaxMin=nTextMaxMin[1:2],
                           cAxisName=paste0(cModMotif, ifelse(lUseCountsPerkbp, " counts per kbp", " counts")),
                           cColor = "red3")
  .DrawPolygonRange(nBreaksV=nBreaksVBase, nTextMaxMin=nTextMaxMin[3:4],
                           cAxisName=paste0(cBaseMotif, ifelse(lUseCountsPerkbp, " counts per kbp", " counts")),
                           cColor = "grey")
  par(mar=opar$mar)
  layout(mat = matrix(1:1, nrow = 1), heights = c(1,1))
  mtext(paste0(cParamFullName, " per ",
               cModMotif," (or ",cBaseMotif,") count ",
               ifelse(lUseCountsPerkbp, "per kbp",""),
               " levels"), side=3, cex=1.5, line = 2)
}

#' DrawPolygonRange Function (ModAnnot)
#'
#' Return a polygon representing the range of Mod (or Base) categories of counts (or counts per kbp).
#' @param nBreaksV Numeric vector giving the breaks to be used for the polygon.
#' Should be equal to seq(from = 1, to= c(number of categories), by=1).
#' @param nTextMaxMin Numeric vector giving the minimum and maximum values to be plotted on the polygon.
#' Should be c(min1, max1, min2, max2).
#' @param cAxisName The name of the axis to draw.
#' @param cColor The color of the polygon to draw.
#' @keywords internal
.DrawPolygonRange <- function(nBreaksV,
                                     nTextMaxMin,
                                     cAxisName,
                                     cColor) {
  plot(nBreaksV,seq(0,1,length.out = length(nBreaksV)),
       type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="",
       xlim=c(nBreaksV[1]-0.5,length(nBreaksV)+0.5), ylim=c(0,1))
  polygon(x = c(nBreaksV, length(nBreaksV)),
          y = c(seq(0, 0.30, length.out=length(nBreaksV)), 0),
          col=cColor)
  par(xpd=TRUE)
  text(x = c(nBreaksV[1], tail(nBreaksV,1)),
       y = c(0.55,0.55),
       nTextMaxMin,
       cex = 1.5)
  text(x = mean(c(nBreaksV[1], tail(nBreaksV,1))),
       y = 0.9,
       cAxisName,
       cex = 1.5)
  par(xpd=FALSE)
}

#' GetModBaseCountsWithinFeature Function (ModAnnot)
#'
#' Cut the Annotation provided into x windows of relative size and return, for each window,
#' the counts (and counts per KiloBase pairs (kbp)) of
#' the base modified (Mod) and the base letter of the modified base (Base).
#' Example: for Mod="6mA", Base="A"; for Mod="5mC", Base="C".
#' @param grangesAnnotations A GRanges object containing the annotation for the genome assembly
#' corresponding to the grangesModPos and gposModTargetBasePos provided. The Genomic features categories must be in a column named "type".
#' @param nWindowsNb The number of output windows by feature. The annotation provided (with input ranges)
#' will be cut into this number of output windows.
#' Each output window will have the same size as the other output windows from the same input range.
#' Defaults to 20.
#' @param grangesModPos A GRanges object containing Modifications Positions data to be counted.
#' @param gposModTargetBasePos A GRanges or GPos object containing Base Positions (which can be targeted by the modification) to be counted.
#' @param lIgnoreStrand If TRUE, Mod and Base will be counted independently of the strand of each feature. If FALSE, only Mod and Base
#' on the same strand as the feature will be counted. Defaults to FALSE.
#' @return A GRanges object based on grangesAnnotations with the counts:
#' \itemize{
#'   \item Modcount:     The   number of "Mod"  within this window  of this feature.
#'   \item ModcountSum:  Total number of "Mod"  within all  windows of this feature.
#'   \item Modprop:      (Modcount / ModcountSum)*100
#'   \item Basecount:    The   number of "Base" within this window  of this feature.
#'   \item BasecountSum: Total number of "Base" within all  windows of this feature.
#'   \item Baseprop:     (Basecount/BasecountSum)*100
#' }
#' @keywords GetModBaseCountsWithinFeature
#' @export
#' @examples
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#'
#' #loading annotation
#' library(rtracklayer)
#' myAnnotations <- readGFFAsGRanges(system.file(package="DNAModAnnot", "extdata",
#'                                           "ptetraurelia_mac_51_annotation_v2.0_sca171819.gff3"))
#'
#' #Preparing a grangesPacBioGFF and a grangesPacBioCSV datasets
#' myGrangesPacBioGFF <-
#'    ImportPacBioGFF(cPacBioGFFPath = system.file(package="DNAModAnnot", "extdata",
#'                                    "ptetraurelia.modifications.sca171819.gff"),
#'                    cNameModToExtract = "m6A",
#'                    cModNameInOutput = "6mA",
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <-
#'    ImportPacBioCSV(cPacBioCSVPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.bases.sca171819.csv"),
#'                    cSelectColumnsToExtract = c("refName", "tpl", "strand", "base",
#'                                                "score", "ipdRatio", "coverage"),
#'                    lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <- myGposPacBioCSV[myGposPacBioCSV$base == "A"]
#'
#' #Retrieve annotations with "Mod" and "Base" counts (and counts per kbp)
#' myAnn_ModBase_counts_by_window <-
#'    GetModBaseCountsWithinFeature(grangesAnnotations=myAnnotations[myAnnotations$type == "gene",],
#'                                  grangesModPos=myGrangesPacBioGFF,
#'                                  gposModTargetBasePos=myGposPacBioCSV,
#'                                  nWindowsNb = 20)
#' myAnn_ModBase_counts_by_window
GetModBaseCountsWithinFeature <- function(grangesAnnotations,
                                              nWindowsNb=20,
                                              grangesModPos,
                                              gposModTargetBasePos,
                                              lIgnoreStrand = FALSE){
  gr_a <- grangesAnnotations[width(grangesAnnotations) >= nWindowsNb]
  ans <- unlist(tile(gr_a, n=nWindowsNb))
  win_id_pos <- rep(1:nWindowsNb, length.out = length(ans))
  win_id_neg <- rep(nWindowsNb:1, length.out = length(ans))
  ans$window_id <- ifelse(strand(ans) == "+", win_id_pos, win_id_neg)
  ans$window_id2 <- rep(1:length(gr_a), each = nWindowsNb)

  gr_b <- grangesModPos
  hits <- findOverlaps(ans, gr_b, ignore.strand = lIgnoreStrand)

  ans$Modcount <- countQueryHits(hits)
  modcountSum <- aggregate(ans$Modcount, by=list(ans$window_id2), sum)
  ans$ModcountSum <- rep(modcountSum$x, each=nWindowsNb)
  ans$Modprop <- 100*ans$Modcount/ans$ModcountSum
  if(length(which(is.nan(ans$Modprop)))>0){
    ans[is.nan(ans$Modprop)]$Modprop <- 0
  }

  gr_b <- gposModTargetBasePos
  hits <- findOverlaps(ans, gr_b, ignore.strand = lIgnoreStrand)

  ans$Basecount <- countQueryHits(hits)
  basecountSum <- aggregate(ans$Basecount, by=list(ans$window_id2), sum)
  ans$BasecountSum <- rep(basecountSum$x, each=nWindowsNb)
  ans$Baseprop <- 100*ans$Basecount/ans$BasecountSum
  if(length(which(is.nan(ans$Baseprop)))>0){
    ans[is.nan(ans$Baseprop)]$Baseprop <- 0
  }

  return(ans)
}

#' DrawModBaseCountsWithinFeature Function (ModAnnot)
#'
#' Return a plot describing the proportion of the base modified (Mod) and the base letter
#' of the modified base (Base) between windows (of relative size) of feature provided.
#' Example: for Mod="6mA", Base="A"; for Mod="5mC", Base="C".
#' @param grangesAnnotationsWithCountsByWindow A GRanges object containing feature annotation with the counts:
#' \itemize{
#'   \item Modcount:     The   number of "Mod"  within this window  of this feature.
#'   \item ModcountSum:  Total number of "Mod"  within all  windows of this feature.
#'   \item Modprop:      (Modcount / ModcountSum)*100
#'   \item Basecount:    The   number of "Base" within this window  of this feature.
#'   \item BasecountSum: Total number of "Base" within all  windows of this feature.
#'   \item Baseprop:     (Basecount/BasecountSum)*100
#' }
#' The Genomic features categories must be in a column named "type".
#' @param cFeatureName Name of the feature (which is cut into windows) to be displayed.
#' @param cBaseMotif The name of the motif with the base letter of the modified base.
#' @param cModMotif The name of the motif with the modification in the output.
#' @keywords DrawModBaseCountsWithinFeature
#' @export
#' @examples
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#'
#' #loading annotation
#' library(rtracklayer)
#' myAnnotations <- readGFFAsGRanges(system.file(package="DNAModAnnot", "extdata",
#'                                           "ptetraurelia_mac_51_annotation_v2.0_sca171819.gff3"))
#'
#' #Preparing a grangesPacBioGFF and a grangesPacBioCSV datasets
#' myGrangesPacBioGFF <-
#'    ImportPacBioGFF(cPacBioGFFPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.modifications.sca171819.gff"),
#'                    cNameModToExtract = "m6A",
#'                    cModNameInOutput = "6mA",
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <-
#'    ImportPacBioCSV(cPacBioCSVPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.bases.sca171819.csv"),
#'                    cSelectColumnsToExtract = c("refName", "tpl", "strand", "base",
#'                                                "score", "ipdRatio", "coverage"),
#'                    lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <- myGposPacBioCSV[myGposPacBioCSV$base == "A"]
#'
#' #Retrieve annotations with "Mod" and "Base" counts (and counts per kbp)
#' myAnn_ModBase_counts_by_window <-
#'    GetModBaseCountsWithinFeature(grangesAnnotations=myAnnotations[myAnnotations$type == "gene",],
#'                                 grangesModPos=myGrangesPacBioGFF,
#'                                 gposModTargetBasePos=myGposPacBioCSV,
#'                                 nWindowsNb = 20)
#'
#' DrawModBaseCountsWithinFeature(myAnn_ModBase_counts_by_window,
#'                                       cFeatureName="gene",
#'                                       cBaseMotif="A",
#'                                       cModMotif="6mA")
DrawModBaseCountsWithinFeature <- function(grangesAnnotationsWithCountsByWindow,
                                           cFeatureName,
                                           cBaseMotif,
                                           cModMotif){
  dataToPlot <- aggregate(grangesAnnotationsWithCountsByWindow$Modprop,
                          by=list(grangesAnnotationsWithCountsByWindow$window_id),
                          function(x) 100*sum(x)/sum(grangesAnnotationsWithCountsByWindow$Modprop))
  dataToPlot2 <- aggregate(grangesAnnotationsWithCountsByWindow$Baseprop,
                           by=list(grangesAnnotationsWithCountsByWindow$window_id),
                           function(x) 100*sum(x)/sum(grangesAnnotationsWithCountsByWindow$Baseprop))
  ylimits <- c(0, max(dataToPlot$x, dataToPlot2$x))

  layout(mat=matrix(1:2, ncol = 2), widths = c(8,2))
  opar <- par()
  par(mar=c(7.1,5.1,2.1,0))

  bar_vect <- barplot(dataToPlot$x, col = "red3",
                      width = 1, space = 1.5,
                      main = paste0("Cumulated ", cModMotif," (or ", cBaseMotif,
                                    ") proportion within ", cFeatureName),
                      ylab = "Cumulated proportion (%)",
                      xlab= "",
                      cex.lab=1.5, cex.names = 1.5, cex.axis=1.5, yaxs="i", ylim = ylimits)
  title(xlab=paste0("Relative position within ",cFeatureName, " (%)"), line=5.5, cex.lab=1.5)
  space_vect <- c(2.5, rep(1.5, length(dataToPlot2$x)-1))
  bar_vect2 <- barplot(dataToPlot2$x, col="grey", beside = TRUE, add = TRUE,
                       width=1, space = space_vect, cex.axis=1.5)

  t_length <- sapply(1:(length(bar_vect)),
                     function(i) {mean(c(bar_vect[i],  bar_vect2[i]))
                     })
  t_text <- as.character( ( 100/length(bar_vect) ) * (1:(length(bar_vect))) )
  t_text2 <- as.numeric(t_text) - (as.numeric(t_text)[2] - as.numeric(t_text)[1])
  t_text <- paste0("[",t_text2, ">", t_text,"]")
  axis(side = 1, labels = rep("", length.out=length(t_length)),
       at = t_length, cex.axis=1.25, lwd=2, las=3, line = 0)
  text(t_length, 0, labels = t_text, srt = 60, xpd = TRUE, adj = c(1.2,0.75), cex=1.25)

  par( mar=c(0,0,0,0))
  plot(1:10, col="transparent", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
  legend("left",
         legend = paste0(c(cModMotif, cBaseMotif),"%"),
         col=c("red3","grey"),
         pch=15, cex=1.25, bty="n"
  )

  par( mar=opar$mar)
  layout(mat=matrix(1:1, ncol = 2), widths = c(1,1))

}

#' GetDistFromFeaturePos Function (ModAnnot)
#'
#' Return, in GRanges objects via a GRangesList, the distance between the "Positions" provided with "grangesData" argument
#' (e.g. position of some target sites) and the feature positions provided with "grangesAnnotations".
#' If the ranges in "grangesAnnotations" are not 1-bp positions,
#' the positions of the boundaries will be used as the feature positions: in this case, 2 GRanges ("Position" vs featureStart;
#' "Position" vs featureEnd) will be exported in the output instead of 1 Granges ("Position" vs featureStart).
#' This function can also directly compute counts or proportion of the provided "Positions"
#' at each nucleotide position around provided features (see "lGetGRangesInsteadOfListCounts" argument).
#' @param grangesAnnotations A GRanges object containing the annotation for the genome assembly
#' corresponding to the grangesModPos and gposModTargetBasePos provided. The Genomic features categories must be in a column named "type".
#' @param cSelectFeature The name of the feature from the annotation to be analysed.
#' Defaults to NULL (all ranges from grangesAnnotations will be kept).
#' @param grangesData A GRanges object containing "Positions" to be counted around features positions.
#' @param lGetGRangesInsteadOfListCounts If FALSE, return, in dataframes via a list, the counts
#' (or proportion if "lGetPropInsteadOfCounts" is TRUE) of "Positions" at each base position
#' around feature positions. Defaults to FALSE.
#' @param lGetPropInsteadOfCounts If "lGetGRangesInsteadOfListCounts" and "lGetPropInsteadOfCounts" are TRUE,
#' return the proportion of "Positions" near feature positions: counts / sum of counts.
#' If the ranges in "grangesAnnotations" are not 1-bp positions,
#' the proportion of "Positions" is calculated near both feature positions (Start and End): counts / (sum of counts near feature1 +
#' sum of counts near feature2). Defaults to TRUE.
#' @param nWindowSizeAroundFeaturePos Size, in base pairs, of the viewing window around the feature positions.
#' @param cWhichStrandVsFeaturePos A character value describing if distance comparison must be made between "Mod"
#' (or "Base") and the feature positions...
#' \itemize{
#'   \item "same":    ...if these positions are on the same strand  only.
#'   \item "opposite":...if these positions are on opposite strands only.
#'   \item "both":    ...for all of these positions: same and opposite strands.
#' }
#' @param lAddCorrectedDistFrom5pTo3p If TRUE, the distance will be corrected to reflect 5' to 3' direction
#' and will be stored in a new column (dist_5to3). Defaults to TRUE.
#' @param cFeaturePosNames A character vector returning the names of the feature positions provided.
#' Defaults to c("Start","End").
#' \itemize{
#'   \item If the width of ranges is equal to 1, the name of the feature will be the first element of the vector.
#'   \item If the width of ranges is above 1, the names of the feature borders will be the first element then the second element.
#' }
#' @return A GRangesList with 1 or 2 GRanges objects containing ranges of "Positions" with their distance to feature positions:
#' \itemize{
#'   \item If the width of annotation ranges is equal to 1, 1 GRanges is  provided ("Position" vs featureStart).
#'   \item If the width of annotation ranges is above    1, 2 GRanges are provided ("Position" vs featureStart; "Position" vs featureEnd).
#' }
#' If "lGetGRangesInsteadOfListCounts" is FALSE, retrieve instead a list with 1 or 2 dataframe(s)
#' containing "Positions" counts by distance to feature positions:
#' \itemize{
#'   \item If the width of annotation ranges is equal to 1, 1 dataframe  is  provided ("Position" vs featureStart).
#'   \item If the width of annotation ranges is above    1, 2 dataframes are provided ("Position" vs featureStart; "Position" vs featureEnd).
#' }
#' If a "Position" is within "nWindowSizeAroundFeaturePos" base pairs of x different feature positions: this "Position"
#' will then reported x times with the distance to each feature position around this "Position".
#' @keywords GetDistFromFeaturePos
#' @importFrom GenomicRanges resize
#' @export
#' @examples
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#'
#' #loading annotation
#' library(rtracklayer)
#' myAnnotations <- readGFFAsGRanges(system.file(package="DNAModAnnot", "extdata",
#'                                           "ptetraurelia_mac_51_annotation_v2.0_sca171819.gff3"))
#'
#' #Preparing a grangesPacBioGFF and a grangesPacBioCSV datasets
#' myGrangesPacBioGFF <-
#'    ImportPacBioGFF(cPacBioGFFPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.modifications.sca171819.gff"),
#'                    cNameModToExtract = "m6A",
#'                    cModNameInOutput = "6mA",
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <-
#'    ImportPacBioCSV(cPacBioCSVPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.bases.sca171819.csv"),
#'                    cSelectColumnsToExtract = c("refName", "tpl", "strand", "base",
#'                                                "score", "ipdRatio", "coverage"),
#'                    lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <- myGposPacBioCSV[myGposPacBioCSV$base == "A"]
#'
#' #Retrieve, in a list, dataframes of Mod counts per Distance values from feature positions
#' myModDistCountsList <- GetDistFromFeaturePos(
#'   grangesAnnotations = myAnnotations,
#'   cSelectFeature = "gene",
#'   grangesData = myGrangesPacBioGFF,
#'   lGetGRangesInsteadOfListCounts = FALSE,
#'   lGetPropInsteadOfCounts = TRUE,
#'   cWhichStrandVsFeaturePos = "both", nWindowSizeAroundFeaturePos = 600,
#'   lAddCorrectedDistFrom5pTo3p = TRUE,
#'   cFeaturePosNames = c("TSS", "TTS")
#' )
#' myModDistCountsList
GetDistFromFeaturePos <- function(grangesAnnotations,
                               cSelectFeature = NULL,
                               grangesData,
                               lGetGRangesInsteadOfListCounts = FALSE,
                               lGetPropInsteadOfCounts = TRUE,
                               nWindowSizeAroundFeaturePos,
                               cWhichStrandVsFeaturePos="both",
                               lAddCorrectedDistFrom5pTo3p=TRUE,
                               cFeaturePosNames=c("Start","End")){
  if(!is.null(cSelectFeature)){grangesAnnotations <- subset(grangesAnnotations, type == cSelectFeature)}

  featureIsPos <- all(width(grangesAnnotations)==1)
  if(featureIsPos){
    print("Features provided are 1bp positions: returning positions...")
    grangesAnnotationsPos1 <- resize(grangesAnnotations, width = 1, fix = "start") + nWindowSizeAroundFeaturePos
  } else {
    print("Features provided are not 1bp positions: returning features boundaries as positions...")
    grangesAnnotationsPos1 <- resize(grangesAnnotations, width = 1, fix = "start") + nWindowSizeAroundFeaturePos
    grangesAnnotationsPos2 <- resize(grangesAnnotations, width = 1, fix = "end")   + nWindowSizeAroundFeaturePos
  }

  print("Returning distance from feature positions...")
  ansStart <- .GetPosDistFromFeaturePos(grangesData=grangesData, grangesAnnotationsPos=grangesAnnotationsPos1,
                                        cWhichStrandVsFeaturePos=cWhichStrandVsFeaturePos, nWindowSizeAroundFeaturePos=nWindowSizeAroundFeaturePos,
                                        lAddCorrectedDistFrom5pTo3p=lAddCorrectedDistFrom5pTo3p)

  if(!featureIsPos){
    print("Returning distance from feature end positions...")
    ansEnd <- .GetPosDistFromFeaturePos(grangesData=grangesData, grangesAnnotationsPos=grangesAnnotationsPos2,
                                        cWhichStrandVsFeaturePos=cWhichStrandVsFeaturePos, nWindowSizeAroundFeaturePos=nWindowSizeAroundFeaturePos,
                                        lAddCorrectedDistFrom5pTo3p=lAddCorrectedDistFrom5pTo3p)
    ans_list <- GRangesList(ansStart, ansEnd)
    names(ans_list) <- paste0(rep(cFeaturePosNames,   each=1), "GR")

    if(!lGetGRangesInsteadOfListCounts) {
      ans_list <- GetListCountsByDist(listGRangesDist = ans_list,
                                      lGetPropInsteadOfCounts = lGetPropInsteadOfCounts,
                                      lAddCorrectedDistFrom5pTo3p = lAddCorrectedDistFrom5pTo3p)
    }
  } else {
    ans_list <- GRangesList(ansStart)
    names(ans_list) <- paste0(cFeaturePosNames[1], "GR")

    if(!lGetGRangesInsteadOfListCounts) {
      ans_list <- GetListCountsByDist(listGRangesDist = ans_list,
                                      lGetPropInsteadOfCounts = lGetPropInsteadOfCounts,
                                      lAddCorrectedDistFrom5pTo3p = lAddCorrectedDistFrom5pTo3p)
    }
  }

  return(ans_list)
}

#' .GetPosDistFromFeaturePos Function (ModAnnot)
#'
#' Return, in a GRanges object, the distance between "Mod" (or "Base") positions (provided with grangesData)
#' and feature positions provided with grangesAnnotationsPos.
#' "Mod": the base modified. "Base": the base letter of the modified base.
#' Example: for Mod="6mA", Base="A"; for Mod="5mC", Base="C".
#' @param grangesData A GRanges object containing the Mod (or Base) positions to be compared with
#' the grangesAnnotationsPos positions.
#' @param grangesAnnotationsPos A GRanges object containing the annotation positions to be compared with
#' the grangesData positions. The Genomic features categories must be in a column named "type".
#' @param nWindowSizeAroundFeaturePos Size, in base pairs, of the viewing window around the feature positions.
#' @param cWhichStrandVsFeaturePos A character value describing if distance comparison must be made between "Mod"
#' (or "Base") and the feature positions...
#' \itemize{
#'   \item "same":    ...if these positions are on the same strand  only.
#'   \item "opposite":...if these positions are on opposite strands only.
#'   \item "both":    ...for all of these positions: same and opposite strands.
#' }
#' @param lAddCorrectedDistFrom5pTo3p If TRUE, the distance will be corrected to reflect 5' to 3' direction
#' and will be stored in a new column (dist_5to3). Defaults to TRUE.
#' @keywords internal
#' @importFrom S4Vectors Pairs
.GetPosDistFromFeaturePos <- function(grangesData, grangesAnnotationsPos, cWhichStrandVsFeaturePos,
                                      nWindowSizeAroundFeaturePos, lAddCorrectedDistFrom5pTo3p){
  if (cWhichStrandVsFeaturePos=="both") {
    hits <- findOverlaps(grangesData, grangesAnnotationsPos, ignore.strand = TRUE)
  } else if (cWhichStrandVsFeaturePos=="same") {
    hits <- findOverlaps(grangesData, grangesAnnotationsPos, ignore.strand = FALSE)
  } else if (cWhichStrandVsFeaturePos=="opposite") {
    hits <- findOverlaps(grangesData, invertStrand(grangesAnnotationsPos), ignore.strand = FALSE)
  }
  pairs <- Pairs(grangesData, grangesAnnotationsPos, hits = hits)

  pairs@first$dist_v <- (start(pairs@first) - start(pairs@second)) - nWindowSizeAroundFeaturePos
  ans <- pairs@first
  if(lAddCorrectedDistFrom5pTo3p){
    ansf <- ans[which(strand(pairs@second) == "+")]
    ansr <- ans[which(strand(pairs@second) == "-")]
    ansf$dist_5to3 <- ansf$dist_v
    ansr$dist_5to3 <- ansr$dist_v * -1
    ans <- c(ansf, ansr)
  }
  return(sort(ans))
}

#' GetListCountsByDist Function (ModAnnot)
#'
#' Return, in dataframes via a list, the counts (or proportion) of provided "Positions" by distance from feature
#' positions. If the input list contains 2 GRanges, 2 dataframes ("Position" vs featureStart; "Position" vs featureEnd)
#' will be exported in the output instead of 1 dataframe ("Position" vs featureStart).
#' @param listGRangesDist A GRangesList with 1 or 2 GRanges objects containing ranges of given "Positions"
#' with their distance to feature positions.
#' @param lAddCorrectedDistFrom5pTo3p If TRUE, the distance will be corrected to reflect 5' to 3' direction
#' and will be stored in a new column (dist_5to3). Defaults to TRUE.
#' @param lGetPropInsteadOfCounts If TRUE, return the proportion of given "Positions" near feature
#' position: counts / sum of counts. If listGRangesDist contains 4 GRanges, the proportion
#' of given "Positions" is calculated near both feature positions: counts / (sum of counts near feature1 +
#' sum of counts near feature2). Defaults to TRUE.
#' @return A list with 1 or 2 dataframe(s) containing "Positions" counts by distance to feature positions:
#' \itemize{
#'   \item If 1 GRanges are provided in listGRangesDist, 1 dataframe  is  provided ("Position" vs featureStart).
#'   \item If 2 GRanges are provided in listGRangesDist, 2 dataframes are provided ("Position" vs featureStart; "Position" vs featureEnd).
#' }
#' If a given "Position" is within nWindowSizeAroundFeaturePos base pairs of x different feature positions: this given "Position"
#' will then reported x times with the distance to each feature position.
#' @keywords GetListCountsByDist
#' @export
#' @examples
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#'
#' #loading annotation
#' library(rtracklayer)
#' myAnnotations <- readGFFAsGRanges(system.file(package="DNAModAnnot", "extdata",
#'                                           "ptetraurelia_mac_51_annotation_v2.0_sca171819.gff3"))
#'
#' #Preparing a grangesPacBioGFF and a grangesPacBioCSV datasets
#' myGrangesPacBioGFF <-
#'    ImportPacBioGFF(cPacBioGFFPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.modifications.sca171819.gff"),
#'                    cNameModToExtract = "m6A",
#'                    cModNameInOutput = "6mA",
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <-
#'    ImportPacBioCSV(cPacBioCSVPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.bases.sca171819.csv"),
#'                    cSelectColumnsToExtract = c("refName", "tpl", "strand", "base",
#'                                                "score", "ipdRatio", "coverage"),
#'                    lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <- myGposPacBioCSV[myGposPacBioCSV$base == "A"]
#'
#' #Retrieve, in a list, dataframes of ModBase counts per Distance values from feature positions
#' myModDistGRangesList <- GetDistFromFeaturePos(
#'   grangesAnnotations = myAnnotations,
#'   cSelectFeature = "gene",
#'   grangesData = myGrangesPacBioGFF,
#'   lGetGRangesInsteadOfListCounts = TRUE,
#'   cWhichStrandVsFeaturePos = "both", nWindowSizeAroundFeaturePos = 600,
#'   lAddCorrectedDistFrom5pTo3p = TRUE,
#'   cFeaturePosNames = c("TSS", "TTS")
#' )
#' myModDistCountsList <- GetListCountsByDist(listGRangesDist = myModDistGRangesList,
#'                                            lAddCorrectedDistFrom5pTo3p = TRUE,
#'                                            lGetPropInsteadOfCounts = TRUE)
#' myModDistCountsList
GetListCountsByDist <- function(listGRangesDist, lAddCorrectedDistFrom5pTo3p = TRUE,
                                lGetPropInsteadOfCounts = TRUE){
  twoPositions <- (length(listGRangesDist)==2)
  returnPropByfPos <- ifelse(twoPositions, FALSE, TRUE)

  table1 <- .GetCountsByDist(grangesDist=listGRangesDist[[1]],
                             lAddCorrectedDistFrom5pTo3p=lAddCorrectedDistFrom5pTo3p)
  if(twoPositions){
    table2 <- .GetCountsByDist(grangesDist=listGRangesDist[[2]],
                               lAddCorrectedDistFrom5pTo3p=lAddCorrectedDistFrom5pTo3p)

    if(lGetPropInsteadOfCounts){
      table1_freq <- 100*table1$Freq / sum(c(table1$Freq, table2$Freq))
      table2$Freq <- 100*table2$Freq / sum(c(table1$Freq, table2$Freq))
      table1$Freq <- table1_freq
    }
    table <- list(table1, table2)
  } else {
    if(lGetPropInsteadOfCounts){ table1$Freq <- 100*table1$Freq / sum(table1$Freq)  }
    table <- list(table1)
  }
  names(table) <- names(listGRangesDist)
  return(table)
}

#' GetCountsByDist Function (ModAnnot)
#'
#' Return, in a table, the counts (or proportion) of given "Positions" by distance from feature
#' positions.
#' @param grangesDist A GRanges objects containing ranges of given "Positions"
#' with their distance to feature positions.
#' @param lAddCorrectedDistFrom5pTo3p If TRUE, the distance will be corrected to reflect 5' to 3' direction
#' and will be stored in a new column (dist_5to3). Defaults to TRUE.
#' @return A table objects containing "Mod" (or "Base") counts by distance to feature positions.
#' If a "Mod" (or "Base") is within nWindowSizeAroundFeaturePos base pairs of x different feature
#' positions: this "Mod" (or "Base")
#' will then reported x times with the distance to each feature position.
#' @keywords internal
.GetCountsByDist <- function(grangesDist, lAddCorrectedDistFrom5pTo3p=TRUE){
  tableCounts <- as.data.frame(table(
    mcols(grangesDist)[[ifelse(lAddCorrectedDistFrom5pTo3p,"dist_5to3","dist_v")]]
  ))
  tableCounts$Var1 <- as.numeric(as.character(tableCounts$Var1))
  return(tableCounts)
}

#' DrawModBasePropDistFromFeature Function (ModAnnot)
#'
#' Return, in dataframes via a list, the counts (or proportion) of "Mod" (or "Base") positions by distance from feature
#' positions. If the input list contains 4 GRanges, 4 dataframes ("Mod" vs featureStart; "Mod" vs featureEnd;
#' "Base" vs featureStart; "Base" vs featureEnd) will be exported in the output instead of 2 dataframes
#' ("Mod" vs featureStart; "Base" vs featureStart).
#' "Mod": the base modified. "Base": the base letter of the modified base.
#' Example: for Mod="6mA", Base="A"; for Mod="5mC", Base="C".
#' @param listModCountsDistDataframe A list with 1 or 2 dataframe(s) containing "Mod" counts by distance to feature positions:
#' \itemize{
#'   \item If 1 dataframe  is  provided in listModCountsDistDataframe, 1 position  will be plotted (featureStart).
#'   \item If 2 dataframes are provided in listModCountsDistDataframe, 2 positions will be plotted (featureStart; featureEnd).
#' }
#' Must have the same length as the list provided with "listBaseCountsDistDataframe".
#' @param listBaseCountsDistDataframe A list with 1 or 2 dataframe(s) containing "Base" counts by distance to feature positions:
#' \itemize{
#'   \item If 1 dataframe  is  provided in listModCountsDistDataframe, 1 position  will be plotted (featureStart).
#'   \item If 2 dataframes are provided in listModCountsDistDataframe, 2 positions will be plotted (featureStart; featureEnd).
#' }
#' Must have the same length as the list provided with "listModCountsDistDataframe".
#' @param cFeaturePosNames A character vector returning the names of the feature positions provided.
#' Defaults to c("Start","End").
#' \itemize{
#'   \item If 2 dataframes are provided in listModBaseCountsDistDataframe, the name of the feature will be the first element of the vector.
#'   \item If 4 dataframes are provided in listModBaseCountsDistDataframe, the names of the feature borders will be the first element then the second element.
#' }
#' @param cBaseMotif The name of the motif with the base letter of the modified base.
#' @param cModMotif The name of the motif with the modification in the output.
#' @param nDensityBaseMotif Numeric vector giving the density of the polygon made with the "Base" counts by distance to feature positions.
#' @keywords DrawModBasePropDistFromFeature
#' @export
#' @examples
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#'
#' #loading annotation
#' library(rtracklayer)
#' myAnnotations <- readGFFAsGRanges(system.file(package="DNAModAnnot", "extdata",
#'                                           "ptetraurelia_mac_51_annotation_v2.0_sca171819.gff3"))
#'
#' #Preparing a grangesPacBioGFF and a grangesPacBioCSV datasets
#' myGrangesPacBioGFF <-
#'    ImportPacBioGFF(cPacBioGFFPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.modifications.sca171819.gff"),
#'                    cNameModToExtract = "m6A",
#'                    cModNameInOutput = "6mA",
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <-
#'    ImportPacBioCSV(cPacBioCSVPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.bases.sca171819.csv"),
#'                    cSelectColumnsToExtract = c("refName", "tpl", "strand", "base",
#'                                                "score", "ipdRatio", "coverage"),
#'                    lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <- myGposPacBioCSV[myGposPacBioCSV$base == "A"]
#'
#' #Retrieve, in a list, dataframes of ModBase counts per Distance values from feature positions
#' myModDistCountsList <- GetDistFromFeaturePos(
#'   grangesAnnotations = myAnnotations,
#'   cSelectFeature = "gene",
#'   grangesData = myGrangesPacBioGFF,
#'   lGetGRangesInsteadOfListCounts = FALSE,
#'   lGetPropInsteadOfCounts = TRUE,
#'   cWhichStrandVsFeaturePos = "both", nWindowSizeAroundFeaturePos = 600,
#'   lAddCorrectedDistFrom5pTo3p = TRUE,
#'   cFeaturePosNames = c("TSS", "TTS")
#' )
#' myBaseDistCountsList <- GetDistFromFeaturePos(
#'   grangesAnnotations = myAnnotations,
#'   cSelectFeature = "gene",
#'   grangesData = myGposPacBioCSV,
#'   lGetGRangesInsteadOfListCounts = FALSE,
#'   lGetPropInsteadOfCounts = TRUE,
#'   cWhichStrandVsFeaturePos = "both", nWindowSizeAroundFeaturePos = 600,
#'   lAddCorrectedDistFrom5pTo3p = TRUE,
#'   cFeaturePosNames = c("TSS", "TTS")
#' )
#' DrawModBasePropDistFromFeature(
#'   listModCountsDistDataframe = myModDistCountsList,
#'   listBaseCountsDistDataframe = myBaseDistCountsList,
#'   cFeaturePosNames = c("TSS", "TTS"),
#'   cBaseMotif = "A",
#'   cModMotif = "6mA"
#' )
DrawModBasePropDistFromFeature <- function(listModCountsDistDataframe,
                                           listBaseCountsDistDataframe,
                                           cFeaturePosNames=c("Start","End"),
                                           cBaseMotif,
                                           cModMotif,
                                           nDensityBaseMotif=50){
  if(length(listModCountsDistDataframe) == length(listBaseCountsDistDataframe)){
    twoPositions <- (length(listModCountsDistDataframe)==2)

    dataToPlot1 <- listModCountsDistDataframe[[1]]
    dataToPlot2 <- listBaseCountsDistDataframe[[1]]

    opar <- par()
    if(twoPositions){
      dataToPlot3 <- listModCountsDistDataframe[[2]]
      dataToPlot4 <- listBaseCountsDistDataframe[[2]]
      ylimits <- c(0, max(dataToPlot1$Freq, dataToPlot2$Freq, dataToPlot3$Freq, dataToPlot4$Freq))
      layout(mat = matrix(1:2, ncol = 2))
      par(mar=c(5.1,10.2,4.1,1.1))
    } else {
      ylimits <- c(0, max(dataToPlot1$Freq, dataToPlot2$Freq))
      par(mar=c(5.1,10.2,4.1,10.2))
    }

    plot(x=dataToPlot1$Var1, y=dataToPlot1$Freq, ylim=ylimits, col="white",
         xaxs="i", yaxs="i",yaxt="n", bty=ifelse(twoPositions,"C","o"),
         xlab=paste0("Distance from ", cFeaturePosNames[1]),
         ylab="", cex.axis=1.25, cex.lab=2)
    polygon(x=c(min(dataToPlot1$Var1), dataToPlot1$Var1, max(dataToPlot1$Var1)) ,
            y=c( 0, dataToPlot1$Freq, 0),
            col='red2', border="red2")
    polygon(x=c(min(dataToPlot2$Var1), dataToPlot2$Var1, max(dataToPlot2$Var1)) ,
            y=c( 0, dataToPlot2$Freq, 0),
            col='grey50',density = nDensityBaseMotif)
    abline(v=0, lty=2, lwd=2, col=rgb(0,0,0,0.3))

    axis(side = 2, lwd = 3, col = "red2", cex.axis=1.25 )
    mtext(paste0(cModMotif,  " proportion"), col = "red2", side = 2, line = 2.5, cex=2)

    if(!twoPositions){
      axis(side = 4, lwd = 3, col = "grey50", cex.axis=1.25 )
      mtext(paste0(cBaseMotif, " proportion"), col = "grey50", side = 4, line = 2.5, cex=2)
    } else {
      par(mar=c(5.1,1.1,4.1,10.2))
      plot(x=dataToPlot3$Var1, y=dataToPlot3$Freq, ylim=ylimits, col="white",
           xaxs="i", yaxs="i",yaxt="n", bty="]",
           xlab=paste0("Distance from ", cFeaturePosNames[2]),
           ylab="", cex.axis=1.25, cex.lab=2)
      polygon(x=c(min(dataToPlot3$Var1), dataToPlot3$Var1, max(dataToPlot3$Var1)) ,
              y=c( 0, dataToPlot3$Freq, 0),
              col='red2', border="red2")
      polygon(x=c(min(dataToPlot4$Var1), dataToPlot4$Var1, max(dataToPlot4$Var1)) ,
              y=c( 0, dataToPlot4$Freq, 0),
              col='grey50',density = nDensityBaseMotif)
      abline(v=0, lty=2, lwd=2, col=rgb(0,0,0,0.3))

      axis(side = 4, lwd = 3, col = "grey50", cex.axis=1.25 )
      mtext(paste0(cBaseMotif, " proportion"), col = "grey50", side = 4, line = 2.5, cex=2)
      par(xpd=NA)
      segments(x0=0, x1=min(dataToPlot3$Var1)*2, y0=0,            lty=2, lwd=1, col="black")
      segments(x0=0, x1=min(dataToPlot3$Var1)*2, y0=max(ylimits), lty=2, lwd=1, col="black")
      par(xpd=FALSE)
    }
    par(mar=opar$mar)
    layout(mat = matrix(1:1, ncol = 1))
    if(twoPositions){
      title(paste0(cModMotif, " and ", cBaseMotif, " distance from ", cFeaturePosNames[1],
                   " and ", cFeaturePosNames[2]), line = 2.5 )
    } else {
      title(paste0(cModMotif, " and ", cBaseMotif, " distance from ", cFeaturePosNames[1]), line = 2.5 )
    }

  } else {
    print("Error: listModCountsDistDataframe and listBaseCountsDistDataframe must have a same length.")
  }

}

#' AddToModBasePropDistFromFeaturePlot Function (ModAnnot)
#'
#' Add an additional axis with data on a plot generated by DrawModBasePropDistFromFeature function.
#' @param dPosCountsDistFeatureStart A dataframe containing the data counts by distance to feature positions.
#' @param dPosCountsDistFeatureEnd A second dataframe containing the data counts
#' by distance to the second feature positions.
#' If NULL, only 1 position  will be plotted (featureStart). Defaults to NULL.
#' @param cSubtitleContent A character vector giving the content of the subtitle to add below the title of the plot.
#' @param cParamYLabel A character vector giving the content of the label on the new Y axis to add on the plot.
#' @param cParamColor The color of the line and new axis to be added. Defaults to "cyan3".
#' @param cParamType The type of plot to be added. See "type" argument plot base function. Defaults to "l".
#' @param cParamLwd The width of the line/points to be added. Defaults to 3.
#' @param cParamLty If a line is plotted, change the type of the line to be added. Defaults to 3.
#' @param lAddAxisOnLeftSide If TRUE, add the new axis on the left side of the plot.
#' If FALSE, add the new axis on the right side of the plot. Defaults to TRUE.
#' @param nYLimits Numeric vector giving the limits for the new Y axis. Defaults to NULL (will use the minimum and
#' the maximum of both data provided (dPosCountsDistFeatureStart and dPosCountsDistFeatureEnd)).
#' @keywords AddToModBasePropDistFromFeaturePlot
#' @export
#' @examples
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#'
#' #loading annotation
#' myAnnotations <- rtracklayer::readGFFAsGRanges(system.file(package="DNAModAnnot", "extdata",
#'                                           "ptetraurelia_mac_51_annotation_v2.0_sca171819.gff3"))
#'
#' #Preparing a grangesPacBioGFF and a grangesPacBioCSV datasets
#' myGrangesPacBioGFF <-
#'    ImportPacBioGFF(cPacBioGFFPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.modifications.sca171819.gff"),
#'                    cNameModToExtract = "m6A",
#'                    cModNameInOutput = "6mA",
#'                    cContigToBeAnalyzed = names(myGenome))
#' myGposPacBioCSV <-
#'    ImportPacBioCSV(cPacBioCSVPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.bases.sca171819.csv"),
#'                    cSelectColumnsToExtract = c("refName", "tpl", "strand", "base",
#'                                                "score", "ipdRatio", "coverage"),
#'                    lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'                    cContigToBeAnalyzed = names(myGenome))
#'
#' #Retrieve, in a list, dataframes of ModBase counts per Distance values from feature positions
#' myModDistCountsList <- GetDistFromFeaturePos(
#'   grangesAnnotations = myAnnotations,
#'   cSelectFeature = "gene",
#'   grangesData = myGrangesPacBioGFF,
#'   lGetGRangesInsteadOfListCounts = FALSE,
#'   lGetPropInsteadOfCounts = TRUE,
#'   cWhichStrandVsFeaturePos = "both", nWindowSizeAroundFeaturePos = 600,
#'   lAddCorrectedDistFrom5pTo3p = TRUE,
#'   cFeaturePosNames = c("TSS", "TTS")
#' )
#' myBaseDistCountsList <- GetDistFromFeaturePos(
#'   grangesAnnotations = myAnnotations,
#'   cSelectFeature = "gene",
#'   grangesData = myGposPacBioCSV,
#'   lGetGRangesInsteadOfListCounts = FALSE,
#'   lGetPropInsteadOfCounts = TRUE,
#'   cWhichStrandVsFeaturePos = "both", nWindowSizeAroundFeaturePos = 600,
#'   lAddCorrectedDistFrom5pTo3p = TRUE,
#'   cFeaturePosNames = c("TSS", "TTS")
#' )
#' DrawModBasePropDistFromFeature(
#'   listModCountsDistDataframe = myModDistCountsList,
#'   listBaseCountsDistDataframe = myBaseDistCountsList,
#'   cFeaturePosNames = c("TSS", "TTS"),
#'   cBaseMotif = "A",
#'   cModMotif = "6mA"
#' )
#'
#' #Loading Bam data
#' myBamfile <- Rsamtools::BamFile(file = system.file(package="DNAModAnnot", "extdata",
#'                      "PTET_MonoNuc_3-2new.pe.sca171819.sorted.bam"))
#' myBam_GRanges <- as(GenomicAlignments::readGAlignments(myBamfile), "GRanges")
#' myBam_GRanges <- GetGposCenterFromGRanges(grangesData=myBam_GRanges)
#'
#' #Retrieve dataframes of Read center counts per Distance values in a list
#' myCountsDist_List_bamfile <- GetDistFromFeaturePos(grangesAnnotations=myAnnotations,
#'                                                    cSelectFeature="gene",
#'                                                    lGetGRangesInsteadOfListCounts =  FALSE,
#'                                                    lGetPropInsteadOfCounts = FALSE,
#'                                                    grangesData = myBam_GRanges,
#'                                                    cWhichStrandVsFeaturePos="both",
#'                                                    nWindowSizeAroundFeaturePos=600,
#'                                                    lAddCorrectedDistFrom5pTo3p=TRUE,
#'                                                    cFeaturePosNames=c("TSS","TTS"))
#'
#' #adding new axis to plot from DrawModBasePropDistFromFeature function
#' AddToModBasePropDistFromFeaturePlot(
#'    dPosCountsDistFeatureStart=myCountsDist_List_bamfile[[1]],
#'    dPosCountsDistFeatureEnd=myCountsDist_List_bamfile[[2]],
#'    cSubtitleContent="Along with nucleosome center distance (MonoNuc_3-2newreplicate)",
#'    cParamYLabel="Nucleosome center count (MonoNuc_3-2newreplicate)",
#'    cParamColor="cyan3",
#'    lAddAxisOnLeftSide=TRUE)
AddToModBasePropDistFromFeaturePlot <- function(dPosCountsDistFeatureStart,
                                                       dPosCountsDistFeatureEnd=NULL,
                                                       cSubtitleContent,
                                                       cParamYLabel,
                                                       cParamColor="cyan3", cParamType="l",
                                                       cParamLwd=3, cParamLty=3,
                                                       lAddAxisOnLeftSide=TRUE,
                                                       nYLimits=NULL){
  twoPositions <- (!is.null(dPosCountsDistFeatureEnd))

  dataToPlot1 <- dPosCountsDistFeatureStart

  opar <- par()
  if(twoPositions){
    dataToPlot2 <- dPosCountsDistFeatureEnd
    if (is.null(nYLimits)) {
      nYLimits <- c(min(dataToPlot1$Freq, dataToPlot2$Freq), max(dataToPlot1$Freq, dataToPlot2$Freq))
    }
    par(new=TRUE)
    layout(mat = matrix(1:2, ncol = 2))
    par(mar=c(5.1,10.2,4.1,1.1))
  } else {
    if (is.null(nYLimits)) {nYLimits <- c(min(dataToPlot1$Freq), max(dataToPlot1$Freq))}
    par(new=TRUE)
    par(mar=c(5.1,10.2,4.1,10.2))
  }

  par(mfg=c(1,1))
  plot(x=dataToPlot1$Var1,
       y=dataToPlot1$Freq,
       ylim=nYLimits, col=cParamColor,
       xaxs="i", yaxs="i", yaxt="n", xaxt="n", bty="n",
       xlab="", ylab="", type=cParamType, lwd=cParamLwd, lty=cParamLty)
  lines(x=dataToPlot1$Var1, y=dataToPlot1$Freq, col=cParamColor, lwd=0.5, lty=1)

  if(twoPositions){
    par(mar=c(5.1,1.1,4.1,10.2))
    par(mfg=c(1,2))
    plot(x=dataToPlot2$Var1,
         y=dataToPlot2$Freq,
         ylim=nYLimits, col=cParamColor,
         xaxs="i", yaxs="i", yaxt="n", xaxt="n", bty="n",
         xlab="", ylab="", type=cParamType, lwd=cParamLwd, lty=cParamLty)
    lines(x=dataToPlot2$Var1, y=dataToPlot2$Freq, col=cParamColor, lwd=0.5, lty=1)
  }
  if(twoPositions){
    if (lAddAxisOnLeftSide) {
      par(mfg=c(1,1))
      par(mar=c(5.1,4.1,4.1,1.1))
    } else {
      par(mfg=c(1,2))
      par(mar=c(5.1,1.1,4.1,4.1))
    }
  } else {
    par(mfg=c(1,1))
    par(mar=c(5.1,4.1,4.1,4.1))
  }
  axis(side = ifelse(lAddAxisOnLeftSide, 2, 4), lwd = 3, col = cParamColor, cex.axis=1.25,
       at=pretty(range(dataToPlot1$Freq),n = 10), line=-1,
       labels = pretty(range(dataToPlot1$Freq),n = 10))

  mtext(cParamYLabel, col = cParamColor, side = ifelse(lAddAxisOnLeftSide, 2, 4), line = 1.5, cex=2)

  layout(mat = matrix(1, ncol = 1))
  par(mfg=c(1,1))
  par(mar=opar$mar)
  title(cSubtitleContent, line = ifelse(lAddAxisOnLeftSide, 1.5, 0.5))

  par(new=FALSE)
}



#' AdaptedIdeogramTrackWithoutBandsData Function (ModAnnot)
#'
#' Subfunction to return an empty IdeogramTrack without bands data.
#' @param grangesGenome A GRanges object containing the width of each contig.
#' @param cContigToViz A character vector describing which contigs this track list should be prepared for.
#' @param cOrgAssemblyName The name of the genome assembly provided.
#' @importFrom Gviz IdeogramTrack
#' @export
#' @keywords internal
AdaptedIdeogramTrackWithoutBandsData <- function(grangesGenome,
                                                 cContigToViz,
                                                 cOrgAssemblyName){
  tmp_start= as.vector(sapply(cContigToViz, function(x){round(seq(from=start(grangesGenome[seqnames(grangesGenome)==x & strand(grangesGenome)=="+"]),
                                                                  to=end(grangesGenome[seqnames(grangesGenome)==x][1]),
                                                                  length.out = 5), 0)[1:4]}))

  tmp_end= as.vector(sapply(cContigToViz, function(x){1 + round(seq(from=start(grangesGenome[seqnames(grangesGenome)==x & strand(grangesGenome)=="+"]),
                                                                    to=end(grangesGenome[seqnames(grangesGenome)==x][1]),
                                                                    length.out = 5), 0)[2:5]}))

  dCytoBandsDataframe <- data.frame(chrom=rep(cContigToViz, each=4),
                                    chromStart=tmp_start,
                                    chromEnd=tmp_end,
                                    name=rep(cContigToViz, each=4),
                                    gieStain="gneg")

  return(IdeogramTrack(genome=cOrgAssemblyName,chromosome=cContigToViz,name=cOrgAssemblyName,bands=dCytoBandsDataframe, cex=1))
}



#' ExportFilesForGViz Function (ModAnnot)
#'
#' Export data as files that can be directly used with Gviz Package (making tracks + displaying).
#' Except for gff3 format, all available format allows streaming while making and displaying tracks with Gviz package.
#' Multiple objects can be exported at the same time. All arguments (except "dnastringsetGenome" argument) must have the same length.
#' @param dnastringsetGenome A DNAStringSet object containing the sequence for each contig.
#' @param cFileNames A character vector giving the names of the files to be exported.
#' @param listObjects A list of objects to be exported.
#' Required format of the objects is described in "cFileFormats" argument documentation.
#' @param cFileFormats A character vector giving the format for each file to be exported.
#' The following exportation formats are supported:
#' \itemize{
#'   \item Fasta  using DNAStrinSet  object : for sequenceTracks (streaming)
#'   \item Bam    using GRanges-like object : for alignmentTracks (streaming) (only the positions of the ranges will be retrieved) or
#'   for dataTracks (streaming) (only coverage of provided ranges will be displayed) or
#'   for annotationTracks (streaming) (positions of the ranges will be retrieved with the names of each range (e.g. gene name))
#'   \item bw (bigwig) using GRanges-like object : for dataTracks (streaming) (only the chosen numeric parameter will be
#'   retrieved for each range provided)
#'   \item gff3   using GRanges-like object : for annotationTracks or geneRegionTracks
#'   (not in streaming: more memory required for displaying)
#' }
#' @param cBigwigParameters A character vector describing the name of the parameter to be stored in bigwig files.
#' Must correspond to the name of a column in the provided Granges object.
#' Use NA as value in the vector if the associated file is not to be exported as a bigwig.
#' Defaults to "rep(NA, length(cFileNames)))"
#' @param lBigwigParametersByStrand A logical vector: if TRUE, the bigwig parameter will be negative
#' for each range that is located on the reverse strand in the GRanges object provided.
#' Use NA as value in the vector if the associated file is not to be exported as a bigwig.
#' Defaults to "rep(NA, length(cFileNames)))"
#' @param cBamXaParameters A character vector describing a parameter to be stored as a "Xa" optional field
#' in the exported bam file.
#' Use NA as value in the vector if the associated file is not to be exported as a bam.
#' If the exported file is a bam and if the value is NA or NULL, the bam will be exported without the "Xa" optional field.
#' Defaults to "rep(NA, length(cFileNames)))"
#' @export
#' @keywords ExportFilesForGViz
#' @examples
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#'
#' #Preparing a grangesPacBioCSV dataset
#' myGposPacBioCSV <-
#'    ImportPacBioCSV(cPacBioCSVPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.bases.sca171819.csv"),
#'                    cSelectColumnsToExtract = c("refName", "tpl", "strand", "base",
#'                                                "score", "ipdRatio", "coverage"),
#'                    lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'                    cContigToBeAnalyzed = names(myGenome))
#'
#' ##NOT RUN!
#' ##Export files for Gviz
#' #ExportFilesForGViz(dnastringsetGenome = myGenome,
#' #                   cFileNames = c("ipdRatio_for_each_A.bw",
#' #                                  "score_for_all_bases.bw",
#' #                                  "newFastaOnlyscaffold51_17.fa"),
#' #                   listGRangesObjects = c(myGposPacBioCSV[myGposPacBioCSV$base == "A"],
#' #                                          myGposPacBioCSV,
#' #                                          myGenome["scaffold51_17"]),
#' #                   cFileFormats = c("bw", "bw", "fa"),
#' #                   cBigwigParameters = c("ipdRatio", "score", NA),
#' #                   lBigwigParametersByStrand = c(TRUE, TRUE, NA),
#' #                   cBamXaParameters = c(NA, NA, NA))
#'
#' #loading annotation
#' library(rtracklayer)
#' myAnnotations <- readGFFAsGRanges(system.file(package="DNAModAnnot", "extdata",
#'                                           "ptetraurelia_mac_51_annotation_v2.0_sca171819.gff3"))
#'
#' ##NOT RUN!
#' ##Export files for Gviz
#' #ExportFilesForGViz(dnastringsetGenome = myGenome,
#' #                   cFileNames = c("genes.bam"),
#' #                   listGRangesObjects = list(myAnnotations[myAnnotations$type == "gene"]),
#' #                   cFileFormats = c("bam"),
#' #                   cBamXaParameters = c("Name") )
ExportFilesForGViz <- function(dnastringsetGenome,
                               cFileNames,
                               listObjects,
                               cFileFormats,
                               cBigwigParameters = rep(NA, length(cFileNames)),
                               lBigwigParametersByStrand = rep(NA, length(cFileNames)),
                               cBamXaParameters = rep(NA, length(cFileNames))){
  if(is.list(listObjects)){
    nCheckLength <- unique(c(length(cFileNames), length(listObjects), length(cFileFormats),
                             length(cBigwigParameters), length(lBigwigParametersByStrand),
                             length(cBamXaParameters)))
    if(length(nCheckLength)==1){


      sapply(1:length(cFileNames), function(x){
        switch (cFileFormats[x],
                bw = .ExportBigwigForGviz(cFileName = cFileNames[x],
                                          grangesObject = listObjects[[x]],
                                          cBigwigParameter = cBigwigParameters[x],
                                          lBigwigParametersByStrand = lBigwigParametersByStrand[x],
                                          dnastringsetGenome = dnastringsetGenome),
                bam = .ExportBamForGviz(cFileName = cFileNames[x],
                                        grangesObject = listObjects[[x]],
                                        cBamXaParameter = cBamXaParameters[[x]],
                                        dnastringsetGenome = dnastringsetGenome),
                fa = .ExportFastaForGviz(cFileName = cFileNames[x],
                                         dnastringsetObject = listObjects[[x]]),
                gff3 = rtracklayer::export.gff3(con = cFileNames[x],
                                                object = listObjects[[x]])
        )
        print(paste0(cFileNames[x], " Finished."))

      })


    } else {
      print("Error: all parameters provided must be of same length (except for dnastringsetGenome).")
    }
  } else {
    print("Error: listObjects parameter must be a list.")
  }
}

#' .ExportFastaForGviz Function (ModAnnot)
#'
#' Subfunction to export DNAStringset object as fasta file.
#' @param cFileName Name of the file to be exported.
#' @param dnastringsetObject A DNAStringSet object to be exported.
#' @export
#' @keywords internal
.ExportFastaForGviz <- function(cFileName,
                                dnastringsetObject){
  Biostrings::writeXStringSet(filepath = cFileName,
                              x = dnastringsetObject)
  write.table(x = Biostrings::fasta.index(filepath = cFileName, seqtype = "DNA"),
              file = paste0(cFileName,".fai"),
              quote = FALSE, sep = "\t")
}

#' .ExportBigwigForGviz Function (ModAnnot)
#'
#' Subfunction to export a GRanges object as a bw (bigwig) file.
#' @param cFileName Name of the file to be exported.
#' @param grangesObject A GRanges object to be exported.
#' @param cBigwigParameter The name of the parameter to be exported within the bigwig file.
#' @param lBigwigParametersByStrand If TRUE, the parameter is negative if associated range is on the reverse strand.
#' @param dnastringsetGenome A DNAStringSet object containing the sequence for each contig.
#' @export
#' @keywords internal
.ExportBigwigForGviz <- function(cFileName,
                                 grangesObject,
                                 cBigwigParameter,
                                 lBigwigParametersByStrand,
                                 dnastringsetGenome){
  grangesObjectTmp <- grangesObject
  GenomicRanges::mcols(grangesObjectTmp) <- NULL
  grangesObjectTmp$score <- GenomicRanges::mcols(grangesObject)[[cBigwigParameter]]
  if(lBigwigParametersByStrand){
    grangesObjectTmp[GenomicRanges::strand(grangesObject)=="-"]$score <-
      grangesObjectTmp[GenomicRanges::strand(grangesObject)=="-"]$score * -1
  }
  GenomicRanges::strand(grangesObjectTmp) <- "*"
  BSgenome::seqinfo(grangesObjectTmp) <-
    BSgenome::seqinfo(dnastringsetGenome)[BSgenome::seqnames(BSgenome::seqinfo(grangesObjectTmp))]
  rtracklayer::export.bw(grangesObjectTmp, con = cFileName)
}

#' .ExportBamForGviz Function (ModAnnot)
#'
#' Subfunction to export a GRanges object as a bam file.
#' @param cFileName Name of the file to be exported.
#' @param grangesObject A GRanges object to be exported.
#' @param cBamXaParameter The name of the parameter to be exported as a "Xa" optional field within the bam file.
#' @param dnastringsetGenome A DNAStringSet object containing the sequence for each contig.
#' @export
#' @keywords internal
.ExportBamForGviz <- function(cFileName,
                              grangesObject,
                              cBamXaParameter,
                              dnastringsetGenome){
  grangesObject <- as(grangesObject,"GAlignments")
  names(grangesObject) <- 1:length(grangesObject)

  if( (!is.na(cBamXaParameter)) | !is.null(cBamXaParameter)) {
    GenomicRanges::mcols(grangesObject)$Xa <- GenomicRanges::mcols(grangesObject)[[cBamXaParameter]]
  }

  BSgenome::seqinfo(grangesObject) <-
    BSgenome::seqinfo(dnastringsetGenome)[BSgenome::seqnames(BSgenome::seqinfo(grangesObject))]

  rtracklayer::export(grangesObject,
                      Rsamtools::BamFile(cFileName),
                      format = "bam")
}

#' ImportBamExtendedAnnotationTrack Function (ModAnnot)
#'
#' Subfunction to import a bam file for displaying using Gviz package
#' while retrieving the content of the "Xa" optional field (will be retrieved as "tag" field).
#' (See Gviz package documentation to see how to use it while making tracks)
#' After making the AnnotationTrack,
#' in order to allow the names of the genomic features to be displayed, the "mapping" sub-list of
#' the generated annotation track must be completed with the new "id" and "group" values chosen while making the annotationTrack.
#' Example: trackAnnotation@mapping <- list(id="tag", group="tag")
#' @export
#' @keywords internal
ImportBamExtendedAnnotationTrack <- function(file, selection) {
  if (!file.exists(paste(file, "bai", sep = ".")) &&
      !file.exists(paste(paste(head(strsplit("xxx.bam", ".", fixed = TRUE)[[1]], -1), collapse = "."), "bai", sep = "."))) {
    stop(
      "Unable to find index for BAM file '", file, "'. You can build an index using the following command:\n\t",
      "library(Rsamtools)\n\tindexBam(\"", file, "\")"
    )
  }
  sinfo <- Rsamtools::scanBamHeader(file)[[1]]

  res <- if (!as.character(GenomicRanges::seqnames(selection)[1]) %in% names(sinfo$targets)) {
    GenomicRanges::mcols(selection) <- S4Vectors::DataFrame(id = "NA", group = "NA")
    selection[0]
  } else {
    param <- Rsamtools::ScanBamParam(what = c("pos", "qwidth", "strand", "qname"), tag="Xa",
                                     which = selection, flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
    x <- Rsamtools::scanBam(file, param = param)[[1]]
    GenomicRanges::GRanges(
      seqnames = GenomicRanges::seqnames(selection), ranges = IRanges::IRanges(start = x[["pos"]], width = x[["qwidth"]]), strand = x[["strand"]],
      id = make.unique(x[["qname"]]), group = x[["qname"]]
    )
  }

  x <- Rsamtools::scanBam(file, param = Rsamtools::ScanBamParam(which = selection, tag="Xa"))[[1]]$tag
  res$tag <- x[[which(names(x) == "Xa")]]

  return(res)
}
