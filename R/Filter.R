### =========================================================================
### Functions from Filter category
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
### FiltContig
### FiltParam
### FiltFdrBased
### FiltPacBio
### FiltDeepSignal
###
### -------------------------------------------------------------------------

#' FiltContig Function (Filter)
#'
#' Filter out data from contigs that do not reach criterias of selection.
#' @param gposModBasePos An UnStitched GPos object containing PacBio CSV data to be filtered.
#' @param grangesModPos A GRanges object containing PacBio GFF data to be filtered.
#' @param cContigToBeRemoved Names of contigs for which the data will be removed. Defaults to NULL.
#' @param dnastringsetGenome A DNAStringSet object containing the sequence for each contig.
#' @param nContigMinSize Minimum size for contigs to keep. Contigs with a size below this value will be removed. Defaults to -1 (= no filter).
#' @param listPctSeqByContig List containing, for each strand, the percentage of sequencing for each contig.
#' This list must be composed of 2 dataframes (one by strand) called f_strand and r_strand.
#' In each dataframe, "refName" column returning names of contigs and "seqPct" column returning percentage of sequencing.
#' @param nContigMinPctOfSeq Minimum percentage of sequencing for contigs to keep.
#' Contigs with a percentage below this value will be removed. Defaults to 95.
#' @param listMeanCovByContig List containing, for each strand, the mean of coverage for each contig.
#' This list must be composed of 2 dataframes (one by strand) called f_strand and r_strand.
#' In each dataframe, "refName" column returning names of contigs and "mean_coverage" column returning mean of coverage.
#' @param nContigMinCoverage Minimum mean coverage for contigs to keep.
#' Contigs with a mean coverage below this value will be removed. Defaults to 20.
#' @return A list with filtered gposModBasePos and filtered grangesModPos.
#' @keywords FiltContig
#' @importFrom BiocGenerics which
#' @export
#' @examples
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' #Preparing a gposPacBioCSV and a grangesPacBioGFF datasets
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
#'
#' #Preparing ParamByStrand Lists
#' myPct_seq_csv <- GetSeqPctByContig(myGposPacBioCSV, grangesGenome = myGrangesGenome)
#' myMean_cov_list <- GetMeanParamByContig(grangesData = myGposPacBioCSV,
#'                                         dnastringsetGenome = myGenome,
#'                                         cParamName = "coverage")
#'
#' #Filtering
#' myFiltered_data <- FiltContig(myGposPacBioCSV, myGrangesPacBioGFF, cContigToBeRemoved = NULL,
#'                               dnastringsetGenome = myGenome, nContigMinSize = 1000,
#'                               listPctSeqByContig = myPct_seq_csv, nContigMinPctOfSeq = 95,
#'                               listMeanCovByContig = myMean_cov_list, nContigMinCoverage = 20)
#' myFiltered_data$csv
#' myFiltered_data$gff
FiltContig <- function(gposModBasePos,
                         grangesModPos,
                         cContigToBeRemoved = NULL,
                         dnastringsetGenome,
                         nContigMinSize = -1,
                         listPctSeqByContig,
                         nContigMinPctOfSeq = 95,
                         listMeanCovByContig,
                         nContigMinCoverage = 20) {
  if(!is.null(cContigToBeRemoved)) {
    #remove contigs manually according to a vector of contigs to be removed
    gposModBasePos <- gposModBasePos[which(!seqnames(gposModBasePos) %in% cContigToBeRemoved),]
    grangesModPos <- grangesModPos[which(!seqnames(grangesModPos) %in% cContigToBeRemoved),]
  }

  if(nContigMinSize > 0) {
    #remove contigs below x size (bp)
    sca_to_remove <- dnastringsetGenome@ranges@NAMES[which(dnastringsetGenome@ranges@width < nContigMinSize)]
    gposModBasePos <- gposModBasePos[which(!seqnames(gposModBasePos) %in% sca_to_remove),]
    grangesModPos <- grangesModPos[which(!seqnames(grangesModPos) %in% sca_to_remove),]
  }

  if(nContigMinPctOfSeq > 0) {
    #remove contigs below x % of sequencing for each strand
    sca_to_remove <- unique( as.character(
      listPctSeqByContig$f_strand[which(listPctSeqByContig$f_strand$seqPct < nContigMinPctOfSeq),"refName"],
      listPctSeqByContig$r_strand[which(listPctSeqByContig$r_strand$seqPct < nContigMinPctOfSeq),"refName"]
    ) )
    gposModBasePos <- gposModBasePos[which(!seqnames(gposModBasePos) %in% sca_to_remove),]
    grangesModPos <- grangesModPos[which(!seqnames(grangesModPos) %in% sca_to_remove),]
  }

  if(nContigMinCoverage > 0) {
    #remove contigs below x mean coverage for each strand
    sca_to_remove <- unique( as.character(
      listMeanCovByContig$f_strand[which(listMeanCovByContig$f_strand$mean_coverage < nContigMinCoverage),"refName"],
      listMeanCovByContig$r_strand[which(listMeanCovByContig$r_strand$mean_coverage < nContigMinCoverage),"refName"]
    ) )
    gposModBasePos <- gposModBasePos[which(!seqnames(gposModBasePos) %in% sca_to_remove),]
    grangesModPos <- grangesModPos[which(!seqnames(grangesModPos) %in% sca_to_remove),]
  }

  return(list(ModBase=gposModBasePos,
              Mod=grangesModPos))
}

#' FiltParam Function (Filter)
#'
#' Filter out modifications which have a chosen parameter that do not reach criterias of selection.
#' @param grangesModPos A GRanges object containing Modifications positions data to be filtered.
#' @param cParamNameForFilter A character vector giving the name of the parameter to be filtered.
#' Must correspond to the name of one column in the object provided with grangesModPos.
#' @param lFiltParam If TRUE, remove modifications which have values of the given parameter that are not included
#' in the intervals provided with "nFiltParamLoBoundaries" and "nFiltParamUpBoundaries". Defaults to FALSE.
#' @param nFiltParamLoBoundaries A numeric vector returning the lower boundaries of intervals.
#' Must have the same length as "nFiltParamUpBoundaries". Defaults to NULL.
#' @param nFiltParamUpBoundaries A numeric vector returning the upper boundaries of intervals.
#' Must have the same length as "nFiltParamLoBoundaries". Defaults to NULL.
#' @param cFiltParamBoundariesToInclude A character vector describing which interval boundaries
#' must be included in the intervals provided. Can be "upperOnly" (only upper boundaries), "lowerOnly" (only lower boundaries),
#' "both" (both upper and lower boundaries) or "none" (do not include upper and lower boundaries).
#' If NULL, both upper and lower boundaries will be included (= "both"). Defaults to NULL.
#' cFiltParamBoundariesToInclude = NULL #can be "upperOnly","lowerOnly","both", "none' (NULL = "both" for all)
#' @param listMeanParamByContig List containing, for each strand, the mean of a given parameter for each contig.
#' This list must be composed of 2 dataframes (one by strand) called f_strand and r_strand.
#' In each dataframe, "refName" column returning names of contigs and "mean_"[parameter name] column returning the mean of the given parameter.
#' If not NULL, remove contigs that are too far away from the mean of the Parameter of all contigs (which are not included in
#' the interval centered on the mean) in the list provided. Defaults to NULL.
#' @param nContigFiltParamLoBound A numeric value to be removed from the mean of the given parameter of all contigs (calculates the
#' lower bound of the interval centered on the mean). Defaults to NULL.
#' @param nContigFiltParamUpBound A numeric value to be added to the mean of the given parameter of all contigs (calculates the
#' upper bound of the interval centered on the mean). Defaults to NULL.
#' @return A filtered grangesModPos object.
#' @keywords FiltParam
#' @importFrom data.table between
#' @export
#' @examples
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#'
#' #Preparing a grangesPacBioGFF dataset
#' myGrangesPacBioGFF <-
#'    ImportPacBioGFF(cPacBioGFFPath = system.file(package="DNAModAnnot", "extdata",
#'                                    "ptetraurelia.modifications.sca171819.gff"),
#'                    cNameModToExtract = "m6A",
#'                    cModNameInOutput = "6mA",
#'                    cContigToBeAnalyzed = c("scaffold51_17", "scaffold51_18", "scaffold51_19"))
#'
#' #Preparing ParamByStrand List
#' myMean_fra_list <- GetMeanParamByContig(grangesData = myGrangesPacBioGFF,
#'                                         dnastringsetGenome = myGenome,
#'                                         cParamName = "frac")
#'
#' #Filtering Out Modif with Frac 20% above/below the mean Frac
#' myGRangesPacBioGFF_filt <- FiltParam(grangesModPos = myGrangesPacBioGFF,
#'                                      cParamNameForFilter = "frac",
#'                                         listMeanParamByContig = myMean_fra_list,
#'                                         nContigFiltParamLoBound = 0.2,
#'                                         nContigFiltParamUpBound = 0.2)
#' myGRangesPacBioGFF_filt
#'
#' #Filtering Out Modif with Frac below 5%
#' myGRangesPacBioGFF_filt <- FiltParam(grangesModPos = myGrangesPacBioGFF,
#'                                      cParamNameForFilter = "frac",
#'                                         listMeanParamByContig = myMean_fra_list,
#'                                         lFiltParam = TRUE,
#'                                         nFiltParamLoBoundaries = 0.05,
#'                                         nFiltParamUpBoundaries = 1.00)
#' myGRangesPacBioGFF_filt
#'
#' #Keeping Modif with Frac between 40% and 60%;
#' #or between 90% and 100%  (adapted to genome of diploid organisms)
#' myGRangesPacBioGFF_filt <- FiltParam(grangesModPos = myGrangesPacBioGFF,
#'                                      cParamNameForFilter = "frac",
#'                                         lFiltParam = TRUE,
#'                                         nFiltParamLoBoundaries = c(0.40, 0.90),
#'                                         nFiltParamUpBoundaries = c(0.60, 1.00))
FiltParam <- function(grangesModPos,
                      cParamNameForFilter,
                      lFiltParam = FALSE,
                      nFiltParamLoBoundaries = NULL,
                      nFiltParamUpBoundaries = NULL,
                      cFiltParamBoundariesToInclude = NULL,
                      listMeanParamByContig = NULL,
                      nContigFiltParamLoBound = NULL,
                      nContigFiltParamUpBound = NULL) {
  cParName <- cParamNameForFilter

  if (!is.null(listMeanParamByContig)) {
    #remove Modifs in contigs that have a mean parameter too far from the mean parameter of all contigs
    low_b = mean(mcols(grangesModPos)[, cParName]) - nContigFiltParamLoBound
    up_b  = mean(mcols(grangesModPos)[, cParName]) + nContigFiltParamUpBound

    sca_to_remove <- unique( c(
      which(!between(listMeanParamByContig$f_strand[[paste0("mean_", cParName)]], lower = low_b, upper = up_b) ),
      which(!between(listMeanParamByContig$r_strand[[paste0("mean_", cParName)]], lower = low_b, upper = up_b) )
    ) )
    sca_to_remove <- listMeanParamByContig$f_strand[sca_to_remove, "refName"]
    grangesModPos <- grangesModPos[which(!seqnames(grangesModPos) %in% sca_to_remove)]
  }

  if (lFiltParam) {
    #keep 6mA that are included in the intervals of parameter provided
    if(length(nFiltParamLoBoundaries)!=length(nFiltParamUpBoundaries)){
      print("ERROR: there must be an equal number of lower and upper bounds.")
    } else {
      mod_to_keep <- unique( unlist(
        sapply(1:length(nFiltParamLoBoundaries),
               function(i){
                 if(cFiltParamBoundariesToInclude[i] == "both" || is.null(cFiltParamBoundariesToInclude)){
                   which(mcols(grangesModPos)[, cParName] >= nFiltParamLoBoundaries[i] & mcols(grangesModPos)[, cParName] <= nFiltParamUpBoundaries[i])
                 } else if(cFiltParamBoundariesToInclude[i]=="upperOnly"){
                   which(mcols(grangesModPos)[, cParName] >  nFiltParamLoBoundaries[i] & mcols(grangesModPos)[, cParName] <= nFiltParamUpBoundaries[i])
                 } else if(cFiltParamBoundariesToInclude[i]=="lowerOnly"){
                   which(mcols(grangesModPos)[, cParName] >= nFiltParamLoBoundaries[i] & mcols(grangesModPos)[, cParName] < nFiltParamUpBoundaries[i])
                 } else if(cFiltParamBoundariesToInclude[i]=="none"){
                   which(mcols(grangesModPos)[, cParName] >  nFiltParamLoBoundaries[i] & mcols(grangesModPos)[, cParName] < nFiltParamUpBoundaries[i])
                 } else {
                   print("ERROR: cFiltParamBoundariesToInclude value not recognized.")
                 }
               } )
      ))
      grangesModPos <- grangesModPos[as.numeric(mod_to_keep), ]
    }
  }

  return(grangesModPos)
}



#' FiltFdrBased Function (Filter)
#'
#' Filter out modifications which have a parameter (tested with FDR estimations) that do not reach criterias of selection.
#' @param grangesModPosWithSeq A List of GRanges object, with the sequence for each motif associated to the modification,
#' and containing PacBio GFF data to be filtered.
#' @param listFdrEstByThrIpdRatio A list of thresholds on ipdRatio for each motif associated to the modification.
#' Defaults to NULL.
#' @param listFdrEstByThrScore A list of thresholds on score for each motif associated to the modification.
#' Defaults to NULL.
#' @return A filtered grangesModPosWithSeq.
#' @keywords FiltFdrBased
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
#'                                    "ptetraurelia.modifications.sca171819.gff"),
#'                    cNameModToExtract = "m6A",
#'                    cModNameInOutput = "6mA",
#'                    cContigToBeAnalyzed = c("scaffold51_17", "scaffold51_18", "scaffold51_19"))
#'
#' myMotif_pct_and_GRangesList <- ExtractListModPosByModMotif(grangesModPos=myGrangesPacBioGFF,
#'                                                         grangesGenome = myGrangesGenome,
#'                                                         dnastringsetGenome = myGenome,
#'                                                         nUpstreamBpToAdd = 0,
#'                                                         nDownstreamBpToAdd = 1,
#'                                                         nModMotifMinProp = 0.05,
#'                                                         cBaseLetterForMod = "A",
#'                                                         cModNameInOutput = "6mA")
#' #Filtering
#' myPosMod_GRangesbyMotif_filt <-
#'    FiltFdrBased(grangesModPosWithSeq = myMotif_pct_and_GRangesList$GRangesbyMotif,
#'                 listFdrEstByThrIpdRatio = list(2.1,2.1, 2.1),
#'                 listFdrEstByThrScore= list(42,42, 42))
#' myPosMod_GRangesbyMotif_filt
FiltFdrBased <- function(grangesModPosWithSeq,
                         listFdrEstByThrIpdRatio=NULL,
                         listFdrEstByThrScore=NULL) {
  if (!is.null(listFdrEstByThrIpdRatio)) {
    #remove Modifs that have an ipdRatio below FDR-based limit for each motif
    names_v <- names(grangesModPosWithSeq)
    grangesModPosWithSeq <- lapply(1:length(grangesModPosWithSeq), function(i) {
      return( subset(grangesModPosWithSeq[[i]],
                     grangesModPosWithSeq[[i]]$ipdRatio >=  listFdrEstByThrIpdRatio[[i]]) )
    } )
    names(grangesModPosWithSeq) <- names_v
  }

  if (!is.null(listFdrEstByThrScore)) {
    #remove Modifs that have an score below FDR-based limit for each motif
    names_v <- names(grangesModPosWithSeq)
    grangesModPosWithSeq <- lapply(1:length(grangesModPosWithSeq), function(i) {
      return( subset(grangesModPosWithSeq[[i]],
                     grangesModPosWithSeq[[i]]$score >=  listFdrEstByThrScore[[i]]) )
    } )
    names(grangesModPosWithSeq) <- names_v
  }

  return(grangesModPosWithSeq)
}

#' FiltPacBio Function (Filter)
#'
#' Filter out data from contigs or Modifications that do not reach criterias of selection.
#' @param grangesPacBioGFF A GRanges object containing PacBio GFF data to be filtered OR A List of GRanges object, with the sequence for each motif associated to the modification,
#' and containing PacBio GFF data to be filtered.
#' @param gposPacBioCSV An UnStitched GPos object containing PacBio CSV data to be filtered. Defaulst to NULL.
#' @param cContigToBeRemoved Names of contigs for which the data will be removed.
#' gposPacBioCSV must be provided if using this argument. Defaults to NULL.
#' @param dnastringsetGenome A DNAStringSet object containing the sequence for each contig.
#' @param nContigMinSize Minimum size for contigs to keep.
#' Contigs with a size below this value will be removed.
#' gposPacBioCSV must be provided if using this argument. Defaults to -1 (= no filter).
#' @param listPctSeqByContig List containing, for each strand, the percentage of sequencing for each contig.
#' This list must be composed of 2 dataframes (one by strand) called f_strand and r_strand.
#' In each dataframe, "refName" column returning names of contigs and "seqPct" column returning percentage of sequencing.
#' gposPacBioCSV must be provided if using this argument.
#' @param nContigMinPctOfSeq Minimum percentage of sequencing for contigs to keep.
#' Contigs with a percentage below this value will be removed.
#' gposPacBioCSV must be provided if using this argument. Defaults to 95.
#' @param listMeanCovByContig List containing, for each strand, the mean of coverage for each contig.
#' This list must be composed of 2 dataframes (one by strand) called f_strand and r_strand.
#' In each dataframe, "refName" column returning names of contigs and "mean_coverage" column returning mean of coverage.
#' gposPacBioCSV must be provided if using this argument.
#' @param nContigMinCoverage Minimum mean coverage for contigs to keep.
#' Contigs with a mean coverage below this value will be removed.
#' gposPacBioCSV must be provided if using this argument. Defaults to 20.
#' @param cParamNameForFilter A character vector giving the name of the parameter to be filtered.
#' Must correspond to the name of one column in the object provided with grangesModPos.
#' @param lFiltParam If TRUE, remove modifications which have values of the given parameter that are not included
#' in the intervals provided with "nFiltParamLoBoundaries" and "nFiltParamUpBoundaries". Defaults to FALSE.
#' @param nFiltParamLoBoundaries A numeric vector returning the lower boundaries of intervals.
#' Must have the same length as "nFiltParamUpBoundaries". Defaults to NULL.
#' @param nFiltParamUpBoundaries A numeric vector returning the upper boundaries of intervals.
#' Must have the same length as "nFiltParamLoBoundaries". Defaults to NULL.
#' @param cFiltParamBoundariesToInclude A character vector describing which interval boundaries
#' must be included in the intervals provided. Can be "upperOnly" (only upper boundaries), "lowerOnly" (only lower boundaries),
#' "both" (both upper and lower boundaries) or "none" (do not include upper and lower boundaries).
#' If NULL, both upper and lower boundaries will be included (= "both"). Defaults to NULL.
#' cFiltParamBoundariesToInclude = NULL #can be "upperOnly","lowerOnly","both", "none' (NULL = "both" for all)
#' @param listMeanParamByContig List containing, for each strand, the mean of a given parameter for each contig.
#' This list must be composed of 2 dataframes (one by strand) called f_strand and r_strand.
#' In each dataframe, "refName" column returning names of contigs and "mean_"[parameter name] column returning the mean of the given parameter.
#' If not NULL, remove contigs that are too far away from the mean of the Parameter of all contigs (which are not included in
#' the interval centered on the mean) in the list provided. Defaults to NULL.
#' @param nContigFiltParamLoBound A numeric value to be removed from the mean of the given parameter of all contigs (calculates the
#' lower bound of the interval centered on the mean). Defaults to NULL.
#' @param nContigFiltParamUpBound A numeric value to be added to the mean of the given parameter of all contigs (calculates the
#' upper bound of the interval centered on the mean). Defaults to NULL.
#' @param nModMinIpdRatio Minimum ipdRatio for all Modifications to be kept.
#' Modifications with an ipdRatio below this value will be removed. Defaults to NULL (no filter).
#' @param nModMinScore Minimum score for all Modifications to be kept.
#' Modifications with a score below this value will be removed. Defaults to NULL (no filter).
#' @param nModMinCoverage Minimum coverage for all Modifications to be kept.
#' Modifications with a coverage below this value will be removed. Defaults to NULL (no filter).
#' @param listFdrEstByThrIpdRatio A list of thresholds on ipdRatio for each motif associated to the modification.
#' @param listFdrEstByThrScore A list of thresholds on score for each motif associated to the modification.
#' @return A list with filtered gposPacBioCSV and filtered gposPacBioGFF.
#' @keywords FiltPacBio
#' @importFrom BiocGenerics which
#' @importFrom data.table between
#' @export
#' @examples
#' #loading genome
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' #Preparing a gposPacBioCSV and a grangesPacBioGFF datasets
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
#'
#' #Preparing ParamByStrand Lists
#' myPct_seq_csv <- GetSeqPctByContig(myGposPacBioCSV, grangesGenome = myGrangesGenome)
#' myMean_cov_list <- GetMeanParamByContig(grangesData = myGposPacBioCSV,
#'                                         dnastringsetGenome = myGenome,
#'                                         cParamName = "coverage")
#'
#' #Filtering
#' myFiltered_data <- FiltPacBio(grangesPacBioGFF=myGrangesPacBioGFF,
#'                               gposPacBioCSV=myGposPacBioCSV, cContigToBeRemoved = NULL,
#'                               dnastringsetGenome = myGenome, nContigMinSize = 1000,
#'                               listPctSeqByContig = myPct_seq_csv, nContigMinPctOfSeq = 95,
#'                               listMeanCovByContig = myMean_cov_list, nContigMinCoverage = 20)
#' myFiltered_data$csv
#' myFiltered_data$gff
FiltPacBio <- function(grangesPacBioGFF, #PacBioGFF: granges or grangesList
                       gposPacBioCSV = NULL,
                       cContigToBeRemoved = NULL,
                       dnastringsetGenome,
                       nContigMinSize = -1,
                       listPctSeqByContig,
                       nContigMinPctOfSeq = -1,
                       listMeanCovByContig,
                       nContigMinCoverage = -1,
                       cParamNameForFilter = NULL,
                       listMeanParamByContig = NULL,
                       nContigFiltParamLoBound = NULL,
                       nContigFiltParamUpBound = NULL,
                       lFiltParam = FALSE,
                       nFiltParamLoBoundaries = NULL,
                       nFiltParamUpBoundaries = NULL,
                       cFiltParamBoundariesToInclude = NULL,
                       nModMinIpdRatio = NULL,
                       nModMinScore = NULL,
                       nModMinCoverage = NULL,
                       listFdrEstByThrIpdRatio = NULL, #by motif
                       listFdrEstByThrScore = NULL #by motif
                       ) {

  if(class(grangesPacBioGFF) == "GRanges"){
    if(is.null(gposPacBioCSV) && any(c(!is.null(cContigToBeRemoved),
                                       nContigMinSize > 0,
                                       nContigMinPctOfSeq > 0,
                                       nContigMinCoverage > 0))) {
      print("Error: gposPacBioCSV argument must be provided with the chosen filters.")

    } else {
      listFiltData <- FiltContig(gposModBasePos = gposPacBioCSV,
                                 grangesModPos = grangesPacBioGFF,
                                 cContigToBeRemoved = cContigToBeRemoved,
                                 dnastringsetGenome = dnastringsetGenome,
                                 nContigMinSize = nContigMinSize,
                                 listPctSeqByContig = listPctSeqByContig,
                                 nContigMinPctOfSeq = nContigMinPctOfSeq,
                                 listMeanCovByContig = listMeanCovByContig,
                                 nContigMinCoverage = nContigMinCoverage)
      gposPacBioCSV <- listFiltData$ModBase
      grangesPacBioGFF <- listFiltData$Mod

      if(!is.null(cParamNameForFilter)){
        grangesPacBioGFF <- FiltParam(grangesModPos = grangesPacBioGFF,
                                      cParamNameForFilter = cParamNameForFilter,
                                      listMeanParamByContig = listMeanParamByContig,
                                      nContigFiltParamLoBound = nContigFiltParamLoBound,
                                      nContigFiltParamUpBound = nContigFiltParamUpBound,
                                      lFiltParam = lFiltParam,
                                      nFiltParamLoBoundaries = nFiltParamLoBoundaries,
                                      nFiltParamUpBoundaries = nFiltParamUpBoundaries,
                                      cFiltParamBoundariesToInclude = cFiltParamBoundariesToInclude)
      }

      if(!is.null(nModMinIpdRatio)){
        grangesPacBioGFF <- grangesPacBioGFF[which(grangesPacBioGFF$ipdRatio >= nModMinIpdRatio),]
      }
      if(!is.null(nModMinScore)){
        grangesPacBioGFF <- grangesPacBioGFF[which(grangesPacBioGFF$score    >= nModMinScore),]
      }
      if(!is.null(nModMinCoverage)){
        grangesPacBioGFF <- grangesPacBioGFF[which(grangesPacBioGFF$coverage >= nModMinCoverage),]
      }


    }

    if(!is.null(gposPacBioCSV)) { seqlevels(gposPacBioCSV) <- seqlevelsInUse(gposPacBioCSV) }
    seqlevels(grangesPacBioGFF) <- seqlevelsInUse(grangesPacBioGFF)

  } else if(class(grangesPacBioGFF) %in% c("CompressedGRangesList", "GRangesList")) {
    grangesPacBioGFF <- FiltFdrBased(grangesModPosWithSeq = grangesPacBioGFF,
                                     listFdrEstByThrIpdRatio = listFdrEstByThrIpdRatio,
                                     listFdrEstByThrScore = listFdrEstByThrScore)

    for(x in 1:length(grangesPacBioGFF)){
      seqlevels(grangesPacBioGFF[[x]]) <- seqlevelsInUse(grangesPacBioGFF[[x]])
    }

  } else {
    print("Error: grangesPacBioGFF argument must be a GRanges/CompressedGranges or GRangesList/CompressedGRangesList object.")
  }

  return(list(csv=gposPacBioCSV,
              gff=grangesPacBioGFF))
}

#' FiltModDeepSignal Function (Filter)
#'
#' Filter out data from contigs or Modifications that do not reach criterias of selection.
#' Can also be used to obtain a gposDeepSignalMod object by simply filtering target sites which have a fraction above 0.
#' @param gposDeepSignalModBase An UnStitched GPos object containing DeepSignal modification target sites data to be filtered. Defaults to NULL.
#' @param gposDeepSignalMod An UnStitched GPos object containing DeepSignal modified sites data to be filtered. Defaults to NULL.
#' @param cContigToBeRemoved Names of contigs for which the data will be removed.
#' gposPacBioCSV must be provided if using this argument. Defaults to NULL.
#' @param dnastringsetGenome A DNAStringSet object containing the sequence for each contig.
#' @param nContigMinSize Minimum size for contigs to keep.
#' Contigs with a size below this value will be removed.
#' gposPacBioCSV must be provided if using this argument. Defaults to -1 (= no filter).
#' @param listPctSeqByContig List containing, for each strand, the percentage of sequencing for each contig.
#' This list must be composed of 2 dataframes (one by strand) called f_strand and r_strand.
#' In each dataframe, "refName" column returning names of contigs and "seqPct" column returning percentage of sequencing.
#' gposPacBioCSV must be provided if using this argument.
#' @param nContigMinPctOfSeq Minimum percentage of sequencing for contigs to keep.
#' Contigs with a percentage below this value will be removed.
#' gposPacBioCSV must be provided if using this argument. Defaults to 95.
#' @param listMeanCovByContig List containing, for each strand, the mean of coverage for each contig.
#' This list must be composed of 2 dataframes (one by strand) called f_strand and r_strand.
#' In each dataframe, "refName" column returning names of contigs and "mean_coverage" column returning mean of coverage.
#' gposPacBioCSV must be provided if using this argument.
#' @param nContigMinCoverage Minimum mean coverage for contigs to keep.
#' Contigs with a mean coverage below this value will be removed.
#' gposPacBioCSV must be provided if using this argument. Defaults to 20.
#' @param cParamNameForFilter A character vector giving the name of the parameter to be filtered.
#' Must correspond to the name of one column in the object provided with grangesModPos.
#' @param lFiltParam If TRUE, remove modifications which have values of the given parameter that are not included
#' in the intervals provided with "nFiltParamLoBoundaries" and "nFiltParamUpBoundaries". Defaults to FALSE.
#' @param nFiltParamLoBoundaries A numeric vector returning the lower boundaries of intervals.
#' Must have the same length as "nFiltParamUpBoundaries". Defaults to NULL.
#' @param nFiltParamUpBoundaries A numeric vector returning the upper boundaries of intervals.
#' Must have the same length as "nFiltParamLoBoundaries". Defaults to NULL.
#' @param cFiltParamBoundariesToInclude A character vector describing which interval boundaries
#' must be included in the intervals provided. Can be "upperOnly" (only upper boundaries), "lowerOnly" (only lower boundaries),
#' "both" (both upper and lower boundaries) or "none" (do not include upper and lower boundaries).
#' If NULL, both upper and lower boundaries will be included (= "both"). Defaults to NULL.
#' cFiltParamBoundariesToInclude = NULL #can be "upperOnly","lowerOnly","both", "none' (NULL = "both" for all)
#' @param listMeanParamByContig List containing, for each strand, the mean of a given parameter for each contig.
#' This list must be composed of 2 dataframes (one by strand) called f_strand and r_strand.
#' In each dataframe, "refName" column returning names of contigs and "mean_"[parameter name] column returning the mean of the given parameter.
#' If not NULL, remove contigs that are too far away from the mean of the Parameter of all contigs (which are not included in
#' the interval centered on the mean) in the list provided. Defaults to NULL.
#' @param nContigFiltParamLoBound A numeric value to be removed from the mean of the given parameter of all contigs (calculates the
#' lower bound of the interval centered on the mean). Defaults to NULL.
#' @param nContigFiltParamUpBound A numeric value to be added to the mean of the given parameter of all contigs (calculates the
#' upper bound of the interval centered on the mean). Defaults to NULL.
#' @param nModMinCoverage Minimum coverage for all Modifications to be kept.
#' Modifications with a coverage below this value will be removed. Defaults to NULL (no filter).
#' @keywords FiltModDeepSignal
#' @export
#' @examples
#' #Loading Nanopore data
#' myDeepSignalModPath <- system.file(
#'   package="DNAModAnnot", "extdata",
#'   "FAB39088-288418386-Chr1.CpG.call_mods.frequency.tsv")
#' mygposDeepSignalModBase <- ImportDeepSignalModFrequency(cDeepSignalModPath=myDeepSignalModPath,
#'                                                         lSortGPos=TRUE,
#'                                                         cContigToBeAnalyzed = "all")
#' mygposDeepSignalModBase
#'
#' #Filtering
#' mygposDeepSignalMod <- FiltDeepSignal(gposDeepSignalModBase = mygposDeepSignalModBase,
#'                                       cParamNameForFilter = "frac",
#'                                       lFiltParam = TRUE,
#'                                       nFiltParamLoBoundaries = 0,
#'                                       nFiltParamUpBoundaries = 1,
#'                                       cFiltParamBoundariesToInclude = "upperOnly")$Mod
#' mygposDeepSignalMod
FiltDeepSignal <- function(gposDeepSignalModBase = NULL,
                              gposDeepSignalMod=NULL,
                              cContigToBeRemoved = NULL,
                              dnastringsetGenome,
                              nContigMinSize = -1,
                              listPctSeqByContig,
                              nContigMinPctOfSeq = -1,
                              listMeanCovByContig,
                              nContigMinCoverage = -1,
                           cParamNameForFilter = NULL,
                           lFiltParam = FALSE,
                           nFiltParamLoBoundaries = NULL,
                           nFiltParamUpBoundaries = NULL,
                           cFiltParamBoundariesToInclude = NULL,
                           listMeanParamByContig = NULL,
                           nContigFiltParamLoBound = NULL,
                           nContigFiltParamUpBound = NULL,
                              nModMinCoverage = NULL
) {

  if(is.null(gposDeepSignalMod)) { gposDeepSignalMod <- gposDeepSignalModBase }


  listFiltData <- FiltContig(gposModBasePos = gposDeepSignalModBase,
                             grangesModPos = gposDeepSignalMod,
                             cContigToBeRemoved = cContigToBeRemoved,
                             dnastringsetGenome = dnastringsetGenome,
                             nContigMinSize = nContigMinSize,
                             listPctSeqByContig = listPctSeqByContig,
                             nContigMinPctOfSeq = nContigMinPctOfSeq,
                             listMeanCovByContig = listMeanCovByContig,
                             nContigMinCoverage = nContigMinCoverage)
  gposDeepSignalModBase <- listFiltData$ModBase
  gposDeepSignalMod <- listFiltData$Mod

  if(!is.null(cParamNameForFilter)){
    gposDeepSignalMod <- FiltParam(grangesModPos = gposDeepSignalMod,
                                  cParamNameForFilter = cParamNameForFilter,
                                  listMeanParamByContig = listMeanParamByContig,
                                  nContigFiltParamLoBound = nContigFiltParamLoBound,
                                  nContigFiltParamUpBound = nContigFiltParamUpBound,
                                  lFiltParam = lFiltParam,
                                  nFiltParamLoBoundaries = nFiltParamLoBoundaries,
                                  nFiltParamUpBoundaries = nFiltParamUpBoundaries,
                                  cFiltParamBoundariesToInclude = cFiltParamBoundariesToInclude)
  }

  if(!is.null(nModMinCoverage)){
    gposDeepSignalMod <- gposDeepSignalMod[which(gposDeepSignalMod$coverage >= nModMinCoverage),]
  }

  if(!is.null(gposDeepSignalModBase)) { seqlevels(gposDeepSignalModBase) <- seqlevelsInUse(gposDeepSignalModBase) }
  seqlevels(gposDeepSignalMod) <- seqlevelsInUse(gposDeepSignalMod)

  return(list(ModBase=gposDeepSignalModBase,
              Mod=gposDeepSignalMod))
}

