### =========================================================================
### Functions from DaLoad category
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
### GetGposCenterFromGRanges
### PredictMissingAnnotation
### ImportPacBioGFF
### ImportPacBioCSV
### GetGenomeGRanges
### ImportDeepSignalModFrequency
###
### -------------------------------------------------------------------------

#' @importFrom graphics abline axis barplot boxplot grid hist layout legend lines mtext par plot polygon segments text title
#' @importFrom grDevices dev.print pdf png rgb
#' @importFrom methods as
#' @importFrom stats aggregate density quantile
#' @importFrom utils write.csv write.table
NULL

#' GetGposCenterFromGRanges Function (DaLoad)
#'
#' Retrieve, in a GPos object, the positions of the center of the ranges from a GRanges object.
#' @param grangesData A GRanges object to be converted as a GPos object. The center of each range will be processed in the GPos object output.
#' @keywords GetGposCenterFromGRanges
#' @export
#' @examples
#' myRanges <- as(c("chrI:300-500:+", "chrI:308-680:+", "chrII:30-550:-"), "GRanges")
#' myCenterPos <- GetGposCenterFromGRanges(grangesData = myRanges)
#' myCenterPos
GetGposCenterFromGRanges <- function(grangesData) {
  grangesData$center <- start(grangesData) + ceiling((end(grangesData) - start(grangesData)) / 2)
  grangesData <- GPos(
    seqnames = seqnames(grangesData),
    pos = grangesData$center,
    strand = strand(grangesData),
    stitch = FALSE
  )
  return(grangesData)
}

#' PredictMissingAnnotation Function (DaLoad)
#'
#' Complete annotation with features, such as "intergenic", "antisense_strand_of_gene" or "exon"|"intron", using available features in the annotation.
#' @param grangesAnnotations A GRanges object with the annotation to be completed.
#' @param grangesGenome A GRanges object with number and width of contigs (both strands).
#' @param cFeaturesColName The name of the column containing feature type annotation ("gene", "exon", "mRNA"...). Defaults to "type".
#' @param cGeneCategories The name of the categories considered as genes in the column containing feature type annotation. Defaults to c("gene").
#' @param lAddIntronRangesUsingExon If TRUE, uses "exon" and "mRNA" categories to add "intron" if "intron" is missing,
#' or uses "intron" and "mRNA" categories to add "exon" if "exon" is missing.
#' This will return an error if 2 categories among "mRNA", "exon" and "intron" are missing. Defaults to FALSE.
#' @keywords PredictMissingAnnotation
#' @importFrom IRanges extractList findOverlapPairs
#' @importFrom GenomeInfoDb sortSeqlevels
#' @export
#' @examples
#' # Loading genome data
#' myGenome <- Biostrings::readDNAStringSet(system.file(
#'   package = "DNAModAnnot", "extdata",
#'   "ptetraurelia_mac_51_sca171819.fa"
#' ))
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' # Loading annotation data
#' myAnnotations <-
#'   rtracklayer::readGFFAsGRanges(system.file(
#'     package = "DNAModAnnot", "extdata",
#'     "ptetraurelia_mac_51_annotation_v2.0_sca171819.gff3"
#'   ))
#'
#' # Completing annotation data
#' levels(myAnnotations$type)
#' myAnnotations <- PredictMissingAnnotation(
#'   grangesAnnotations = myAnnotations,
#'   grangesGenome = myGrangesGenome,
#'   cFeaturesColName = "type",
#'   cGeneCategories = c("gene"),
#'   lAddIntronRangesUsingExon = TRUE
#' )
#' levels(myAnnotations$type)
PredictMissingAnnotation <- function(grangesAnnotations, grangesGenome,
                                     cFeaturesColName = "type",
                                     cGeneCategories = c("gene"),
                                     lAddIntronRangesUsingExon = FALSE) {
  colnames(mcols(grangesAnnotations))[which(colnames(mcols(grangesAnnotations)) == cFeaturesColName)] <- "type"

  if (!cGeneCategories == "gene") {
    grangesAnnotations$gene_type <- as.character(NA)
    type_to_modif <- which(grangesAnnotations$type %in% cGeneCategories)
    mcols(grangesAnnotations)[type_to_modif, "gene_type"] <- as.character(mcols(grangesAnnotations)[type_to_modif, cFeaturesColName])
    mcols(grangesAnnotations)[type_to_modif, cFeaturesColName] <- "gene"
  }
  gene_check <- "gene" %in% levels(grangesAnnotations$type)
  if (!gene_check) {
    print("No \"gene\" category recognized. Please, correct the annotation or provide names of sub-categories of genes from this genome annotation.")
  } else {
    # adding intergenic Granges
    intergenic_check <- "intergenic" %in% levels(grangesAnnotations$type)
    if (!intergenic_check) {
      gr_a <- grangesGenome

      gr_b <- subset(grangesAnnotations, type == "gene")
      hits <- findOverlaps(gr_a, gr_b, ignore.strand = TRUE)
      toSubtract <- reduce(extractList(gr_b, as(hits, "List")),
        ignore.strand = TRUE
      )
      ans <- psetdiff(gr_a, toSubtract, ignore.strand = TRUE)
      ans <- unlist(ans)
      ans <- subset(ans, width(ans) > 0L)
      mcols(ans)[cFeaturesColName] <- as.factor("intergenic")

      if (length(grangesAnnotations) >= length(ans)) {
        grangesAnnotations <- c(grangesAnnotations, ans)
      } else {
        grangesAnnotations <- c(ans, grangesAnnotations)
      }

      gr_a <- grangesGenome

      gr_b <- subset(grangesAnnotations, type == "gene")
      pairs <- findOverlapPairs(gr_a, invertStrand(gr_b), ignore.strand = FALSE)
      ans <- pintersect(pairs, ignore.strand = FALSE)
      mcols(ans)[cFeaturesColName] <- as.factor("antisense_strand_of_gene")

      if (length(grangesAnnotations) >= length(ans)) {
        grangesAnnotations <- c(grangesAnnotations, ans)
      } else {
        grangesAnnotations <- c(ans, grangesAnnotations)
      }
    }

    if (lAddIntronRangesUsingExon) {
      # adding exon/intron Granges
      mRNA_check <- "mRNA" %in% levels(grangesAnnotations$type)
      exon_check <- "exon" %in% levels(grangesAnnotations$type)
      intron_check <- "intron" %in% levels(grangesAnnotations$type)

      if (mRNA_check & (exon_check | intron_check)) {
        if (!exon_check) {
          gr_a <- subset(grangesAnnotations, type == "mRNA")
          gr_b <- subset(grangesAnnotations, type == "intron")
          hits <- findOverlaps(gr_a, gr_b, ignore.strand = FALSE)
          toSubtract <- reduce(extractList(gr_b, as(hits, "List")),
            ignore.strand = TRUE
          )
          ans <- psetdiff(gr_a, toSubtract, ignore.strand = TRUE)
          ans <- unlist(ans)
          ans <- subset(ans, width > 0L)
          mcols(ans)[cFeaturesColName] <- as.factor("exon")

          if (length(grangesAnnotations) >= length(ans)) {
            grangesAnnotations <- c(grangesAnnotations, ans)
          } else {
            grangesAnnotations <- c(ans, grangesAnnotations)
          }
        } else if (!intron_check) {
          gr_a <- subset(grangesAnnotations, type == "mRNA")
          gr_b <- subset(grangesAnnotations, type == "exon")
          hits <- findOverlaps(gr_a, gr_b, ignore.strand = FALSE)
          toSubtract <- reduce(extractList(gr_b, as(hits, "List")),
            ignore.strand = TRUE
          )
          ans <- psetdiff(gr_a, toSubtract, ignore.strand = TRUE)
          ans <- unlist(ans)
          ans <- subset(ans, width > 0L)
          mcols(ans)[cFeaturesColName] <- as.factor("intron")

          if (length(grangesAnnotations) >= length(ans)) {
            grangesAnnotations <- c(grangesAnnotations, ans)
          } else {
            grangesAnnotations <- c(ans, grangesAnnotations)
          }
        } else {
          print("\"Exon\"/\"Intron\" categories already detected: no need to fill exon/intron positions.")
        }
      } else {
        print("\"Exon\"/\"Intron\" categories cannot be predicted because of missing annotation (mRNA and (exon or intron)).")
      }
    }
  }
  grangesAnnotations <- sortSeqlevels(grangesAnnotations)
  grangesAnnotations <- sort(grangesAnnotations)
  return(grangesAnnotations)
}

#' ImportPacBioGFF Function (DaLoad)
#'
#' Import PacBio GFF file, extract one modification, rename this modification and convert it
#' as a GRanges object with new colnames similar to PacBio CSV file containing data from all bases sequenced.
#' @param cPacBioGFFPath Path to a PacBio GFF file containing modification detection data.
#' @param cNameModToExtract Name of modification to be extracted.
#' @param cModNameInOutput Name for the extracted modification in the output.
#' @param cContigToBeAnalyzed Names of contigs for which the data will be kept.
#' If NULL, data from all contigs available will be imported. Defaults to NULL.
#' @param lKeepSequence If TRUE, the sequence of the base will be retained in one column.
#' Otherwise, it will be discarded to reduce object size. Defaults to TRUE.
#' @keywords ImportPacBioGFF
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges seqnames mcols mcols<-
#' @importFrom GenomeInfoDb seqlevels seqlevelsInUse seqlevels<-
#' @export
#' @examples
#' # Loading genome data
#' myGenome <- Biostrings::readDNAStringSet(system.file(
#'   package = "DNAModAnnot", "extdata",
#'   "ptetraurelia_mac_51_sca171819.fa"
#' ))
#' names(myGenome)
#'
#' # Loading PacBio data
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
#' myGrangesPacBioGFF
#'
#' # Loading PacBio data for 2 scaffolds only
#' myGrangesPacBioGFF <-
#'   ImportPacBioGFF(
#'     cPacBioGFFPath = system.file(
#'       package = "DNAModAnnot", "extdata",
#'       "ptetraurelia.modifications.sca171819.gff"
#'     ),
#'     cNameModToExtract = "m6A",
#'     cModNameInOutput = "6mA",
#'     cContigToBeAnalyzed = c("scaffold51_18", "scaffold51_19")
#'   )
#' myGrangesPacBioGFF
ImportPacBioGFF <- function(cPacBioGFFPath,
                            cNameModToExtract,
                            cModNameInOutput,
                            cContigToBeAnalyzed = NULL,
                            lKeepSequence = TRUE) {
  grangesPacBioGFF <- import(cPacBioGFFPath)
  if (!is.null(cContigToBeAnalyzed)) {
    grangesPacBioGFF <- subset(grangesPacBioGFF, seqnames(grangesPacBioGFF) %in% cContigToBeAnalyzed)
    seqlevels(grangesPacBioGFF) <- seqlevelsInUse(grangesPacBioGFF)
  }
  grangesPacBioGFF <- subset(grangesPacBioGFF, type == cNameModToExtract)
  grangesPacBioGFF$type <- cModNameInOutput

  grangesPacBioGFF <- subset(grangesPacBioGFF,
    select = c(
      "type",
      "coverage", "score", "IPDRatio",
      "frac", "fracLow", "fracUp",
      "identificationQv"
    )
  )

  names(mcols(grangesPacBioGFF))[1] <- "base"
  names(mcols(grangesPacBioGFF))[4] <- "ipdRatio"

  if (lKeepSequence) {
    grangesPacBioGFF$base <- as.character(grangesPacBioGFF$base)
  } else {
    grangesPacBioGFF$base <- NULL
  }
  grangesPacBioGFF$coverage <- as.integer(grangesPacBioGFF$coverage)
  grangesPacBioGFF$score <- as.integer(grangesPacBioGFF$score)
  grangesPacBioGFF$ipdRatio <- as.numeric(grangesPacBioGFF$ipdRatio)
  grangesPacBioGFF$frac <- as.numeric(grangesPacBioGFF$frac)
  grangesPacBioGFF$fracLow <- as.numeric(grangesPacBioGFF$fracLow)
  grangesPacBioGFF$fracUp <- as.numeric(grangesPacBioGFF$fracUp)
  grangesPacBioGFF$identificationQv <- as.integer(grangesPacBioGFF$identificationQv)

  return(grangesPacBioGFF)
}

#' ImportPacBioCSV Function (DaLoad)
#'
#' Import PacBio CSV file and convert it as an UnStitched GPos object.
#' @param cPacBioCSVPath Path to a PacBio CSV file containing data from all bases sequenced.
#' @param cSelectColumnsToExtract Names of columns to extract from PacBio CSV file. Less there are columns, faster the file will be loaded.
#' The columns "refName", "tpl" and "strand" are mandatory to convert to a GPos object.
#' Defaults to c("refName", "tpl", "strand", "base", "score", "ipdRatio", "coverage")
#' @param lKeepExtraColumnsInGPos If FALSE, only the contig names, start/end positions and strand will be displayed in the resulting GPos object.
#' Defaults to TRUE.
#' @param lSortGPos If TRUE, the GPos object will be sorted before being returned: the function will take a longer time
#' to proceed but the GPos Object will require less memory.
#' @param cContigToBeAnalyzed Names of contigs for which the data will be kept.
#' If NULL, data from all contigs available will be imported. Defaults to NULL.
#' @param lKeepSequence If TRUE, the sequence of the base will be retained in one column.
#' Otherwise, it will be discarded to reduce object size. Defaults to TRUE.
#' @keywords ImportPacBioCSV
#' @importFrom data.table fread setDF
#' @importFrom GenomicRanges GPos
#' @importFrom IRanges IPos
#' @importFrom BiocGenerics sort
#' @export
#' @examples
#' # Loading genome data
#' myGenome <- Biostrings::readDNAStringSet(system.file(
#'   package = "DNAModAnnot", "extdata",
#'   "ptetraurelia_mac_51_sca171819.fa"
#' ))
#' names(myGenome)
#'
#' # Loading PacBio data
#' myGrangesPacBioCSV <-
#'   ImportPacBioCSV(
#'     cPacBioCSVPath = system.file(
#'       package = "DNAModAnnot", "extdata",
#'       "ptetraurelia.bases.sca171819.csv"
#'     ),
#'     cSelectColumnsToExtract = c("refName", "tpl", "strand", "base",
#'                                 "score", "ipdRatio", "coverage"),
#'     lKeepExtraColumnsInGPos = TRUE,
#'     lSortGPos = TRUE,
#'     cContigToBeAnalyzed = names(myGenome)
#'   )
#' myGrangesPacBioCSV
#'
#' # Loading PacBio data for 2 scaffolds only
#' myGrangesPacBioCSV <-
#'   ImportPacBioCSV(
#'     cPacBioCSVPath = system.file(
#'       package = "DNAModAnnot", "extdata",
#'       "ptetraurelia.bases.sca171819.csv"
#'     ),
#'     cSelectColumnsToExtract = c(
#'       "refName", "tpl", "strand", "base",
#'       "score", "ipdRatio", "coverage"
#'     ),
#'     lKeepExtraColumnsInGPos = TRUE,
#'     lSortGPos = TRUE,
#'     cContigToBeAnalyzed = c("scaffold51_18", "scaffold51_19")
#'   )
#' myGrangesPacBioCSV
ImportPacBioCSV <- function(cPacBioCSVPath,
                            cSelectColumnsToExtract = c("refName", "tpl", "strand", "base", "score", "ipdRatio", "coverage"),
                            lKeepExtraColumnsInGPos = TRUE,
                            lSortGPos = TRUE,
                            cContigToBeAnalyzed = NULL,
                            lKeepSequence = TRUE) {
  print("Loading csv file...")
  gposPacBioCSV <- fread(cPacBioCSVPath,
    header = TRUE,
    sep = ",",
    select = cSelectColumnsToExtract
  )

  if (!is.null(cContigToBeAnalyzed)) {
    print("Filtering contigs to be analysed...")
    gposPacBioCSV <- gposPacBioCSV[gposPacBioCSV$refName %in% cContigToBeAnalyzed]
  }

  print("Adjusting strand column...")
  gposPacBioCSV[, strand := as.character(strand)][strand == "0", strand := "+"]
  gposPacBioCSV[, strand := as.character(strand)][strand == "1", strand := "-"]
  gposPacBioCSV[, strand := as.character(strand)][!strand %in% c("+", "-"), strand := "*"]

  if (lKeepExtraColumnsInGPos) {
    print("Creating GPos object with extra columns...")
    gposPacBioCSV <- GPos(
      seqnames = gposPacBioCSV$refName,
      pos = IPos(pos = gposPacBioCSV$tpl),
      strand = gposPacBioCSV$strand,
      setDF(gposPacBioCSV[, -1:-3]),
      stitch = FALSE
    )
    if (!lKeepSequence) {
      gposPacBioCSV$base <- NULL
    }
  } else {
    gposPacBioCSV <- GPos(
      seqnames = gposPacBioCSV$refName,
      pos = IPos(pos = gposPacBioCSV$tpl),
      strand = gposPacBioCSV$strand,
      stitch = FALSE
    )
  }

  if (lSortGPos) {
    print("Sorting GPos Object...")
    gposPacBioCSV <- sort(gposPacBioCSV)
  }

  return(gposPacBioCSV)
}

#' GetGenomeGRanges Function (DaLoad)
#'
#' Return a GRanges object from a DNAStringSet object with ranges from 1 to the sequences width (for both strands).
#' @param dnastringsetGenome A DNAStringSet object to convert as a GRanges object.
#' If sequences have no names, these will be named with "seq" and a number
#' (corresponding to the order in the DNAStringSet object) in the GRanges object.
#' @keywords GetGenomeGRanges
#' @import Biostrings
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics sort
#' @importFrom data.table :=
#' @export
#' @examples
#' mySeqs <- Biostrings::DNAStringSet(c("ACCATTGATTAT", "AATATCGACTA", "GACTAT"))
#' myRanges <- GetGenomeGRanges(mySeqs)
#' myRanges
GetGenomeGRanges <- function(dnastringsetGenome) {
  if (is.null(names(dnastringsetGenome))) {
    names(dnastringsetGenome) <- paste0("seq", 1:length(dnastringsetGenome))
  }
  grangesGenome <- c(
    GRanges(dnastringsetGenome@ranges@NAMES,
      IRanges(start = 1, end = dnastringsetGenome@ranges@width),
      strand = "+"
    ),
    GRanges(dnastringsetGenome@ranges@NAMES,
      IRanges(start = 1, end = dnastringsetGenome@ranges@width),
      strand = "-"
    )
  )
  return(sort(grangesGenome))
}

#' ImportDeepSignalModFrequency Function (DaLoad)
#'
#' Import DeepSignal call_modification_frequency.py output file and convert it as an UnStitched GPos object.
#' @param cDeepSignalModPath Path to a DeepSignal call_modification_frequency.py output file containing data from all target sites.
#' @param cColumnNames Names for each column in the DeepSignal call_modification_frequency.py output file.
#' Should not be changed unless some columns are missing in the file to be imported.
#' Defaults to c("chrom", "pos", "strand", "pos_in_strand",
#' "prob_0_sum", "prob_1_sum", "count_modified", "count_unmodified",
#' "coverage", "modification_frequency", "k_mer")
#' @param cSelectColumnsToExtract Names of columns to extract from DeepSignal call_modification_frequency.py output file.
#' Less there are columns, faster the file will be loaded.
#' The columns "chrom", "pos" and "strand" are mandatory to convert to a GPos object.
#' Defaults to c("chrom", "pos", "strand", "prob_0_sum", "prob_1_sum", "count_modified", "count_unmodified",
#' "coverage", "modification_frequency", "k_mer")
#' @param lSortGPos If TRUE, the GPos object will be sorted before being returned: the function will take a longer time
#' to proceed but the GPos Object will require less memory.
#' @param cContigToBeAnalyzed Names of contigs for which the data will be kept.
#' If NULL, data from all contigs available will be imported. Defaults to NULL.
#' @param lKeepSequence If TRUE, the sequence of the base will be retained in one column.
#' Otherwise, it will be discarded to reduce object size. Defaults to TRUE.
#' @keywords ImportDeepSignalModFrequency
#' @export
#' @examples
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
#' mygposDeepSignalModBase
ImportDeepSignalModFrequency <- function(cDeepSignalModPath,
                                         cColumnNames = c(
                                           "chrom", "pos", "strand", "pos_in_strand",
                                           "prob_0_sum", "prob_1_sum", "count_modified", "count_unmodified",
                                           "coverage", "modification_frequency", "k_mer"
                                         ),
                                         cSelectColumnsToExtract = c(
                                           "chrom", "pos", "strand",
                                           "prob_0_sum", "prob_1_sum", "count_modified", "count_unmodified",
                                           "coverage", "modification_frequency", "k_mer"
                                         ),
                                         lSortGPos = TRUE,
                                         cContigToBeAnalyzed,
                                         lKeepSequence = TRUE) {
  cColumnNames <- gsub(cColumnNames,
    pattern = "modification_frequency", replacement = "frac"
  )
  cSelectColumnsToExtract <- gsub(cSelectColumnsToExtract,
    pattern = "modification_frequency", replacement = "frac"
  )

  print("Loading tsv file...")
  gposDeepSignalMod <- fread(cDeepSignalModPath,
    header = FALSE,
    sep = "\t"
  )
  colnames(gposDeepSignalMod) <- cColumnNames
  gposDeepSignalMod <- gposDeepSignalMod[, colnames(gposDeepSignalMod) %in% cSelectColumnsToExtract, with = FALSE]

  if (cContigToBeAnalyzed != "all") {
    print("Filtering contigs to be analysed...")
    gposDeepSignalMod <- gposDeepSignalMod[gposDeepSignalMod$chrom %in% cContigToBeAnalyzed]
  }

  print("Creating GPos object with extra columns...")
  gposDeepSignalMod <- GPos(
    seqnames = gposDeepSignalMod$chrom,
    pos = IPos(pos = gposDeepSignalMod$pos + 1),
    strand = gposDeepSignalMod$strand,
    setDF(gposDeepSignalMod[, -1:-3]),
    stitch = FALSE
  )
  seqlevels(gposDeepSignalMod) <- seqlevelsInUse(gposDeepSignalMod)

  if (!lKeepSequence) {
    gposDeepSignalMod$k_mer <- NULL
  }

  if (lSortGPos) {
    print("Sorting GPos Object...")
    gposDeepSignalMod <- sort(gposDeepSignalMod)
  }

  return(gposDeepSignalMod)
}
