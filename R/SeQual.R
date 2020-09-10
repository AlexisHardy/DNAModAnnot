### =========================================================================
### Functions from SeQual category
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
### GetAssemblyReport
### .CalculateNLx
### GetContigCumulLength
### DrawContigCumulLength
### GetSeqPctByContig
###
### -------------------------------------------------------------------------

#' GetAssemblyReport Function (SeQual)
#' importFrom base data.frame max sum return
#' importFrom Biostrings alphabetFrequency letterFrequency
#'
#' Return a report with global characteristics of genome assembly.
#' @param dnastringsetGenome A DNAStringSet object containing the sequence for each contig.
#' @param cOrgAssemblyName The name of the genome assembly provided.
#' @keywords GetAssemblyReport
#' @import Biostrings
#' @export
#' @examples
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#'
#' myReport <- GetAssemblyReport(dnastringsetGenome = myGenome,
#'                                   cOrgAssemblyName = "ptetraurelia_mac_51")
#'
#' myReport
GetAssemblyReport <- function(dnastringsetGenome,
                                  cOrgAssemblyName) {
  freq_v <- alphabetFrequency(dnastringsetGenome, as.prob = TRUE,
                              baseOnly=TRUE, collapse=TRUE)[1:4]
  dTable <- data.frame(
    row.names = cOrgAssemblyName,
    nbcontigs = length(dnastringsetGenome),
    nbcontigs_ab0 = length(which(dnastringsetGenome@ranges@width >= 0)),
    nbcontigs_ab1000 = length(which(dnastringsetGenome@ranges@width >= 1000)),
    nbcontigs_ab5000 = length(which(dnastringsetGenome@ranges@width >= 5000)),
    nbcontigs_ab10000 = length(which(dnastringsetGenome@ranges@width >= 10000)),
    nbcontigs_ab25000 = length(which(dnastringsetGenome@ranges@width >= 25000)),
    nbcontigs_ab50000 = length(which(dnastringsetGenome@ranges@width >= 50000)),
    largest_contig_size = max(dnastringsetGenome@ranges@width),
    total_length = sum(dnastringsetGenome@ranges@width),
    total_length_ab0 = sum(dnastringsetGenome[dnastringsetGenome@ranges@width >= 0, ]@ranges@width),
    total_length_ab1000 = sum(dnastringsetGenome[dnastringsetGenome@ranges@width >= 1000, ]@ranges@width),
    total_length_ab5000 = sum(dnastringsetGenome[dnastringsetGenome@ranges@width >= 5000, ]@ranges@width),
    total_length_ab10000 = sum(dnastringsetGenome[dnastringsetGenome@ranges@width >= 10000, ]@ranges@width),
    total_length_ab25000 = sum(dnastringsetGenome[dnastringsetGenome@ranges@width >= 25000, ]@ranges@width),
    total_length_ab50000 = sum(dnastringsetGenome[dnastringsetGenome@ranges@width >= 50000, ]@ranges@width),
    n50 = .CalculateNLx(dnastringsetGenome, 50)[1],
    n75 = .CalculateNLx(dnastringsetGenome, 75)[1],
    l50 = .CalculateNLx(dnastringsetGenome, 50)[2],
    l75 = .CalculateNLx(dnastringsetGenome, 75)[2],
    gc_pct = sum( (freq_v / sum(freq_v))[c("G", "C")] )*100,
    nb_N = letterFrequency(dnastringsetGenome, letters="N", collapse=TRUE)
  )
  return(t(dTable))
}

#' CalculateNLx Function (SeQual)
#'
#' Calculate Nx and Lx.
#' @param dnastringsetGenome A DNAStringSet object containing the sequence for each contig.
#' @param nNLx A number for the Nx and Lx to calculate. Defaults to 50 (so N50 and L50 calculation).
#' @keywords internal
.CalculateNLx <- function(dnastringsetGenome, nNLx = 50) {
  w_sort_v <- sort(dnastringsetGenome@ranges@width, decreasing = TRUE)
  vector_o <- c(NA, NA)
  names(vector_o) <- paste(c("N", "L"), nNLx, sep="")
  vector_o[1] <- head( w_sort_v[ cumsum(w_sort_v) >= (sum(dnastringsetGenome@ranges@width) * (nNLx/100)) ], 1)
  vector_o[2] <- head( which(    cumsum(w_sort_v) >= (sum(dnastringsetGenome@ranges@width) * (nNLx/100)) ), 1)
  return(vector_o)
}

#' GetContigCumulLength Function (SeQual)
#'
#' Return a dataframe with the length (and cumulative length) of contigs (ordered from largest to smallest contig) provided in genome assembly.
#' @param dnastringsetGenome A DNAStringSet object containing the sequence for each contig.
#' @return A dataframe:
#' \itemize{
#'   \item contig_names: The names of each contig.
#'   \item Mbp_length: The size of each contig (in Megabase pairs (Mbp)).
#'   \item cumsum_Mbp_length: The cumulative size (from largest to smallest contig) for each contig (in Megabase pairs (Mbp)).
#' }
#' @keywords GetContigCumulLength
#' @export
#' @examples
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#'
#' myContig_cumul_len_t <- GetContigCumulLength(dnastringsetGenome = myGenome)
#'
#' myContig_cumul_len_t
GetContigCumulLength <- function(dnastringsetGenome) {
  sca_s <- dnastringsetGenome[order(dnastringsetGenome@ranges@width, decreasing = TRUE)]
  mbp_length <- sca_s@ranges@width / 1000000
  data_table <- data.frame(
    contig_names = sca_s@ranges@NAMES,
    Mbp_length = mbp_length,
    cumsum_Mbp_length = cumsum(mbp_length)
  )
  return(data_table)
}

#' DrawContigCumulLength Function (SeQual)
#'
#' Return a line-plot describing the cumulative length of contigs (ordered from largest to smallest contig).
#' @param nContigCumsumLength A numeric vector containing the cumulative length of contigs (ordered from largest to smallest contig).
#' @param cOrgAssemblyName The name of the genome assembly corresponding to these contigs.
#' @param lGridInBackground If TRUE, add a grid in the background. Defaults to FALSE.
#' @param cLineColor The color of the line. Defaults to "red".
#' @param nLineWidth The width of the line. Defaults to 2.
#' @param nLineType The type of the line. Defaults to 1.
#' @keywords DrawContigCumulLength
#' @export
#' @examples
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#'
#' myContig_cumul_len_t <- GetContigCumulLength(dnastringsetGenome = myGenome)
#'
#' DrawContigCumulLength(nContigCumsumLength = myContig_cumul_len_t$cumsum_Mbp_length,
#'                              cOrgAssemblyName="ptetraurelia_mac_51")
DrawContigCumulLength <- function(nContigCumsumLength,
                                         cOrgAssemblyName,
                                         lGridInBackground = FALSE,
                                         cLineColor = "red",
                                         nLineWidth = 2,
                                         nLineType = 1) {
  opar <- par()
  par(mar = c(5.5, 5, 2.5, 2))
  plot(nContigCumsumLength,
       ylim = c(0, tail(nContigCumsumLength, 1)),

       type = "s", lty = nLineType, lwd = nLineWidth, col = cLineColor,

       xlab = "Contigs (ordered from largest to smallest)",
       ylab = "Cumulative length (Mbp)",
       main = paste("Cumulative length of ", cOrgAssemblyName, " genome assembly", sep=""),

       las = 1, lab = c(10, 10, 7),
       xaxs = "i", yaxs = "i", xaxt="n"
  )
  axis(1,at=1:length(nContigCumsumLength),labels=1:length(nContigCumsumLength))
  if(lGridInBackground) {
    grid()
  }
  par(mar=opar$mar)
}

#' GetSeqPctByContig Function (SeQual)
#'
#' Return a list with the percentage of sequencing by strand for all scaffolds of genome assembly provided.
#' This function is not adapted for data from DeepSignal.
#' @param gposPacBioCSV An UnStitched GPos object containing PacBio CSV data to be analysed.
#' @param grangesGenome A GRanges object containing the width of each contig.
#' @return A list composed of 3 dataframes: 1 dataframe by strand and 1 dataframe with both strands. In each dataframe:
#' \itemize{
#'   \item refName: The names of each contig.
#'   \item strand: The strand of each contig.
#'   \item width: The width of each contig.
#'   \item nb_sequenced: The number of bases sequenced by strand for each contig.
#'   \item seqPct: The percentage of bases sequenced for each strand for each contig (percentage of sequencing).
#' }
#' @keywords GetSeqPctByContig
#' @importFrom data.table data.table
#' @importFrom GenomicRanges strand
#' @export
#' @examples
#' myGenome <- Biostrings::readDNAStringSet(system.file(package="DNAModAnnot", "extdata",
#'                                          "ptetraurelia_mac_51_sca171819.fa"))
#' myGrangesGenome <- GetGenomeGRanges(myGenome)
#'
#' #Preparing a gposPacBioCSV dataset
#' myGposPacBioCSV <-
#'    ImportPacBioCSV(cPacBioCSVPath = system.file(package="DNAModAnnot", "extdata",
#'                                  "ptetraurelia.bases.sca171819.csv"),
#'                    cSelectColumnsToExtract = c("refName", "tpl", "strand", "base",
#'                                                "score", "ipdRatio", "coverage"),
#'                    lKeepExtraColumnsInGPos = TRUE, lSortGPos = TRUE,
#'                    cContigToBeAnalyzed = names(myGenome))
#'
#' myPct_seq_csv <- GetSeqPctByContig(myGposPacBioCSV, grangesGenome = myGrangesGenome)
#' myPct_seq_csv
GetSeqPctByContig <- function(gposPacBioCSV, grangesGenome) {
  table_count <- table(  data.table(seqnames=as.factor(seqnames(gposPacBioCSV)),
                                    strand=as.factor(strand(gposPacBioCSV)))  )
  table_count <- table_count[,c("+","-")]
  table_count <- table_count[rownames(table_count) %in% seqnames(grangesGenome),]
  table_count <- data.frame(refName=rep(rownames(table_count),2),
                            strand=c( rep("+",nrow(table_count)), rep("-",nrow(table_count))),
                            nb_sequenced=as.numeric(table_count))

  table_contig <- data.frame(refName=seqnames(grangesGenome),
                             strand=strand(grangesGenome),
                             width=width(grangesGenome))

  table_contig <- merge(table_contig, table_count, by=c("refName", "strand"), all.x = TRUE)

  table_contig[is.na(table_contig$nb_sequenced),"nb_sequenced"] <- 0
  table_contig$seqPct <- 100 * table_contig$nb_sequenced / table_contig$width

  table_contig <- table_contig[order(table_contig$width, decreasing = TRUE),]
  f_strand <- subset(table_contig, strand == "+")
  r_strand <- subset(table_contig, strand == "-")

  return(list(both_strand=table_contig,
              f_strand=f_strand,
              r_strand=r_strand))
}
