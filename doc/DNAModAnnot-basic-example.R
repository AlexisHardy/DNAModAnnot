## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install(c('Biostrings', 'BSgenome', 'Gviz', 'Logolas'))
#  
#  # Installation from tar.gz file
#  setwd("path/to/package/file/")
#  
#  install.packages("DNAModAnnot_0.0.0.9014.tar.gz", repos = NULL, type = 'source')

## ----setup--------------------------------------------------------------------
library(DNAModAnnot)

## -----------------------------------------------------------------------------
# Loading genome assembly (fasta file or DNAStringSet)
ptetraurelia_genome_fa <- system.file(
  package = "DNAModAnnot", "extdata",
  "ptetraurelia_mac_51_sca171819.fa"
)
ptetraurelia_genome <- Biostrings::readDNAStringSet(ptetraurelia_genome_fa)
ptetraurelia_genome_range <- GetGenomeGRanges(ptetraurelia_genome)

## -----------------------------------------------------------------------------
# Statistics summary of genome
report_assembly <- GetAssemblyReport(
  dnastringsetGenome = ptetraurelia_genome,
  cOrgAssemblyName = "ptetraurelia_mac_51"
)

## ----plot1, fig.width = 7, fig.asp = .62, fig.align = "center"----------------
# Cumulative length of genome assembly
contig_cumulative_length <- GetContigCumulLength(ptetraurelia_genome)
DrawContigCumulLength(
  nContigCumsumLength = contig_cumulative_length$cumsum_Mbp_length,
  cOrgAssemblyName = "ptetraurelia_mac_51",
  lGridInBackground = FALSE
)

## ----warning=FALSE, results="hide"--------------------------------------------
# loading SMRT-seq data (PacBio gff and csv files)
PacBioGFF_path <- system.file(package = "DNAModAnnot", "extdata", "ptetraurelia.modifications.sca171819.gff")
PacBioGFF_granges <- ImportPacBioGFF(
  cPacBioGFFPath = PacBioGFF_path,
  cNameModToExtract = "m6A",
  cModNameInOutput = "6mA",
  cContigToBeAnalyzed = names(ptetraurelia_genome)
)
PacBioGFF_path <- system.file(package = "DNAModAnnot", "extdata", "ptetraurelia.bases.sca171819.csv")
PacBioCSV_gpos <- ImportPacBioCSV(
  cPacBioCSVPath = PacBioGFF_path,
  cSelectColumnsToExtract = c(
    "refName", "tpl", "strand",
    "base", "score",
    "ipdRatio", "coverage"
  ),
  lKeepExtraColumnsInGPos = TRUE,
  lSortGPos = TRUE,
  cContigToBeAnalyzed = names(ptetraurelia_genome)
)

## ----barplotfig1, fig.width = 7, fig.asp = .62, fig.align = "center"----------
# Percentage of sequencing by scaffold and by strand
contig_percentage_sequencing <- GetSeqPctByContig(PacBioCSV_gpos,
  grangesGenome = ptetraurelia_genome_range
)
DrawBarplotBothStrands(
  nParamByContigForward = contig_percentage_sequencing$f_strand$seqPct,
  nParamByContigReverse = contig_percentage_sequencing$r_strand$seqPct,
  cContigNames = contig_percentage_sequencing$f_strand$refName,
  cGraphName = "Percentage of sequencing per contig"
)

## ----barplotfig2, fig.width = 7, fig.asp = .62, fig.align = "center"----------
# Mean coverage by scaffold and by strand
contig_mean_coverage <- GetMeanParamByContig(
  grangesData = PacBioCSV_gpos,
  dnastringsetGenome = ptetraurelia_genome,
  cParamName = "coverage"
)
DrawBarplotBothStrands(
  nParamByContigForward = contig_mean_coverage$f_strand$mean_coverage,
  nParamByContigReverse = contig_mean_coverage$r_strand$mean_coverage,
  cContigNames = contig_mean_coverage$f_strand$refName,
  cGraphName = "Mean Coverage per contig"
)

## ----histfig, fig.width = 10, fig.asp = .62, fig.align = "center"-------------
# coverage distribution of all bases sequenced
DrawDistriHistBox(PacBioCSV_gpos$coverage,
  cGraphName = "Coverage distribution of all bases sequenced",
  cParamName = "Coverage",
  lTrimOutliers = FALSE
)

## -----------------------------------------------------------------------------
# First filter: contigs
PacBio_filtered_data <- FiltPacBio(
  gposPacBioCSV = PacBioCSV_gpos,
  grangesPacBioGFF = PacBioGFF_granges,
  cContigToBeRemoved = NULL,
  dnastringsetGenome = ptetraurelia_genome,
  nContigMinSize = 1000,
  listPctSeqByContig = contig_percentage_sequencing,
  nContigMinPctOfSeq = 1,
  listMeanCovByContig = contig_mean_coverage,
  nContigMinCoverage = 20
)
PacBioCSV_gpos_filt1 <- PacBio_filtered_data$csv
PacBioGFF_granges_filt1 <- PacBio_filtered_data$gff

## ----warning=FALSE, results="hide"--------------------------------------------
# report basic Mod for remaining scaffolds
report_modifications <- GetModReportPacBio(
  grangesGenome = ptetraurelia_genome_range,
  grangesPacBioGFF = PacBioGFF_granges_filt1,
  gposPacBioCSV = PacBioCSV_gpos_filt1,
  cOrgAssemblyName = "ptetraurelia_mac_51",
  dnastringsetGenome = ptetraurelia_genome,
  cBaseLetterForMod = "A",
  cModNameInOutput = "6mA"
)

## ----barplotfig3, fig.width = 7, fig.asp = .62, fig.align = "center"----------
# Modif ratio by scaffold and by strand
contig_modification_ratio <- GetModRatioByContig(PacBioGFF_granges_filt1,
                                                 PacBioCSV_gpos_filt1[PacBioCSV_gpos_filt1$base == "A"],
                                                 dnastringsetGenome = ptetraurelia_genome,
                                                 cBaseLetterForMod = "A"
)
DrawBarplotBothStrands(
  nParamByContigForward = contig_modification_ratio$f_strand$Mod_ratio,
  nParamByContigReverse = contig_modification_ratio$r_strand$Mod_ratio,
  cContigNames = contig_modification_ratio$f_strand$refName,
  cGraphName = "Modif/Base ratio per contig (Sequenced sites only)"
)

## ----logofig, fig.width = 5, fig.asp = .62, fig.align = "center", warning=FALSE, results="hide"----
# Sequence logo associated to Modif detected
PacBioGFF_granges_with_sequence <- GetGRangesWindowSeqandParam(PacBioGFF_granges_filt1,
  ptetraurelia_genome_range,
  dnastringsetGenome = ptetraurelia_genome,
  nUpstreamBpToAdd = 5,
  nDownstreamBpToAdd = 5
)
DrawModLogo(
  dnastringsetSeqAroundMod = as(PacBioGFF_granges_with_sequence$sequence, "DNAStringSet"),
  cLogoType = "EDLogo",
  nGenomicBgACGT = c(0.35, 0.15, 0.15, 0.35)
)

## -----------------------------------------------------------------------------
# Second filter: fraction
PacBioGFF_granges_filt2 <- FiltPacBio(
  grangesPacBioGFF = PacBioGFF_granges_filt1,
  cParamNameForFilter = "frac", 
  lFiltParam = TRUE,
  nFiltParamLoBoundaries = 0.05, nFiltParamUpBoundaries = 1, 
  cFiltParamBoundariesToInclude = "upperOnly"
)$gff 

## ----warning=FALSE, results="hide"--------------------------------------------
# Extract Mod Data by motif over-represented (at least % of motifs)
motif_pct_and_PacBioGFF_grangeslist <- ExtractListModPosByModMotif(
  grangesModPos = PacBioGFF_granges_filt2,
  grangesGenome = ptetraurelia_genome_range,
  dnastringsetGenome = ptetraurelia_genome,
  nUpstreamBpToAdd = 0, nDownstreamBpToAdd = 1,
  nModMotifMinProp = 0.05,
  cBaseLetterForMod = "A",
  cModNameInOutput = "6mA"
)

## ----fdrplotfig, fig.width = 10, fig.asp = .62, fig.align = "center", warning=FALSE, results="hide"----
# Motif FDR estimation for each motif over-represented (at least % of motifs)
BaseCSV_granges_filt1 <- as(PacBioCSV_gpos_filt1[PacBioCSV_gpos_filt1$base == "A"], "GRanges")
BaseCSV_granges_with_sequence <- GetGRangesWindowSeqandParam(
  grangesData = BaseCSV_granges_filt1,
  grangesGenome = ptetraurelia_genome_range,
  dnastringsetGenome = ptetraurelia_genome,
  nUpstreamBpToAdd = 0,
  nDownstreamBpToAdd = 1
)
score_fdr_by_motif_list <- GetFdrEstListByThresh(
  grangesDataWithSeq = BaseCSV_granges_with_sequence,
  grangesDataWithSeqControl = NULL,
  cNameParamToTest = "score",
  nRoundDigits = 1,
  cModMotifsAsForeground = motif_pct_and_PacBioGFF_grangeslist$motifs_to_analyse
)
score_fdr_by_motif_limit <- GetFdrBasedThreshLimit(score_fdr_by_motif_list,
  nFdrPropForFilt = 0.05,
  lUseBestThrIfNoFdrThr = TRUE
)
DrawFdrEstList(
  listFdrEstByThr = score_fdr_by_motif_list,
  cNameParamToTest = "score",
  nFdrPropForFilt = 0.05
)

## -----------------------------------------------------------------------------
PacBioGFF_grangeslist_filt <- FiltPacBio(
  grangesPacBioGFF = motif_pct_and_PacBioGFF_grangeslist$GRangesbyMotif,
  listFdrEstByThrIpdRatio = NULL,
  listFdrEstByThrScore = score_fdr_by_motif_limit
)$gff

## -----------------------------------------------------------------------------
# Select one motif for the analysis by motif
motifs_base <- motif_pct_and_PacBioGFF_grangeslist$motifs_to_analyse[motif_pct_and_PacBioGFF_grangeslist$motifs_to_analyse == "AT"]
motifs_modifications <- motif_pct_and_PacBioGFF_grangeslist$mod_motif[motifs_base == motif_pct_and_PacBioGFF_grangeslist$motifs_to_analyse]
PacBioGFF_granges_filtAT <- PacBioGFF_grangeslist_filt[[motifs_base]]
BaseCSV_granges_filtAT <- BaseCSV_granges_with_sequence[BaseCSV_granges_with_sequence$sequence == motifs_base, ]

## -----------------------------------------------------------------------------
# Loading GFF annotation file (Minimum annotation required: gene (seqnames, start, end, strand))
annotations_path <- system.file(package = "DNAModAnnot", "extdata", "ptetraurelia_mac_51_annotation_v2.0_sca171819.gff3")
annotations_range <- rtracklayer::readGFFAsGRanges(annotations_path)
annotations_range <- PredictMissingAnnotation(
  grangesAnnotations = annotations_range,
  grangesGenome = ptetraurelia_genome_range,
  cFeaturesColName = "type",
  cGeneCategories = c("gene"),
  lAddIntronRangesUsingExon = TRUE
)

## ----annbyffig, fig.width = 8, fig.asp = .62, fig.align = "center"------------
# Mod annotation by feature
annotations_range_ModBase_counts <- GetModBaseCountsByFeature(
  grangesAnnotations = annotations_range,
  grangesModPos = PacBioGFF_granges_filtAT,
  gposModTargetBasePos = BaseCSV_granges_filtAT,
  lIgnoreStrand = FALSE
)
DrawModBasePropByFeature(
  grangesAnnotationsWithCounts = annotations_range_ModBase_counts,
  cFeaturesToCompare = c("gene", "intergenic"),
  lUseCountsPerkbp = TRUE,
  cBaseMotif = motifs_base,
  cModMotif = motifs_modifications
)

## -----------------------------------------------------------------------------
# Loading additional Data: RNA-seq dataset
expression_file_path <- system.file(package = "DNAModAnnot", "extdata", "ptetraurelia.gene_expression.sca171819.tsv")
expression_dataframe <- read.table(
  file = expression_file_path,
  header = TRUE, sep = "\t", dec = ",",
  quote = "\"", fill = FALSE
)
expression_dataframe <- expression_dataframe[, c("ID", "T30")]

## -----------------------------------------------------------------------------
# merge RNA-seq data with Annotation GRanges object
genes_range_ModBase_counts_param <- annotations_range_ModBase_counts[annotations_range_ModBase_counts$type == "gene"]
GenomicRanges::mcols(genes_range_ModBase_counts_param) <- merge(
  x = GenomicRanges::mcols(genes_range_ModBase_counts_param),
  by.x = "Name",
  y = expression_dataframe,
  by.y = "ID"
)

## ----annQfig, fig.width = 9, fig.asp = .62, fig.align = "center"--------------
# Comparison of quantitative parameter with Mod annotation
DrawParamPerModBaseCategories(
  grangesAnnotationsWithCounts = genes_range_ModBase_counts_param,
  cParamColname = "T30",
  cParamFullName = "Gene expression at T30",
  cParamYLabel = "Normalized RNA-seq read counts (T30)",
  cSelectFeature = "gene",
  lUseCountsPerkbp = TRUE,
  cBaseMotif = motifs_base,
  cModMotif = motifs_modifications,
  lBoxPropToCount = FALSE
)

## ----annwithinffig, fig.width = 8, fig.asp = .62, fig.align = "center"--------
# Mod annotation within feature
gene_annotation_range <- annotations_range[annotations_range$type == "gene", ]
gene_annotation_range <- GetModBaseCountsWithinFeature(
  grangesAnnotations = gene_annotation_range,
  grangesModPos = PacBioGFF_granges_filtAT,
  gposModTargetBasePos = BaseCSV_granges_filtAT,
  lIgnoreStrand = FALSE,
  nWindowsNb = 20
)
DrawModBaseCountsWithinFeature(
  grangesAnnotationsWithCountsByWindow = gene_annotation_range,
  cFeatureName = "gene",
  cBaseMotif = motifs_base,
  cModMotif = motifs_modifications
)

## ----anndistfig, fig.width = 10, fig.asp = 1, warning=FALSE, results="hide"----
# ModBase distance from feature/feature limit
Mod_distance_feature_countslist <- GetDistFromFeaturePos(
  grangesAnnotations = annotations_range,
  cSelectFeature = "gene",
  grangesData = PacBioGFF_granges_filtAT,
  lGetGRangesInsteadOfListCounts = FALSE,
  lGetPropInsteadOfCounts = TRUE,
  cWhichStrandVsFeaturePos = "both", nWindowSizeAroundFeaturePos = 600,
  lAddCorrectedDistFrom5pTo3p = TRUE,
  cFeaturePosNames = c("TSS", "TTS")
)
Base_distance_feature_countslist <- GetDistFromFeaturePos(
  grangesAnnotations = annotations_range,
  cSelectFeature = "gene",
  grangesData = BaseCSV_granges_filtAT,
  lGetGRangesInsteadOfListCounts = FALSE,
  lGetPropInsteadOfCounts = TRUE,
  cWhichStrandVsFeaturePos = "both", nWindowSizeAroundFeaturePos = 600,
  lAddCorrectedDistFrom5pTo3p = TRUE,
  cFeaturePosNames = c("TSS", "TTS")
)
DrawModBasePropDistFromFeature(
  listModCountsDistDataframe = Mod_distance_feature_countslist,
  listBaseCountsDistDataframe = Base_distance_feature_countslist,
  cFeaturePosNames = c("TSS", "TTS"),
  cBaseMotif = motifs_base,
  cModMotif = motifs_modifications
)

# ModBase and Reads center (from Bam file) distance from feature/feature limit
bamfile_path <- system.file(package = "DNAModAnnot", "extdata", "PTET_MonoNuc_3-2new.pe.sca171819.sorted.bam")
bamfile_object <- Rsamtools::BamFile(file = bamfile_path)
bamfile_ranges <- as(GenomicAlignments::readGAlignments(bamfile_object), "GRanges")
bamfile_ranges <- GetGposCenterFromGRanges(grangesData = bamfile_ranges)
bamfile_distance_feature_countslist <- GetDistFromFeaturePos(
  grangesAnnotations = annotations_range,
  cSelectFeature = "gene",
  grangesData = bamfile_ranges,
  lGetGRangesInsteadOfListCounts = FALSE,
  lGetPropInsteadOfCounts = FALSE,
  cWhichStrandVsFeaturePos = "both", nWindowSizeAroundFeaturePos = 600,
  lAddCorrectedDistFrom5pTo3p = TRUE,
  cFeaturePosNames = c("TSS", "TTS")
)
AddToModBasePropDistFromFeaturePlot(
  dPosCountsDistFeatureStart = bamfile_distance_feature_countslist[[1]],
  dPosCountsDistFeatureEnd = bamfile_distance_feature_countslist[[2]],
  cSubtitleContent = "Along with nucleosome center distance (MonoNuc_3-2newreplicate)",
  cParamYLabel = "Nucleosome center count (MonoNuc_3-2newreplicate)",
  cParamColor = "cyan3",
  lAddAxisOnLeftSide = TRUE
)

## -----------------------------------------------------------------------------
ipdRatio6mABigwig <- system.file(package = "DNAModAnnot", "extdata", "ipdRatio6mATsites.bw")
ptet51GenesBam <- system.file(package = "DNAModAnnot", "extdata", "ptet51Genes.bam")

## ----eval=FALSE---------------------------------------------------------------
#  # #Gviz file making for vignette (not run)
#  ExportFilesForGViz(cFileNames = c(ipdRatio6mABigwig,
#                                    ptet51GenesBam),
#                     listObjects = list(PacBioGFF_granges_filtAT,
#                                               annotations_range[annotations_range$type == "gene"]),
#                     cFileFormats = c("bw", "bam"),
#                     lBigwigParametersByStrand = c(TRUE, NA),
#                     cBigwigParameters = c("ipdRatio",NA),
#                     cBamXaParameters = c(NA, "Name"),
#                     dnastringsetGenome = ptetraurelia_genome)

## -----------------------------------------------------------------------------
cContigToViz <- unique(GenomicRanges::seqnames(ptetraurelia_genome_range))

#Generating Ideogram TRACK--------
options(ucscChromosomeNames=FALSE)
trackIdeogram <- AdaptedIdeogramTrackWithoutBandsData(grangesGenome = ptetraurelia_genome_range,
                                                      cContigToViz = cContigToViz,
                                                      cOrgAssemblyName = "ptetraurelia_mac_51")

## -----------------------------------------------------------------------------
#Generating classic Gviz tracks

#GenomeAxis TRACK------
trackGenomeAxis <- Gviz::GenomeAxisTrack(cex=1)

#Sequence TRACK--------
trackSequence <- Gviz::SequenceTrack(ptetraurelia_genome_fa, 
                                     chromosome=cContigToViz, 
                                     add53=TRUE,
                                     complement=FALSE, cex=0.8, stream = TRUE)

#DATATRACK--------
trackData6mATipdRatio <- Gviz::DataTrack(ipdRatio6mABigwig, 
                                         stream = TRUE, name = "6mAT\nipdRatio", type="histogram", 
                                         col.histogram=c("red"), fill = "red",
                                         background.title = "darkred", col = NULL
)

trackDataNuclCoverage <- Gviz::DataTrack(bamfile_path,
                                         stream = TRUE, name = "Nucleosome\ncoverage\n(MNase-seq)", 
                                         type="histogram", 
                                         col.histogram=c("blue"), fill="blue",
                                         background.title = "darkblue", col = NULL
)

## -----------------------------------------------------------------------------
#Generating Annotation TRACK with streaming--------
trackAnnotation <- Gviz::AnnotationTrack(ptet51GenesBam, 
                                   name = "Gene", stacking = "squish", 
                                   stream = TRUE, importFunction = ImportBamExtendedAnnotationTrack,
                                   group="tag", groupAnnotation = "tag", 
                                   just.group="below",
                                   fontcolor.group = "black", fontsize.group=18,
                                   fill="lightblue", shape="fixedArrow", 
                                   arrowHeadWidth = 50, lwd = 3,
                                   background.title = "darkblue"
)

## -----------------------------------------------------------------------------
trackAnnotation@mapping <- list(id="tag", group="tag") 

## ----annlocalfig, fig.width=7, fig.height = 8, warning=FALSE, results="hide"----
#Plotting Tracks--------
Gviz::plotTracks(trackList = list(trackIdeogram, trackGenomeAxis, 
                                  trackData6mATipdRatio, trackSequence, 
                                  trackDataNuclCoverage, trackAnnotation),
                 chromosome = "scaffold51_17", 
                 from = 361000, to = 365000)

## -----------------------------------------------------------------------------
sessionInfo()

