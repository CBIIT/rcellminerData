library(rcellminer)
library(stringr)

#--------------------------------------------------------------------------------------------------
# HELPER FUNCTIONS
#--------------------------------------------------------------------------------------------------

validateDataTab <- function(dataTab, keyCol = "Gene name", featureDataColNums = 1:9){
  stopifnot(keyCol %in% colnames(dataTab))

  for (j in featureDataColNums){
    dataTab[, j] <- stringr::str_trim(dataTab[, j])
  }

  # Make sure expected numeric data is numeric (symbols to be read in as NAs are
  # properly handled).
  stopifnot(all(c(lapply(dataTab[, -featureDataColNums], is.numeric), recursive = TRUE)))

  # Set row names of data table to names in key column after checks.
  stopifnot(all(!is.na(dataTab[, keyCol])))
  stopifnot(all(dataTab[, keyCol] != ""))
  stopifnot(all(dataTab[, keyCol] != "1-Mar"))    # No Excel conversion to dates.
  stopifnot(all(!duplicated(dataTab[, keyCol])))
  rownames(dataTab) <- dataTab[, keyCol]

  return(dataTab)
}

#--------------------------------------------------------------------------------------------------
# LOAD DATA: MRNA EXPRESSION.
#--------------------------------------------------------------------------------------------------
# http://discovery.nci.nih.gov/cellminerint/loadDownload.do
# Select: [RNA: 5 Platform Gene Transcript, select: Average z score]
# Processing of data file (RNA__5_Platform_Gene_Transcript_Average_z_scores.xls):
# --- Save as text file
# --- Delete first 10 rows to get to table.
# --- Clean up column names (delete superscripts for footnotes)

#----[z score data]------------------------------------------------------------
filePath <- "inst/extdata/cellminer_2_0/RNA__5_Platform_Gene_Transcript_Average_z_scores.txt"
expTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="-")

featureDataCols <- 1:9
expTabOrig <- validateDataTab(expTabOrig, keyCol = "Gene name",
                              featureDataColNums = featureDataCols)

expData <- ExpressionSet(as.matrix(expTabOrig[, -featureDataCols]))
featureData(expData) <- new("AnnotatedDataFrame", data=expTabOrig[, featureDataCols])

###################################################################################################
# Set CellMiner NCI-60 Cell Line Names
cmNci60Names <- colnames(exprs(expData))
stopifnot(identical(cmNci60Names, stringr::str_trim(cmNci60Names)))

# Make CellMiner 1.6, 2.0 cell line name match tab.
cmNci60Names_1_6 <- colnames(rcellminer::getAllFeatureData(rcellminerData::molData)[["exp"]])
CellMinerNci60LineTab <- data.frame(CellMiner_1_6 = cmNci60Names_1_6,
                                    CellMiner_2_0 = cmNci60Names,
                                    stringsAsFactors = FALSE)
save(CellMinerNci60LineTab, file = "inst/extdata/CellMinerNci60LineTab.Rdata")
###################################################################################################

#----[average log2 intensity data]---------------------------------------------
# http://discovery.nci.nih.gov/cellminerint/loadDownload.do
# Select: [RNA: 5 Platform Gene Transcript, select: Averaged intensities]
# Processing of data file (RNA__5_Platform_Gene_Transcript_Averaged_intensities.xls):
# --- Save as text file
# --- Delete first 10 rows to get to table.
# --- Clean up column names (delete superscripts for footnotes)

filePath <- "inst/extdata/cellminer_2_0/RNA__5_Platform_Gene_Transcript_Averaged_intensities.txt"
xaiTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="-")

featureDataCols <- 1:9
xaiTabOrig <- validateDataTab(xaiTabOrig, keyCol = "Gene name",
                              featureDataColNums = featureDataCols)

xaiData <- ExpressionSet(as.matrix(xaiTabOrig[, -featureDataCols]))
featureData(xaiData) <- new("AnnotatedDataFrame", data=xaiTabOrig[, featureDataCols])

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(xaiData)), cmNci60Names))

#--------------------------------------------------------------------------------------------------
# LOAD DATA: GENE COPY.
#--------------------------------------------------------------------------------------------------
# http://discovery.nci.nih.gov/cellminerint/loadDownload.do
# Select: [DNA: Combined aCGH, select: Gene summary]
# Processing of data file (DNA__Combined_aCGH_Gene_summary.xls):
# --- Save as text file
# --- Delete first 10 rows to get to table.
# --- Clean up column names (delete superscripts for footnotes)

filePath <- "inst/extdata/cellminer_2_0/DNA__Combined_aCGH_Gene_summary.txt"
copTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="-")

featureDataCols <- 1:9
copTabOrig <- validateDataTab(copTabOrig, keyCol = "Probe id",
                              featureDataColNums = featureDataCols)

copData <- ExpressionSet(as.matrix(copTabOrig[, -featureDataCols]))
featureData(copData) <- new("AnnotatedDataFrame", data=copTabOrig[, featureDataCols])

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(copData)), cmNci60Names))

#--------------------------------------------------------------------------------------------------
# LOAD DATA: GENE METHYLATION
#--------------------------------------------------------------------------------------------------
# http://discovery.nci.nih.gov/cellminerint/loadDownload.do
# Select: [DNA: Illumina 450K methylation, select: Gene average]
# Processing of data file (DNA__Illumina_450K_methylation_Gene_average.xls):
# --- Save as text file
# --- Delete first 10 rows to get to table.
# --- Clean up column names (delete superscripts for footnotes)

filePath <- "inst/extdata/cellminer_2_0/DNA__Illumina_450K_methylation_Gene_average.txt"
metTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="-")

featureDataCols <- 1:9
metTabOrig <- validateDataTab(metTabOrig, keyCol = "Probe id",
                              featureDataColNums = featureDataCols)

metData <- ExpressionSet(as.matrix(metTabOrig[, -featureDataCols]))
featureData(metData) <- new("AnnotatedDataFrame", data=metTabOrig[, featureDataCols])

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(metData)), cmNci60Names))
#--------------------------------------------------------------------------------------------------
# LOAD DATA: EXOME/MUTATION.
#--------------------------------------------------------------------------------------------------

#----[gene level function altering mutations]-------------------------------------------------
# http://discovery.nci.nih.gov/cellminerint/loadDownload.do
# Select: [DNA: Exome Seq, select: protein function affecting]
# Processing of data file (DNA__Exome_Seq_Protein_function_affecting.xls):
# --- Save as text file
# --- Delete first 10 rows to get to table.
# --- Clean up column names (delete superscripts for footnotes)

filePath <- "inst/extdata/cellminer_2_0/DNA__Exome_Seq_Protein_function_affecting.txt"
mutTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="-")

featureDataCols <- 1:9
mutTabOrig <- validateDataTab(mutTabOrig, keyCol = "Gene name",
                              featureDataColNums = featureDataCols)

# NAs indicate that not enough reads were available to determine variant allele percent
# conversion (see CellMiner spreadsheet footnotes); treated as zeros for analyses.
for (cLine in colnames(mutTabOrig[, -featureDataCols])){
  naIndexSet <- which(is.na(mutTabOrig[, cLine]))
  mutTabOrig[naIndexSet, cLine] <- 0
}

mutData <- ExpressionSet(as.matrix(mutTabOrig[, -featureDataCols]))
featureData(mutData) <- new("AnnotatedDataFrame", data=mutTabOrig[, featureDataCols])

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(mutData)), cmNci60Names))

#----[variant level exome sequencing data]--------------------------------------------------

filePath <- "inst/extdata/cellminer_2_0/DNA__Exome_Seq_none.txt"
exoTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="-")

featureDataCols <- 1:18
exoTabOrig <- validateDataTab(exoTabOrig, keyCol = "Probe id",
                              featureDataColNums = featureDataCols)

exoData <- ExpressionSet(as.matrix(exoTabOrig[, -featureDataCols]))
featureData(exoData) <- new("AnnotatedDataFrame", data=exoTabOrig[, featureDataCols])

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(exoData)), cmNci60Names))

#--------------------------------------------------------------------------------------------------
# LOAD DATA: PROTEIN EXPRESSION (RPLA)
#--------------------------------------------------------------------------------------------------
# http://discovery.nci.nih.gov/cellminerint/loadDownload.do
# Select: [Protein: Lysate Array, select: log2].

filePath <- "inst/extdata/cellminer_2_0/"
proTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="-")

#--------------------------------------------------------------------------------------------------
# LOAD DATA: PROTEIN EXPRESSION (SWATH-MS)
#--------------------------------------------------------------------------------------------------

filePath <- "inst/extdata/cellminer_2_0/"
swaTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="-")

#--------------------------------------------------------------------------------------------------
# LOAD DATA: MICRORNA EXPRESSION.
#--------------------------------------------------------------------------------------------------

filePath <- "inst/extdata/cellminer_2_0/"
mirTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="-")

#--------------------------------------------------------------------------------------------------
# LOAD DATA: CELL LINE  METADATA.
#--------------------------------------------------------------------------------------------------

filePath <- "inst/extdata/cellminer_2_0/"
mdaTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="-")

#--------------------------------------------------------------------------------------------------
# LOAD DATA: DRUG ACTIVITY.
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
# Make NCI-60 sample info (shared by molData and drugData objects to be constructed).
#--------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
# Make NCI-60 MolData object.
#--------------------------------------------------------------------------------------------------

nci60ESetList <- list()

nci60ESetList[["exp"]] <- expData
nci60ESetList[["xai"]] <- xaiData

nci60ESetList[["cop"]] <- copData
nci60ESetList[["met"]] <- metData

nci60ESetList[["mut"]] <- mutData
nci60ESetList[["exo"]] <- exoData

nci60ESetList[["pro"]] <- proData
nci60ESetList[["swa"]] <- swaData

nci60ESetList[["mir"]] <- mirData
nci60ESetList[["mda"]] <- mdaData

molData <- new("MolData", eSetList = nci60ESetList, sampleData = nci60Miame)

save(molData, file = "data/molData.RData")

#--------------------------------------------------------------------------------------------------
# Make NCI-60 DrugData object.
#--------------------------------------------------------------------------------------------------

drugData <- new("DrugData", act = actData, repeatAct = repeatActData, sampleData = nci60Miame)

save(drugData, file = "data/drugData.RData")

#--------------------------------------------------------------------------------------------------
