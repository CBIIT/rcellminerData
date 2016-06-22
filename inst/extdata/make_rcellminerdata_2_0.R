library(rcellminer)
library(stringr)

#--------------------------------------------------------------------------------------------------
# HELPER FUNCTIONS
#--------------------------------------------------------------------------------------------------

validateDataTab <- function(dataTab, keyCol = "Gene name", featureDataColNums = 1:9){
  stopifnot(keyCol %in% colnames(dataTab))

  # Remove any whitespace in feature data columns.
  for (j in featureDataColNums){
    if (is.character(dataTab[, j])){
      dataTab[, j] <- stringr::str_trim(dataTab[, j])
    }
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
# cmNci60Names_1_6 <- colnames(rcellminer::getAllFeatureData(rcellminerData::molData)[["exp"]])
# CellMinerNci60LineTab <- data.frame(CellMiner_1_6 = cmNci60Names_1_6,
#                                     CellMiner_2_0 = cmNci60Names,
#                                     stringsAsFactors = FALSE)
# save(CellMinerNci60LineTab, file = "inst/extdata/CellMinerNci60LineTab.Rdata")
load("inst/extdata/CellMinerNci60LineTab.Rdata")
stopifnot(identical(cmNci60Names, CellMinerNci60LineTab$CellMiner_2_0))
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

filePath <- "inst/extdata/cellminer_2_0/Protein__Lysate_Array_log2.txt"
proTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="-")

featureDataCols <- 1:9
proTabOrig <- validateDataTab(proTabOrig, keyCol = "Probe id",
                              featureDataColNums = featureDataCols)

proData <- ExpressionSet(as.matrix(proTabOrig[, -featureDataCols]))
featureData(proData) <- new("AnnotatedDataFrame", data=proTabOrig[, featureDataCols])

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(proData)), cmNci60Names))

#--------------------------------------------------------------------------------------------------
# LOAD DATA: PROTEIN EXPRESSION (SWATH-MS)
#--------------------------------------------------------------------------------------------------

# filePath <- "inst/extdata/cellminer_2_0/"
# swaTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
#                          check.names = FALSE, comment.char="", quote="", na.strings="-")

swaExpMat <- rcellminer::getAllFeatureData(nci60imsb::molData)[["swathms_avg"]]
stopifnot(identical(colnames(swaExpMat), CellMinerNci60LineTab$CellMiner_1_6))
colnames(swaExpMat) <- CellMinerNci60LineTab$CellMiner_2_0

swaAnnot <- rcellminer::getFeatureAnnot(nci60imsb::molData)[["swathms_avg"]]
stopifnot(identical(rownames(swaAnnot), rownames(swaExpMat)))

swaData <- ExpressionSet(swaExpMat)
featureData(swaData) <- new("AnnotatedDataFrame", data=swaAnnot)

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(swaData)), cmNci60Names))

#--------------------------------------------------------------------------------------------------
# LOAD DATA: MICRORNA EXPRESSION.
#--------------------------------------------------------------------------------------------------
# http://discovery.nci.nih.gov/cellminerint/loadDownload.do
# Select: [RNA: Agilent Human microRNA (V2)].

filePath <- "inst/extdata/cellminer_2_0/RNA__Agilent_Human_microRNA_(V2)_GeneSpringGX.txt"
mirTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="-")

featureDataCols <- 1:11
mirTabOrig <- validateDataTab(mirTabOrig, keyCol = "Probe id",
                              featureDataColNums = featureDataCols)

mirData <- ExpressionSet(as.matrix(mirTabOrig[, -featureDataCols]))
featureData(mirData) <- new("AnnotatedDataFrame", data=mirTabOrig[, featureDataCols])

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(mirData)), cmNci60Names))

#--------------------------------------------------------------------------------------------------
# LOAD DATA: CELL LINE  METADATA.
#--------------------------------------------------------------------------------------------------

filePath <- "inst/extdata/cellminer_2_0/CELLMINER_CELL_LINE_METADATA.txt"
mdaTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings=c("", "NA", "?"))

if (identical(mdaTabOrig$`Cell Line Name`, CellMinerNci60LineTab$CellMiner_1_6)){
  mdaTabOrig$`Cell Line Name` <- CellMinerNci60LineTab$CellMiner_2_0
}
rownames(mdaTabOrig) <- mdaTabOrig$`Cell Line Name`

quantFeatures <- c("age", "Epithelial", "p53", "mdr", "doubling time")
mdaQuantTab <- mdaTabOrig[, quantFeatures]
colnames(mdaQuantTab) <- c("AGE", "IS_EPITHELIAL", "IS_P53_MUT", "MDR", "DOUBLING_TIME")

mdaQuantTab$AGE <- as.integer(mdaQuantTab$AGE)

mdaQuantTab$IS_EPITHELIAL[str_trim(mdaQuantTab$IS_EPITHELIAL) == "yes"] <- 1
mdaQuantTab$IS_EPITHELIAL[str_trim(mdaQuantTab$IS_EPITHELIAL) == "no"] <- 0
mdaQuantTab$IS_EPITHELIAL <- as.integer(mdaQuantTab$IS_EPITHELIAL)

mdaQuantTab$IS_P53_MUT[str_trim(mdaQuantTab$IS_P53_MUT) == "MT"] <- 1
mdaQuantTab$IS_P53_MUT[str_trim(mdaQuantTab$IS_P53_MUT) == "WT"] <- 0
mdaQuantTab$IS_P53_MUT <- as.integer(mdaQuantTab$IS_P53_MUT)

mdaQuantTab$MDR <- as.numeric(mdaQuantTab$MDR)

mdaQuantTab$DOUBLING_TIME <- as.numeric(mdaQuantTab$DOUBLING_TIME)

mdaTabSampleInfo <- mdaTabOrig[, setdiff(colnames(mdaTabOrig), quantFeatures)]

stopifnot(all(c(lapply(mdaQuantTab, is.numeric), recursive = TRUE)))
mdaData <- ExpressionSet(t(mdaQuantTab))
stopifnot(is.numeric(exprs(mdaData)))

mdaAnnot <- data.frame(Name = rownames(exprs(mdaData)), Footnote = NA, stringsAsFactors = FALSE)
rownames(mdaAnnot) <- mdaAnnot$Name
mdaAnnot["AGE", "Footnote"] <- "Information from Stinson, et al., (Anticancer Res. 1992 Jul-Aug;12(4):1035-53), DTP, ATCC, and other sources."
mdaAnnot["IS_EPITHELIAL", "Footnote"] <- ""
mdaAnnot["IS_P53_MUT", "Footnote"] <- "p53 status as determined by yeast growth functional assay: PM O'Conner, et al. (Cancer Res. 1997 Oct 1;57(19):4285-300)."
mdaAnnot["MDR", "Footnote"] <- "MDR Function: from DTP site (Lee JS etal., Mol Pharmacol. 1994 Oct;46(4):627-38)."
mdaAnnot["DOUBLING_TIME", "Footnote"] <- "Doubling times described at NCI/DTP site."

featureData(mdaData) <- new("AnnotatedDataFrame", data=mdaAnnot)

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(mdaData)), cmNci60Names))

#--------------------------------------------------------------------------------------------------
# LOAD DATA: DRUG ACTIVITY.
#--------------------------------------------------------------------------------------------------

# activity data -------------------------------------------------------------------------
filePath <- "inst/extdata/cellminer_2_0/DTP_NCI60_ZSCORE.txt"
actTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="",
                         na.strings=c("", "na", "-"))
actTabOrig <- actTabOrig[, c(1:6, 67, 68, 7:66)]

featureDataCols <- 1:8
actTabOrig <- validateDataTab(actTabOrig, keyCol = "NSC #",
                              featureDataColNums = featureDataCols)

drugInfoTab <- actTabOrig[, featureDataCols]
colnames(drugInfoTab) <- c("NSC", "NAME", "FDA_STATUS", "MOA",
                           "PUBCHEM_ID", "SMILES", "TOTAL_EXPS", "TOTAL_EXPS_AFTER_QC")
drugInfoTab$NSC <- as.character(drugInfoTab$NSC)
drugInfoTab$PUBCHEM_ID <- as.integer(drugInfoTab$PUBCHEM_ID)
drugInfoTab$TOTAL_EXPS <- as.integer(drugInfoTab$TOTAL_EXPS)
drugInfoTab$TOTAL_EXPS_AFTER_QC <- as.integer(drugInfoTab$TOTAL_EXPS_AFTER_QC)

actData <- ExpressionSet(as.matrix(actTabOrig[, -featureDataCols]))
stopifnot(identical(rownames(exprs(actData)), rownames(drugInfoTab)))
featureData(actData) <- new("AnnotatedDataFrame", data=drugInfoTab)

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(actData)), cmNci60Names))

# repeat activity data ------------------------------------------------------------------
filePath <- "inst/extdata/cellminer_2_0/DTP_NCI60_EXPS_USED_FOR_ZSCORE_ACT.txt"
usedInZscoreAct <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                              check.names = FALSE, comment.char="", quote="",
                              na.strings=c("", "na", "-"))
usedInZscoreAct$NSC_EXP_NAME <- paste0(usedInZscoreAct[, "NSC #"], "_",
                                       usedInZscoreAct[, "Experiment name"])

filePath <- "inst/extdata/cellminer_2_0/DTP_NCI60_RAW.txt"
rawActTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                            check.names = FALSE, comment.char="", quote="",
                            na.strings=c("", "na", "-"))
rawActTabOrig[, "NSC #"] <- as.character(rawActTabOrig[, "NSC #"])
rawActTabOrig$NSC_EXP_NAME <- paste0(rawActTabOrig[, "NSC #"], "_",
                                     rawActTabOrig[, "Experiment name"])
rawActTabOrig$used_in_zscore <- FALSE
rawActTabOrig <- rawActTabOrig[, c(70, 71, 1:9, 10:69)]


featureDataCols <- 1:11
rawActTabOrig <- validateDataTab(rawActTabOrig, keyCol = "NSC_EXP_NAME",
                                 featureDataColNums = featureDataCols)

for (i in seq_len(nrow(rawActTabOrig))){
  nscExpId <- rawActTabOrig[i, "NSC_EXP_NAME"]
  if (nscExpId %in% usedInZscoreAct$NSC_EXP_NAME){
    rawActTabOrig[i, "used_in_zscore"] <- TRUE
  }
  cat(i, "\n")
}

nscSet <- unique(rawActTabOrig[, "NSC #"])
i <- 0
for (nscStr in nscSet){
  i <- i + 1
  nscRawActData <- rawActTabOrig[(rawActTabOrig[, "NSC #"] == nscStr), ]
  expectedNumUsedInZscore <- unique(nscRawActData[, "Total after quality control"])
  stopifnot(sum(nscRawActData$used_in_zscore) == expectedNumUsedInZscore)
  cat(i, "\n")
}

drugRepeatInfoTab <- rawActTabOrig[, c("NSC_EXP_NAME", "NSC #", "Experiment name",
                                       "used_in_zscore")]
colnames(drugRepeatInfoTab) <- c("NSC_EXP_NAME", "nsc", "experiment", "used_in_zscore")
stopifnot(is.character(drugRepeatInfoTab$nsc))

repeatActData <- ExpressionSet(as.matrix(rawActTabOrig[, -featureDataCols]))
stopifnot(identical(rownames(exprs(repeatActData)), rownames(drugRepeatInfoTab)))
featureData(repeatActData) <- new("AnnotatedDataFrame", data=drugRepeatInfoTab)

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(repeatActData)), cmNci60Names))

#--------------------------------------------------------------------------------------------------
# Make NCI-60 sample info (shared by molData and drugData objects to be constructed).
#--------------------------------------------------------------------------------------------------
cellLineInfo <- loadNciColorSet(returnDf = TRUE)
if (identical(cellLineInfo$abbrCellLines, CellMinerNci60LineTab$CellMiner_1_6)){
  cellLineInfo$abbrCellLines <- CellMinerNci60LineTab$CellMiner_2_0
}
stopifnot(identical(cellLineInfo$abbrCellLines, cmNci60Names))

cellLineOncoTreeTab <- read.table(file="inst/extdata/CellLineToOncoTree.txt",
                                  header=TRUE, sep="\t", stringsAsFactors=FALSE,
                                  check.names = FALSE, comment.char="", quote="",
                                  na.strings=c("", "NA"))
# update "inst/extdata/CellLineToOncoTree.txt" if necessary -----------------------------
if (identical(cellLineOncoTreeTab[1:60, "Name"], CellMinerNci60LineTab$CellMiner_1_6)){
  cellLineOncoTreeTab[1:60, "Name"] <- CellMinerNci60LineTab$CellMiner_2_0
  write.table(cellLineOncoTreeTab, file="inst/extdata/CellLineToOncoTree.txt",
              quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, na="NA")

}
# ---------------------------------------------------------------------------------------
cellLineOncoTreeTab <- cellLineOncoTreeTab[(cellLineOncoTreeTab$DataSource == "NCI-60"), ]
stopifnot(identical(cellLineOncoTreeTab$Name, cmNci60Names))

stopifnot(identical(mdaTabSampleInfo$`Cell Line Name`, cmNci60Names))

nci60Miame <- new("MIAME", name="CellMiner 2.0", lab="NCI/DTB",
                  samples=list(Name = cmNci60Names,
                               TissueType = cellLineInfo$tissues,
                               OncoTree1 = cellLineOncoTreeTab$OncoTree1,
                               OncoTree2 = cellLineOncoTreeTab$OncoTree2,
                               OncoTree3 = cellLineOncoTreeTab$OncoTree3,
                               OncoTree4 = cellLineOncoTreeTab$OncoTree4,
                               Gender = mdaTabSampleInfo[, "sex"],
                               PriorTreatment = mdaTabSampleInfo[, "prior treatment"],
                               Histology = mdaTabSampleInfo[, "histology"],
                               Source = mdaTabSampleInfo[, "source"],
                               Ploidy = mdaTabSampleInfo[, "ploidy"],
                               Institute = mdaTabSampleInfo[, "Institute"],
                               Contributor = mdaTabSampleInfo[, "Contributor"],
                               Reference = mdaTabSampleInfo[, "Reference"]))

#--------------------------------------------------------------------------------------------------
# Make NCI-60 MolData object.
#--------------------------------------------------------------------------------------------------

nci60ESetList <- list()

nci60ESetList[["exp"]] <- expData
nci60ESetList[["xai"]] <- xaiData

nci60ESetList[["cop"]] <- copData
nci60ESetList[["met"]] <- metData
nci60ESetList[["mir"]] <- mirData

nci60ESetList[["mut"]] <- mutData
nci60ESetList[["exo"]] <- exoData

nci60ESetList[["pro"]] <- proData
nci60ESetList[["swa"]] <- swaData

nci60ESetList[["mda"]] <- mdaData

molData <- new("MolData", eSetList = nci60ESetList, sampleData = nci60Miame)

save(molData, file = "data/molData.RData")

#--------------------------------------------------------------------------------------------------
# Make NCI-60 DrugData object.
#--------------------------------------------------------------------------------------------------

drugData <- new("DrugData", act = actData, repeatAct = repeatActData, sampleData = nci60Miame)

save(drugData, file = "data/drugData.RData")

#--------------------------------------------------------------------------------------------------
# UPDATE DRUG DATA TO EXCLUDE ERRONEOUSLY ADDED COMPOUNDS
#--------------------------------------------------------------------------------------------------
library(rcellminer)
library(stringr)
excludedNscSet <- as.character(read.table("~/TMP/excludedNscSet.txt", header = FALSE)[, 1])
excludedNscSet <- stringr::str_trim(excludedNscSet)

# load original data ----------------------------------------------------------
nci60Miame <- rcellminerData::drugData@sampleData

actDataOrigMat   <- exprs(getAct(rcellminerData::drugData))
actDataOrigAnnot <- rcellminer::getFeatureAnnot(rcellminerData::drugData)[["drug"]]
stopifnot(identical(rownames(actDataOrigMat), rownames(actDataOrigAnnot)))

repActDataOrigMat   <- exprs(getRepeatAct(rcellminerData::drugData))
repActDataOrigAnnot <- rcellminer::getFeatureAnnot(rcellminerData::drugData)[["drugRepeat"]]
stopifnot(identical(rownames(repActDataOrigMat), rownames(repActDataOrigAnnot)))
# -----------------------------------------------------------------------------

# determine proper compound set -----------------------------------------------
# length(excludedNscSet)
# dim(actDataOrigMat)
# length(intersect(excludedNscSet, rownames(actDataOrigMat)))
properNscSet <- setdiff(rownames(actDataOrigMat), excludedNscSet)
# -----------------------------------------------------------------------------

# make actData ----------------------------------------------------------------
actDataMat   <- actDataOrigMat[properNscSet, ]
actDataAnnot <- actDataOrigAnnot[properNscSet, ]

stopifnot(identical(rownames(actDataMat), rownames(actDataAnnot)))

####################################
actData <- ExpressionSet(actDataMat)
featureData(actData) <- new("AnnotatedDataFrame", data=actDataAnnot)
# -----------------------------------------------------------------------------

# make repeatActData ----------------------------------------------------------
stopifnot(is.character(repActDataOrigAnnot$nsc))
properNscIndexSet <- which(repActDataOrigAnnot$nsc %in% properNscSet)

repActDataMat   <- repActDataOrigMat[properNscIndexSet, ]
repActDataAnnot <- repActDataOrigAnnot[properNscIndexSet, ]

stopifnot(identical(rownames(repActDataMat), rownames(repActDataAnnot)))
stopifnot(all(repActDataAnnot$nsc %in% rownames(exprs(actData))))

#############################################
repeatActData <- ExpressionSet(repActDataMat)
featureData(repeatActData) <- new("AnnotatedDataFrame", data=repActDataAnnot)
#------------------------------------------------------------------------------

drugData <- new("DrugData", act = actData, repeatAct = repeatActData, sampleData = nci60Miame)

save(drugData, file = "data/drugData.RData")
#--------------------------------------------------------------------------------------------------
