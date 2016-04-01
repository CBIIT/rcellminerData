library(rcellminer)
library(xlsx)
library(stringr)

#--------------------------------------------------------------------------------------------------
# LOAD DATA: GENE COPY.
#--------------------------------------------------------------------------------------------------
# http://discovery.nci.nih.gov/cellminerint/loadDownload.do
# Select: [DNA: Combined aCGH, select: gene summary].
# NOTE: EMPTY ROW BETWEEN HEADER AND DATA REMOVED.

filePath <- "inst/extdata/cellminer_1_6/NCI60_DNA__Combined_aCGH_gene_summary.txt"
copTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                     check.names = FALSE, comment.char="", quote="", na.strings="null")
copTabOrig[, "Probe Name"] <- str_trim(copTabOrig[, "Probe Name"])

stopifnot(all(copTabOrig[, "Probe Name"] != ""))
stopifnot(all(copTabOrig[, "Probe Name"] != "1-Mar"))   # No Excel gene to date conversion.
stopifnot(all(!duplicated(copTabOrig[, "Probe Name"])))

rownames(copTabOrig) <- copTabOrig[, "Probe Name"]

featureDataCols <- 1:4
stopifnot(all(c(lapply(copTabOrig[, -featureDataCols], is.numeric), recursive = TRUE)))

copData <- ExpressionSet(as.matrix(copTabOrig[, -featureDataCols]))
featureData(copData) <- new("AnnotatedDataFrame", data=copTabOrig[, featureDataCols])

#--------------------------------------------------------------------------------------------------
# LOAD DATA: MRNA EXPRESSION.
#--------------------------------------------------------------------------------------------------
# http://discovery.nci.nih.gov/cellminerint/loadDownload.do
# Select: [RNA: 5 Platform Gene Transcript, select: Average z score].

#----[z score data]------------------------------------------------------------
filePath <- "inst/extdata/cellminer_1_6/GENE_ZSCORES.txt"
expTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE, skip = 8,
                         check.names = FALSE, comment.char="", quote="", na.strings="na")

# remove empty column that causes problems.
expTabOrig <- expTabOrig[, -67]
# expTabOrig[(nrow(expTabOrig)-3):nrow(expTabOrig), 1:5]

expTabOrig[, "Gene symbol"] <- str_trim(expTabOrig[, "Gene symbol"])

stopifnot(all(expTabOrig[, "Gene symbol"] != ""))
stopifnot(all(expTabOrig[, "Gene symbol"] != "1-Mar"))    # No Excel gene to date conversion.
stopifnot(all(!duplicated(expTabOrig[, "Gene symbol"])))

rownames(expTabOrig) <- expTabOrig[, "Gene symbol"]

featureDataCols <- c(1:4, 65, 66)
stopifnot(all(c(lapply(expTabOrig[, -featureDataCols], is.numeric), recursive = TRUE)))

expData <- ExpressionSet(as.matrix(expTabOrig[, -featureDataCols]))
featureData(expData) <- new("AnnotatedDataFrame", data=expTabOrig[, featureDataCols])

#----[average log2 intensity data]---------------------------------------------
filePath <- "inst/extdata/cellminer_1_6/GENE_QA_AVERAGES.txt"
xaiTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE, skip = 7,
                         check.names = FALSE, comment.char="", quote="", na.strings="na")

# remove empty rows/columns that causes problems.
xaiTabOrig <- xaiTabOrig[, -(64:ncol(xaiTabOrig))]
# xaiTabOrig[(nrow(xaiTabOrig)-3):nrow(xaiTabOrig), 1:5]
xaiTabOrig <- xaiTabOrig[-nrow(xaiTabOrig), ]

xaiTabOrig[, "Gene symbol"] <- str_trim(xaiTabOrig[, "Gene symbol"])

stopifnot(all(xaiTabOrig[, "Gene symbol"] != ""))
stopifnot(all(xaiTabOrig[, "Gene symbol"] != "1-Mar"))    # No Excel gene to date conversion.
stopifnot(all(!duplicated(xaiTabOrig[, "Gene symbol"])))

rownames(xaiTabOrig) <- xaiTabOrig[, "Gene symbol"]

featureDataCols <- c(1, 62, 63)
stopifnot(all(c(lapply(xaiTabOrig[, -featureDataCols], is.numeric), recursive = TRUE)))

xaiData <- ExpressionSet(as.matrix(xaiTabOrig[, -featureDataCols]))
featureData(xaiData) <- new("AnnotatedDataFrame", data=xaiTabOrig[, featureDataCols])

#--------------------------------------------------------------------------------------------------
# LOAD DATA: EXOME/MUTATION.
#--------------------------------------------------------------------------------------------------
# http://discovery.nci.nih.gov/cellminerint/loadDownload.do
# Select: [DNA: Exome Seq, select: none].
# NOTE: EMPTY ROW BETWEEN HEADER AND DATA REMOVED.

filePath <- "inst/extdata/cellminer_1_6/NCI60_DNA__Exome_Seq_none.txt"
exoTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings=c("-", "", " "))
stopifnot(all(!duplicated(exoTabOrig[, "Probe Name"])))
rownames(exoTabOrig) <- exoTabOrig[, "Probe Name"]
exoTabOrig <- exoTabOrig[sort(rownames(exoTabOrig)), ]

filePath <- "inst/extdata/cellminer_1_6/exome_gene_probeIds_mapping.txt"
geneMapTab <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings=c("-", "", " "))
stopifnot(all(!duplicated(geneMapTab[, "probe.ids"])))
rownames(geneMapTab) <- geneMapTab[, "probe.ids"]
geneMapTab <- geneMapTab[sort(rownames(geneMapTab)), ]

stopifnot(identical(exoTabOrig[, "Probe Name"], geneMapTab[, "probe.ids"]))

exoInfoCols <- c(1, 3:14)
exoInfoTab <- cbind(geneMapTab, exoTabOrig[, exoInfoCols])
colnames(exoInfoTab) <- str_replace_all(colnames(exoInfoTab), pattern = " ", replacement = "_")

exoDataCols <- c(15:74)
exoDataMat <- as.matrix(exoTabOrig[, exoDataCols])
stopifnot(all(apply(exoDataMat, MARGIN = 2, is.numeric)))
exoDataMat[which(is.na(exoDataMat))] <- 0

binMutData <- getBinaryMutationData(exoInfoTab, exoDataMat)

exoData <- ExpressionSet(exoDataMat)
featureData(exoData) <- new("AnnotatedDataFrame", data=exoInfoTab)

mutData <- ExpressionSet(binMutData)

#--------------------------------------------------------------------------------------------------
# LOAD DATA: MICRORNA EXPRESSION.
#--------------------------------------------------------------------------------------------------
# http://discovery.nci.nih.gov/cellminerint/loadDownload.do
# Select: [RNA: Agilent Human microRNA (V2)].
# NOTE: EMPTY ROW BETWEEN HEADER AND DATA REMOVED.

filePath <- "inst/extdata/cellminer_1_6/NCI60_RNA__Agilent_Human_microRNA_(V2)_GeneSpringGX.txt"
mirTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="na")
mirTabOrig[, "Probe Name"] <- str_trim(mirTabOrig[, "Probe Name"])

stopifnot(all(mirTabOrig[, "Probe Name"] != ""))
stopifnot(all(mirTabOrig[, "Probe Name"] != "1-Mar"))   # No Excel gene to date conversion.
stopifnot(all(!duplicated(mirTabOrig[, "Probe Name"])))

rownames(mirTabOrig) <- mirTabOrig[, "Probe Name"]

featureDataCols <- 1:4
stopifnot(all(c(lapply(mirTabOrig[, -featureDataCols], is.numeric), recursive = TRUE)))

mirData <- ExpressionSet(as.matrix(mirTabOrig[, -featureDataCols]))
featureData(mirData) <- new("AnnotatedDataFrame", data=mirTabOrig[, featureDataCols])
#--------------------------------------------------------------------------------------------------
# LOAD DATA: PROTEIN EXPRESSION.
#--------------------------------------------------------------------------------------------------
# http://discovery.nci.nih.gov/cellminerint/loadDownload.do
# Select: [Protein: Lysate Array, select: log2].
# NOTE: EMPTY ROW BETWEEN HEADER AND DATA REMOVED.

filePath <- "inst/extdata/cellminer_1_6/nci60_Protein__Lysate_Array_log2.txt"
proTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                         check.names = FALSE, comment.char="", quote="", na.strings="na")

proTabOrig[, "Probe Name"] <- str_trim(proTabOrig[, "Probe Name"])

stopifnot(all(proTabOrig[, "Probe Name"] != ""))
stopifnot(all(proTabOrig[, "Probe Name"] != "1-Mar"))   # No Excel gene to date conversion.
stopifnot(all(!duplicated(proTabOrig[, "Probe Name"])))

rownames(proTabOrig) <- proTabOrig[, "Probe Name"]

featureDataCols <- 1:4
stopifnot(all(c(lapply(proTabOrig[, -featureDataCols], is.numeric), recursive = TRUE)))

proData <- ExpressionSet(as.matrix(proTabOrig[, -featureDataCols]))
featureData(proData) <- new("AnnotatedDataFrame", data=proTabOrig[, featureDataCols])

#--------------------------------------------------------------------------------------------------
# LOAD DATA: CELL LINE  METADATA.
#--------------------------------------------------------------------------------------------------

filePath <- "inst/extdata/cellminer_1_6/nci60_cellline_metadata.xls"
mdaTabOrig <- read.xlsx2(file = filePath, sheetIndex = 1, stringsAsFactors = FALSE)
stopifnot(all(!duplicated(mdaTabOrig$Cell.Line.Name)))
rownames(mdaTabOrig) <- mdaTabOrig$Cell.Line.Name

quantFeatures <- c("age", "Epithelial", "p53", "mdr", "doublingtime")
mdaQuantTab <- mdaTabOrig[, quantFeatures]
colnames(mdaQuantTab) <- c("age", "is_epithelial", "is_p53_mut", "mdr", "doublingtime")

mdaQuantTab$age <- as.integer(mdaQuantTab$age)

mdaQuantTab$is_epithelial[str_trim(mdaQuantTab$is_epithelial) == "yes"] <- 1
mdaQuantTab$is_epithelial[str_trim(mdaQuantTab$is_epithelial) == "no"] <- 0
mdaQuantTab$is_epithelial <- as.integer(mdaQuantTab$is_epithelial)

mdaQuantTab$is_p53_mut[str_trim(mdaQuantTab$is_p53_mut) == "?"] <- NA
mdaQuantTab$is_p53_mut[str_trim(mdaQuantTab$is_p53_mut) == "MT"] <- 1
mdaQuantTab$is_p53_mut[str_trim(mdaQuantTab$is_p53_mut) == "WT"] <- 0
mdaQuantTab$is_p53_mut <- as.integer(mdaQuantTab$is_p53_mut)

mdaQuantTab$mdr <- as.numeric(mdaQuantTab$mdr)

mdaQuantTab$doublingtime <- as.numeric(mdaQuantTab$doublingtime)

mdaTabSampleInfo <- mdaTabOrig[, setdiff(colnames(mdaTabOrig), quantFeatures)]

# RENAME COLUMN NAMED 'Institute.'
colnames(mdaTabSampleInfo) <- c(colnames(mdaTabSampleInfo)[1:7], "Institute",
                                colnames(mdaTabSampleInfo)[9:10])

mdaData <- ExpressionSet(t(mdaQuantTab))

#--------------------------------------------------------------------------------------------------
# LOAD DATA: DRUG ACTIVITY.
#--------------------------------------------------------------------------------------------------
# http://discovery.nci.nih.gov/cellminerint/loadDownload.do
# Select: [Compound activity: DTP NCI-60].

filePath <- "inst/extdata/cellminer_1_6/DTP_NCI60.txt"
actTabOrig <- read.table(file=filePath, header=TRUE, sep="\t", stringsAsFactors=FALSE, skip = 8,
                         check.names = FALSE, comment.char="", quote="", na.strings="na")

# remove empty rows/columns that causes problems.
actTabOrig <- actTabOrig[, -68]

actTabOrig[, "NSC #"] <- str_trim(actTabOrig[, "NSC #"])

stopifnot(all(actTabOrig[, "NSC #"] != ""))
stopifnot(all(!duplicated(actTabOrig[, "NSC #"])))

rownames(actTabOrig) <- actTabOrig[, "NSC #"]

featureDataCols <- c(1:5, 66, 67)
stopifnot(all(c(lapply(actTabOrig[, -featureDataCols], is.numeric), recursive = TRUE)))

actData <- ExpressionSet(as.matrix(actTabOrig[, -featureDataCols]))
drugInfoTab <- actTabOrig[, featureDataCols]
colnames(drugInfoTab) <- c("NSC", "NAME", "FDA_STATUS", "MOA", "PUBCHEM_ID", "TOTAL_EXPS", "TOTAL_EXPS_AFTER_QC")

naChars <- c("-", " ", "")
drugInfoTab$NAME[drugInfoTab$NAME %in% naChars] <- NA
drugInfoTab$FDA_STATUS[drugInfoTab$FDA_STATUS %in% naChars] <- NA
drugInfoTab$MOA[drugInfoTab$MOA %in% naChars] <- NA

stopifnot(!exists("annot"))
stopifnot(!exists("experiment.info"))
stopifnot(!exists("logGI50"))
load("inst/extdata/cellminer_1_6/GI50.public.int1_6.Rdata")
stopifnot(identical(colnames(logGI50), colnames(exprs(actData))))

nscToSmiles <- unlist(by(annot, INDICES = as.character(annot$nsc), FUN = function(x) x$smiles[1]))
stopifnot(is.character(drugInfoTab$NSC))
stopifnot(all(drugInfoTab$NSC %in% names(nscToSmiles)))
drugInfoTab$SMILES <- nscToSmiles[drugInfoTab$NSC]

featureData(actData) <- new("AnnotatedDataFrame", data=drugInfoTab)

repeatActData <- ExpressionSet(logGI50)
featureData(repeatActData) <- new("AnnotatedDataFrame",
  data=cbind(nsc = as.character(annot$nsc), experiment.info, stringsAsFactors = FALSE))

#--------------------------------------------------------------------------------------------------
# Make NCI-60 sample info (shared by molData and drugData objects to be constructed).
#--------------------------------------------------------------------------------------------------
cellLineInfo <- loadNciColorSet(returnDf = TRUE)
stopifnot(identical(cellLineInfo$abbrCellLines, colnames(exprs(actData))))

cellLineOncoTreeTab <- read.table(file="inst/extdata/CellLineToOncoTree.txt",
                                  header=TRUE, sep="\t", stringsAsFactors=FALSE,
                                  check.names = FALSE, comment.char="", quote="", na.strings="")
cellLineOncoTreeTab <- cellLineOncoTreeTab[(cellLineOncoTreeTab$DataSource == "NCI-60"), ]
stopifnot(identical(cellLineInfo$abbrCellLines, cellLineOncoTreeTab$Name))

nci60Miame <- new("MIAME", name="CellMiner", lab="NCILMP",
                  samples=list(Name = colnames(exprs(actData)),
                               TissueType = cellLineInfo$tissues,
                               OncoTree1 = cellLineOncoTreeTab$OncoTree1,
                               OncoTree2 = cellLineOncoTreeTab$OncoTree2,
                               OncoTree3 = cellLineOncoTreeTab$OncoTree3,
                               OncoTree4 = cellLineOncoTreeTab$OncoTree4,
                               Gender = mdaTabSampleInfo$sex,
                               PriorTreatment = mdaTabSampleInfo$prior.treatment,
                               Histology = mdaTabSampleInfo$histology,
                               Source = mdaTabSampleInfo$source,
                               Ploidy = mdaTabSampleInfo$ploidy,
                               Institute = mdaTabSampleInfo$Institute,
                               Contributor = mdaTabSampleInfo$Contributor,
                               Reference = mdaTabSampleInfo$Reference))

#--------------------------------------------------------------------------------------------------
# Make NCI-60 MolData object.
#--------------------------------------------------------------------------------------------------

nci60ESetList <- list()
nci60ESetList[["cop"]] <- copData
nci60ESetList[["exp"]] <- expData
nci60ESetList[["xai"]] <- xaiData
nci60ESetList[["exo"]] <- exoData
nci60ESetList[["mut"]] <- mutData
nci60ESetList[["mir"]] <- mirData
nci60ESetList[["pro"]] <- proData
nci60ESetList[["mda"]] <- mdaData

molData <- new("MolData", eSetList = nci60ESetList, sampleData = nci60Miame)

save(molData, file = "data/molData.RData")

#--------------------------------------------------------------------------------------------------
# Make NCI-60 DrugData object.
#--------------------------------------------------------------------------------------------------

drugData <- new("DrugData", act = actData, repeatAct = repeatActData, sampleData = nci60Miame)

save(drugData, file = "data/drugData.RData")

#--------------------------------------------------------------------------------------------------
