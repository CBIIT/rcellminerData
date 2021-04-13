#--------------------------------------------------------------------------------------------------
library(rcellminer)
library(stringr)
tmpEnv <- new.env()
## check package version :  here downloaded from Discover /data/cellminercdb
load("data/molData.RData", envir = tmpEnv)
# load("data/drugData.RData", envir = tmpEnv)

nciMiame <- tmpEnv$molData@sampleData
cells <- nciMiame@samples$Name
# To check this MIAME below !!!


## nciSclcESetList <- tmpEnv$molData@eSetList
mdaExp <- rcellminer::getAllFeatureData(tmpEnv$molData)[["mda"]]
mda_info = tmpEnv$molData@eSetList$mda@featureData@data

dim(mdaExp); dim(mda_info)
# 41 - 60 ; 41 -  4
#--------------------------------------------------------------------------------------------------
# Imported new data from cellminer 2.5
#--------------------------------------------------------------------------------------------------
# KEEP REP STRESS

mdadat = read.delim("inst/extdata/phenodata_70.txt", stringsAsFactors = F, check.names = F, row.names = 1)
mdainf = read.delim("inst/extdata/phenoinfo_70_correct.txt", stringsAsFactors = F, check.names = F, row.names = 1)
dim(mdadat); dim(mdainf)
# 70 60 - 70 4
stopifnot(is.numeric(matrix(mdadat))) # F
mdadat2 = apply(mdadat, 2, as.numeric)
class(mdadat2)
stopifnot(is.numeric(mdadat2))
rownames(mdadat2) = rownames(mdadat)

stopifnot(identical(rownames(mdadat2), rownames(mdainf)))
stopifnot(identical(colnames(mdadat2), cells))
##
dim(mdadat2); dim(mdainf)
#
mdaData <- ExpressionSet(as.matrix(mdadat2))
featureData(mdaData) <- new("AnnotatedDataFrame", data=mdainf)


# saving now
nciSetList=tmpEnv$molData@eSetList

nciSetList[["mda"]] <- mdaData

molData <- new("MolData", eSetList = nciSetList, sampleData = nciMiame)

save(molData, file = "data/molData.RData")

## END
