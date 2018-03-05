#--------------------------------------------------------------------------------------------------
library(rcellminer)
library(stringr)
tmpEnv <- new.env()
## check package version :  here downloaded from Discover /data/cellminercdb
load("data/molData.RData", envir = tmpEnv)
load("data/drugData.RData", envir = tmpEnv)

nciSclcMiame <- tmpEnv$molData@sampleData

# To check this MIAME below !!!


## nciSclcESetList <- tmpEnv$molData@eSetList
mdaExp <- rcellminer::getAllFeatureData(tmpEnv$molData)[["mda"]]
mda_info = tmpEnv$molData@eSetList$mda@featureData@data

cells <- colnames(mdaExp)
drug1 = exprs(tmpEnv$drugData@act)
actdrug1 = exprs(tmpEnv$drugData@repeatAct)
cellsdrug=colnames(drug1)
cellsactdrug=colnames(actdrug1)
# checking columns
length(which(cells==cellsdrug)) ## 60
length(which(cells==cellsactdrug)) ## 60
#--------------------------------------------------------------------------------------------------
# LOAD AND PREPARE METADATA CELL LINE DATA
#--------------------------------------------------------------------------------------------------
minfo= read.delim("./inst/extdata/extra-metadata-nci60.txt",header = F)[1:2,-1]

expDataTab=read.delim("./inst/extdata/extra-metadata-nci60.txt",skip=2,row.names = 1, na.strings = c("na","NA"))
#expDataTab <- read.csv("inst/extdata/METHYL.csv", check.names = FALSE, stringsAsFactors = FALSE)
dim(expDataTab)
# 60 - 33
# checking order cell lines
length(which(cells==rownames(expDataTab))) ## 60 OK !!!
length(unique(colnames(expDataTab))) ## 33 OK !!!

colnames(minfo)=colnames(expDataTab)
rownames(minfo)=c("Name","Footnote")
## add new data to current mda dataframe

new_mda=rbind(mdaExp,t(expDataTab))
dim(new_mda) #39-60

new_info=rbind(mda_info, t(minfo))
dim(new_info) # 39-2
##
mdaData <- ExpressionSet(new_mda)
stopifnot(is.numeric(exprs(mdaData)))


featureData(mdaData) <- new("AnnotatedDataFrame", data=new_info)

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(mdaData)), cells))
stopifnot(identical(rownames(exprs(mdaData)), rownames(new_info)))

# saving now
nciSetList=tmpEnv$molData@eSetList

nciSetList[["mda"]] <- mdaData

molData <- new("MolData", eSetList = nciSetList, sampleData = nciSclcMiame)

save(molData, file = "data/molData.RData")

## END
