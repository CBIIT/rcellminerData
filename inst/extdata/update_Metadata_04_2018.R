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


#--------------------------------------------------------------------------------------------------
# LOAD AND PREPARE METADATA CELL LINE DATA
#--------------------------------------------------------------------------------------------------
minfo=read.delim("./inst/extdata/updated_metadata_042018.txt",row.names = 1)
length(which(rownames(minfo)==rownames(mda_info)))
## OK 39

## add new data to current mda dataframe

##
mdaData <- ExpressionSet(mdaExp)
stopifnot(is.numeric(exprs(mdaData)))


featureData(mdaData) <- new("AnnotatedDataFrame", data=minfo)


# saving now
nciSetList=tmpEnv$molData@eSetList

nciSetList[["mda"]] <- mdaData

molData <- new("MolData", eSetList = nciSetList, sampleData = nciSclcMiame)

save(molData, file = "data/molData.RData")

## END
