#--------------------------------------------------------------------------------------------------
library(rcellminer)
library(stringr)
tmpEnv <- new.env()
## check package version :  here downloaded from Discover /data/cellminercdb
load("data/molData.RData", envir = tmpEnv)
load("data/drugData.RData", envir = tmpEnv)

nciMiame <- tmpEnv$molData@sampleData

# To check this MIAME below !!!


## nciSclcESetList <- tmpEnv$molData@eSetList
mdaExp <- rcellminer::getAllFeatureData(tmpEnv$molData)[["mda"]]
mda_info = tmpEnv$molData@eSetList$mda@featureData@data

dim(mdaExp); dim(mda_info)
# 39 - 60 ; 39 -  4
#--------------------------------------------------------------------------------------------------
# Remove Is_Epithelial feature and rename KOHN_EMT_PC1 TO EMT
#--------------------------------------------------------------------------------------------------
# REMOVE is_epithelial in row 2
rownames(mdaExp)[2]; rownames(mda_info)[2] # "IS_EPITHELIAL"

mdaExp=mdaExp[-2,]
mda_info=mda_info[-2,]
dim(mdaExp); dim(mda_info)
#38 - 60 ; 38 -  4
stopifnot(identical(rownames(mdaExp),rownames(mda_info)))

## rename KOHN_EMT_PC1 TO EMT
rownames(mdaExp)[5]; rownames(mda_info)[5] # "KOHN_EMT_PC1"

rownames(mdaExp)[5]= "EMT"; rownames(mda_info)[5]= "EMT"
mda_info$ID = as.character(mda_info$ID)
mda_info$ID[5]="EMT"
stopifnot(identical(rownames(mdaExp),rownames(mda_info)))
##
mdaData <- ExpressionSet(mdaExp)
stopifnot(is.numeric(exprs(mdaData)))


featureData(mdaData) <- new("AnnotatedDataFrame", data=mda_info)


# saving now
nciSetList=tmpEnv$molData@eSetList

nciSetList[["mda"]] <- mdaData

molData <- new("MolData", eSetList = nciSetList, sampleData = nciMiame)

save(molData, file = "data/molData.RData")

## END
