# add RepStress anew metadata
#--------------------------------------------------------------------------------------------------
library(rcellminer)
library(stringr)
##
tmpEnv <- new.env()
#
load("data/molData.RData", envir = tmpEnv)
# load("data/drugData.RData", envir = tmpEnv)

nci60SetList=tmpEnv$molData@eSetList


nci60Miame <- tmpEnv$molData@sampleData
cells <- nci60Miame@samples$Name

## read new data

## 1 repstress -----------------------------------------------------------------------------------------

mdaExp <- exprs(tmpEnv$molData@eSetList$mda)
mda_info = tmpEnv$molData@eSetList$mda@featureData@data
dim(mdaExp); dim(mda_info)
# 70 60 - 70 -4
repstress = read.csv("inst/extdata/nci60_RepStress.csv", row.names = 1)
dim(repstress)
stopifnot(identical(rownames(repstress), colnames(mdaExp)))
colnames(repstress) = "RepStress"

mdaExp = rbind(mdaExp,t(repstress))
dim(mdaExp) # 71 - 60

mda_info = rbind(mda_info, c("RepStress","Replication Stress transcript expression signature score based on first principal component weights of 18 genes","Replication stress signature","Jo et al, Mol Cancer Therapeutics, 2021. PMID: 34045232"))
dim(mda_info) # 71 - 4
rownames(mda_info)[71] = "RepStress"
stopifnot(identical(rownames(mda_info), rownames(mdaExp)))

mdaData <- ExpressionSet(mdaExp)
featureData(mdaData) <- new("AnnotatedDataFrame", data=mda_info)

# saving now

nci60SetList[["mda"]] <- mdaData

molData <- new("MolData", eSetList = nci60SetList, sampleData = nci60Miame)

save(molData, file = "data/molData.RData")

## END
