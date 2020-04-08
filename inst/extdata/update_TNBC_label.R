# ----------------------------------------------------------------------------------------
# update label from "TNBC" to "Triple Negative" status for Breast cell lines (curated information from Sudhir Varma. December 2019)
#
#--------------------------------------------------------------------------------------------------
library(rcellminer)
library(stringr)

##
tmpEnv <- new.env()
load("data/molData.RData", envir = tmpEnv)
load("data/drugData.RData", envir = tmpEnv)

nci60Miame <- tmpEnv$molData@sampleData

## update label

tn = nci60Miame@samples$TNBC
table(tn)
tn[which(tn=="TNBC")] = "Breast Triple Negative"
table(tn)

# update moldata and drugdata objects
##
nci60Miame@samples$TNBC = tn

nci60ESetList=tmpEnv$molData@eSetList
act  = tmpEnv$drugData@act
repact  = tmpEnv$drugData@repeatAct


## take old miame

molData <- new("MolData", eSetList = nci60ESetList, sampleData = nci60Miame)

save(molData, file = "data/molData.RData")

drugData <- new("DrugData", act = act, repeatAct = repact,
                sampleData = nci60Miame)

save(drugData, file = "data/drugData.RData")


### STOP Here
### --------------------------------------------------------------
