# --------------------------------------------------------------------------------------
# Update drug name, MOA and FDA status from SandboxDB done on August 14, 2020
# --------------------------------------------------------------------------------------
library(rcellminer)
library(stringr)
tmpEnv <- new.env()
##
load("data/molData.RData", envir = tmpEnv)
load("data/drugData.RData", envir = tmpEnv)

nci60Miame <- tmpEnv$molData@sampleData
nci60ESetList=tmpEnv$molData@eSetList

drugact = exprs(tmpEnv$drugData@act)
druginfo = tmpEnv$drugData@act@featureData@data

repAct = tmpEnv$drugData@repeatAct

# drugRepact = exprs(tmpEnv$drugData@repeatAct)
# drugRepinfo = tmpEnv$drugData@repeatAct@featureData@data

dim(drugact); dim(druginfo)        # 22257    60 ; 22257    8
## dim(drugRepact); dim(drugRepinfo)  # 34623    60 ; 34623     4
### ------------------------------------------
## read new data
## -------------------------------------------

drugann = read.delim("./inst/extdata/NCI60_druginfo4_22257_Sandbox_08142020.txt",stringsAsFactors = F, check.names = F)
dim(drugann) # 22257 - 4
rownames(drugann) = drugann$nsc

stopifnot(identical(rownames(drugann),rownames(druginfo)))

## --------------------------------------------
## update name, MOA and testing status
## --------------------------------------------
druginfo$NAME = drugann$name
druginfo$FDA_STATUS = drugann$testing_status
druginfo$MOA = drugann$mechanism

dim(druginfo)
##
actData <- ExpressionSet(drugact)
featureData(actData) <- new("AnnotatedDataFrame", data=druginfo)


drugData <- new("DrugData", act = actData, repeatAct = repAct, sampleData = nci60Miame)

save(drugData, file = "data/drugData.RData")

