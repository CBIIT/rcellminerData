# ----------------------------------------------------------------------------------------
# add Triple Negative status for Breast cell lines (curated information from Sudhir Varma. December 2019)
#
#--------------------------------------------------------------------------------------------------
library(rcellminer)
library(stringr)
## read TNBC status file
tnbc = read.delim("./inst/extdata/TNBC_cellminerCDB_curated.txt",stringsAsFactors = F)
## select NCI60 cell lines -------------
i=which(tnbc$source=="NCI60" & !is.na(tnbc$Combined) & tnbc$Combined!=""); length(i) # 3
tnbc = tnbc[i,]; dim(tnbc) # 3 -11
##
tmpEnv <- new.env()
load("data/molData.RData", envir = tmpEnv)
load("data/drugData.RData", envir = tmpEnv)

nci60Miame <- tmpEnv$molData@sampleData

cellname= nci60Miame@samples$Name
nb = length(cellname) # 60
tn = replicate(nb, NA)
# update for all cell lines
j= match(tnbc$Name, cellname); length(j)
tn[j]= "TNBC"

## add tn annotation

nci60Miame@samples$TNBC=tn

# update moldata and drugdata objects
##
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
