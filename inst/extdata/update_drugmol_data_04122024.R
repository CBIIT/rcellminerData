#--------------------------------------------------------------------------------------------------
#  import public drug activity and good experiments from cellminer database 2.12
#--------------------------------------------------------------------------------------------------

library(rcellminer)
library(stringr)


tmpEnv <- new.env()
#
load("data/molData.RData", envir = tmpEnv)
load("data/drugData.RData", envir = tmpEnv)

nci60Miame <- tmpEnv$drugData@sampleData
celllines = nci60Miame@samples$Name
nciSetList=tmpEnv$molData@eSetList
##

curdruginfo = tmpEnv$drugData@act@featureData@data
curRepdruginfo = tmpEnv$drugData@repeatAct@featureData@data
dim(curdruginfo); dim(curRepdruginfo)
# 24980       8  ; 38448     4
#--------------------------------------------------------------------------------------------------
# read new data
#--------------------------------------------------------------------------------------------------

drugact = read.delim("inst/extdata/public_drug_act_v2_12.txt", row.names = 1, stringsAsFactors = F,check.names = F)
druginfo = read.table("inst/extdata/public_drug_info_v2_12.txt", encoding = "UTF-8", header = T,stringsAsFactors = F, row.names = 1)
dim(drugact); dim(druginfo)
# [1] 25147    60
# [1] 25147    11 (not 10)
stopifnot(identical(rownames(drugact), rownames(druginfo))) # TRUE
stopifnot(identical(colnames(drugact), celllines ))
## iconv(druginfo$NAME, to = "ASCII", sub="?")

druginfo = druginfo[, colnames(curdruginfo)]
colnames(druginfo); dim(druginfo) # 25147    8

stopifnot(identical(rownames(drugact), rownames(druginfo))) # TRUE
stopifnot(identical(colnames(drugact), celllines)) # TRUE compare to current cell lines


drugRepact = read.delim("inst/extdata/public_drug_Rep_act_v2_12.txt", row.names = 1, stringsAsFactors = F,check.names = F)
drugRepinfo = read.delim("inst/extdata/public_drug_Rep_info_v2_12.txt", row.names = 1, stringsAsFactors = F,check.names = F)
dim(drugRepact); dim(drugRepinfo)
# [1] 38722    60
# [1] 38722    4
stopifnot(identical(rownames(drugRepact), rownames(drugRepinfo))) # TRUE
stopifnot(identical(colnames(drugRepact), celllines)) # TRUE

colnames(drugRepinfo) = toupper(colnames(drugRepinfo))
stopifnot(identical(colnames((curRepdruginfo)), colnames(drugRepinfo)))
##
actData <- ExpressionSet(as.matrix(drugact))
featureData(actData) <- new("AnnotatedDataFrame", data=druginfo)

repeatActData <- ExpressionSet(as.matrix(drugRepact))
featureData(repeatActData) <- new("AnnotatedDataFrame", data=drugRepinfo)
# saving drug now

drugData <- new("DrugData", act = actData, repeatAct = repeatActData, sampleData = nci60Miame)

save(drugData, file = "data/drugData.RData")

## Add Cell Surface data -----------------------------------------

surfdat = read.delim("inst/extdata/NCI60_cellsurface_329_log2.txt", stringsAsFactors = F, check.names = F)
dim(surfdat) # 329 - 66

surfinfo = surfdat[, 1:6]; dim(surfinfo) # 329 - 6
surfdat = surfdat[, 7:66]; dim(surfdat)  # 329 - 60
rownames(surfinfo) = surfinfo$receptor_esc
rownames(surfdat) = surfinfo$receptor_esc
colnames(surfinfo)[1] = "receptor_id"

stopifnot(identical(colnames(surfdat), celllines))
range(surfdat, na.rm=T) ## 0.00000 11.95674

surData <- ExpressionSet(as.matrix(surfdat))
featureData(surData) <- new("AnnotatedDataFrame", data=surfinfo)

# saving  surface receptors data (flow cytometry)

nciSetList[["sur"]] <- surData

molData <- new("MolData", eSetList = nciSetList, sampleData = nci60Miame)

save(molData, file = "data/molData.RData")


## END
