#--------------------------------------------------------------------------------------------------
#  import public drug activity and good experiments from cellminer database 2.5
#--------------------------------------------------------------------------------------------------

library(rcellminer)
library(stringr)


tmpEnv <- new.env()
#
# load("data/molData.RData", envir = tmpEnv)
load("data/drugData.RData", envir = tmpEnv)

nci60Miame <- tmpEnv$drugData@sampleData
celllines = nci60Miame@samples$Name

curdruginfo = tmpEnv$drugData@act@featureData@data
curRepdruginfo = tmpEnv$drugData@repeatAct@featureData@data
#--------------------------------------------------------------------------------------------------
# read new data
#--------------------------------------------------------------------------------------------------

drugact = read.delim("inst/extdata/public_drug_act_v2.5_correct.txt", row.names = 1, stringsAsFactors = F,check.names = F)
druginfo = read.table("inst/extdata/public_drug_info_v2.5_correct.txt", encoding = "UTF-8", header = T,stringsAsFactors = F, row.names = 1)
dim(drugact); dim(druginfo)
# [1] 23765    60
# [1] 23765   10
stopifnot(identical(rownames(drugact), rownames(druginfo))) # TRUE

## iconv(druginfo$NAME, to = "ASCII", sub="?")

druginfo = druginfo[, colnames(curdruginfo)]
colnames(druginfo); dim(druginfo) # 23765     8

stopifnot(identical(rownames(drugact), rownames(druginfo))) # TRUE
stopifnot(identical(colnames(drugact), celllines)) # TRUE


drugRepact = read.delim("inst/extdata/public_drug_Rep_act_v2.5_correct.txt", row.names = 1, stringsAsFactors = F,check.names = F)
drugRepinfo = read.delim("inst/extdata/public_drug_Rep_info_v2.5_correct.txt", row.names = 1, stringsAsFactors = F,check.names = F)
dim(drugRepact); dim(drugRepinfo)
# [1] 36668    60
# [1] 36668     4
stopifnot(identical(rownames(drugRepact), rownames(drugRepinfo))) # TRUE
stopifnot(identical(colnames(drugRepact), celllines)) # TRUE

colnames(drugRepinfo) = toupper(colnames(drugRepinfo))

##
actData <- ExpressionSet(as.matrix(drugact))
featureData(actData) <- new("AnnotatedDataFrame", data=druginfo)

repeatActData <- ExpressionSet(as.matrix(drugRepact))
featureData(repeatActData) <- new("AnnotatedDataFrame", data=drugRepinfo)

##


# saving now -------------------------------

drugData <- new("DrugData", act = actData, repeatAct = repeatActData, sampleData = nci60Miame)

save(drugData, file = "data/drugData.RData")
## END
