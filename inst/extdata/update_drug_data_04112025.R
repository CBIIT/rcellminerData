#--------------------------------------------------------------------------------------------------
#  import public drug activity and good experiments from cellminer database 2.13
#--------------------------------------------------------------------------------------------------

library(rcellminer)
library(stringr)


tmpEnv <- new.env()
#
load("data/drugData.RData", envir = tmpEnv)

nci60Miame <- tmpEnv$drugData@sampleData
celllines = nci60Miame@samples$Name
# nciSetList=tmpEnv$molData@eSetList
##

curdruginfo = tmpEnv$drugData@act@featureData@data
curRepdruginfo = tmpEnv$drugData@repeatAct@featureData@data
dim(curdruginfo); dim(curRepdruginfo)
# 25293     8  ;  38932     4
#--------------------------------------------------------------------------------------------------
# read new data
#--------------------------------------------------------------------------------------------------

drugact = read.delim("inst/extdata/public_drug_act_2_14.txt", row.names = 1, stringsAsFactors = F,check.names = F)
druginfo = read.table("inst/extdata/public_drug_info_2_14.txt", encoding = "UTF-8", header = T,stringsAsFactors = F, row.names = 1)
dim(drugact); dim(druginfo)
# [1] 25374    60
# [1] 25374    11 (not 10)
stopifnot(identical(rownames(drugact), rownames(druginfo))) # TRUE
stopifnot(identical(colnames(drugact), celllines ))
## iconv(druginfo$NAME, to = "ASCII", sub="?")

colnames(curdruginfo)[5] = "PUBCHEM_SID" # replaced pubchem_id with pubchem_sid
length(intersect(colnames(druginfo), colnames(curdruginfo))) # 8

druginfo = druginfo[, colnames(curdruginfo)]
colnames(druginfo); dim(druginfo) # 25374   8

stopifnot(identical(rownames(drugact), rownames(druginfo))) # TRUE
stopifnot(identical(colnames(drugact), celllines)) # TRUE compare to current cell lines


drugRepact = read.delim("inst/extdata/public_drug_Rep_act_2_14.txt", row.names = 1, stringsAsFactors = F,check.names = F)
drugRepinfo = read.delim("inst/extdata/public_drug_Rep_info_2_14.txt", row.names = 1, stringsAsFactors = F,check.names = F)
dim(drugRepact); dim(drugRepinfo)
# [1] 39137    60
# [1] 39137   4
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

## END
