# adding 60 drugs that became public imported from cellminercdb_private-----------------------------------
setwd("./inst/extdata")

## activity 60 drugs ---------------------
#
#--------------------------------------------------------------------------------------------------
# MAKE AND SAVE DrugData OBJECT
#--------------------------------------------------------------------------------------------------


library(rcellminer)

# current ones
nci60Miame <- rcellminerData::drugData@sampleData

actnow <- exprs(rcellminerData::drugData@act)
featactnow <- rcellminerData::drugData@act@featureData@data
dim(actnow); dim(featactnow) # 21183 - 60 -8

repeatactnow <- exprs(rcellminerData::drugData@repeatAct)
featrepeatactnow <- rcellminerData::drugData@repeatAct@featureData@data
dim(repeatactnow); dim(featrepeatactnow) # 43126 - 60 -4
## --------- new staff -------------
## remove drug "762064"
#
indexdrug=which(rownames(actnow)=="762064")
# check same index for metadata
stopifnot(identical(rownames(actnow)[indexdrug],rownames(featactnow)[indexdrug]))
#OK
nsc=apply(array(rownames(repeatactnow)),1,function(x) unlist(strsplit(x,"_"))[1])
# experiment=apply(array(rownames(repeatactnow)),1,function(x) unlist(strsplit(x,"_"))[2])
indexrep=which(nsc=="762064")
# check same index for metadata

stopifnot(identical(rownames(repeatactnow)[indexrep],rownames(featrepeatactnow)[indexrep]))
# OK
#
actMix=actnow[-indexdrug,]
annotMix=featactnow[-indexdrug,]
dim(actMix); dim(annotMix) # 21182 - 60 -8
stopifnot(identical(rownames(actMix),rownames(annotMix)))

actData <- ExpressionSet(actMix)
featureData(actData) <- new("AnnotatedDataFrame", data=annotMix)

RactMix=repeatactnow[-indexrep,]
RannotMix=featrepeatactnow[-indexrep,]
dim(RactMix); dim(RannotMix) #43124 - 60 4
stopifnot(identical(rownames(RactMix),rownames(RannotMix)))

RactData <- ExpressionSet(RactMix)
featureData(RactData) <- new("AnnotatedDataFrame", data=RannotMix)

# finalize !
drugData <- new("DrugData", act = actData, repeatAct = RactData, sampleData = nci60Miame)

save(drugData, file = "../../data/drugData.RData")

# now update package rcellminerUtilsCDB
