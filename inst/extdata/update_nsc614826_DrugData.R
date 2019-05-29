setwd("./inst/extdata")

library(rcellminer)

# current ones
nci60Miame <- rcellminerData::drugData@sampleData

actnow <- exprs(rcellminerData::drugData@act)
featactnow <- rcellminerData::drugData@act@featureData@data
dim(actnow); dim(featactnow) # 21182 - 60 -8

repeatactnow <- exprs(rcellminerData::drugData@repeatAct)
featrepeatactnow <- rcellminerData::drugData@repeatAct@featureData@data
dim(repeatactnow); dim(featrepeatactnow) # 43124 - 60 -4
## --------- new staff -------------
## update drug data for nsc  "614826"
#
indexdrug=which(rownames(actnow)=="614826")
# check same index for metadata
stopifnot(identical(rownames(actnow)[indexdrug],rownames(featactnow)[indexdrug]))
#OK
featactnow[indexdrug,"NAME"]; featactnow[indexdrug,"PUBCHEM_ID"]
#  "Acetalax" 99355751
featactnow[indexdrug,"NAME"] = "Bisacodyl"
featactnow[indexdrug,"PUBCHEM_ID"] = 487150
# OK

actData <- ExpressionSet(actnow)
featureData(actData) <- new("AnnotatedDataFrame", data=featactnow)


RactData <- ExpressionSet(repeatactnow)
featureData(RactData) <- new("AnnotatedDataFrame", data=featrepeatactnow)

# finalize !
drugData <- new("DrugData", act = actData, repeatAct = RactData, sampleData = nci60Miame)

save(drugData, file = "../../data/drugData.RData")

# now update package rcellminerUtilsCDB
