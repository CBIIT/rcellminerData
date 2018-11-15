# Read BMN673 data -----------------------------------
setwd("./inst/extdata")

## activity bmn-673 ---------------------
vtitle.act=read.delim("bmn-673.txt",header=F,stringsAsFactors=F)[1,-1]
activity=read.delim("bmn-673.txt",stringsAsFactors=F,row.names=1)
colnames(activity)=as.character(vtitle.act[1,])
# convert to numeric matrix
mydata=data.matrix(activity)
## or
## m=apply(activity,2,as.numeric)

## experiments for bmn-673 -----------------------------
vtitle.rep=read.delim("bmn-673-repAct.txt",header=F,stringsAsFactors=F)[1,-1]
repact=read.delim("bmn-673-repAct.txt",stringsAsFactors=F,row.names=1)
colnames(repact)=as.character(vtitle.rep[1,])
mydata.rep=data.matrix(repact)

## features for bmn-673 -----------------------------
metadata = read.delim("bmn-673-meta-ready.txt",stringsAsFactors=F) # for activity
rownames(metadata)=metadata[,1]

metadata.rep = read.delim("bmn-673-meta-repAct.txt",stringsAsFactors=F,row.names=1) # for rep activity

dim(mydata); dim(mydata.rep); dim(metadata); dim(metadata.rep)

stopifnot(identical(colnames(mydata), colnames(mydata.rep)))
#--------------------------------------------------------------------------------------------------
# MAKE AND SAVE DrugData OBJECT
#--------------------------------------------------------------------------------------------------


library(rcellminer)

stopifnot(identical(colnames(mydata), colnames(exprs(rcellminerData::drugData@act))))

stopifnot(identical(colnames(mydata), colnames(exprs(getAct(rcellminerData::drugData)))))
# current ones
nci60Miame <- rcellminerData::drugData@sampleData

actnow <- exprs(rcellminerData::drugData@act)
featactnow <- rcellminerData::drugData@act@featureData@data
dim(actnow); dim(featactnow) # 21121 - 60 -8

repeatactnow <- exprs(rcellminerData::drugData@repeatAct)
featrepeatactnow <- rcellminerData::drugData@repeatAct@featureData@data
dim(repeatactnow); dim(featrepeatactnow) # 43007 - 60 -4
## --------- new staff -------------
colnames(metadata)=colnames(featactnow)
colnames(metadata.rep)=colnames(featrepeatactnow)

actMix=rbind(actnow,mydata)
annotMix=rbind(featactnow,metadata)
dim(actMix); dim(annotMix) # 21123 - 60 -8

actData <- ExpressionSet(actMix)
featureData(actData) <- new("AnnotatedDataFrame", data=annotMix)

RactMix=rbind(repeatactnow,mydata.rep)
RannotMix=rbind(featrepeatactnow,metadata.rep)
dim(RactMix); dim(RannotMix) #43012 - 60 4

RactData <- ExpressionSet(RactMix)
featureData(RactData) <- new("AnnotatedDataFrame", data=RannotMix)


drugData <- new("DrugData", act = actData, repeatAct = RactData, sampleData = nci60Miame)

save(drugData, file = "../../data/drugData.RData")

