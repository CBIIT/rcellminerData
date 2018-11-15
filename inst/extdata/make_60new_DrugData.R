# adding 60 drugs that became public imported from cellminercdb_private-----------------------------------
setwd("./inst/extdata")

## activity 60 drugs ---------------------
vtitle.act=read.delim("act60.txt",header=F,stringsAsFactors=F)[1,-1]
activity=read.delim("act60.txt",stringsAsFactors=F,row.names=1)
colnames(activity)=as.character(vtitle.act[1,])
# convert to numeric matrix
mydata=data.matrix(activity)

## experiments for 60 drugs -----------------------------
vtitle.rep=read.delim("rep.act60.txt",header=F,stringsAsFactors=F)[1,-1]
repact=read.delim("rep.act60.txt",stringsAsFactors=F,row.names=1)
colnames(repact)=as.character(vtitle.rep[1,])
mydata.rep=data.matrix(repact)

## features for 60 drugs -----------------------------
metadata = read.delim("act60meta_complete.txt",stringsAsFactors=F,row.names = 1) # for activity
# imported from database

metadata.rep = read.delim("rep.act60meta.txt",stringsAsFactors=F,row.names=1) # for rep activity

dim(mydata); dim(mydata.rep); dim(metadata); dim(metadata.rep)

stopifnot(identical(colnames(mydata), colnames(mydata.rep)))

stopifnot(identical(rownames(mydata), rownames(metadata)))
stopifnot(identical(rownames(mydata.rep), metadata.rep[,1]))
#
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
dim(actnow); dim(featactnow) # 21123 - 60 -8

repeatactnow <- exprs(rcellminerData::drugData@repeatAct)
featrepeatactnow <- rcellminerData::drugData@repeatAct@featureData@data
dim(repeatactnow); dim(featrepeatactnow) # 43012 - 60 -4
## --------- new staff -------------
metadata$pubchem_id = NA
colnames(metadata)=colnames(featactnow)
# ok now make metadat for rep activity

nsc=apply(array(metadata.rep[,1]),1,function(x) unlist(strsplit(x,"_"))[1])
experiment=apply(array(metadata.rep[,1]),1,function(x) unlist(strsplit(x,"_"))[2])
metadata.rep$nsc=nsc
metadata.rep$experiment=experiment
metadata.rep$zs=TRUE # all used experiments
rownames(metadata.rep)=metadata.rep[,1]

colnames(metadata.rep)=colnames(featrepeatactnow)
#ok
# last check new 60 are not in current public
length(intersect(rownames(actnow),rownames(mydata)))
# zero OK!
#
actMix=rbind(actnow,mydata)
annotMix=rbind(featactnow,metadata)
dim(actMix); dim(annotMix) # 21183 - 60 -8

actData <- ExpressionSet(actMix)
featureData(actData) <- new("AnnotatedDataFrame", data=annotMix)

RactMix=rbind(repeatactnow,mydata.rep)
RannotMix=rbind(featrepeatactnow,metadata.rep)
dim(RactMix); dim(RannotMix) #43126 - 60 4

RactData <- ExpressionSet(RactMix)
featureData(RactData) <- new("AnnotatedDataFrame", data=RannotMix)


drugData <- new("DrugData", act = actData, repeatAct = RactData, sampleData = nci60Miame)

save(drugData, file = "../../data/drugData.RData")

# now update package rcellminerUtilsCDB
