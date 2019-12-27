# --------------------------------------------------------------------------------------
# adding 1030 drugs that became public
# update after checking (by Jeffrey) the DTP website for public drug
# and if some exist on the internal rcellminerDataInt then move them to public one
# --------------------------------------------------------------------------------------
setwd("./inst/extdata")

nscs = as.character(read.delim("publicNSC_Oct2019.txt",check.names = F, stringsAsFactors = F,header = F)[,1])
length(nscs) # 53279
publics = rcellminerData::drugData@act@featureData@data
privates = rcellminerDataInt::drugData@act@featureData@data

dim(publics); dim(privates)
# 21182     8
# 59224     5

# checking public ones
length(intersect(publics$NSC,nscs)) #  21182 versus website

length(intersect(publics$NSC,privates$NSC)) #  21182 versus private website

# new ones ??
length(intersect(privates$NSC,nscs)) # 22212 public with activities stored
temp = intersect(privates$NSC,nscs)
newones= setdiff(temp,publics$NSC)
length(newones) # 1030 from private to public with activities stored

##   new to explore later from DTP for activity
# temp2 = setdiff(nscs,privates$NSC)
# length(temp2) # 31067 without any activity in our internal package
##              #  31067  did not pass QC ??
# -----------------------------------------------------------------------------
#  move these 1030 privates drugs to the public domain: activity + experiments
# -----------------------------------------------------------------------------
library(rcellminer)
vindex= match(newones,privates$NSC)

drugdat = exprs(rcellminerDataInt::drugData@act)
stopifnot(identical(rownames(drugdat),privates$NSC))

drugdat1030 = drugdat[vindex,]
dim(drugdat1030) # 1030 - 60
stopifnot(identical(rownames(drugdat1030),newones))

druginfo1030 = privates[vindex,]
dim(druginfo1030) # 1030 - 5
stopifnot(identical(rownames(druginfo1030),newones))
# now experiments or repAct ------
expdat = exprs(rcellminerDataInt::drugData@repeatAct)
expinfo = rcellminerDataInt::drugData@repeatAct@featureData@data
dim(expdat); dim(expinfo)
# [1] 80739    60
# [1] 80739     1

stopifnot(identical(rownames(expdat),expinfo[,1]))


lnsc =  apply(as.array(rownames(expdat)),1,function(x) unlist(strsplit(x,"_"))[1])
length(lnsc) # 80739

vindex2 = which(lnsc %in% newones)
length(vindex2) # 1603


length(unique(lnsc[vindex2])) # 1030
length(intersect(newones,unique(lnsc[vindex2]))) # 1030

expdat1603 = expdat[vindex2,]
expinfo1603 = expinfo[vindex2,1,drop=F]
dim(expdat1603); dim(expinfo1603)
# [1] 1603   60
# [1] 1603    1

stopifnot(identical(rownames(expdat1603),expinfo1603[,1]))

#--------------------------------------------------------------------------------------------------
# MAKE AND SAVE DrugData OBJECT
#--------------------------------------------------------------------------------------------------
stopifnot(identical(colnames(drugdat1030), colnames(exprs(rcellminerData::drugData@act))))
# current ones
nci60Miame <- rcellminerData::drugData@sampleData
actnow <- exprs(rcellminerData::drugData@act)
featactnow <- rcellminerData::drugData@act@featureData@data
dim(actnow); dim(featactnow) # 21182 - 60 -8

repeatactnow <- exprs(rcellminerData::drugData@repeatAct)
featrepeatactnow <- rcellminerData::drugData@repeatAct@featureData@data
dim(repeatactnow); dim(featrepeatactnow) # 43124 - 60 -4
## --------- new staff -------------
dim(drugdat1030); dim(druginfo1030) # 1030 - 60 - 5
colnames(druginfo1030)

druginfo1030[,c("FDA_STATUS","PUBCHEM_ID","TOTAL_EXPS","TOTAL_EXPS_AFTER_QC")] = NA ## !!!!!
druginfo1030$CONFIDENTIAL_FLAG = NULL

druginfo1030 =druginfo1030[,colnames(featactnow)]
dim(druginfo1030) # 1030 - 8

# ok now make metadat for rep activity

dim(expdat1603); dim(expinfo1603) # 1603 - 60 -1
dim(featrepeatactnow)  # 43124     4



nsc=apply(array(expinfo1603[,1]),1,function(x) unlist(strsplit(x,"_"))[1])
experiment=apply(array(expinfo1603[,1]),1,function(x) unlist(strsplit(x,"_"))[2])
expinfo1603$nsc=nsc
expinfo1603$experiment=experiment
expinfo1603$used_in_zscore = NA # all used experiments !!!!!
rownames(expinfo1603)=expinfo1603[,1]

colnames(expinfo1603)[1]="NSC_EXP_NAME"
#ok
# last check new 60 are not in current public
length(intersect(rownames(actnow),rownames(drugdat1030)))
# zero OK!
#
actMix=rbind(actnow,drugdat1030)
annotMix=rbind(featactnow,druginfo1030)
dim(actMix); dim(annotMix) # 22212 - 60 -8

actData <- ExpressionSet(actMix)
featureData(actData) <- new("AnnotatedDataFrame", data=annotMix)
###

RactMix=rbind(repeatactnow,expdat1603)
RannotMix=rbind(featrepeatactnow,expinfo1603)
dim(RactMix); dim(RannotMix) #44727 - 60 4

RactData <- ExpressionSet(RactMix)
featureData(RactData) <- new("AnnotatedDataFrame", data=RannotMix)


drugData <- new("DrugData", act = actData, repeatAct = RactData, sampleData = nci60Miame)

save(drugData, file = "../../data/drugData.RData")

# now update package rcellminerUtilsCDB for new drug matching ??
