# --------------------------------------------------------------------------------------
# removing the 1030 drugs data taken from private database that became public
# Taking fom the internal db v2.4.2 the 1030 drugs + new 45 ones
# and later update MDA description as well as SWATH data
# --------------------------------------------------------------------------------------
library(rcellminer)
library(stringr)
tmpEnv <- new.env()
##
load("data/molData.RData", envir = tmpEnv)
load("data/drugData.RData", envir = tmpEnv)

nci60Miame <- tmpEnv$molData@sampleData
nci60ESetList=tmpEnv$molData@eSetList

cellname= nci60Miame@samples$Name
length(cellname) # 60

drugact = exprs(tmpEnv$drugData@act)
druginfo = tmpEnv$drugData@act@featureData@data

drugRepact = exprs(tmpEnv$drugData@repeatAct)
drugRepinfo = tmpEnv$drugData@repeatAct@featureData@data

dim(drugact); dim(druginfo)        # 22212    60 ; 22212    8
dim(drugRepact); dim(drugRepinfo)  # 44727    60 ; 44727    4
### ------------------------------------------
## read new data
## -------------------------------------------

list1030 = as.character(read.csv("./inst/extdata/list1030drugs.csv")[,2])
act1075  = read.csv("./inst/extdata/drug_activity_zscore_for_1075drugs.csv",check.names = F,stringsAsFactors = F)
info1075 = read.csv("./inst/extdata/drug_1075_public_v2_4_2.csv",check.names = F,stringsAsFactors = F)
repact1075 = read.csv("./inst/extdata/drug_goodexperiments_for_1075drugs.csv",check.names = F,stringsAsFactors = F)

length(list1030)             # 1030
dim(act1075); dim(info1075)  # 1075 63 ; 1075 10
dim(repact1075)              # 1690 63

## --------------------------------------------
## Remove 1030 drug act, 1030 rep act and all BAD experiments
## --------------------------------------------
# .. start with removing current  1030 drugs .................................
i = match (list1030, rownames(drugact))
dim(drugact[-i,])
j = which(is.na(druginfo$TOTAL_EXPS))
stopifnot(identical(i,j))

drugact = drugact[-i,]
druginfo = druginfo[-i,]
dim(drugact); dim(druginfo) # 21182    60 ; 21182    8
stopifnot(identical(rownames(drugact), rownames(druginfo))) # OK

## .. now remove rep act for 1030 drugs
stopifnot(identical(rownames(drugRepact),rownames(drugRepinfo)))
r = which(is.na(drugRepinfo$used_in_zscore))
length(r) # 1603
nsc1030 = unique(drugRepinfo$nsc[r])
# length(intersect(list1030,nsc1030)) # 1030
stopifnot(identical(list1030,nsc1030)) ## OK

drugRepact = drugRepact[-r,]
drugRepinfo = drugRepinfo[-r,]
dim(drugRepact); dim(drugRepinfo) # 43124    60 ; 43124    4
stopifnot(identical(rownames(drugRepact), rownames(drugRepinfo))) # OK

# ..Removing the BAD experiments  ...............................................
b = which(drugRepinfo$used_in_zscore==FALSE)
length(b) # 10191
drugRepact = drugRepact[-b,]
drugRepinfo = drugRepinfo[-b,]
dim(drugRepact); dim(drugRepinfo) # 32933  60 ; 32933    4
stopifnot(identical(rownames(drugRepact), rownames(drugRepinfo))) # OK

# .. Adding now the new 1075 drugs ..............................................

stopifnot(identical(act1075$nsc,info1075$nsc)) # FALSE
info1075 = info1075[order(info1075$nsc),]
stopifnot(identical(act1075$nsc,info1075$nsc)) # TRUE
info1075$pubchem_id[which(info1075$pubchem_id==-1)]=NA

# drug activity
drugdat1075 = act1075[,4:63]
rownames(drugdat1075) = act1075$nsc
dim(drugdat1075) #  1075   60
length(intersect(rownames(drugact),rownames(drugdat1075)))
# zero OK!
stopifnot(identical(colnames(drugact),colnames(drugdat1075))) # FALSE
drugdat1075 = drugdat1075[,cellname]
stopifnot(identical(colnames(drugact),colnames(drugdat1075))) # TRUE !!
# rounding
drugdat1075 = round(drugdat1075,2)
#
actMix=rbind(drugact,drugdat1075)
dim(actMix) # 22257 - 60
#
# Info 1075 now

druginfo1075 = info1075[,c("nsc","name","testing_status","mechanism","pubchem_id","smiles","total_probes","total_good_probes")]
rownames(druginfo1075) = druginfo1075$nsc
colnames(druginfo1075) = colnames(druginfo)
annotMix=rbind(druginfo,druginfo1075)
dim(annotMix) # 22257 - 8

stopifnot(identical(rownames(actMix),rownames(annotMix))) # TRUE !!

## change drug name for nsc 38721 from Dacarbazine to Mitotane ...............
annotMix[which(annotMix$NSC=="38721"),2] = "Mitotane"
## ...........................................................................

actData <- ExpressionSet(as.matrix(actMix))
featureData(actData) <- new("AnnotatedDataFrame", data=annotMix)

### NOW add new experiments for 1075 drugs

dim(repact1075) # 1690 - 63
drugRepact1075 = repact1075[,4:63]
rownames(drugRepact1075) = paste0(repact1075$nsc,"_",repact1075$name)
length(intersect(rownames(drugRepact1075),rownames(drugRepact))) ## OK zero !

drugRepact1075 = drugRepact1075[,cellname]
stopifnot(identical(colnames(drugRepact),colnames(drugRepact1075))) # TRUE !!
# rounding
drugRepact1075 = round(drugRepact1075,2)
#
# colnames(drugRepinfo)
drugRepinfo1075 =data.frame(cbind(NSC_EXP_NAME=paste0(repact1075$nsc,"_",repact1075$name),nsc=repact1075$nsc,experiment=repact1075$name),stringsAsFactors = F)
drugRepinfo1075$used_in_zscore = TRUE
rownames(drugRepinfo1075) = drugRepinfo1075[,1]

stopifnot(identical(rownames(drugRepinfo1075),rownames(drugRepact1075)))
stopifnot(identical(colnames(drugRepinfo1075),colnames(drugRepinfo)))

# ......

RactMix=rbind(drugRepact,drugRepact1075)
RannotMix=rbind(drugRepinfo,drugRepinfo1075)
dim(RactMix); dim(RannotMix) # 34623    60 ; 34623    4
stopifnot(identical(rownames(RactMix),rownames(RannotMix)))


RactData <- ExpressionSet(as.matrix(RactMix))
featureData(RactData) <- new("AnnotatedDataFrame", data=RannotMix)


drugData <- new("DrugData", act = actData, repeatAct = RactData, sampleData = nci60Miame)

save(drugData, file = "data/drugData.RData")

# now update package rcellminerUtilsCDB for new drug matching ??
