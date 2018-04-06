# Read Almanac data -----------------------------------
setwd("./inst/extdata")

vtitle=read.delim("DTB_NCI60_cell_lines_Almanac_paired_drug_combo_scores.txt",header=F,stringsAsFactors=F)[1,-1]
almanac=read.delim("DTB_NCI60_cell_lines_Almanac_paired_drug_combo_scores.txt",stringsAsFactors=F,row.names=1)
# be sure to put right cell lines
colnames(almanac)=as.character(vtitle[1,])
# convert to numeric matrix
mydata=data.matrix(almanac)
## or
## m=apply(almanac,2,as.numeric)

# extract unique NSC information from paired drug data
nsc1=apply(array(rownames(mydata)),1,function(x) {unlist(strsplit(x,"_"))[1]})
nsc2=apply(array(rownames(mydata)),1,function(x) {unlist(strsplit(x,"_"))[2]})
nsc=union(nsc1,nsc2); length(nsc) # 105
#--------------------------------------------------------------------------------------------------
# GET DRUG ANNOTATION TABLE.
#--------------------------------------------------------------------------------------------------
if (!require(RMySQL))
  stop("Package RMySQL must be installed")

pswd="****"
con = dbConnect(dbDriver("MySQL"),
                host = "discovery.nci.nih.gov",
                user = "****",
                port = 1111,
                password = pswd)

rs <- dbSendQuery(con, statement = paste("use", "nci60private2_2"))

query2 <- "Select  nsc, name, testing_status,  smiles, confidential_flag, mechanism from drug where prefix='S'"
rs <- dbSendQuery(con, statement = query2)

drugAnnotPublic <- fetch(rs, n = -1)
dbDisconnect(con)
n=dim(drugAnnotPublic)[1]

rownames(drugAnnotPublic)=drugAnnotPublic$nsc

drugAnnotPublic105=drugAnnotPublic[nsc,]; dim(drugAnnotPublic105)
#
### save(drugAnnotPublic105, file = "drugAnnotPublic105.RData")
#

#----------------------------------------------------------------------
# ANNOTATE PAIRED DRUG
#--------------------------------------------------------------------------------------------------

# head(rcellminerData::drugData@act@featureData@data)
vlab=colnames(rcellminerData::drugData@act@featureData@data)
# [1] "NSC"                 "NAME"                "FDA_STATUS"          "MOA"
# [5] "PUBCHEM_ID"          "SMILES"              "TOTAL_EXPS"          "TOTAL_EXPS_AFTER_QC"

## to do: for each paired drug concatenate ONLY name, moa with FDA approved status.
## put rest with NA


NSC=rownames(mydata)


nb=dim(mydata)[1]
vmoa=replicate(nb,NA); vname=replicate(nb,NA)


for (i in 1:nb) {
vmoa[i]=paste0(drugAnnotPublic105[nsc1[i],"mechanism"]," - ",drugAnnotPublic105[nsc2[i],"mechanism"])
vname[i]=paste0(drugAnnotPublic105[nsc1[i],"name"]," - ",drugAnnotPublic105[nsc2[i],"name"])
cat(i, " \n")
}

pairedAnnot=cbind(NSC,NAME = vname,FDA_STATUS = "FDA approved both",MOA = vmoa,PUBCHEM_ID = NA,SMILES = NA, TOTAL_EXPS = NA,TOTAL_EXPS_AFTER_QC = NA)
rownames(pairedAnnot)=NSC
### save(pairedAnnot, file = "pairedAnnot.RData")

#--------------------------------------------------------------------------------------------------
# MAKE AND SAVE DrugData OBJECT
#--------------------------------------------------------------------------------------------------


library(rcellminer)

stopifnot(identical(rownames(mydata), rownames(pairedAnnot)))
stopifnot(identical(colnames(mydata), colnames(exprs(getAct(rcellminerData::drugData)))))
# current ones
nci60Miame <- rcellminerData::drugData@sampleData

actnow <- exprs(rcellminerData::drugData@act)
featactnow <- rcellminerData::drugData@act@featureData@data

repeatactnow <- rcellminerData::drugData@repeatAct
## --------- new stuff -------------
actMix=rbind(actnow,mydata)
annotMix=rbind(featactnow,pairedAnnot)

actData <- ExpressionSet(actMix)
featureData(actData) <- new("AnnotatedDataFrame", data=annotMix)


drugData <- new("DrugData", act = actData, repeatAct = repeatactnow, sampleData = nci60Miame)

save(drugData, file = "../../data/drugData.RData")

