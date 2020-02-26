# ----------------------------------------------------------------------------------------
# add APM score as metadata based on GSVA score
#
#--------------------------------------------------------------------------------------------------
library(rcellminer)
library(stringr)
tmpEnv <- new.env()
##
load("data/molData.RData", envir = tmpEnv)

nci60Miame <- tmpEnv$molData@sampleData
cellname= nci60Miame@samples$Name
nb = length(cellname) # 77
##
nci60ESetList=tmpEnv$molData@eSetList
## update metadata with APM score

mda.data = exprs(nci60ESetList$mda)
mda.info = nci60ESetList$mda@featureData@data

mda.info = data.frame(apply(mda.info,2,as.character),stringsAsFactors = F)

dim(mda.data); dim(mda.info)
# [1]   38 60
# [1]    38 4
apms = read.delim("./inst/extdata/NCI60_scores_SVGA_APM_60samples_9012genes_filteredIQR.txt",stringsAsFactors = F,row.names=1)
dim(apms) # 60   1

apm_scores = replicate(nb,NA)
names(apm_scores)=cellname

# checking
length(intersect(rownames(apms),cellname)) # 60

apm_scores[rownames(apms)] = apms[,1]
length(which(is.na(apm_scores))) # 0
stopifnot(identical(names(apm_scores),colnames(mda.data)))

mda.data = rbind(mda.data,apm_scores)
dim(mda.data) # 39 60
rownames(mda.data)[nrow(mda.data)]= "APM"

mda.info = rbind(mda.info,c("APM","APM Gsva enrichment score","","For details see PMID: 31767055"))
rownames(mda.info)= mda.info$ID
dim(mda.info)
# update mda in elist to do
mdaData <- ExpressionSet(mda.data)
stopifnot(identical(rownames(mda.data), rownames(mda.info)))
featureData(mdaData) <- new("AnnotatedDataFrame", data=mda.info)

nci60ESetList[["mda"]] <- mdaData


## take old miame

molData <- new("MolData", eSetList = nci60ESetList, sampleData = nci60Miame)

save(molData, file = "data/molData.RData")


### STOP Here
### --------------------------------------------------------------
