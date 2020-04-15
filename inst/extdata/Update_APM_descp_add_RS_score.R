# ----------------------------------------------------------------------------------------
# update APM description and add Replication stress score as new metadata
#
#--------------------------------------------------------------------------------------------------
library(rcellminer)
library(stringr)
tmpEnv <- new.env()
##
load("data/molData.RData", envir = tmpEnv)
#load("data/drugData.RData", envir = tmpEnv)

nci60Miame <- tmpEnv$molData@sampleData
nci60ESetList=tmpEnv$molData@eSetList

cellname= nci60Miame@samples$Name
length(cellname) # 60

## update metadata with new RS score

mda.data = exprs(nci60ESetList$mda)
mda.info = nci60ESetList$mda@featureData@data

dim(mda.data); dim(mda.info)
# [1]   39 60
# [1]   39 4
rs = read.delim("./inst/extdata/NCI60_repstress_exp_scores.txt",stringsAsFactors = F,row.names=1)
dim(rs) # 60    1

nbcol = ncol(mda.data)
rs_scores = replicate(nbcol,NA)
names(rs_scores)=cellname

# checking
length(intersect(rownames(rs),cellname)) # 60

rs_scores[rownames(rs)] = rs[,1]
length(which(is.na(rs_scores))) # 0
stopifnot(identical(names(rs_scores),colnames(mda.data)))

mda.data = rbind(mda.data,rs_scores)
dim(mda.data) # 40 60
nb = nrow(mda.data)
rownames(mda.data)[nb]= "RepStress"

mda.info = rbind(mda.info,c("RepStress","Replication Stress signature score","","For details see help section"))
rownames(mda.info)[nb]= "RepStress"
dim(mda.info) # 40 4
## update APM description
mda.info["APM",2] = "Antigen Presentation Machinery (APM) signature enrichment score using GSVA"


# update mda
mdaData <- ExpressionSet(mda.data)
stopifnot(identical(rownames(mda.data), rownames(mda.info)))
featureData(mdaData) <- new("AnnotatedDataFrame", data=mda.info)

nci60ESetList[["mda"]] <- mdaData

##

## take old miame

molData <- new("MolData", eSetList = nci60ESetList, sampleData = nci60Miame)

save(molData, file = "data/molData.RData")



### STOP Here
### --------------------------------------------------------------
