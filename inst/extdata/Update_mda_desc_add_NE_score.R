# ----------------------------------------------------------------------------------------
# update MDA description and add NE (neuroendocrine) score as new metadata
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

## update metadata with new ne score

mda.data = exprs(nci60ESetList$mda)
mda.info = nci60ESetList$mda@featureData@data

dim(mda.data); dim(mda.info)
# [1]   40 60
# [1]    40 4
ne = read.delim("./inst/extdata/NCI60-NEscores-49genes-60samples-04142020.txt",stringsAsFactors = F,row.names=1)
dim(ne) #60    3

nbcol = ncol(mda.data)
ne_scores = replicate(nbcol,NA)
names(ne_scores)=cellname

# checking
length(intersect(rownames(ne),cellname)) # 60

ne_scores[rownames(ne)] = ne$score
length(which(is.na(ne_scores))) # 0
stopifnot(identical(names(ne_scores),colnames(mda.data)))

mda.data = rbind(mda.data,ne_scores)
dim(mda.data) # 41 60
nb = nrow(mda.data)
rownames(mda.data)[nb]= "NE"

mda.info = rbind(mda.info,c("NE","Neuro-endocrine transcript expression signature score, based on the identification of the 25 genes most associated and 25 genes least associated with Neuro-Endocrine (NE) status as determined by inverted microscope examination of living cells and cytological morphology","","Gazdar et al, Translational Lung Cancer Research, 2018. PMID:29535911"))
rownames(mda.info)[nb]= "NE"

## update APM and rep Stress descriptions
mda.info["APM",2] = "Antigen Presentation Machinery (APM) transcript expression signature score based on a Gene Set Variation analysis (GSVA) enrichment score of 18 genes"
mda.info["APM",4] = "Wang et al, eLIFE, 2019. PMID: 31767055"

mda.info["RepStress",2] = "Replication Stress transcript expression signature score based on first principal component weights of 18 genes"
mda.info["RepStress",4] = "Signature under development by Drs Anish Thomas and Vinodh Rajapakse"

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
