#--------------------------------------------------------------------------------------------------
library(rcellminer)
library(stringr)
## process rnaseq data
load("./inst/extdata/RNAseq.geneMay14.Rdata")

dim(gene.annot)
#[1] 23826  5
rownames(gene.annot)=gene.annot$gene
colnames(gene.annot)[1]="ID"

rnaseq=gene.data
rownames(rnaseq)=gene.annot$ID
dim(rnaseq) # 23826 x 60

rnaseq=log2(rnaseq+1)
range(rnaseq) #  0.00000 24.45673

### ---------------------------------

## add it to package rcellminerData
tmpEnv <- new.env()
#
load("data/molData.RData", envir = tmpEnv)
# load("data/drugData.RData", envir = tmpEnv)

nciSclcMiame <- tmpEnv$molData@sampleData

# To check this MIAME below !!!

mdaExp <- rcellminer::getAllFeatureData(tmpEnv$molData)[["mda"]]
cells <- colnames(mdaExp)

# checking columns
length(which(cells==colnames(rnaseq))) ## 60

#--------------------------------------------------------------------------------------------------
# PREPARE METADATA and DATA
#--------------------------------------------------------------------------------------------------

##
xsqData <- ExpressionSet(rnaseq)
stopifnot(is.numeric(exprs(xsqData)))


featureData(xsqData) <- new("AnnotatedDataFrame", data=gene.annot)

# Column (NCI-60 cell line) consistency check.
stopifnot(identical(colnames(exprs(xsqData)), cells))
stopifnot(identical(rownames(exprs(xsqData)), rownames(gene.annot)))

# saving now
nciSetList=tmpEnv$molData@eSetList

nciSetList[["xsq"]] <- xsqData

molData <- new("MolData", eSetList = nciSetList, sampleData = nciSclcMiame)

save(molData, file = "data/molData.RData")

## END
