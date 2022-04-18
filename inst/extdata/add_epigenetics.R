# add epigenetics data  H3K27ac and  H3K4me3. Data received from Lorinc. in log2(tmm-fpkm)
#--------------------------------------------------------------------------------------------------
library(rcellminer)
library(stringr)
##
tmpEnv <- new.env()
#
load("data/molData.RData", envir = tmpEnv)
# load("data/drugData.RData", envir = tmpEnv)

nci60SetList=tmpEnv$molData@eSetList


nci60Miame <- tmpEnv$molData@sampleData
cells <- nci60Miame@samples$Name

## read new data
##  histone data --------------------------------------------------------------------------------------

h3k27 = read.delim("./inst/extdata/H3K27ac.gene_promoter_signal.tsv",stringsAsFactors = F, check.names = F, row.names = 1)
h3k4m = read.delim("./inst/extdata/H3K4me3.gene_promoter_signal.tsv",stringsAsFactors = F, check.names = F, row.names = 1)

dim(h3k27) # 22073    60
dim(h3k4m) # 19615    60

stopifnot(identical(colnames(h3k27), colnames(h3k4m))) ## T

cnames = colnames(h3k27)
cnames= gsub("_",":",cnames)
### ---------------------------------

## add it to package rcellminerData

length(intersect(cnames, cells)) # 52

  cnames[which(cnames=="BR:HS:578T")] = "BR:HS 578T"
  cnames[which(cnames=="CO:COLO:205")] = "CO:COLO 205"
  cnames[which(cnames=="LE:HL-60")] = "LE:HL-60(TB)"
  cnames[which(cnames=="ME:LOX:IMVI")] = "ME:LOX IMVI"
  cnames[which(cnames=="LC:A549")] = "LC:A549/ATCC"
  cnames[which(cnames=="OV:NCI:ADR-RES")] = "OV:NCI/ADR-RES"
  cnames[which(cnames=="RE:RXF:393")] = "RE:RXF 393"

  # not found "BR:MDA-MB-468" in NCI60
  # not found "ME:MDA-N" in epigenetics data
  length(intersect(cnames, cells)) # 59

setdiff(cells,cnames)
 "ME:MDA-N"

setdiff(cnames,cells)
# "BR:MDA-MB-468"

colnames(h3k27) = cnames
colnames(h3k4m) = cnames

h3k27$"ME:MDA-N" = NA
h3k4m$"ME:MDA-N" = NA

dim(h3k27) # 22073    61
dim(h3k4m) # 19615    61

length(intersect(colnames(h3k27), cells)) # 60
length(intersect(colnames(h3k4m), cells)) # 60

h3k27 = h3k27[,cells]
h3k4m = h3k4m[,cells]
dim(h3k27) # 22073    60
dim(h3k4m) # 19615    60

stopifnot(identical(colnames(h3k27), colnames(h3k4m))) ## T
stopifnot(identical(colnames(h3k27), cells)) ## T

# -----------------------------------------------------------------
# add data to molecular list
# -----------------------------------------------------------------

# next gene annotation >> double column ID gene

h3k27.info = read.delim("./inst/extdata/H3K27ac.gene_to_peak_assignments.tsv",stringsAsFactors = F, check.names = F)
h3k4m.info = read.delim("./inst/extdata/H3K4me3.gene_to_peak_assignments.tsv",stringsAsFactors = F, check.names = F)
dim(h3k27.info); dim(h3k4m.info)
# 22073     4 ; 19615     4
stopifnot(identical(rownames(h3k27), h3k27.info$HUGO_symbol )) ## T
stopifnot(identical(rownames(h3k4m), h3k4m.info$HUGO_symbol)) ## T

colnames(h3k27.info)[1] = "ID"
colnames(h3k4m.info)[1] = "ID"

rownames(h3k27.info) = h3k27.info$ID
rownames(h3k4m.info) = h3k4m.info$ID

stopifnot(is.numeric(as.matrix(h3k27)))
stopifnot(is.numeric(as.matrix(h3k4m)))


h3k27Data <- ExpressionSet(as.matrix(h3k27))
featureData(h3k27Data) <- new("AnnotatedDataFrame", data=h3k27.info)

h3k4mData <- ExpressionSet(as.matrix(h3k4m))
featureData(h3k4mData) <- new("AnnotatedDataFrame", data=h3k4m.info)


nci60SetList[["his"]] <- h3k27Data
nci60SetList[["hs4"]] <- h3k4mData

#--------------------------------------------------------------------------------------------------
# PREPARE METADATA and DATA
#--------------------------------------------------------------------------------------------------

molData <- new("MolData", eSetList = nci60SetList, sampleData = nci60Miame)

save(molData, file = "data/molData.RData")

## END
