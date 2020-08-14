# --------------------------------------------------------------------------------------
#  update MDA description as well as SWATH data
# --------------------------------------------------------------------------------------
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

swadata = exprs(nci60ESetList$swa)
sawinfo = nci60ESetList$swa@featureData@data
dim(swadata); dim(sawinfo) # 3142 60; 3142 4
### ----------------------------------------------------------------
## read new SWATH data downloaded from cellminer db internal v.2.4.2
## -----------------------------------------------------------------
newswa = read.delim("./inst/extdata/Protein_SWATH_int2.4.2_05192020.txt",sep="\t",stringsAsFactors = F, check.names = F)
dim(newswa) # 3172 77
ibad = which(newswa$`Entrez gene id e`==0 | is.na(newswa$`Entrez gene id e`))
newswa= newswa[-ibad,1:69]
dim(newswa) # 3167 - 69
#
length(unique(newswa$`Gene name d`)) # 3162
# Read protein access number ........................

uniprot = read.csv("./inst/extdata/SWATH_uniprot_id.csv",stringsAsFactors = F)
dim(uniprot) # 3167 - 3

finswa = merge(uniprot,newswa,by.x=1,by.y=1)
dim(finswa)# 3167 - 71
stopifnot(identical(finswa$symbol, finswa$`Gene name d`))

# duplicated genes
ind= which(duplicated(finswa$symbol))
dn = finswa$symbol[ind]

ind2= which(duplicated(finswa$`Gene name d`))
dn2 = finswa$`Gene name d`[ind2]

ind3 = which(finswa$`Gene name d` %in% dn2)
View(finswa[ind3,])

ID = finswa$symbol
ID[ind3] = paste0(ID[ind3],"_",finswa$uniprot_id[ind3])
rownames(finswa) = ID

newswaData = finswa[,12:71]

dim(newswaData) # 3167   60
newswaData = apply(newswaData, 2, as.numeric) # now a matrix
rownames(newswaData) = ID
stopifnot(identical(cellname,colnames(newswaData)))
#
newswaInfo = cbind(ID=ID,finswa[,c(1,3,4:9,11)])
dim(newswaInfo) # 3167   10
newswaInfo$ID = as.character(newswaInfo$ID)

stopifnot(identical(rownames(newswaData),rownames(newswaInfo))) ##TRUE
colnames(newswaInfo)

# [1] "ID"               "probe_nm"         "uniprot_id"       "Gene name d"      "Entrez gene id e" "Chromosome f"
# [7] "Start f"          "End f"            "Cytoband f"       "M score h"

colnames(newswaInfo)[4:10] = c("Gene_name","Entrez_gene","Chromosome","Start_pos","End_pos","Cytoband","M_score")

# to compare current to new 3142 vs 3171
length(intersect(rownames(newswaData),rownames(swadata))) # 3134 only

# prepare new SWATH
swaData <- ExpressionSet(as.matrix(newswaData))
featureData(swaData) <- new("AnnotatedDataFrame", data=newswaInfo)

nci60ESetList[["swa"]] = swaData


## --------------------------------------------
## update metadata description
## --------------------------------------------
##

mda = exprs(nci60ESetList$mda)
mdainfo = nci60ESetList$mda@featureData@data

dim(mda); dim(mdainfo) # 41 60; 41 4

mymda = c("IS_P53_MUT","EMT","gammaH2AX","Total_H2AX","Ratio_gX_TotX","APM","RepStress","NE")
mypar = c("TP53 mutational status","Epithelial-Mesenchymal signature","DNA damage marker","DNA damage marker","DNA damage marker","Antigen presentation machinery signature","Replication stress signature","Neuro-endocrine signature")
mdainfo[mymda,"Parameter"] =  mypar
mdainfo[mymda[3:5],"Reference"] = c("PMID 28158293","PMID 28158293","PMID 28158293")
View(mdainfo[mymda,])


mdaData <- ExpressionSet(mda)
featureData(mdaData) <- new("AnnotatedDataFrame", data=mdainfo)


# saving now

nci60ESetList[["mda"]] <- mdaData

molData <- new("MolData", eSetList = nci60ESetList, sampleData = nci60Miame)

save(molData, file = "data/molData.RData")


