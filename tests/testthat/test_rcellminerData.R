context("rcellminerData validation")

test_that("NCI-60 cell line names are consistent across molecular and drug data", {
  molDB <- rcellminer::getAllFeatureData(rcellminerData::molData)
  drugAct <- exprs(rcellminer::getAct(rcellminerData::drugData))
  drugRepAct <- exprs(rcellminer::getRepeatAct(rcellminerData::drugData))

  # CellMiner 2.0
  nci60Lines <- c("BR:MCF7", "BR:MDA-MB-231", "BR:HS 578T", "BR:BT-549", "BR:T-47D",
                  "CNS:SF-268", "CNS:SF-295", "CNS:SF-539", "CNS:SNB-19", "CNS:SNB-75", "CNS:U251",
                  "CO:COLO 205", "CO:HCC-2998", "CO:HCT-116", "CO:HCT-15", "CO:HT29", "CO:KM12", "CO:SW-620",
                  "LE:CCRF-CEM", "LE:HL-60(TB)", "LE:K-562", "LE:MOLT-4", "LE:RPMI-8226", "LE:SR",
                  "ME:LOX IMVI", "ME:MALME-3M", "ME:M14", "ME:SK-MEL-2", "ME:SK-MEL-28", "ME:SK-MEL-5", "ME:UACC-257",
                    "ME:UACC-62", "ME:MDA-MB-435", "ME:MDA-N", "LC:A549/ATCC", "LC:EKVX",
                  "LC:HOP-62", "LC:HOP-92", "LC:NCI-H226", "LC:NCI-H23", "LC:NCI-H322M", "LC:NCI-H460", "LC:NCI-H522",
                  "OV:IGROV1", "OV:OVCAR-3", "OV:OVCAR-4", "OV:OVCAR-5", "OV:OVCAR-8", "OV:SK-OV-3", "OV:NCI/ADR-RES",
                  "PR:PC-3", "PR:DU-145",
                  "RE:786-0", "RE:A498", "RE:ACHN", "RE:CAKI-1", "RE:RXF 393", "RE:SN12C", "RE:TK-10", "RE:UO-31")

  expect_identical(colnames(drugAct), nci60Lines)
  expect_identical(colnames(drugRepAct), nci60Lines)

  for (molDataType in names(molDB)){
    expect_identical(colnames(molDB[[molDataType]]), nci60Lines)
  }

  cellLineInfo_molData <- rcellminer::getSampleData(rcellminerData::molData)
  cellLineInfo_drugData <- rcellminer::getSampleData(rcellminerData::drugData)
  expect_identical(cellLineInfo_molData, cellLineInfo_drugData)
  expect_identical(cellLineInfo_molData$Name, nci60Lines)
})

test_that("feature annotation data is matched to data matrix rows", {
  molDB_dat <- rcellminer::getAllFeatureData(rcellminerData::molData)
  molDB_annot <- rcellminer::getFeatureAnnot(rcellminerData::molData)

  expect_identical(names(molDB_dat), names(molDB_annot))
  for (molDataType in names(molDB_dat)){
    expect_true(!is.null(rownames(molDB_dat[[molDataType]])))
    expect_identical(rownames(molDB_dat[[molDataType]]),
                     rownames(molDB_annot[[molDataType]]))
  }

  drugAct <- exprs(rcellminer::getAct(rcellminerData::drugData))
  drugAnnot <- rcellminer::getFeatureAnnot(rcellminerData::drugData)[["drug"]]
  expect_true(!is.null(rownames(drugAct)))
  expect_identical(rownames(drugAct), rownames(drugAnnot))

  drugRepAct <- exprs(rcellminer::getRepeatAct(rcellminerData::drugData))
  drugRepAnnot <- rcellminer::getFeatureAnnot(rcellminerData::drugData)[["drugRepeat"]]
  expect_true(!is.null(rownames(drugRepAct)))
  expect_identical(rownames(drugRepAct), rownames(drugRepAnnot))
})

