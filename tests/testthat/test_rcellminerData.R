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

test_that("feature annotation data is matched to numeric data matrix rows", {
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

test_that("numeric data is properly loaded into data matrices", {
  molDB_dat <- rcellminer::getAllFeatureData(rcellminerData::molData)
  drugAct <- exprs(rcellminer::getAct(rcellminerData::drugData))
  drugRepAct <- exprs(rcellminer::getRepeatAct(rcellminerData::drugData))

  # molDB_dat[["exp"]]-------------------------------------------------------------------
  expect_identical(rownames(molDB_dat[["exp"]])[1], "WASH7P")
  expect_identical(rownames(molDB_dat[["exp"]])[nrow(molDB_dat[["exp"]])], "HLA-DRB3")
  expect_identical(molDB_dat[["exp"]]["WASH7P", "BR:MCF7"], 0.499)
  expect_identical(molDB_dat[["exp"]]["WASH7P", "RE:UO-31"], -0.315)
  expect_identical(molDB_dat[["exp"]]["HLA-DRB3", "BR:MCF7"], -0.761)
  expect_identical(molDB_dat[["exp"]]["HLA-DRB3", "RE:UO-31"], 0.136)

  # molDB_dat[["xai"]]-------------------------------------------------------------------
  expect_identical(rownames(molDB_dat[["xai"]])[1], "LOC729737")
  expect_identical(rownames(molDB_dat[["xai"]])[nrow(molDB_dat[["xai"]])], "HLA-DRB3")
  expect_identical(molDB_dat[["xai"]]["LOC729737", "BR:MCF7"], 9.412)
  expect_identical(molDB_dat[["xai"]]["LOC729737", "RE:UO-31"], 8.594)
  expect_identical(molDB_dat[["xai"]]["HLA-DRB3", "BR:MCF7"], 3.913)
  expect_identical(molDB_dat[["xai"]]["HLA-DRB3", "RE:UO-31"], 4.93)

  # molDB_dat[["cop"]]-------------------------------------------------------------------
  expect_identical(rownames(molDB_dat[["cop"]])[1], "C1orf222")
  expect_identical(rownames(molDB_dat[["cop"]])[nrow(molDB_dat[["cop"]])], "GOLGA2P3Y")
  expect_identical(molDB_dat[["cop"]]["C1orf222", "BR:MCF7"], -0.229)
  expect_identical(molDB_dat[["cop"]]["C1orf222", "RE:UO-31"], -0.016)
  expect_identical(molDB_dat[["cop"]]["GOLGA2P3Y", "BR:MCF7"], -0.483)
  expect_identical(molDB_dat[["cop"]]["GOLGA2P3Y", "RE:UO-31"], -0.426)

  # molDB_dat[["met"]]-------------------------------------------------------------------
  expect_identical(rownames(molDB_dat[["met"]])[1], "C1orf222")
  expect_identical(rownames(molDB_dat[["met"]])[nrow(molDB_dat[["met"]])], "RBMY2FP")
  expect_identical(molDB_dat[["met"]]["C1orf222", "BR:MCF7"], 0.851)
  expect_identical(molDB_dat[["met"]]["C1orf222", "RE:UO-31"], 0.077)
  expect_identical(molDB_dat[["met"]]["EIF1AY", "BR:MCF7"], 0.362)
  expect_identical(molDB_dat[["met"]]["EIF1AY", "RE:UO-31"], 0.273)

  # molDB_dat[["mir"]]-------------------------------------------------------------------
  expect_identical(rownames(molDB_dat[["mir"]])[1], "hsa-miR-200b")
  expect_identical(rownames(molDB_dat[["mir"]])[nrow(molDB_dat[["mir"]])], "hsa-miR-105(2)")
  expect_identical(molDB_dat[["mir"]]["hsa-miR-200b", "BR:MCF7"], 3.531)
  expect_identical(molDB_dat[["mir"]]["hsa-miR-200b", "RE:UO-31"], 3.11)
  expect_identical(molDB_dat[["mir"]]["hsa-miR-105(2)", "BR:MCF7"], -4.002)
  expect_identical(molDB_dat[["mir"]]["hsa-miR-105(2)", "RE:UO-31"], -3.998)

  # molDB_dat[["mut"]]-------------------------------------------------------------------
  expect_identical(rownames(molDB_dat[["mut"]])[1], "C1orf222")
  expect_identical(rownames(molDB_dat[["mut"]])[nrow(molDB_dat[["mut"]])], "DDX3Y")
  expect_identical(molDB_dat[["mut"]]["SMCP", "BR:MCF7"], 39)
  expect_identical(molDB_dat[["mut"]]["ELF3", "RE:UO-31"], 43)

  # molDB_dat[["exo"]]-------------------------------------------------------------------
  expect_identical(rownames(molDB_dat[["exo"]])[1], "chr1:865628_G_A")
  expect_identical(rownames(molDB_dat[["exo"]])[nrow(molDB_dat[["exo"]])], "chrY:22942897_T_C")
  expect_identical(molDB_dat[["exo"]]["chrY:22942897_T_C", "LE:SR"], 100)
  expect_identical(molDB_dat[["exo"]]["chr1:153954646_C_T", "BR:MCF7"], 52.941)

  # molDB_dat[["pro"]]-------------------------------------------------------------------
  expect_identical(rownames(molDB_dat[["pro"]])[1], "JAK1_22")
  expect_identical(rownames(molDB_dat[["pro"]])[nrow(molDB_dat[["pro"]])], "MSN_9")
  expect_identical(molDB_dat[["pro"]]["JAK1_22", "BR:MCF7"], 0.809)
  expect_identical(molDB_dat[["pro"]]["JAK1_22", "RE:UO-31"], 0.877)
  expect_identical(molDB_dat[["pro"]]["MSN_9", "BR:MCF7"], -1.673)
  expect_identical(molDB_dat[["pro"]]["MSN_9", "RE:UO-31"], 2.097)

  # molDB_dat[["swa"]]-------------------------------------------------------------------

  # molDB_dat[["mda"]]-------------------------------------------------------------------

  # drugAct -----------------------------------------------------------------------------
  expect_identical(rownames(drugAct)[1], "1")

  # No longer applicable with corrected NSC set.
  #expect_identical(rownames(drugAct)[nrow(drugAct)], "782785")
  expect_identical(drugAct["1", "BR:MCF7"], -0.27)
  expect_identical(drugAct["1", "RE:UO-31"], 1.06)
  #expect_identical(drugAct["782785", "BR:MCF7"], 1.15)
  #expect_identical(drugAct["782785", "RE:UO-31"], 0.03)

  # drugRepAct --------------------------------------------------------------------------
  expect_identical(rownames(drugRepAct)[1], "1_0809RS22")

  # No longer applicable with corrected NSC set.
  #expect_identical(rownames(drugRepAct)[nrow(drugRepAct)], "782763_1503RS50")
  expect_identical(drugRepAct["1_0809RS22", "BR:MCF7"], 4.76)
  expect_identical(drugRepAct["1_0809RS22", "RE:UO-31"], 4.96)
  #expect_identical(drugRepAct["782763_1503RS50", "BR:MCF7"], 7.55)
  #expect_identical(drugRepAct["782763_1503RS50", "RE:UO-31"], 8)
})
