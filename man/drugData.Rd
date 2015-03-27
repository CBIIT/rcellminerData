\name{drugData}
\alias{drugData}
\docType{data}
\title{
CellMiner compound data: NCI-60 activity and metadata
}
\description{
  The drugData object, of rcellminer package class DrugData, organizes data for 20861 compounds
  screened over the NCI-60 cancer cell lines by the NCI Developmental Therapeutics Program.
  These include activity profiles, as well as other information such as
  chemical structure, and where appropriate, name, mechanism of action,
  clinical testing status, etc.
}
\usage{data(drugData)}
\details{
  The activity data has two forms. A standardized 'z-score' data set provides
  a single NCI-60 profile for each compound, possibly averaging data derived
  from multiple experiments (subject to quality standards described in Reinhold et al.).
  A repeat activity data set provides data at the level of individual NCI-60 profiling
  experiments.
}
\source{
  CellMiner Website:
  <URL: http://discover.nci.nih.gov/cellminer/>
}
\value{
  the drugData object that contains Cellminer drug information as a DrugData object
}
\references{
  Reinhold, W.C., et al. (2012) CellMiner: a web-based suite of genomic and
  pharmacologic tools to explore transcript and drug patterns in the NCI-60 cell
  line set, Cancer Research, 72, 3499-3511
}
\concept{
  rcellminerData
}
\keyword{
  datasets
}
