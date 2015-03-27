\name{molData}
\alias{molData}
\docType{data}
\title{
CellMiner molecular profiling data for the NCI-60 cell lines.
}
\description{
  The molData object, of rcellminer package class MolData, organizes CellMiner's
  molecular profiling data for the NCI-60 cancer cell lines. Also provided are cell
  line annotations, including tissue of origin, detailed histology, and source information.
}
\usage{data(molData)}
\details{
Data matrices for the following types of data are provided (see examples
below for access details).

  \itemize{
    \item cop: Copy Number
    \item mut: Mutation
    \item exp: Expression (Z-Score)
    \item xai: Expression (Avg. log2 Intensity)
    \item pro: Protein (RPLA)
    \item mir: MicroRNA
    \item mda: Metadata (Phenotype Data)
    \item swa: Protein (SWATH-MS)
  }
}
\source{
  CellMiner Website:
  <URL: http://discover.nci.nih.gov/cellminer/>
}
\value{
  the molData object that contains Cellminer molecularp profiling information as a MolData object
}
\references{
  Reinhold, W.C., et al. (2012) CellMiner: a web-based suite of genomic and
  pharmacologic tools to explore transcript and drug patterns in the NCI-60 cell
  line set, Cancer research, 72, 3499-3511
}
\concept{
  rcellminerData
}
\keyword{
  datasets
}
